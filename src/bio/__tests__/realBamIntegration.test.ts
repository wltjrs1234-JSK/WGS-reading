import zlib from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { GmodBamSource } from '../bamSource';
import { buildThreeWayAnalysis } from '../consensus';
import { IndexedFastaSource } from '../fastaIndex';
import { DEFAULT_CONSENSUS_OPTIONS } from '../types';

function crc32(buf: Uint8Array): number {
  let c = 0xffffffff;
  for (let i = 0; i < buf.length; i++) {
    c ^= buf[i];
    for (let k = 0; k < 8; k++) {
      c = (c >>> 1) ^ (c & 1 ? 0xedb88320 : 0);
    }
  }
  return (c ^ 0xffffffff) >>> 0;
}

function makeBgzfBlock(uncompressed: Uint8Array): Uint8Array {
  const deflated = zlib.deflateRawSync(uncompressed);
  const blockSize = 18 + deflated.length + 8;
  const out = new Uint8Array(blockSize);
  const dv = new DataView(out.buffer);

  // GZIP + BGZF extra header
  out[0] = 0x1f;
  out[1] = 0x8b;
  out[2] = 0x08;
  out[3] = 0x04; // FEXTRA
  out[9] = 0xff;
  dv.setUint16(10, 6, true); // XLEN=6
  out[12] = 66; // 'B'
  out[13] = 67; // 'C'
  dv.setUint16(14, 2, true); // SLEN=2
  dv.setUint16(16, blockSize - 1, true); // BSIZE

  out.set(deflated, 18);
  const trailerPos = 18 + deflated.length;
  dv.setUint32(trailerPos, crc32(uncompressed), true);
  dv.setUint32(trailerPos + 4, uncompressed.length >>> 0, true);
  return out;
}

const BGZF_EOF = new Uint8Array([
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
  0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
]);

const SEQ_NIBBLE: Record<string, number> = {
  '=': 0,
  A: 1,
  C: 2,
  M: 3,
  G: 4,
  R: 5,
  S: 6,
  V: 7,
  T: 8,
  W: 9,
  Y: 10,
  H: 11,
  K: 12,
  D: 13,
  B: 14,
  N: 15,
};

function encodeBamRecord(opts: {
  refId: number;
  pos: number;
  mapq: number;
  flags: number;
  name: string;
  cigarOps: Array<{ opCode: number; len: number }>;
  seq: string;
}): Uint8Array {
  const nameBytes = new TextEncoder().encode(opts.name + '\0');
  const seqLen = opts.seq.length;
  const seqByteLen = (seqLen + 1) >> 1;
  const bodySize =
    32 + nameBytes.length + opts.cigarOps.length * 4 + seqByteLen + seqLen;

  const buf = new Uint8Array(4 + bodySize);
  const dv = new DataView(buf.buffer);
  dv.setInt32(0, bodySize, true);
  dv.setInt32(4, opts.refId, true);
  dv.setInt32(8, opts.pos, true);
  dv.setUint32(
    12,
    ((4681 << 16) | ((opts.mapq & 0xff) << 8) | (nameBytes.length & 0xff)) >>> 0,
    true,
  );
  dv.setUint32(16, (((opts.flags & 0xffff) << 16) | (opts.cigarOps.length & 0xffff)) >>> 0, true);
  dv.setInt32(20, seqLen, true);
  dv.setInt32(24, -1, true); // next_refID
  dv.setInt32(28, -1, true); // next_pos
  dv.setInt32(32, 0, true);  // tlen

  let p = 36;
  buf.set(nameBytes, p);
  p += nameBytes.length;

  for (const c of opts.cigarOps) {
    dv.setUint32(p, ((c.len << 4) | (c.opCode & 0xf)) >>> 0, true);
    p += 4;
  }

  for (let i = 0; i < seqLen; i += 2) {
    const hi = SEQ_NIBBLE[opts.seq[i]] ?? 15;
    const lo = i + 1 < seqLen ? SEQ_NIBBLE[opts.seq[i + 1]] ?? 15 : 0;
    buf[p++] = (hi << 4) | lo;
  }

  for (let i = 0; i < seqLen; i++) {
    buf[p++] = 35; // Phred Q35
  }

  return buf;
}

function buildBinaryBamAndBai(
  refName: string,
  refLen: number,
  records: Array<{
    name: string;
    pos: number;
    seq: string;
    cigarOps: Array<{ opCode: number; len: number }>;
  }>,
): { bamBytes: Uint8Array; baiBytes: Uint8Array } {
  const samHeader = `@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:${refName}\tLN:${refLen}\n`;
  const headerTextBytes = new TextEncoder().encode(samHeader);
  const refNameBytes = new TextEncoder().encode(refName + '\0');

  const headerBuf = new Uint8Array(
    4 + 4 + headerTextBytes.length + 4 + 4 + refNameBytes.length + 4,
  );
  const hdv = new DataView(headerBuf.buffer);
  headerBuf.set([66, 65, 77, 1], 0); // BAM\x01
  hdv.setInt32(4, headerTextBytes.length, true);
  headerBuf.set(headerTextBytes, 8);
  let hp = 8 + headerTextBytes.length;
  hdv.setInt32(hp, 1, true); // n_ref = 1
  hp += 4;
  hdv.setInt32(hp, refNameBytes.length, true);
  hp += 4;
  headerBuf.set(refNameBytes, hp);
  hp += refNameBytes.length;
  hdv.setInt32(hp, refLen, true);

  const recChunks = records.map((r) =>
    encodeBamRecord({
      refId: 0,
      pos: r.pos,
      mapq: 60,
      flags: 0,
      name: r.name,
      cigarOps: r.cigarOps,
      seq: r.seq,
    }),
  );

  const totalRecLen = recChunks.reduce((acc, b) => acc + b.length, 0);
  const allRecBytes = new Uint8Array(totalRecLen);
  let rp = 0;
  for (const chunk of recChunks) {
    allRecBytes.set(chunk, rp);
    rp += chunk.length;
  }

  const block0 = makeBgzfBlock(headerBuf);
  const block1 = makeBgzfBlock(allRecBytes);

  const bamBytes = new Uint8Array(block0.length + block1.length + BGZF_EOF.length);
  bamBytes.set(block0, 0);
  bamBytes.set(block1, block0.length);
  bamBytes.set(BGZF_EOF, block0.length + block1.length);

  // Virtual offsets: chunkBeg = (block0.length << 16) | 0, chunkEnd = ((block0.length + block1.length) << 16) | 0
  const chunkBeg = BigInt(block0.length) << 16n;
  const chunkEnd = BigInt(block0.length + block1.length) << 16n;

  const baiBytes = new Uint8Array(4 + 4 + 4 + 4 + 4 + 8 + 8 + 4 + 8);
  const bdv = new DataView(baiBytes.buffer);
  baiBytes.set([66, 65, 73, 1], 0); // BAI\x01
  bdv.setInt32(4, 1, true); // n_ref = 1
  bdv.setInt32(8, 1, true); // n_bin = 1
  bdv.setUint32(12, 4681, true); // bin = 4681 (covers [0..16384))
  bdv.setInt32(16, 1, true); // n_chunk = 1
  bdv.setBigUint64(20, chunkBeg, true);
  bdv.setBigUint64(28, chunkEnd, true);
  bdv.setInt32(36, 1, true); // n_intv = 1
  bdv.setBigUint64(40, chunkBeg, true);

  return { bamBytes, baiBytes };
}

describe('Real binary BAM + BAI end-to-end integration (@gmod/bam + BlobFile)', () => {
  it('parses binary BAM/BAI Files via BlobFile and reconstructs 3-way consensus sequences', async () => {
    const refSeq = 'ATGCGATCGATCGATC'; // 16 bp
    const fastaFile = new File([`>chrV\n${refSeq}\n`], 'S288C_chrV.fasta');
    const refSource = await IndexedFastaSource.fromFiles(fastaFile);

    // Strain 1: 6 reads carrying SNP (pos 2 G->A) + 2bp Insertion ('TT' after pos 5) -> 6M2I10M
    const s1Built = buildBinaryBamAndBai(
      'chrV',
      16,
      Array.from({ length: 6 }, (_, i) => ({
        name: `s1_read_${i}`,
        pos: 0,
        seq: 'ATACGATTTCGATCGATC',
        cigarOps: [
          { opCode: 0, len: 6 }, // 6M
          { opCode: 1, len: 2 }, // 2I
          { opCode: 0, len: 10 }, // 10M
        ],
      })),
    );

    // Strain 2 (with alias contig 'V' instead of 'chrV'): 6 reads carrying 3bp Deletion at pos 8..11 -> 8M3D5M
    const s2Built = buildBinaryBamAndBai(
      'V',
      16,
      Array.from({ length: 6 }, (_, i) => ({
        name: `s2_read_${i}`,
        pos: 0,
        seq: 'ATGCGATCCGATC',
        cigarOps: [
          { opCode: 0, len: 8 }, // 8M
          { opCode: 2, len: 3 }, // 3D
          { opCode: 0, len: 5 }, // 5M
        ],
      })),
    );

    const s1BamFile = new File([s1Built.bamBytes], 'strain1.bam');
    const s1BaiFile = new File([s1Built.baiBytes], 'strain1.bam.bai');
    const s2BamFile = new File([s2Built.bamBytes], 'strain2.bam');
    const s2BaiFile = new File([s2Built.baiBytes], 'strain2.bam.bai');

    const strain1 = new GmodBamSource('s1', 'Strain 1', s1BamFile, s1BaiFile);
    const strain2 = new GmodBamSource('s2', 'Strain 2', s2BamFile, s2BaiFile);

    const s1Chrom = await strain1.resolveContigName('chrV');
    const s2Chrom = await strain2.resolveContigName('chrV');
    expect(s1Chrom).toBe('chrV');
    expect(s2Chrom).toBe('V');

    const [loadedRef, s1Reads, s2Reads] = await Promise.all([
      refSource.getSequence('chrV', 0, 16),
      strain1.getRecords(s1Chrom!, 0, 16),
      strain2.getRecords(s2Chrom!, 0, 16),
    ]);

    expect(s1Reads).toHaveLength(6);
    expect(s1Reads[0].cigar).toBe('6M2I10M');
    expect(s2Reads).toHaveLength(6);
    expect(s2Reads[0].cigar).toBe('8M3D5M');

    const analysis = buildThreeWayAnalysis(
      'chrV',
      0,
      16,
      loadedRef,
      s1Reads,
      s2Reads,
      DEFAULT_CONSENSUS_OPTIONS,
    );

    expect(analysis.ungappedRef).toBe('ATGCGATCGATCGATC');
    expect(analysis.ungappedS1).toBe('ATACGATTTCGATCGATC');
    expect(analysis.ungappedS2).toBe('ATGCGATCCGATC');
    expect(analysis.alignedRef).toBe('ATGCGA--TCGATCGATC');
    expect(analysis.alignedS1).toBe('ATACGATTTCGATCGATC');
    expect(analysis.alignedS2).toBe('ATGCGA--TC---CGATC');
  });
});

