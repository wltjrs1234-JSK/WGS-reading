import { MemoryAlignmentSource } from './bamSource';
import { MemoryFastaSource } from './fastaIndex';
import { AnnotationIndex, parseGff3Text } from './gff3Parser';
import type { AlignedRead, ContigInfo } from './types';

function seededRand(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

function makeDeterministicSeq(length: number, seed: number): string {
  const rand = seededRand(seed);
  const bases = ['A', 'C', 'G', 'T'];
  const chars = new Array<string>(length);
  for (let i = 0; i < length; i++) {
    chars[i] = bases[Math.floor(rand() * 4)];
  }
  return chars.join('');
}

function replaceAt(seq: string, pos: number, sub: string): string {
  return seq.slice(0, pos) + sub + seq.slice(pos + sub.length);
}

interface SyntheticEdit {
  pos: number;
  kind: 'SNP' | 'INS' | 'DEL';
  alt: string;   // for SNP: single base, for INS: inserted string, for DEL: unused
  delLen?: number;
}

/**
 * 레퍼런스 구간을 기준으로 실제 CIGAR 문자열과 염기서열을 갖춘 120bp Paired-end 스타일 리드 생성
 */
function generateReadsForContig(
  refSeq: string,
  edits: SyntheticEdit[],
  gapRange: [number, number] | null,
  coverage: number,
  seed: number,
  prefix: string,
): AlignedRead[] {
  const rand = seededRand(seed);
  const readLen = 120;
  const step = Math.max(4, Math.floor((readLen * 1) / coverage));
  const reads: AlignedRead[] = [];

  for (let start = 0; start + 40 < refSeq.length; start += step) {
    const jitter = Math.floor(rand() * 5);
    const rStart = Math.min(refSeq.length - 40, start + jitter);
    const targetRefSpan = Math.min(readLen, refSeq.length - rStart);
    const rEnd = rStart + targetRefSpan;

    // gapRange 내부와 완전히 겹치면 커버리지 0x(유전자 결실) 시뮬레이션
    if (gapRange && rStart >= gapRange[0] && rEnd <= gapRange[1]) {
      continue;
    }

    let refCursor = rStart;
    const seqParts: string[] = [];
    const cigarOps: Array<{ op: string; len: number }> = [];

    const pushOp = (op: string, len: number) => {
      if (len <= 0) return;
      const last = cigarOps[cigarOps.length - 1];
      if (last && last.op === op) {
        last.len += len;
      } else {
        cigarOps.push({ op, len });
      }
    };

    while (refCursor < rEnd) {
      // gap 경계 클리핑
      if (gapRange && refCursor >= gapRange[0] && refCursor < gapRange[1]) {
        break;
      }

      const edit = edits.find((e) => e.pos === refCursor);
      if (edit && edit.kind === 'DEL' && edit.delLen) {
        const dLen = Math.min(edit.delLen, rEnd - refCursor);
        pushOp('D', dLen);
        refCursor += dLen;
        continue;
      }

      let base = refSeq[refCursor];
      if (edit && edit.kind === 'SNP') {
        base = edit.alt;
      } else if (rand() < 0.002) {
        // 현실적인 0.2% 시퀀싱 에러 노이즈
        const altPool = ['A', 'C', 'G', 'T'].filter((b) => b !== base);
        base = altPool[Math.floor(rand() * 3)];
      }
      seqParts.push(base);
      pushOp('M', 1);

      if (edit && edit.kind === 'INS') {
        seqParts.push(edit.alt);
        pushOp('I', edit.alt.length);
      }
      refCursor++;
    }

    if (seqParts.length < 30) continue;

    const actualEnd = refCursor;
    const strand = rand() > 0.5 ? '+' : '-';
    const cigar = cigarOps.map((o) => `${o.len}${o.op}`).join('');
    const seq = seqParts.join('');
    const qual = new Uint8Array(seq.length).fill(35);

    reads.push({
      id: `${prefix}_${reads.length}`,
      name: `${prefix}:read_${reads.length}`,
      start: rStart,
      end: actualEnd,
      strand,
      seq,
      cigar,
      mapq: 60,
      qual,
      flags: strand === '-' ? 0x10 : 0x0,
    });
  }

  return reads;
}

import { YEAST_CHROMOSOMES } from './yeastAliases';

export interface DemoDataset {
  reference: MemoryFastaSource;
  strain1: MemoryAlignmentSource;
  strain2: MemoryAlignmentSource;
  annotation: AnnotationIndex;
  defaultChrom: string;
  defaultStart: number;
  defaultEnd: number;
}

/**
 * S. cerevisiae 데모 데이터셋 생성 (전 염색체 chrI~chrXVI + chrM 포함):
 * - chrV (4,000 bp 모델): URA3 유전자(1,000..1,804, +) 및 주변부
 *   - 균주 1 (BY4741-mut): Missense SNP 2개 + 3bp Insertion
 *   - 균주 2 (Industrial-var): Synonymous SNP + 6bp In-frame Deletion + 3' 부분 결실(Gap)
 * - chrM (2,500 bp 모델): COX1 유전자(500..740, +) — 미토콘드리아 Table 3 검증용
 * - 나머지 chrI~chrXVI: 각 염색체별 대표 데모 구간(1,500 bp) 포함해 전체 17개 염색체 표시
 */
export function createYeastDemoDataset(): DemoDataset {
  let chrVSeq = makeDeterministicSeq(4000, 202611);
  const ura3Core =
    'ATGTCGAAAGCTACATATAAGGAACGTGCTGCTACTCATCCTAGTCCTGTTGCTGCCAAGCTATTTAATATCATGCACGAAAAGCAAACAAACTTGTGTGCTTCATTGGATGTTCGTACCACCAAGGAATTACTGGAGTTAGTTGAAGCATTAGGTCCCAAAATTTGTTTACTAAAAACACATGTGGATATCTTGACTGATTTTTCCATGGAGGGCACAGTTAAGCCGCTAAAGGCATTATCCGCCAAGTACAATTTTTTACTCTTCGAAGACAGAAAATTTGCTGACATTGGTAATACAGTCAAATTGCAGTACTCTGCGGGTGTATACAGAATAGCAGAATGGGCAGACATTACGAATGCACACGGTGTGGTGGGCCCAGGTATTGTTAGCGGTTTGAAGCAGGCGGCAGAAGAAGTAACAAAGGAACCTAGAGGCCTTTTGATGTTAGCAGAATTGTCATGCAAGGGCTCCCTAGCTACTGGAGAATATACTAAGGGTACTGTTGACATTGCGAAGAGCGACAAAGATTTTGTTATCGGCTTTATTGCTCAAAGAGACATGGGTGGAAGAGATGAAGGTTACGATTGGTTGATTATGACACCCGGTGTGGGTTTAGATGACAAGGGAGACGCATTGGGTCAACAGTATAGAACCGTGGATGATGTGGTCTCTACAGGATCTGACATTATTATTGTTGGAAGAGGACTATTTGCAAAGGGAAGGGATGCTAAGGTAGAGGGTGAACGTTACAGAAAAGCAGGCTGGGAAGCATATTTGAGAAGATGCGGCCAGCAAAACTAA';
  chrVSeq = replaceAt(chrVSeq, 1000, ura3Core);

  let chrMSeq = makeDeterministicSeq(2500, 991155);
  const mitoOrf =
    'ATGGTACAATGATTATATTCAACTAATGCTAAAGATATAGGAGTATTATACTTAATATTTGCATTATTTTCTGGAATAGTAGGTACTGCAATATGAAGTATAATTAGAATAGAATTAGCACCAGCTGGTTCACAATACTTAGGTACTAATCCTCAATTATTTAATTCATTAATTGTTACTAGTCATGCAATATTAATAATATTTTTCTTAGTTATACCTGTTATAATTGGTGGATTTGGAAATTGA';
  chrMSeq = replaceAt(chrMSeq, 500, mitoOrf + 'TAA');

  const refMap = new Map<string, string>();
  const s1ReadsMap = new Map<string, AlignedRead[]>();
  const s2ReadsMap = new Map<string, AlignedRead[]>();
  const s1Contigs: ContigInfo[] = [];
  const s2Contigs: ContigInfo[] = [];

  const s1ChrVEdits: SyntheticEdit[] = [
    { pos: 1030, kind: 'SNP', alt: chrVSeq[1030] === 'C' ? 'G' : 'C' },
    { pos: 1110, kind: 'INS', alt: 'GGT' },
    { pos: 1215, kind: 'SNP', alt: chrVSeq[1215] === 'G' ? 'T' : 'G' },
  ];
  const s2ChrVEdits: SyntheticEdit[] = [
    { pos: 1065, kind: 'DEL', alt: '', delLen: 6 },
    { pos: 1215, kind: 'SNP', alt: chrVSeq[1215] === 'G' ? 'T' : 'G' },
    { pos: 1280, kind: 'SNP', alt: chrVSeq[1280] === 'A' ? 'G' : 'A' },
  ];
  const s1ChrMEdits: SyntheticEdit[] = [{ pos: 511, kind: 'SNP', alt: 'G' }];
  const s2ChrMEdits: SyntheticEdit[] = [{ pos: 525, kind: 'SNP', alt: 'G' }];

  for (let i = 0; i < YEAST_CHROMOSOMES.length; i++) {
    const def = YEAST_CHROMOSOMES[i];
    if (def.canonical === '2-micron') continue;

    if (def.canonical === 'chrV') {
      refMap.set('chrV', chrVSeq);
      s1Contigs.push({ name: 'chrV', length: chrVSeq.length });
      s2Contigs.push({ name: 'V', length: chrVSeq.length });
      s1ReadsMap.set('chrV', generateReadsForContig(chrVSeq, s1ChrVEdits, null, 28, 101, 'S1'));
      s2ReadsMap.set('V', generateReadsForContig(chrVSeq, s2ChrVEdits, [1600, 1850], 28, 201, 'S2'));
    } else if (def.canonical === 'chrM') {
      refMap.set('chrM', chrMSeq);
      s1Contigs.push({ name: 'chrM', length: chrMSeq.length });
      s2Contigs.push({ name: 'MT', length: chrMSeq.length });
      s1ReadsMap.set('chrM', generateReadsForContig(chrMSeq, s1ChrMEdits, null, 24, 102, 'S1M'));
      s2ReadsMap.set('MT', generateReadsForContig(chrMSeq, s2ChrMEdits, null, 24, 202, 'S2M'));
    } else {
      const seq = makeDeterministicSeq(1500, 3000 + i * 37);
      const snpPos = 420 + i * 15;
      const altBase = seq[snpPos] === 'A' ? 'G' : 'A';
      refMap.set(def.canonical, seq);
      s1Contigs.push({ name: def.canonical, length: seq.length });
      s2Contigs.push({ name: def.roman, length: seq.length });
      s1ReadsMap.set(
        def.canonical,
        generateReadsForContig(seq, [{ pos: snpPos, kind: 'SNP', alt: altBase }], null, 20, 400 + i, `S1_${def.roman}`),
      );
      s2ReadsMap.set(
        def.roman,
        generateReadsForContig(seq, [], null, 20, 700 + i, `S2_${def.roman}`),
      );
    }
  }

  const reference = new MemoryFastaSource('S288C_Demo_17Chr.fasta', refMap);
  const strain1 = new MemoryAlignmentSource(
    'demo_s1',
    '균주 1 (BY4741-mut.bam)',
    s1Contigs,
    s1ReadsMap,
  );
  const strain2 = new MemoryAlignmentSource(
    'demo_s2',
    '균주 2 (Industrial-var.bam)',
    s2Contigs,
    s2ReadsMap,
  );

  const gffText = `##gff-version 3
chrV\tSGD\tgene\t1001\t1804\t.\t+\t.\tID=YEL021W;Name=YEL021W;gene=URA3;Note=Orotidine-5'-phosphate%20decarboxylase%20(Demo%20locus)
chrV\tSGD\tCDS\t1001\t1804\t.\t+\t0\tParent=YEL021W
chrV\tSGD\tARS\t420\t580\t.\t.\t.\tID=ARS507;Name=ARS507;Note=Autonomously%20Replicating%20Sequence
chrV\tSGD\tgene\t2150\t2950\t.\t-\t.\tID=YEL020C;Name=YEL020C;gene=TIM9;Note=Mitochondrial%20intermembrane%20space%20protein
chrV\tSGD\tCDS\t2150\t2950\t.\t-\t0\tParent=YEL020C
chrM\tSGD\tgene\t501\t740\t.\t+\t.\tID=Q0045;Name=Q0045;gene=COX1;Note=Subunit%20I%20of%20cytochrome%20c%20oxidase%20(Mito%20Table%203)
chrM\tSGD\tCDS\t501\t740\t.\t+\t0\tParent=Q0045
`;

  const annotation = new AnnotationIndex(parseGff3Text(gffText));

  return {
    reference,
    strain1,
    strain2,
    annotation,
    defaultChrom: 'chrV',
    defaultStart: 950,
    defaultEnd: 1380,
  };
}

