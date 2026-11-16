import { resolveMatchingRefName, sortYeastContigs } from './yeastAliases';
import type { ContigInfo, FastaIndexEntry } from './types';

export interface ReferenceSequenceSource {
  readonly label: string;
  getContigs(): ContigInfo[];
  resolveContigName(query: string, manualOverrides?: Record<string, string>): string | null;
  getSequence(chrom: string, start: number, end: number): Promise<string>;
}

interface ContigByteSpan {
  startOffset: number;
  endOffset: number;
  isUniform: boolean;
}

/**
 * .fai 텍스트 포맷 파싱
 */
export function parseFaiText(text: string): FastaIndexEntry[] {
  const entries: FastaIndexEntry[] = [];
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    const cols = line.split('\t');
    if (cols.length < 5) continue;
    entries.push({
      name: cols[0],
      length: Number(cols[1]),
      offset: Number(cols[2]),
      lineBases: Number(cols[3]),
      lineWidth: Number(cols[4]),
    });
  }
  return entries;
}

/**
 * 브라우저 File 객체에서 비압축 FASTA를 청크 단위로 스캔해 인메모리 .fai 인덱스 구축
 * - UTF-8 BOM 자동 무시
 * - 빈 줄/불규칙 줄 길이가 섞인 비표준 FASTA도 감지해 contig 바이트 범위 슬라이싱 폴백 지원
 */
async function buildFastaIndexFromFile(file: File): Promise<{
  entries: FastaIndexEntry[];
  contigSpans: Map<string, ContigByteSpan>;
}> {
  const CHUNK_SIZE = 4 * 1024 * 1024;
  const entries: FastaIndexEntry[] = [];
  const contigSpans = new Map<string, ContigByteSpan>();

  let currentName: string | null = null;
  let currentStartOffset = 0;
  let currentLength = 0;
  let lineBases = 0;
  let lineWidth = 0;
  let isUniform = true;
  let sawLastShortLine = false;
  let sawBlankLineInsideContig = false;

  let inHeader = false;
  let headerBytes: number[] = [];
  let lineBaseCount = 0;
  let lineByteCount = 0;
  let atLineStart = true;

  const decoder = new TextDecoder('ascii');
  let globalOffset = 0;

  const finalizeContig = (endByteOffset: number) => {
    if (currentName !== null) {
      entries.push({
        name: currentName,
        length: currentLength,
        offset: currentStartOffset,
        lineBases: lineBases || currentLength || 1,
        lineWidth: lineWidth || (lineBases || currentLength || 1) + 1,
      });
      contigSpans.set(currentName, {
        startOffset: currentStartOffset,
        endOffset: endByteOffset,
        isUniform: isUniform && !sawBlankLineInsideContig,
      });
    }
  };

  for (let pos = 0; pos < file.size; pos += CHUNK_SIZE) {
    const slice = file.slice(pos, Math.min(pos + CHUNK_SIZE, file.size));
    const buf = new Uint8Array(await slice.arrayBuffer());

    let startI = 0;
    // UTF-8 BOM 건너뛰기
    if (pos === 0 && buf.length >= 3 && buf[0] === 0xef && buf[1] === 0xbb && buf[2] === 0xbf) {
      startI = 3;
    }

    for (let i = startI; i < buf.length; i++) {
      const b = buf[i];
      const absPos = globalOffset + i;

      if (inHeader) {
        if (b === 10) {
          const headerText = decoder.decode(new Uint8Array(headerBytes));
          const token = headerText.trim().split(/\s+/)[0] || `contig_${entries.length + 1}`;
          currentName = token;
          currentStartOffset = absPos + 1;
          currentLength = 0;
          lineBases = 0;
          lineWidth = 0;
          isUniform = true;
          sawLastShortLine = false;
          sawBlankLineInsideContig = false;
          inHeader = false;
          headerBytes = [];
          atLineStart = true;
          lineBaseCount = 0;
          lineByteCount = 0;
        } else if (b !== 13) {
          headerBytes.push(b);
        }
        continue;
      }

      if (atLineStart && b === 62) {
        // '>'
        if (lineBaseCount > 0 && currentName !== null) {
          if (lineBases === 0) {
            lineBases = lineBaseCount;
            lineWidth = lineByteCount;
          }
          currentLength += lineBaseCount;
        }
        finalizeContig(absPos);
        inHeader = true;
        headerBytes = [];
        atLineStart = false;
        continue;
      }

      atLineStart = false;
      lineByteCount++;

      if (b === 10) {
        if (currentName !== null) {
          if (lineBaseCount === 0) {
            if (currentLength > 0) sawBlankLineInsideContig = true;
          } else {
            if (lineBases === 0) {
              lineBases = lineBaseCount;
              lineWidth = lineByteCount;
            } else {
              if (sawLastShortLine) {
                isUniform = false;
              }
              if (lineBaseCount !== lineBases || lineByteCount !== lineWidth) {
                if (lineBaseCount < lineBases) {
                  sawLastShortLine = true;
                } else {
                  isUniform = false;
                }
              }
            }
            currentLength += lineBaseCount;
          }
        }
        lineBaseCount = 0;
        lineByteCount = 0;
        atLineStart = true;
      } else if (b !== 13 && b !== 32 && b !== 9) {
        lineBaseCount++;
      }
    }

    globalOffset += buf.length;
  }

  if (inHeader && headerBytes.length > 0) {
    const headerText = decoder.decode(new Uint8Array(headerBytes));
    currentName = headerText.trim().split(/\s+/)[0] || `contig_${entries.length + 1}`;
    currentStartOffset = globalOffset;
  } else if (lineBaseCount > 0 && currentName !== null) {
    if (lineBases === 0) {
      lineBases = lineBaseCount;
      lineWidth = lineByteCount;
    }
    currentLength += lineBaseCount;
  }
  finalizeContig(globalOffset);

  return { entries, contigSpans };
}

async function verifyFaiMatchesFile(
  file: File,
  entries: FastaIndexEntry[],
): Promise<boolean> {
  if (entries.length === 0) return false;
  for (const e of entries) {
    if (
      !Number.isFinite(e.length) ||
      !Number.isFinite(e.offset) ||
      !Number.isFinite(e.lineBases) ||
      !Number.isFinite(e.lineWidth) ||
      e.length <= 0 ||
      e.offset < 1 ||
      e.lineBases <= 0 ||
      e.lineWidth < e.lineBases ||
      e.offset >= file.size
    ) {
      return false;
    }
  }

  // 첫 컨티그의 선두 바이트와 개행 폭이 실제 FASTA 파일과 일치하는지 확인 (CRLF 변환 등 불일치 감지)
  const first = entries[0];
  if (first.offset > 0) {
    const prevByte = new Uint8Array(await file.slice(first.offset - 1, first.offset).arrayBuffer())[0];
    if (prevByte !== 10 && prevByte !== 13) return false;
  }

  const probeBases = Math.min(first.length, first.lineBases);
  const probeText = await file.slice(first.offset, first.offset + first.lineWidth).text();
  if (probeText.length < probeBases) return false;

  for (let i = 0; i < probeBases; i++) {
    const ch = probeText[i];
    if (ch === '>' || ch === '\n' || ch === '\r') return false;
  }
  if (first.length > first.lineBases) {
    // 첫 줄 다음 문자들이 개행이어야 lineWidth가 파일과 정확히 일치
    const nlSlice = probeText.slice(first.lineBases);
    if (nlSlice.length > 0 && nlSlice[0] !== '\n' && nlSlice[0] !== '\r') {
      return false;
    }
  }
  return true;
}

export class IndexedFastaSource implements ReferenceSequenceSource {
  readonly label: string;
  private file: File;
  private entriesByName = new Map<string, FastaIndexEntry>();
  private contigSpans: Map<string, ContigByteSpan>;
  private irregularCache = new Map<string, string>();
  private contigs: ContigInfo[] = [];

  private constructor(
    file: File,
    entries: FastaIndexEntry[],
    contigSpans: Map<string, ContigByteSpan>,
  ) {
    this.file = file;
    this.label = file.name;
    this.contigSpans = contigSpans;
    const rawContigs: ContigInfo[] = [];
    for (const e of entries) {
      this.entriesByName.set(e.name, e);
      rawContigs.push({ name: e.name, length: e.length });
    }
    this.contigs = sortYeastContigs(rawContigs);
  }

  static async fromFiles(fastaFile: File, faiFile?: File): Promise<ReferenceSequenceSource> {
    if (fastaFile.name.toLowerCase().endsWith('.gz')) {
      return MemoryFastaSource.fromGzipFile(fastaFile);
    }

    let entries: FastaIndexEntry[] = [];
    let contigSpans = new Map<string, ContigByteSpan>();
    let usedExternalFai = false;

    if (faiFile) {
      const faiText = await faiFile.text();
      const parsed = parseFaiText(faiText);
      if (await verifyFaiMatchesFile(fastaFile, parsed)) {
        entries = parsed;
        usedExternalFai = true;
        for (let i = 0; i < entries.length; i++) {
          const e = entries[i];
          const nextOffset = i + 1 < entries.length ? entries[i + 1].offset : fastaFile.size;
          contigSpans.set(e.name, { startOffset: e.offset, endOffset: nextOffset, isUniform: true });
        }
      }
    }

    if (!usedExternalFai) {
      const built = await buildFastaIndexFromFile(fastaFile);
      entries = built.entries;
      contigSpans = built.contigSpans;
    }

    if (entries.length === 0) {
      throw new Error(`유효한 FASTA 서열(>헤더)을 찾지 못했습니다: ${fastaFile.name}`);
    }

    return new IndexedFastaSource(fastaFile, entries, contigSpans);
  }

  getContigs(): ContigInfo[] {
    return this.contigs;
  }

  resolveContigName(query: string, manualOverrides?: Record<string, string>): string | null {
    return resolveMatchingRefName(
      query,
      this.contigs.map((c) => c.name),
      manualOverrides,
      undefined,
      this.contigs,
    );
  }

  async getSequence(chrom: string, start: number, end: number): Promise<string> {
    const entry = this.entriesByName.get(chrom);
    if (!entry) return '';

    const clampedStart = Math.max(0, start);
    const clampedEnd = Math.min(entry.length, end);
    if (clampedEnd <= clampedStart) return '';

    const spanMeta = this.contigSpans.get(chrom);
    if (spanMeta && !spanMeta.isUniform) {
      let cached = this.irregularCache.get(chrom);
      if (cached === undefined) {
        const raw = await this.file.slice(spanMeta.startOffset, spanMeta.endOffset).text();
        cached = raw.replace(/[\r\n\s]/g, '').toUpperCase();
        this.irregularCache.set(chrom, cached);
      }
      return cached.slice(clampedStart, clampedEnd);
    }

    const startLine = Math.floor(clampedStart / entry.lineBases);
    const startCol = clampedStart % entry.lineBases;
    const byteStart = entry.offset + startLine * entry.lineWidth + startCol;

    const endBaseIdx = clampedEnd - 1;
    const endLine = Math.floor(endBaseIdx / entry.lineBases);
    const endCol = endBaseIdx % entry.lineBases;
    const byteEnd = entry.offset + endLine * entry.lineWidth + endCol + 1;

    const rawText = await this.file.slice(byteStart, byteEnd).text();
    return rawText.replace(/[\r\n\s]/g, '').toUpperCase();
  }
}

export class MemoryFastaSource implements ReferenceSequenceSource {
  readonly label: string;
  private seqs: Map<string, string>;
  private contigs: ContigInfo[];

  constructor(label: string, sequences: Map<string, string>) {
    this.label = label;
    this.seqs = new Map();
    const rawContigs: ContigInfo[] = [];
    for (const [name, seq] of sequences.entries()) {
      const upper = seq.toUpperCase();
      this.seqs.set(name, upper);
      rawContigs.push({ name, length: upper.length });
    }
    this.contigs = sortYeastContigs(rawContigs);
  }

  static async fromGzipUrl(label: string, url: string): Promise<MemoryFastaSource> {
    const resp = await fetch(url);
    if (!resp.ok || !resp.body) {
      throw new Error(`내장 지놈(${label}) 로드 실패: HTTP ${resp.status}`);
    }
    const stream = resp.body.pipeThrough(new DecompressionStream('gzip'));
    const text = await new Response(stream).text();
    return MemoryFastaSource.fromFastaText(label, text);
  }

  static async fromGzipFile(file: File): Promise<MemoryFastaSource> {
    const stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
    const text = await new Response(stream).text();
    return MemoryFastaSource.fromFastaText(file.name, text);
  }

  static fromFastaText(label: string, text: string): MemoryFastaSource {
    const clean = text.charCodeAt(0) === 0xfeff ? text.slice(1) : text;
    const map = new Map<string, string>();
    let currentName: string | null = null;
    let chunks: string[] = [];

    for (const rawLine of clean.split(/\r?\n/)) {
      const line = rawLine.trim();
      if (!line) continue;
      if (line.startsWith('>')) {
        if (currentName !== null) {
          map.set(currentName, chunks.join(''));
        }
        currentName = line.slice(1).trim().split(/\s+/)[0] || `contig_${map.size + 1}`;
        chunks = [];
      } else if (currentName !== null) {
        chunks.push(line);
      }
    }
    if (currentName !== null) {
      map.set(currentName, chunks.join(''));
    }
    return new MemoryFastaSource(label, map);
  }

  getContigs(): ContigInfo[] {
    return this.contigs;
  }

  resolveContigName(query: string, manualOverrides?: Record<string, string>): string | null {
    return resolveMatchingRefName(
      query,
      this.contigs.map((c) => c.name),
      manualOverrides,
      undefined,
      this.contigs,
    );
  }

  async getSequence(chrom: string, start: number, end: number): Promise<string> {
    const seq = this.seqs.get(chrom);
    if (!seq) return '';
    return seq.slice(Math.max(0, start), Math.min(seq.length, end));
  }
}
