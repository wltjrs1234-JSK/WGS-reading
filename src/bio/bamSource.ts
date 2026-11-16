import { BamFile } from '@gmod/bam';
import { BlobFile } from 'generic-filehandle2';
import { resolveMatchingRefName } from './yeastAliases';
import type { AlignedRead, ContigInfo } from './types';

export interface AlignmentSource {
  readonly id: string;
  readonly label: string;
  getContigs(): Promise<ContigInfo[]>;
  resolveContigName(
    query: string,
    manualOverrides?: Record<string, string>,
    targetLength?: number,
  ): Promise<string | null>;
  getRecords(chrom: string, start: number, end: number, maxRecords?: number): Promise<AlignedRead[]>;
  getOverviewCoverage(chrom: string, start: number, end: number, numBins: number): Promise<number[]>;
}

/**
 * 로컬 PC의 BAM + BAI(또는 CSI) 파일을 @gmod/bam + BlobFile(File.slice)로 읽는 소스
 * 외부 서버 업로드 없이 브라우저 샌드박스 내에서 바이트 단위 랜덤 액세스로 동작합니다.
 */
const MAX_BAM_CACHE_BYTES = 64 * 1024 * 1024;
const BAM_CACHE_IDLE_MS = 60_000;
const MAX_RECORDS_PER_WINDOW = 25_000;

export class GmodBamSource implements AlignmentSource {
  readonly id: string;
  readonly label: string;
  private bam: BamFile;
  private headerLoaded = false;

  constructor(id: string, label: string, bamFile: File, indexFile: File) {
    this.id = id;
    this.label = label;
    const isCsi = indexFile.name.toLowerCase().endsWith('.csi');
    this.bam = new BamFile({
      bamFilehandle: new BlobFile(bamFile),
      maxCacheBytes: MAX_BAM_CACHE_BYTES,
      cacheIdleTimeoutMs: BAM_CACHE_IDLE_MS,
      ...(isCsi
        ? { csiFilehandle: new BlobFile(indexFile) }
        : { baiFilehandle: new BlobFile(indexFile) }),
    });
  }

  private async ensureHeader(): Promise<void> {
    if (!this.headerLoaded) {
      await this.bam.getHeader();
      this.headerLoaded = true;
    }
  }

  async getContigs(): Promise<ContigInfo[]> {
    await this.ensureHeader();
    return (this.bam.indexToChr ?? []).map((c) => ({
      name: c.refName,
      length: c.length,
    }));
  }

  async resolveContigName(
    query: string,
    manualOverrides?: Record<string, string>,
    targetLength?: number,
  ): Promise<string | null> {
    const contigs = await this.getContigs();
    return resolveMatchingRefName(
      query,
      contigs.map((c) => c.name),
      manualOverrides,
      targetLength,
      contigs,
    );
  }

  async getRecords(
    chrom: string,
    start: number,
    end: number,
    maxRecords = MAX_RECORDS_PER_WINDOW,
  ): Promise<AlignedRead[]> {
    await this.ensureHeader();
    const records = await this.bam.getRecordsForRange(chrom, Math.max(0, start), Math.max(1, end));

    // flags만 먼저 선별해 unmapped/secondary/QC-fail/duplicate를 제외
    const validIndices: number[] = [];
    for (let idx = 0; idx < records.length; idx++) {
      if ((records[idx].flags & 0x704) === 0) {
        validIndices.push(idx);
      }
    }

    // 플라스미드(URA3)·rDNA·chrM 등 초고심도(수만 x) 구간에서 브라우저 OOM/프리즈를 막기 위해 균일 다운샘플링
    const limit = Number.isFinite(maxRecords) && maxRecords > 0 ? maxRecords : validIndices.length;
    const step = validIndices.length > limit ? validIndices.length / limit : 1;
    const targetLen = step > 1 ? limit : validIndices.length;
    const out: AlignedRead[] = new Array(targetLen);

    for (let k = 0; k < targetLen; k++) {
      const idx = validIndices[step > 1 ? Math.min(validIndices.length - 1, Math.floor(k * step)) : k];
      const r = records[idx];
      out[k] = {
        id: `${r.name}_${r.start}_${r.flags}_${idx}`,
        name: r.name,
        start: r.start,
        end: r.end,
        strand: r.strand === -1 ? '-' : '+',
        seq: r.seq,
        cigar: r.CIGAR,
        mapq: r.mq ?? 60,
        qual: r.qual ?? undefined,
        flags: r.flags,
      };
    }
    return out;
  }

  async getOverviewCoverage(
    chrom: string,
    start: number,
    end: number,
    numBins: number,
  ): Promise<number[]> {
    await this.ensureHeader();
    const bins = new Array<number>(numBins).fill(0);
    const span = Math.max(1, end - start);

    try {
      const entries = await this.bam.indexCov(chrom, start, end);
      let hasPositive = false;
      if (entries.length > 0) {
        const counts = new Array<number>(numBins).fill(0);
        for (const e of entries) {
          const score = Number.isFinite(e.score) ? e.score : 0;
          if (score > 0) hasPositive = true;
          const b0 = Math.max(0, Math.floor(((e.start - start) / span) * numBins));
          const b1 = Math.min(numBins - 1, Math.floor(((e.end - start) / span) * numBins));
          for (let b = b0; b <= b1; b++) {
            bins[b] += score;
            counts[b]++;
          }
        }
        if (hasPositive) {
          for (let b = 0; b < numBins; b++) {
            if (counts[b] > 0) bins[b] /= counts[b];
          }
          return bins;
        }
      }
    } catch {
      // fallback below if indexCov is unsupported or lacks pseudo-bin stats
    }

    // 일부 samtools/picard BAI 인덱스에는 pseudo-bin stats(lineCount)가 빠져 있어 indexCov가 0을 반환함.
    // 1차로 이미 파싱된 BAI 선형 인덱스 오프셋(16kb 간격)을 직접 읽어 I/O 없이 즉시 밀도를 산출합니다.
    try {
      const bamAny = this.bam as unknown as {
        getSeqId?: (c: string) => Promise<number | undefined>;
        index?: {
          parse?: () => Promise<{
            indices?: (seqId: number) => { linearBlockPositions?: Float64Array } | undefined;
          }>;
        };
      };
      const seqId = await bamAny.getSeqId?.(chrom);
      const parsedIdx = seqId !== undefined ? await bamAny.index?.parse?.() : undefined;
      const linearPos =
        seqId !== undefined ? parsedIdx?.indices?.(seqId)?.linearBlockPositions : undefined;

      if (linearPos && linearPos.length > 1) {
        const linearInterval = 16384;
        const binWidth = span / numBins;
        for (let b = 0; b < numBins; b++) {
          const bs = start + b * binWidth;
          const be = start + (b + 1) * binWidth;
          const i0 = Math.max(0, Math.min(linearPos.length - 2, Math.floor(bs / linearInterval)));
          const i1 = Math.max(i0 + 1, Math.min(linearPos.length - 1, Math.ceil(be / linearInterval)));
          const delta = Math.max(0, linearPos[i1] - linearPos[i0]);
          const bpSpan = Math.max(1, (i1 - i0) * linearInterval);
          bins[b] = (delta / bpSpan) * 16; // KB-normalized 상대 깊이
        }
        if (bins.some((v) => v > 0)) {
          return bins;
        }
      }

      const binWidth = span / numBins;
      const chunks = await this.bam.blocksForRange(chrom, start, end);
      if (chunks.length > 0) {
        const step = Math.max(1, Math.floor(numBins / 12));
        const buckets: number[] = [];
        for (let b = 0; b < numBins; b += step) {
          buckets.push(b);
        }
        await Promise.all(
          buckets.map(async (b) => {
            const bs = Math.floor(start + b * binWidth);
            const be = Math.min(end, Math.ceil(start + (b + step) * binWidth));
            const subChunks = await this.bam.blocksForRange(chrom, bs, be);
            let byteMass = 0;
            for (const ch of subChunks) {
              const delta = ch.maxv.blockPosition - ch.minv.blockPosition;
              byteMass += Math.max(1, delta / 1024);
            }
            for (let k = b; k < Math.min(numBins, b + step); k++) {
              bins[k] = byteMass;
            }
          }),
        );
      }
    } catch {
      // ignore
    }

    return bins;
  }
}

/**
 * 테스트 및 내장 S. cerevisiae 데모 데이터를 위한 메모리 기반 AlignmentSource
 */
export class MemoryAlignmentSource implements AlignmentSource {
  readonly id: string;
  readonly label: string;
  private contigs: ContigInfo[];
  private readsByChrom: Map<string, AlignedRead[]>;

  constructor(
    id: string,
    label: string,
    contigs: ContigInfo[],
    readsByChrom: Map<string, AlignedRead[]>,
  ) {
    this.id = id;
    this.label = label;
    this.contigs = contigs;
    this.readsByChrom = readsByChrom;
  }

  async getContigs(): Promise<ContigInfo[]> {
    return this.contigs;
  }

  async resolveContigName(
    query: string,
    manualOverrides?: Record<string, string>,
    targetLength?: number,
  ): Promise<string | null> {
    return resolveMatchingRefName(
      query,
      this.contigs.map((c) => c.name),
      manualOverrides,
      targetLength,
      this.contigs,
    );
  }

  async getRecords(chrom: string, start: number, end: number, _maxRecords?: number): Promise<AlignedRead[]> {
    const list = this.readsByChrom.get(chrom) ?? [];
    return list.filter((r) => r.end > start && r.start < end);
  }

  async getOverviewCoverage(
    chrom: string,
    start: number,
    end: number,
    numBins: number,
  ): Promise<number[]> {
    const bins = new Array<number>(numBins).fill(0);
    const list = this.readsByChrom.get(chrom) ?? [];
    const span = Math.max(1, end - start);
    const binWidth = span / numBins;

    for (const r of list) {
      if (r.end <= start || r.start >= end) continue;
      const b0 = Math.max(0, Math.floor((r.start - start) / binWidth));
      const b1 = Math.min(numBins - 1, Math.floor((r.end - 1 - start) / binWidth));
      for (let b = b0; b <= b1; b++) {
        const ovStart = Math.max(r.start, start + b * binWidth);
        const ovEnd = Math.min(r.end, start + (b + 1) * binWidth);
        if (ovEnd > ovStart) {
          bins[b] += (ovEnd - ovStart) / binWidth;
        }
      }
    }
    return bins;
  }
}

