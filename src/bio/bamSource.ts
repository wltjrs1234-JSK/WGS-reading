import { BamFile } from '@gmod/bam';
import {
  readBlobRangeInChunks,
  readBlobSliceSafe,
  SAFE_IO_SLICE_BYTES,
  type MainThreadBlobRelay,
} from '../utils/safeBlobReader';
import { resolveMatchingRefName } from './yeastAliases';
import type { AlignedRead, ContigInfo, SvRawReadSignal } from './types';

export interface GetRecordsOptions {
  includeQual?: boolean;
  onProgress?: (loadedBytes: number, totalBytes: number) => void;
  signal?: AbortSignal;
}

export interface GmodBamSourceOptions {
  onStatus?: (message: string) => void;
  mainThreadBamRelay?: MainThreadBlobRelay;
  mainThreadIndexRelay?: MainThreadBlobRelay;
}

export interface AlignmentSource {
  readonly id: string;
  readonly label: string;
  getContigs(): Promise<ContigInfo[]>;
  resolveContigName(
    query: string,
    manualOverrides?: Record<string, string>,
    targetLength?: number,
  ): Promise<string | null>;
  getRecords(
    chrom: string,
    start: number,
    end: number,
    maxRecords?: number,
    opts?: GetRecordsOptions,
  ): Promise<AlignedRead[]>;
  getOverviewCoverage(chrom: string, start: number, end: number, numBins: number): Promise<number[]>;
  getSvSignals?(chrom: string, start: number, end: number, maxRecords?: number): Promise<SvRawReadSignal[]>;
}

interface ReadOpts {
  signal?: AbortSignal;
}

interface SerializedFileOptions {
  mainThreadRelay?: MainThreadBlobRelay;
  onWholeLoadProgress?: (loaded: number, total: number) => void;
}

const BAM_PAGE_SIZE = SAFE_IO_SLICE_BYTES; // 512 KB read-ahead page
const MAX_CACHED_BAM_PAGES = 32; // ~16 MB RAM page cache per BAM file
const MAX_WHOLE_INDEX_BYTES = 96 * 1024 * 1024; // cache .bai/.csi up to 96 MB cleanly in memory

/**
 * 사내 DLP/DRM 환경 및 대용량 BAM 랜덤 액세스 병목을 막는 직렬화·분할 캐시 Filehandle:
 * 1) 모든 디스크 File.slice() I/O를 전역 동시성 1 큐(runGlobalSerializedIo)로 직렬화
 * 2) 512 KB 단위 분할 읽기 + 타임아웃(8초) + FileReader/메인스레드 폴백으로 중간 멈춤 차단
 * 3) 인덱스(.bai/.csi)는 1회 메모리 적재, BAM은 512 KB Read-Ahead LRU 페이지 캐시 사용
 */
export class SerializedCachedBlobFile {
  readonly source = undefined;
  private readonly file: File;
  private readonly cacheWholeFile: boolean;
  private readonly mainThreadRelay?: MainThreadBlobRelay;
  private onWholeLoadProgress?: (loaded: number, total: number) => void;
  private activeProgressHook: ((loaded: number, total: number) => void) | null = null;
  private wholeBytesPromise: Promise<Uint8Array> | null = null;
  private pageInFlight = new Map<number, Promise<Uint8Array>>();
  private pages = new Map<number, Uint8Array>();

  constructor(file: File, cacheWholeFile = false, opts: SerializedFileOptions = {}) {
    this.file = file;
    this.cacheWholeFile = cacheWholeFile && file.size <= MAX_WHOLE_INDEX_BYTES;
    this.mainThreadRelay = opts.mainThreadRelay;
    this.onWholeLoadProgress = opts.onWholeLoadProgress;
  }

  setWholeLoadProgress(cb?: (loaded: number, total: number) => void): void {
    this.onWholeLoadProgress = cb;
  }

  setActiveProgressHook(cb: ((loaded: number, total: number) => void) | null): void {
    this.activeProgressHook = cb;
  }

  private getOrLoadWholeBytes(signal?: AbortSignal): Promise<Uint8Array> {
    if (!this.wholeBytesPromise) {
      this.wholeBytesPromise = readBlobRangeInChunks(this.file, 0, this.file.size, {
        signal,
        mainThreadRelay: this.mainThreadRelay,
        onProgress: (loaded, total) => this.onWholeLoadProgress?.(loaded, total),
      }).catch((err) => {
        this.wholeBytesPromise = null;
        throw err;
      });
    }
    return this.wholeBytesPromise;
  }

  private async getPage(pageIdx: number, signal?: AbortSignal): Promise<Uint8Array> {
    const hit = this.pages.get(pageIdx);
    if (hit) {
      this.pages.delete(pageIdx);
      this.pages.set(pageIdx, hit);
      return hit;
    }

    let inFlight = this.pageInFlight.get(pageIdx);
    if (!inFlight) {
      const pStart = pageIdx * BAM_PAGE_SIZE;
      const pEnd = Math.min(this.file.size, pStart + BAM_PAGE_SIZE);
      inFlight = readBlobSliceSafe(this.file, pStart, pEnd, {
        signal,
        mainThreadRelay: this.mainThreadRelay,
      })
        .then((loaded) => {
          this.pages.set(pageIdx, loaded);
          if (this.pages.size > MAX_CACHED_BAM_PAGES) {
            const oldestKey = this.pages.keys().next().value;
            if (oldestKey !== undefined) this.pages.delete(oldestKey);
          }
          return loaded;
        })
        .finally(() => {
          this.pageInFlight.delete(pageIdx);
        });
      this.pageInFlight.set(pageIdx, inFlight);
    }
    return inFlight;
  }

  async read(length: number, position = 0, opts: ReadOpts = {}): Promise<Uint8Array> {
    if (length <= 0 || position >= this.file.size) return new Uint8Array(0);
    opts.signal?.throwIfAborted();

    const safePos = Math.max(0, position);
    const end = Math.min(this.file.size, safePos + length);
    const actualLen = end - safePos;
    if (actualLen <= 0) return new Uint8Array(0);

    if (this.cacheWholeFile) {
      const full = await this.getOrLoadWholeBytes(opts.signal);
      opts.signal?.throwIfAborted();
      return full.subarray(safePos, end);
    }

    // 큰 단일 청크(> 1 MB)도 512 KB씩 스트리밍 읽기해 진행률 업데이트 및 OS 응답성을 확보합니다.
    if (actualLen > BAM_PAGE_SIZE * 2) {
      return readBlobRangeInChunks(this.file, safePos, end, {
        signal: opts.signal,
        mainThreadRelay: this.mainThreadRelay,
        onProgress: (loaded, total) => this.activeProgressHook?.(loaded, total),
      });
    }

    const firstPage = Math.floor(safePos / BAM_PAGE_SIZE);
    const lastPage = Math.floor((end - 1) / BAM_PAGE_SIZE);

    if (firstPage === lastPage) {
      const page = await this.getPage(firstPage, opts.signal);
      const off = safePos - firstPage * BAM_PAGE_SIZE;
      return page.subarray(off, off + actualLen);
    }

    const out = new Uint8Array(actualLen);
    let writeOff = 0;
    for (let p = firstPage; p <= lastPage; p++) {
      opts.signal?.throwIfAborted();
      const page = await this.getPage(p, opts.signal);
      const pBase = p * BAM_PAGE_SIZE;
      const sOff = Math.max(0, safePos - pBase);
      const eOff = Math.min(page.length, end - pBase);
      if (eOff > sOff) {
        out.set(page.subarray(sOff, eOff), writeOff);
        writeOff += eOff - sOff;
        this.activeProgressHook?.(writeOff, actualLen);
      }
    }
    return writeOff === actualLen ? out : out.subarray(0, writeOff);
  }

  async readFile(options?: { encoding?: string; signal?: AbortSignal } | string): Promise<Uint8Array | string> {
    const opts = typeof options === 'string' ? { encoding: options } : options ?? {};
    opts.signal?.throwIfAborted();

    const bytes = this.cacheWholeFile
      ? await this.getOrLoadWholeBytes(opts.signal)
      : await readBlobRangeInChunks(this.file, 0, this.file.size, {
          signal: opts.signal,
          mainThreadRelay: this.mainThreadRelay,
          onProgress: (loaded, total) => this.activeProgressHook?.(loaded, total),
        });

    if (opts.encoding === 'utf8' || opts.encoding === 'utf-8') {
      return new TextDecoder('utf8').decode(bytes);
    }
    return bytes;
  }

  async stat(): Promise<{ size: number }> {
    return { size: this.file.size };
  }

  async close(): Promise<void> {
    this.pages.clear();
    this.pageInFlight.clear();
  }
}

interface RawBamChunk {
  minv: { blockPosition: number; dataPosition: number };
  maxv: { blockPosition: number; dataPosition: number };
  fetchedSize(): number;
}

interface RawBamRecord {
  ref_id: number;
  start: number;
  end: number;
  flags: number;
  mq?: number;
  name: string;
  seq: string;
  CIGAR: string;
  qual: Uint8Array | null;
  strand: number;
  next_refid: number;
  next_pos: number;
  template_length: number;
  getTag(name: string): unknown;
  isPaired(): boolean;
  isMateUnmapped(): boolean;
  isReverseComplemented(): boolean;
  isMateReverseComplemented(): boolean;
  isSupplementary(): boolean;
}

function materializeAlignedRead(r: RawBamRecord, idxTag: number, includeQual: boolean): AlignedRead {
  return {
    id: `${r.name}_${r.start}_${r.flags}_${idxTag}`,
    name: r.name,
    start: r.start,
    end: r.end,
    strand: r.strand === -1 ? '-' : '+',
    seq: r.seq,
    cigar: r.CIGAR,
    mapq: r.mq ?? 60,
    qual: includeQual ? r.qual ?? undefined : undefined,
    flags: r.flags,
  };
}

function yieldToWorkerLoop(): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, 0));
}

/**
 * 로컬 PC의 BAM + BAI(또는 CSI) 파일을 @gmod/bam + SerializedCachedBlobFile로 읽는 소스
 * 외부 서버 업로드 없이 브라우저 샌드박스 내에서 바이트 단위 랜덤 액세스로 동작합니다.
 */
const MAX_BAM_CACHE_BYTES = 48 * 1024 * 1024; // 48 MB/트랙 (다중 트랙 OOM 방지)
const BAM_CACHE_IDLE_MS = 45_000;
const MAX_RECORDS_PER_WINDOW = 25_000;

export class GmodBamSource implements AlignmentSource {
  readonly id: string;
  readonly label: string;
  private bam: BamFile;
  private headerLoaded = false;
  private bamHandle: SerializedCachedBlobFile;
  private indexHandle: SerializedCachedBlobFile;
  private onStatus?: (message: string) => void;
  private cachedContigs: ContigInfo[] | null = null;
  private contigResolveCache = new Map<string, string | null>();

  constructor(
    id: string,
    label: string,
    bamFile: File,
    indexFile: File,
    opts: GmodBamSourceOptions = {},
  ) {
    this.id = id;
    this.label = label;
    this.onStatus = opts.onStatus;
    const isCsi = indexFile.name.toLowerCase().endsWith('.csi');
    this.bamHandle = new SerializedCachedBlobFile(bamFile, false, {
      mainThreadRelay: opts.mainThreadBamRelay,
    });
    this.indexHandle = new SerializedCachedBlobFile(indexFile, true, {
      mainThreadRelay: opts.mainThreadIndexRelay,
    });
    this.bam = new BamFile({
      bamFilehandle: this.bamHandle as never,
      maxCacheBytes: MAX_BAM_CACHE_BYTES,
      cacheIdleTimeoutMs: BAM_CACHE_IDLE_MS,
      ...(isCsi
        ? { csiFilehandle: this.indexHandle as never }
        : { baiFilehandle: this.indexHandle as never }),
    });
  }

  private async ensureHeader(signal?: AbortSignal): Promise<void> {
    if (!this.headerLoaded) {
      let lastProg = 0;
      this.indexHandle.setWholeLoadProgress((loaded, total) => {
        const now = performance.now();
        if (now - lastProg >= 180 || loaded >= total) {
          lastProg = now;
          const mbL = (loaded / (1024 * 1024)).toFixed(1);
          const mbT = Math.max(0.1, total / (1024 * 1024)).toFixed(1);
          this.onStatus?.(`'${this.label}' 인덱스 읽는 중 (${mbL}/${mbT} MB)...`);
        }
      });
      try {
        await this.bam.getHeader({ signal });
        this.headerLoaded = true;
      } catch (err) {
        const msg = err instanceof Error ? err.message : String(err);
        if (msg.includes('Not a BAI file')) {
          throw new Error(`'${this.label}' 인덱스(.bai) 형식이 손상되었거나 올바른 BAI 파일이 아닙니다.`);
        }
        if (msg.includes('Not a BAM file')) {
          throw new Error(`'${this.label}' 파일이 올바른 바이너리 BAM(BGZF) 형식이 아닙니다.`);
        }
        throw err;
      } finally {
        this.indexHandle.setWholeLoadProgress(undefined);
      }
    }
  }

  async getContigs(): Promise<ContigInfo[]> {
    await this.ensureHeader();
    if (!this.cachedContigs) {
      this.cachedContigs = (this.bam.indexToChr ?? []).map((c) => ({
        name: c.refName,
        length: c.length,
      }));
    }
    return this.cachedContigs;
  }

  async resolveContigName(
    query: string,
    manualOverrides?: Record<string, string>,
    targetLength?: number,
  ): Promise<string | null> {
    if (manualOverrides && manualOverrides[query]) {
      const contigs = await this.getContigs();
      return resolveMatchingRefName(
        query,
        contigs.map((c) => c.name),
        manualOverrides,
        targetLength,
        contigs,
      );
    }
    const cacheKey = `${query}::${targetLength ?? 0}`;
    if (this.contigResolveCache.has(cacheKey)) {
      return this.contigResolveCache.get(cacheKey) ?? null;
    }
    const contigs = await this.getContigs();
    const resolved = resolveMatchingRefName(
      query,
      contigs.map((c) => c.name),
      undefined,
      targetLength,
      contigs,
    );
    this.contigResolveCache.set(cacheKey, resolved);
    return resolved;
  }

  private async iterateChunkRecordsInRange(
    chrom: string,
    start: number,
    end: number,
    opts: { signal?: AbortSignal; onProgress?: (loadedBytes: number, totalBytes: number) => void },
    onValidRecord: (r: RawBamRecord) => boolean | void,
  ): Promise<void> {
    const signal = opts.signal;
    signal?.throwIfAborted();
    await this.ensureHeader(signal);

    const safeStart = Math.max(0, start);
    const safeEnd = Math.max(safeStart + 1, end);
    const bamInternal = this.bam as unknown as {
      getSeqId(chr: string, o?: { signal?: AbortSignal }): Promise<number | undefined>;
      blocksForRange(chr: string, s: number, e: number, o?: { signal?: AbortSignal }): Promise<RawBamChunk[]>;
      _cachedChunkFeatures(chunk: RawBamChunk, o: { signal?: AbortSignal }): Promise<RawBamRecord[]>;
    };

    const chrId = await bamInternal.getSeqId(chrom, { signal });
    if (chrId === undefined) return;

    const chunks = await bamInternal.blocksForRange(chrom, safeStart, safeEnd, { signal });
    if (chunks.length === 0) {
      opts.onProgress?.(0, 0);
      return;
    }

    let totalBytes = 0;
    for (let i = 0; i < chunks.length; i++) {
      totalBytes += chunks[i].fetchedSize();
    }
    let completedChunkBytes = 0;
    const onProgress = opts.onProgress;
    onProgress?.(0, totalBytes);

    try {
      // 청크를 1개씩 순차 디코딩하여(기존 6개 동시 읽기 대신) 피크 메모리 폭발·OOM을 막습니다.
      let stoppedEarly = false;
      for (let ci = 0; ci < chunks.length && !stoppedEarly; ci++) {
        signal?.throwIfAborted();
        const chunk = chunks[ci];
        const chunkBytes = chunk.fetchedSize();
        const baseDone = completedChunkBytes;

        this.bamHandle.setActiveProgressHook(
          onProgress ? (subDone) => onProgress(baseDone + Math.min(chunkBytes, subDone), totalBytes) : null,
        );

        const chunkStartMs = performance.now();
        const features = await bamInternal._cachedChunkFeatures(chunk, { signal });
        this.bamHandle.setActiveProgressHook(null);

        if (features.length > 0) {
          const first = features[0];
          if (first.ref_id > chrId || (first.ref_id === chrId && first.start >= safeEnd)) {
            break;
          }
        }

        for (let fi = 0; fi < features.length; fi++) {
          const r = features[fi];
          if (r.ref_id > chrId) break;
          if (r.ref_id === chrId) {
            if (r.start >= safeEnd) break;
            const rEnd = r.end > r.start ? r.end : r.start + 1;
            if (rEnd > safeStart && (r.flags & 0x704) === 0) {
              if (onValidRecord(r) === false) {
                stoppedEarly = true;
                break;
              }
            }
          }
        }

        completedChunkBytes += chunkBytes;
        onProgress?.(completedChunkBytes, totalBytes);

        if (!stoppedEarly && ci + 1 < chunks.length && performance.now() - chunkStartMs >= 8) {
          await yieldToWorkerLoop();
        }
      }

      onProgress?.(totalBytes, totalBytes);
    } finally {
      this.bamHandle.setActiveProgressHook(null);
    }
  }

  async getRecords(
    chrom: string,
    start: number,
    end: number,
    maxRecords = MAX_RECORDS_PER_WINDOW,
    opts?: GetRecordsOptions,
  ): Promise<AlignedRead[]> {
    const safeStart = Math.max(0, start);
    const safeEnd = Math.max(safeStart + 1, end);
    const includeQual = opts?.includeQual === true;
    const limit = Number.isFinite(maxRecords) && maxRecords > 0 ? Math.floor(maxRecords) : 0;
    const out: AlignedRead[] = [];
    let totalValidSeen = 0;
    let rngState = ((safeStart * 2654435761) ^ (safeEnd * 40503) ^ 0x9e3779b9) >>> 0;

    // 고심도(URA3 플라스미드 등) 구간은 스트리밍 저수지 샘플링으로 객체 생성을 limit 이내로 제한합니다.
    await this.iterateChunkRecordsInRange(
      chrom,
      safeStart,
      safeEnd,
      { signal: opts?.signal, onProgress: opts?.onProgress },
      (r) => {
        const seenIdx = totalValidSeen++;
        if (limit === 0 || out.length < limit) {
          out.push(materializeAlignedRead(r, seenIdx, includeQual));
        } else {
          rngState = (Math.imul(rngState, 1664525) + 1013904223) >>> 0;
          const slot = Math.floor((rngState / 4294967296) * (seenIdx + 1));
          if (slot < limit) {
            out[slot] = materializeAlignedRead(r, seenIdx, includeQual);
          }
        }
      },
    );

    if (limit > 0 && totalValidSeen > limit) {
      out.sort((a, b) => a.start - b.start || a.end - b.end);
    }
    return out;
  }

  async getSvSignals(
    chrom: string,
    start: number,
    end: number,
    maxRecords = 50_000,
  ): Promise<SvRawReadSignal[]> {
    await this.ensureHeader();
    const indexToChr = this.bam.indexToChr ?? [];
    const signals: SvRawReadSignal[] = [];

    await this.iterateChunkRecordsInRange(chrom, start, end, {}, (r) => {
      const mapq = r.mq ?? 60;
      if (mapq < 15) return;

      const rawSa = r.getTag('SA');
      const saTag = typeof rawSa === 'string' && rawSa.length > 0 ? rawSa : undefined;
      const isPaired = r.isPaired();
      const mateMapped = isPaired && !r.isMateUnmapped();
      const nextRefId = r.next_refid;
      const tlen = r.template_length;
      const sameStrandOrientation =
        mateMapped &&
        nextRefId === r.ref_id &&
        r.isReverseComplemented() === r.isMateReverseComplemented();
      const isInterChrom = mateMapped && nextRefId >= 0 && nextRefId !== r.ref_id;
      const isLargeSpan = mateMapped && nextRefId === r.ref_id && Math.abs(tlen) >= 2_000;

      if (!saTag && !isInterChrom && !isLargeSpan && !sameStrandOrientation) return;

      const mateChrom =
        mateMapped && nextRefId >= 0 && nextRefId < indexToChr.length
          ? indexToChr[nextRefId].refName
          : null;

      signals.push({
        name: r.name,
        chrom,
        start: r.start,
        end: r.end,
        strand: r.strand === -1 ? '-' : '+',
        mapq,
        isPaired,
        mateChrom,
        matePos: Math.max(0, r.next_pos),
        mateStrand: r.isMateReverseComplemented() ? '-' : '+',
        tlen,
        saTag,
        isSupplementary: r.isSupplementary(),
      });
      if (signals.length >= maxRecords) return false;
    });

    return signals;
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
        for (let b = 0; b < numBins; b += step) {
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
        }
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
  private svSignalsByChrom?: Map<string, SvRawReadSignal[]>;

  constructor(
    id: string,
    label: string,
    contigs: ContigInfo[],
    readsByChrom: Map<string, AlignedRead[]>,
    svSignalsByChrom?: Map<string, SvRawReadSignal[]>,
  ) {
    this.id = id;
    this.label = label;
    this.contigs = contigs;
    this.readsByChrom = readsByChrom;
    this.svSignalsByChrom = svSignalsByChrom;
  }

  async getSvSignals(chrom: string, start: number, end: number): Promise<SvRawReadSignal[]> {
    const list = this.svSignalsByChrom?.get(chrom) ?? [];
    return list.filter((r) => r.end > start && r.start < end);
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

