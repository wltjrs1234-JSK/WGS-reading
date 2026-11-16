import { GmodBamSource, type AlignmentSource } from '../bio/bamSource';
import {
  annotateCrossStrainSpecificity,
  buildBarcodeCandidatesForRegion,
  collectBarcodeSignalsFromReads,
  HOTSPOT_LOCUS_NAMES,
  type BarcodeBuildTelemetry,
  type BarcodeCandidate,
  type BarcodeScanOptions,
} from '../bio/barcodeScanner';
import { forEachCigarOp } from '../bio/cigar';
import {
  buildAnalysisFromCalls,
  buildThreeWayAnalysis,
  computeConsensusForRegion,
  computeRegionMappingStats,
  serializeCallsForIntervals,
} from '../bio/consensus';
import { createYeastDemoDataset } from '../bio/demoData';
import {
  deriveAlleleTicksFromCalls,
  MAX_VISIBLE_PILEUP_ROWS,
  packReadsWithHeap,
} from '../bio/readPacking';
import {
  IndexedFastaSource,
  MemoryFastaSource,
  type ReferenceSequenceSource,
} from '../bio/fastaIndex';
import { reverseComplement } from '../bio/geneticCode';
import { AnnotationIndex, parseGff3File } from '../bio/gff3Parser';
import {
  buildBlastReferenceIndex,
  evaluateHitAgainstStrainCalls,
  findForeignMatchesInReads,
  parseBlastQueries,
  searchBlastReference,
  type BlastContigSeq,
  type BlastForeignReadHit,
  type BlastRawHit,
  type BlastReferenceIndex,
} from '../bio/localBlast';
import { createBuiltInS288CAnnotationIndex } from '../bio/s288cBuiltIn';
import { looksLikeYeastGenome } from '../bio/yeastAliases';
import type {
  AlignedRead,
  ConsensusOptions,
  ConsensusPositionCall,
  ContigInfo,
  VariantRecord,
} from '../bio/types';
import type {
  AllChromsOverviewResult,
  AnnotationOrigin,
  BarcodeScanPayload,
  BarcodeScanResult,
  BlastHitWithStrains,
  BlastOverlappingGene,
  BlastQueryResult,
  BlastRunPayload,
  BlastSearchResult,
  ExtractedExtraStrain,
  ExtractionResult,
  StrainWindowData,
  TrackSummaryMeta,
  WindowQueryResult,
  WorkspaceLoadedSummary,
  WorkerRequest,
  WorkerResponse,
} from './wgsProtocol';

const DETAIL_MAX_SPAN = 30_000;
const EXTRACT_CHUNK_BP = 20_000;
const EXTRACT_MAX_RECORDS_PER_CHUNK = 200_000;
const OVERVIEW_BINS = 480;
const DETAIL_COV_BINS = 2400;
const SCAN_CHUNK_BP = 25_000;
const BLAST_MAX_INDEX_BP = 100_000_000;
const BLAST_FOREIGN_CHUNK_BP = 50_000;
const BARCODE_SCAN_CHUNK_BP = 50_000;

let genomeId = 'sacCer3';
let refSource: ReferenceSequenceSource | null = null;
let cachedSacCer3Ref: MemoryFastaSource | null = null;
let cachedBuiltInAnnotation: AnnotationIndex | null = null;
let cachedBlastRefSource: ReferenceSequenceSource | null = null;
let cachedBlastIndex: BlastReferenceIndex | null = null;
let strainSources: AlignmentSource[] = [];
let annotationIndex: AnnotationIndex | null = null;
let annotationOrigin: AnnotationOrigin = 'none';
let isDemoMode = false;
let activeScanToken = 0;
let activeBlastToken = 0;
let activeBarcodeToken = 0;

function getBuiltInAnnotation(): AnnotationIndex {
  if (!cachedBuiltInAnnotation) {
    cachedBuiltInAnnotation = createBuiltInS288CAnnotationIndex();
  }
  return cachedBuiltInAnnotation;
}

function replyOk(reqId: number, data: unknown) {
  const msg: WorkerResponse = { reqId, ok: true, data };
  self.postMessage(msg);
}

function replyErr(reqId: number, err: unknown) {
  const message = err instanceof Error ? err.message : String(err);
  const msg: WorkerResponse = { reqId, ok: false, error: message };
  self.postMessage(msg);
}

function reportLoadProgress(message: string) {
  const ev: WorkerResponse = { type: 'load_progress', message };
  self.postMessage(ev);
}

async function resolveContigForStrain(
  strain: AlignmentSource,
  targetChrom: string,
  targetLength: number | undefined,
  overrides?: Record<string, Record<string, string>>,
): Promise<string | null> {
  return strain.resolveContigName(targetChrom, overrides?.[strain.id], targetLength);
}

async function buildTrackSummaries(
  refContigs: ContigInfo[],
  defaultChrom: string,
): Promise<TrackSummaryMeta[]> {
  // 사내 보안(DLP/DRM) 환경에서 다수 대용량 BAM+BAI 동시 오픈 경합을 피하기 위해 순차 초기화
  const summaries: TrackSummaryMeta[] = [];
  for (let i = 0; i < strainSources.length; i++) {
    const s = strainSources[i];
    reportLoadProgress(
      `'${s.label}' BAM 헤더·BAI 인덱스 읽는 중 (${i + 1}/${strainSources.length})...`,
    );
    const strainContigs = await s.getContigs();
    const mappings: Record<string, string | null> = {};
    for (const rc of refContigs) {
      mappings[rc.name] = await s.resolveContigName(rc.name, undefined, rc.length);
    }
    summaries.push({
      id: s.id,
      label: s.label,
      contigs: strainContigs,
      matchedChrom: mappings[defaultChrom] ?? null,
      mappings,
    });
  }
  return summaries;
}

async function makeWorkspaceSummary(
  defaultChrom: string,
  defaultStart: number,
  defaultEnd: number,
): Promise<WorkspaceLoadedSummary> {
  if (!refSource) throw new Error('레퍼런스 지놈이 없습니다.');
  const contigs = refSource.getContigs();
  const strains = await buildTrackSummaries(contigs, defaultChrom);
  return {
    genomeId,
    referenceLabel: refSource.label,
    contigs,
    gffCount: annotationIndex?.totalCount ?? 0,
    annotationOrigin,
    strains,
    defaultChrom,
    defaultStart,
    defaultEnd,
    isDemo: isDemoMode,
  };
}

function findTargetContig(contigs: ContigInfo[], chromQuery: string): ContigInfo | undefined {
  const exact = contigs.find((c) => c.name === chromQuery);
  if (exact) return exact;
  if (refSource) {
    const resolved = refSource.resolveContigName(chromQuery);
    if (resolved) {
      return contigs.find((c) => c.name === resolved);
    }
  }
  return contigs[0];
}

async function handleQueryAllChroms(
  contigOverrides?: Record<string, Record<string, string>>,
): Promise<AllChromsOverviewResult> {
  if (!refSource) throw new Error('레퍼런스가 로드되지 않았습니다.');
  const contigs = refSource.getContigs();

  const strains: AllChromsOverviewResult['strains'] = [];
  for (let i = 0; i < strainSources.length; i++) {
    const s = strainSources[i];
    reportLoadProgress(
      `'${s.label}' 전체 염색체 커버리지 계산 중 (${i + 1}/${strainSources.length})...`,
    );
    const byChrom: AllChromsOverviewResult['strains'][number]['byChrom'] = [];
    for (const c of contigs) {
      const matchedChrom = await resolveContigForStrain(s, c.name, c.length, contigOverrides);
      if (!matchedChrom) {
        byChrom.push({
          chrom: c.name,
          matchedChrom: null,
          bins: new Array<number>(24).fill(0),
          meanDepth: 0,
        });
        continue;
      }
      const numBins = Math.max(12, Math.min(80, Math.round(c.length / 20_000)));
      const bins = await s.getOverviewCoverage(matchedChrom, 0, c.length, numBins);
      const meanDepth = bins.reduce((a, b) => a + b, 0) / Math.max(1, bins.length);
      byChrom.push({
        chrom: c.name,
        matchedChrom,
        bins,
        meanDepth,
      });
    }
    strains.push({
      id: s.id,
      label: s.label,
      byChrom,
    });
  }

  return { contigs, strains };
}

async function handleQueryWindow(payload: {
  chrom: string;
  start: number;
  end: number;
  viewSpan?: number;
  primaryStrainIndex: number;
  secondaryStrainIndex: number | null;
  options: ConsensusOptions;
  contigOverrides?: Record<string, Record<string, string>>;
}): Promise<WindowQueryResult> {
  if (!refSource) {
    throw new Error('레퍼런스 FASTA가 아직 로드되지 않았습니다.');
  }

  const contigs = refSource.getContigs();
  const contig = findTargetContig(contigs, payload.chrom);
  if (!contig) {
    throw new Error('레퍼런스에 유효한 염색체/컨티그가 없습니다.');
  }

  const start = Math.max(0, Math.floor(payload.start));
  const end = Math.min(contig.length, Math.max(start + 1, Math.floor(payload.end)));
  const span = end - start;
  const effectiveSpan = Math.max(1, Math.floor(payload.viewSpan ?? span));
  const isOverviewOnly = effectiveSpan > DETAIL_MAX_SPAN;

  const features = annotationIndex
    ? annotationIndex.getFeaturesInRange(contig.name, start, end)
    : [];

  if (isOverviewOnly) {
    const strains: StrainWindowData[] = await Promise.all(
      strainSources.map(async (s) => {
        const matchedChrom = await resolveContigForStrain(
          s,
          contig.name,
          contig.length,
          payload.contigOverrides,
        );
        if (!matchedChrom) {
          return {
            id: s.id,
            label: s.label,
            matchedChrom: null,
            reads: [],
            overviewBins: new Array<number>(OVERVIEW_BINS).fill(0),
            meanDepth: 0,
          };
        }
        const overviewBins = await s.getOverviewCoverage(
          matchedChrom,
          start,
          end,
          OVERVIEW_BINS,
        );
        const meanDepth =
          overviewBins.reduce((a, b) => a + b, 0) / Math.max(1, overviewBins.length);
        return {
          id: s.id,
          label: s.label,
          matchedChrom,
          reads: [],
          overviewBins,
          meanDepth,
        };
      }),
    );

    return {
      chrom: contig.name,
      start,
      end,
      isOverviewOnly: true,
      refSeq: '',
      features,
      strains,
      analysis: null,
    };
  }

  const refSeq = await refSource.getSequence(contig.name, start, end);

  // 다수 BAM 파일 동시 디스크/DLP I/O 및 WASM 힙 경합을 막기 위해 균주별로 순차 매핑
  const rawStrainsData: Array<{
    id: string;
    label: string;
    matchedChrom: string | null;
    allReads: AlignedRead[];
    packedReads: ReturnType<typeof packReadsWithHeap>['visiblePacked'];
    totalRows: number;
    overviewBins: number[];
    meanDepth: number;
  }> = [];

  for (let sIdx = 0; sIdx < strainSources.length; sIdx++) {
    const s = strainSources[sIdx];
    const matchedChrom = await resolveContigForStrain(
      s,
      contig.name,
      contig.length,
      payload.contigOverrides,
    );
    if (!matchedChrom) {
      rawStrainsData.push({
        id: s.id,
        label: s.label,
        matchedChrom: null,
        allReads: [],
        packedReads: [],
        totalRows: 0,
        overviewBins: new Array<number>(Math.min(OVERVIEW_BINS, span)).fill(0),
        meanDepth: 0,
      });
      continue;
    }

    reportLoadProgress(
      `'${s.label}' ${contig.name} 구간 리드 읽는 중 (${sIdx + 1}/${strainSources.length})...`,
    );
    const rawReads = await s.getRecords(matchedChrom, start, end);
    const reads = rawReads.filter((r) => (r.flags & 0x704) === 0);
    const numBins = Math.min(DETAIL_COV_BINS, span);
    const overviewBins = new Array<number>(numBins).fill(0);
    const binWidth = span / numBins;
    let totalBases = 0;

    const addCovSpan = (covStart: number, covEnd: number) => {
      const ovStart = covStart > start ? covStart : start;
      const ovEnd = covEnd < end ? covEnd : end;
      if (ovEnd <= ovStart) return;
      totalBases += ovEnd - ovStart;
      const b0 = Math.max(0, Math.floor((ovStart - start) / binWidth));
      const b1 = Math.min(numBins - 1, Math.floor((ovEnd - 1 - start) / binWidth));
      for (let b = b0; b <= b1; b++) {
        const bs = start + b * binWidth;
        const be = start + (b + 1) * binWidth;
        const segOv = (ovEnd < be ? ovEnd : be) - (ovStart > bs ? ovStart : bs);
        if (segOv > 0) overviewBins[b] += segOv / binWidth;
      }
    };

    for (const r of reads) {
      if (r.mapq < payload.options.minMapQ) continue;
      let refPos = r.start;
      let sawOp = false;
      forEachCigarOp(r.cigar, (op, len) => {
        sawOp = true;
        if (op === 'M' || op === '=' || op === 'X') {
          addCovSpan(refPos, refPos + len);
          refPos += len;
        } else if (op === 'D' || op === 'N') {
          refPos += len;
        }
      });
      if (!sawOp) {
        addCovSpan(r.start, r.end);
      }
    }

    const packed = packReadsWithHeap(reads, MAX_VISIBLE_PILEUP_ROWS);

    rawStrainsData.push({
      id: s.id,
      label: s.label,
      matchedChrom,
      allReads: reads,
      packedReads: packed.visiblePacked,
      totalRows: packed.totalRows,
      overviewBins,
      meanDepth: totalBases / Math.max(1, span),
    });
  }

  const hasPrimary =
    strainSources.length > 0 &&
    payload.primaryStrainIndex >= 0 &&
    payload.primaryStrainIndex < rawStrainsData.length;
  const s1Reads: AlignedRead[] = hasPrimary
    ? rawStrainsData[payload.primaryStrainIndex]?.allReads ?? []
    : [];
  const hasSecondary =
    hasPrimary &&
    payload.secondaryStrainIndex !== null &&
    payload.secondaryStrainIndex >= 0 &&
    payload.secondaryStrainIndex < rawStrainsData.length;
  const s2Reads: AlignedRead[] | null = hasSecondary
    ? rawStrainsData[payload.secondaryStrainIndex!]?.allReads ?? []
    : null;

  const analysis = hasPrimary
    ? buildThreeWayAnalysis(
        contig.name,
        start,
        end,
        refSeq,
        s1Reads,
        s2Reads,
        payload.options,
        true,
      )
    : null;

  const s1Ticks = analysis ? deriveAlleleTicksFromCalls(analysis.s1Calls) : undefined;
  const s2Ticks =
    analysis && analysis.hasStrain2 ? deriveAlleleTicksFromCalls(analysis.s2Calls) : undefined;

  const strains: StrainWindowData[] = rawStrainsData.map((sd, idx) => {
    let alleleTicks =
      idx === payload.primaryStrainIndex
        ? s1Ticks
        : hasSecondary && idx === payload.secondaryStrainIndex
          ? s2Ticks
          : undefined;

    if (!alleleTicks && sd.matchedChrom && sd.allReads.length > 0 && refSeq) {
      const extraCalls = computeConsensusForRegion(
        refSeq,
        sd.allReads,
        start,
        end,
        payload.options,
      );
      alleleTicks = deriveAlleleTicksFromCalls(extraCalls);
    }

    return {
      id: sd.id,
      label: sd.label,
      matchedChrom: sd.matchedChrom,
      reads: [], // 가시 행(28줄)만 packedReads로 전송해 워커-메인 스레드 직렬화 부하를 95%+ 절감
      packedReads: sd.packedReads,
      totalRows: sd.totalRows,
      alleleTicks: alleleTicks ?? [],
      overviewBins: sd.overviewBins,
      meanDepth: sd.meanDepth,
    };
  });

  return {
    chrom: contig.name,
    start,
    end,
    isOverviewOnly: false,
    refSeq,
    features,
    strains,
    analysis,
  };
}

async function handleExtractRegion(payload: {
  chrom: string;
  start: number;
  end: number;
  primaryStrainIndex: number;
  secondaryStrainIndex: number | null;
  options: ConsensusOptions;
  exons?: Array<{ start: number; end: number }>;
  contigOverrides?: Record<string, Record<string, string>>;
}): Promise<ExtractionResult> {
  if (!refSource) throw new Error('레퍼런스가 로드되지 않았습니다.');

  const contigs = refSource.getContigs();
  const contig = findTargetContig(contigs, payload.chrom);
  const targetChrom = contig?.name ?? payload.chrom;

  const span = payload.end - payload.start;
  if (span <= 0) throw new Error('유효하지 않은 선택 구간입니다.');
  if (span > 250_000) {
    throw new Error('브라우저 메모리 보호를 위해 1회 서열 추출은 250,000 bp 이하로 제한됩니다.');
  }

  const refSeq = await refSource.getSequence(targetChrom, payload.start, payload.end);

  const matchedChroms = await Promise.all(
    strainSources.map((s) =>
      resolveContigForStrain(s, targetChrom, contig?.length, payload.contigOverrides),
    ),
  );

  // 청크 단위로 순차 읽기 & 200,000 리드 상한을 적용해 대용량·고심도 유전자에서도 리드 누락/다운샘플링을 방지합니다.
  const callsPerStrain: ConsensusPositionCall[][] = strainSources.map(() => []);
  for (let chunkStart = payload.start; chunkStart < payload.end; chunkStart += EXTRACT_CHUNK_BP) {
    const chunkEnd = Math.min(payload.end, chunkStart + EXTRACT_CHUNK_BP);
    const chunkRefSeq = refSeq.slice(chunkStart - payload.start, chunkEnd - payload.start);

    for (let i = 0; i < strainSources.length; i++) {
      const matched = matchedChroms[i];
      const reads = matched
        ? await strainSources[i].getRecords(
            matched,
            chunkStart,
            chunkEnd,
            EXTRACT_MAX_RECORDS_PER_CHUNK,
          )
        : [];
      const chunkCalls = computeConsensusForRegion(
        chunkRefSeq,
        reads,
        chunkStart,
        chunkEnd,
        payload.options,
      );
      for (const c of chunkCalls) {
        callsPerStrain[i].push(c);
      }
    }
  }

  const s1Idx = payload.primaryStrainIndex;
  const s2Idx = payload.secondaryStrainIndex;
  const hasStrain1 = strainSources.length > 0 && s1Idx >= 0 && s1Idx < strainSources.length;
  const hasSecondary = s2Idx !== null && s2Idx >= 0 && s2Idx < strainSources.length;

  const emptyZeroCalls = () =>
    computeConsensusForRegion(refSeq, [], payload.start, payload.end, payload.options);

  const s1Calls = hasStrain1 ? callsPerStrain[s1Idx] : emptyZeroCalls();
  const s2Calls = hasSecondary ? callsPerStrain[s2Idx!] : null;

  const analysis = buildAnalysisFromCalls(
    targetChrom,
    payload.start,
    payload.end,
    refSeq,
    s1Calls,
    s2Calls,
    false,
  );

  const validExons = (payload.exons ?? [])
    .map((ex) => ({
      start: Math.max(payload.start, ex.start),
      end: Math.min(payload.end, ex.end),
    }))
    .filter((ex) => ex.end > ex.start)
    .sort((a, b) => a.start - b.start);
  const splicedRefSlice =
    validExons.length > 0
      ? validExons
          .map((ex) => refSeq.slice(ex.start - payload.start, ex.end - payload.start))
          .join('')
      : '';

  const extraStrains: ExtractedExtraStrain[] = [];
  for (let i = 0; i < strainSources.length; i++) {
    if (i === s1Idx || (hasSecondary && i === s2Idx)) continue;
    const calls = callsPerStrain[i];
    const serialized = serializeCallsForIntervals(calls, payload.start, [
      { start: payload.start, end: payload.end },
    ]);
    const mappingStats = computeRegionMappingStats(calls, refSeq);

    let splicedUngappedSeq: string | undefined;
    let splicedAlignedSeq: string | undefined;
    let splicedMappingStats: ExtractedExtraStrain['splicedMappingStats'];
    if (validExons.length > 0) {
      const sp = serializeCallsForIntervals(calls, payload.start, validExons);
      splicedUngappedSeq = sp.ungappedSeq;
      splicedAlignedSeq = sp.alignedSeq;
      splicedMappingStats = computeRegionMappingStats(sp.subsetCalls, splicedRefSlice);
    }

    extraStrains.push({
      label: strainSources[i].label,
      ungappedSeq: serialized.ungappedSeq,
      alignedSeq: serialized.alignedSeq,
      splicedUngappedSeq,
      splicedAlignedSeq,
      mappingStats,
      splicedMappingStats,
    });
  }

  return {
    analysis,
    hasStrain1,
    strain1Label: hasStrain1 ? strainSources[s1Idx].label : 'Reference Only',
    strain2Label: hasSecondary ? strainSources[s2Idx!]?.label ?? 'Strain 2' : null,
    extraStrains,
  };
}

async function handleScanChromosome(payload: {
  chrom: string;
  primaryStrainIndex: number;
  secondaryStrainIndex: number | null;
  options: ConsensusOptions;
  contigOverrides?: Record<string, Record<string, string>>;
}): Promise<VariantRecord[]> {
  if (!refSource) throw new Error('레퍼런스가 로드되지 않았습니다.');
  if (strainSources.length === 0) return [];
  const contig = findTargetContig(refSource.getContigs(), payload.chrom);
  if (!contig) throw new Error(`염색체 ${payload.chrom}를 찾을 수 없습니다.`);

  const token = ++activeScanToken;
  const s1 = strainSources[payload.primaryStrainIndex];
  const hasSecondary =
    payload.secondaryStrainIndex !== null &&
    payload.secondaryStrainIndex >= 0 &&
    payload.secondaryStrainIndex < strainSources.length;
  const s2 = hasSecondary ? strainSources[payload.secondaryStrainIndex!] : undefined;

  const s1Chrom = s1
    ? await resolveContigForStrain(s1, contig.name, contig.length, payload.contigOverrides)
    : null;
  const s2Chrom = s2
    ? await resolveContigForStrain(s2, contig.name, contig.length, payload.contigOverrides)
    : null;

  const allVariants: VariantRecord[] = [];

  for (let start = 0; start < contig.length; start += SCAN_CHUNK_BP) {
    if (token !== activeScanToken) {
      throw new Error('염색체 스캔이 취소되었습니다.');
    }

    const end = Math.min(contig.length, start + SCAN_CHUNK_BP);
    const [refSeq, s1Reads, s2Reads] = await Promise.all([
      refSource.getSequence(contig.name, start, end),
      s1 && s1Chrom ? s1.getRecords(s1Chrom, start, end) : Promise.resolve([]),
      hasSecondary && s2 && s2Chrom ? s2.getRecords(s2Chrom, start, end) : Promise.resolve(null),
    ]);

    const analysis = buildThreeWayAnalysis(
      contig.name,
      start,
      end,
      refSeq,
      s1Reads,
      hasSecondary ? s2Reads ?? [] : null,
      payload.options,
    );

    allVariants.push(...analysis.variants);

    self.postMessage({
      type: 'scan_progress',
      chrom: contig.name,
      scannedBp: end,
      totalBp: contig.length,
      variantsSoFar: allVariants.length,
    } satisfies WorkerResponse);
  }

  return allVariants;
}

function reportBlastProgress(message: string, percent: number) {
  const ev: WorkerResponse = { type: 'blast_progress', message, percent };
  self.postMessage(ev);
}

async function getOrBuildBlastIndex(
  token: number,
  caller: 'blast' | 'barcode' = 'blast',
): Promise<BlastReferenceIndex> {
  if (!refSource) throw new Error('레퍼런스 지놈이 아직 로드되지 않았습니다.');
  if (cachedBlastRefSource === refSource && cachedBlastIndex) {
    return cachedBlastIndex;
  }

  const contigs = refSource.getContigs();
  const totalLen = contigs.reduce((acc, c) => acc + c.length, 0);
  if (totalLen > BLAST_MAX_INDEX_BP) {
    throw new Error(
      `현재 레퍼런스 크기(${totalLen.toLocaleString()} bp)가 브라우저 로컬 BLAST 인덱스 허용치(100 Mb)를 초과합니다.`,
    );
  }

  if (caller === 'blast') {
    reportBlastProgress('Reference 지놈 k-mer 인덱스 구축 중...', 8);
  } else {
    reportBarcodeProgress('sacCer3 내인성 서열 비교를 위한 지놈 인덱스 구축 중...', 96, 0);
  }
  const contigSeqs: BlastContigSeq[] = [];
  for (let i = 0; i < contigs.length; i++) {
    if (caller === 'blast' && token !== activeBlastToken) {
      throw new Error('BLAST 검색이 취소되었습니다.');
    }
    if (caller === 'barcode' && token !== activeBarcodeToken) {
      throw new Error('바코드 스캔이 취소되었습니다.');
    }
    const c = contigs[i];
    const seq = (await refSource.getSequence(c.name, 0, c.length)).toUpperCase();
    contigSeqs.push({ name: c.name, seq });
  }

  const built = buildBlastReferenceIndex(contigSeqs);
  cachedBlastRefSource = refSource;
  cachedBlastIndex = built;
  return built;
}

async function handleRunBlast(payload: BlastRunPayload): Promise<BlastSearchResult> {
  if (!refSource) throw new Error('레퍼런스 지놈이 아직 로드되지 않았습니다.');
  const t0 = performance.now();
  const token = ++activeBlastToken;
  const queries = parseBlastQueries(payload.rawInput);
  const index = await getOrBuildBlastIndex(token);

  const contigs = refSource.getContigs();
  const contigLenByName = new Map(contigs.map((c) => [c.name, c.length]));

  const matchedChromsPerStrain = await Promise.all(
    strainSources.map(async (s) => {
      const m = new Map<string, string | null>();
      for (const c of contigs) {
        m.set(c.name, await resolveContigForStrain(s, c.name, c.length, payload.contigOverrides));
      }
      return m;
    }),
  );

  const results: BlastQueryResult[] = [];

  for (let qi = 0; qi < queries.length; qi++) {
    if (token !== activeBlastToken) throw new Error('BLAST 검색이 취소되었습니다.');
    const q = queries[qi];
    const basePct = 15 + Math.round((qi / queries.length) * 55);
    reportBlastProgress(
      `Reference 검색 및 균주 BAM 대조 중 (${qi + 1}/${queries.length}: ${q.id})...`,
      basePct,
    );

    const rawHits = searchBlastReference(index, q, {
      maxHitsPerQuery: payload.maxHitsPerQuery,
      minIdentityPct: payload.minIdentityPct,
      shortPrimerMaxMismatches: payload.shortPrimerMaxMismatches,
    });

    const enrichedHits: BlastHitWithStrains[] = [];
    for (const hit of rawHits) {
      if (token !== activeBlastToken) throw new Error('BLAST 검색이 취소되었습니다.');

      const overlappingFeatures = annotationIndex
        ? annotationIndex.getFeaturesInRange(hit.chrom, hit.refStart, hit.refEnd)
        : [];
      const overlappingGenes: BlastOverlappingGene[] = overlappingFeatures.map((f) => {
        let relation: BlastOverlappingGene['relation'] = 'partial_overlap';
        if (hit.refStart >= f.start && hit.refEnd <= f.end) relation = 'inside_gene';
        else if (hit.refStart <= f.start && hit.refEnd >= f.end) relation = 'exact_or_contains';
        return {
          id: f.id,
          name: f.name,
          systematicName: f.systematicName,
          strand: f.strand,
          start: f.start,
          end: f.end,
          relation,
        };
      });

      const refSubSeq = await refSource.getSequence(hit.chrom, hit.refStart, hit.refEnd);
      const strainComparisons = [];

      for (let si = 0; si < strainSources.length; si++) {
        const s = strainSources[si];
        const matchedChrom = matchedChromsPerStrain[si].get(hit.chrom) ?? null;
        const reads = matchedChrom
          ? await s.getRecords(matchedChrom, hit.refStart, hit.refEnd, 100_000)
          : [];
        const calls = matchedChrom
          ? computeConsensusForRegion(refSubSeq, reads, hit.refStart, hit.refEnd, payload.options)
          : [];
        strainComparisons.push(
          evaluateHitAgainstStrainCalls(hit, s.id, s.label, matchedChrom, calls, reads),
        );
      }

      enrichedHits.push({
        ...hit,
        alignedQueryDisplay:
          hit.strand === '-' ? reverseComplement(hit.alignedQueryPlus) : hit.alignedQueryPlus,
        alignedRefDisplay:
          hit.strand === '-' ? reverseComplement(hit.alignedRefPlus) : hit.alignedRefPlus,
        overlappingGenes,
        strains: strainComparisons,
      });
    }

    results.push({
      queryId: q.id,
      queryLength: q.seq.length,
      hits: enrichedHits,
      foreignHits: [],
    });
  }

  // Optional Soft-clip / Insertion Foreign Read Scanning in BAM tracks
  const foreignScope = payload.foreignBamScope ?? 'off';
  if (foreignScope !== 'off' && strainSources.length > 0) {
    const targetContigs =
      foreignScope === 'current_chrom' &&
      payload.currentChrom &&
      payload.currentChrom !== 'all'
        ? contigs.filter((c) => c.name === payload.currentChrom)
        : contigs;

    const totalTasks = Math.max(1, targetContigs.length * strainSources.length);
    let doneTasks = 0;

    for (let si = 0; si < strainSources.length; si++) {
      const s = strainSources[si];
      const byQuery = queries.map(() => [] as BlastForeignReadHit[]);

      for (const contig of targetContigs) {
        if (token !== activeBlastToken) throw new Error('BLAST 검색이 취소되었습니다.');
        const matchedChrom = matchedChromsPerStrain[si].get(contig.name) ?? null;
        doneTasks++;
        reportBlastProgress(
          `BAM Soft-clip 외래 서열 스캔 중: '${s.label}' · ${contig.name} (${doneTasks}/${totalTasks})`,
          70 + Math.round((doneTasks / totalTasks) * 28),
        );
        if (!matchedChrom) continue;

        const chromLen = contigLenByName.get(contig.name) ?? contig.length;
        for (let cs = 0; cs < chromLen; cs += BLAST_FOREIGN_CHUNK_BP) {
          if (token !== activeBlastToken) throw new Error('BLAST 검색이 취소되었습니다.');
          const ce = Math.min(chromLen, cs + BLAST_FOREIGN_CHUNK_BP);
          const reads = await s.getRecords(matchedChrom, cs, ce, 30_000);
          if (reads.length === 0) continue;

          for (let qi = 0; qi < queries.length; qi++) {
            const fHits = findForeignMatchesInReads(
              queries[qi],
              s.id,
              s.label,
              contig.name,
              reads,
            );
            for (const fh of fHits) {
              const dup = byQuery[qi].find(
                (ex) =>
                  ex.strainId === fh.strainId &&
                  ex.chrom === fh.chrom &&
                  Math.abs(ex.anchorPos - fh.anchorPos) <= 50 &&
                  ex.clipSide === fh.clipSide,
              );
              if (dup) {
                dup.supportingReads += fh.supportingReads;
              } else {
                byQuery[qi].push(fh);
              }
            }
          }
        }
      }

      for (let qi = 0; qi < queries.length; qi++) {
        results[qi].foreignHits.push(...byQuery[qi]);
      }
    }

    for (const r of results) {
      r.foreignHits.sort((a, b) => b.supportingReads - a.supportingReads || b.matchLen - a.matchLen);
      r.foreignHits = r.foreignHits.slice(0, 30);
    }
  }

  reportBlastProgress('완료', 100);
  return {
    referenceLabel: refSource.label,
    totalRefLength: index.totalLength,
    elapsedMs: Math.round(performance.now() - t0),
    queries: results,
    strainLabels: strainSources.map((s) => s.label),
  };
}

function reportBarcodeProgress(message: string, percent: number, foundSoFar: number) {
  const ev: WorkerResponse = { type: 'barcode_progress', message, percent, foundSoFar };
  self.postMessage(ev);
}

interface ScanInterval {
  chrom: string;
  start: number;
  end: number;
  label: string;
}

function buildBarcodeIntervals(
  scope: BarcodeScanPayload['scope'],
  contigs: ContigInfo[],
  currentChrom?: string,
): ScanInterval[] {
  const rawIntervals: ScanInterval[] = [];

  if (scope === 'hotspots') {
    const seenFeatureIds = new Set<string>();
    if (annotationIndex) {
      for (const locus of HOTSPOT_LOCUS_NAMES) {
        const hits = annotationIndex.searchGenes(locus, 3);
        for (const f of hits) {
          const matchExact =
            f.name.toUpperCase() === locus.toUpperCase() ||
            f.systematicName?.toUpperCase() === locus.toUpperCase();
          if (!matchExact || seenFeatureIds.has(f.id)) continue;
          seenFeatureIds.add(f.id);
          const resolvedChr = refSource?.resolveContigName(f.chrom) ?? f.chrom;
          const contig = contigs.find((c) => c.name === resolvedChr);
          if (!contig) continue;
          rawIntervals.push({
            chrom: contig.name,
            start: Math.max(0, f.start - 2_500),
            end: Math.min(contig.length, f.end + 2_500),
            label: f.name,
          });
        }
      }
    }
    if (rawIntervals.length === 0) {
      const targetContig = contigs.find((c) => c.name === currentChrom) ?? contigs[0];
      if (targetContig) {
        rawIntervals.push({
          chrom: targetContig.name,
          start: 0,
          end: Math.min(targetContig.length, 60_000),
          label: targetContig.name,
        });
      }
    }
  } else {
    const targetContigs =
      scope === 'current_chrom' && currentChrom && currentChrom !== 'all'
        ? contigs.filter((c) => c.name === currentChrom)
        : contigs;
    for (const c of targetContigs) {
      for (let s = 0; s < c.length; s += BARCODE_SCAN_CHUNK_BP) {
        rawIntervals.push({
          chrom: c.name,
          start: s,
          end: Math.min(c.length, s + BARCODE_SCAN_CHUNK_BP),
          label: `${c.name}:${(s + 1).toLocaleString()}`,
        });
      }
    }
  }

  // 동일 염색체 겹치는 핫스팟 구간 병합해 중복 읽기 제거
  rawIntervals.sort((a, b) => a.chrom.localeCompare(b.chrom) || a.start - b.start);
  const merged: ScanInterval[] = [];
  for (const iv of rawIntervals) {
    const prev = merged[merged.length - 1];
    if (prev && prev.chrom === iv.chrom && iv.start <= prev.end) {
      prev.end = Math.max(prev.end, iv.end);
      prev.label = `${prev.label}/${iv.label}`;
    } else {
      merged.push({ ...iv });
    }
  }
  return merged;
}

function isEndogenousSelfMatch(queryLen: number, hit: BlastRawHit): boolean {
  const globalIdentityRatio = hit.identities / Math.max(1, queryLen);
  if (queryLen < 25) {
    // 12~24 bp 짧은 바코드는 효모 지놈(12.16 Mb) 내 우연한 16 bp 부분 일치(2~3 불일치)로
    // 오인 제외되지 않도록 전체 길이 대비 무갭 최대 1 mismatch 수준만 내인성 서열로 판정
    return (
      hit.gaps === 0 &&
      hit.queryCoveragePct >= 90 &&
      hit.identities >= queryLen - 1 &&
      globalIdentityRatio >= 0.9
    );
  }
  return (
    hit.queryCoveragePct >= 80 &&
    hit.identityPct >= 88 &&
    globalIdentityRatio >= 0.78 &&
    hit.eValue <= 0.05
  );
}

async function handleScanBarcodes(payload: BarcodeScanPayload): Promise<BarcodeScanResult> {
  if (!refSource) throw new Error('레퍼런스 지놈이 로드되지 않았습니다.');
  if (strainSources.length === 0) {
    throw new Error('먼저 BAM + BAI 파일을 1개 이상 불러와주세요.');
  }

  const t0 = performance.now();
  const token = ++activeBarcodeToken;
  const contigs = refSource.getContigs();
  const contigLenByName = new Map(contigs.map((c) => [c.name, c.length]));
  const intervals = buildBarcodeIntervals(payload.scope, contigs, payload.currentChrom);

  const rawOpts = payload.options ?? {};
  const scanOpts: BarcodeScanOptions = {
    ...rawOpts,
    minMapQ: rawOpts.minMapQ ?? 10,
    excludeRepeatLoci: rawOpts.excludeRepeatLoci ?? true,
    excludeEndogenousRefMatches: rawOpts.excludeEndogenousRefMatches ?? true,
  };

  const targetIndices =
    payload.targetStrainIndex !== undefined &&
    payload.targetStrainIndex !== null &&
    payload.targetStrainIndex >= 0 &&
    payload.targetStrainIndex < strainSources.length
      ? [payload.targetStrainIndex]
      : strainSources.map((_, i) => i);

  const totalSteps = Math.max(1, intervals.length * targetIndices.length);
  let doneSteps = 0;
  let scannedBpTotal = 0;
  const rawCandidates: BarcodeCandidate[] = [];
  const telemetry: BarcodeBuildTelemetry = {
    repeatLociFiltered: 0,
    tandemOrArtifactFiltered: 0,
  };

  for (const sIdx of targetIndices) {
    const strain = strainSources[sIdx];
    const matchedChromMap = new Map<string, string | null>();
    for (const iv of intervals) {
      if (!matchedChromMap.has(iv.chrom)) {
        const len = contigLenByName.get(iv.chrom);
        matchedChromMap.set(
          iv.chrom,
          await resolveContigForStrain(strain, iv.chrom, len, payload.contigOverrides),
        );
      }
    }

    const strainFound: BarcodeCandidate[] = [];

    for (const iv of intervals) {
      if (token !== activeBarcodeToken) throw new Error('바코드 스캔이 취소되었습니다.');
      doneSteps++;
      const matchedChrom = matchedChromMap.get(iv.chrom) ?? null;
      reportBarcodeProgress(
        `'${strain.label}' ${iv.label} 구간 분석 중 (${doneSteps}/${totalSteps})`,
        Math.round((doneSteps / totalSteps) * 92),
        rawCandidates.length + strainFound.length,
      );
      if (!matchedChrom) continue;

      if (sIdx === targetIndices[0]) {
        scannedBpTotal += iv.end - iv.start;
      }

      const contigLen = contigLenByName.get(iv.chrom) ?? iv.end;
      const reads = await strain.getRecords(
        matchedChrom,
        Math.max(0, iv.start - 5),
        Math.min(contigLen, iv.end + 5),
        80_000,
      );
      if (reads.length === 0) continue;

      const signals = collectBarcodeSignalsFromReads(iv.chrom, reads, iv.start, iv.end, scanOpts);
      if (signals.length === 0) continue;

      const features = annotationIndex
        ? annotationIndex.getFeaturesInRange(iv.chrom, iv.start - 50, iv.end + 50)
        : [];
      const refCtxStart = Math.max(0, iv.start - 200);
      const refCtxEnd = Math.min(contigLen, iv.end + 200);
      const refCtxSeq = (await refSource.getSequence(iv.chrom, refCtxStart, refCtxEnd)).toUpperCase();

      const found = buildBarcodeCandidatesForRegion(
        strain.id,
        strain.label,
        signals,
        reads,
        features,
        scanOpts,
        { start: refCtxStart, seq: refCtxSeq },
        telemetry,
      );

      for (const cand of found) {
        const dup = strainFound.find(
          (ex) =>
            ex.chrom === cand.chrom &&
            ex.signalType === cand.signalType &&
            Math.abs(ex.anchorPos - cand.anchorPos) <= 12,
        );
        if (!dup) {
          strainFound.push(cand);
        } else if (cand.supportReads > dup.supportReads) {
          Object.assign(dup, cand);
        }
      }
    }

    rawCandidates.push(...strainFound);
  }

  let endogenousRefFiltered = 0;
  let keptCandidates = rawCandidates;
  if (scanOpts.excludeEndogenousRefMatches !== false && rawCandidates.length > 0) {
    const needSelfCheck = rawCandidates.filter(
      (c) =>
        c.motifHits.length === 0 &&
        !c.extractedBarseq20mer &&
        c.assembledSeq.length >= 12,
    );
    if (needSelfCheck.length > 0) {
      reportBarcodeProgress(
        `sacCer3 내인성 서열 중복 검증 중 (후보 ${needSelfCheck.length}건 대조)...`,
        95,
        rawCandidates.length,
      );
      const blastIdx = await getOrBuildBlastIndex(token, 'barcode');
      const dropSet = new Set<BarcodeCandidate>();
      for (const cand of needSelfCheck) {
        if (token !== activeBarcodeToken) throw new Error('바코드 스캔이 취소되었습니다.');
        const selfHits = searchBlastReference(
          blastIdx,
          { id: cand.id, seq: cand.assembledSeq },
          {
            maxHitsPerQuery: 3,
            minIdentityPct: 86,
            maxEValue: 5,
            shortPrimerMaxMismatches: cand.assembledSeq.length < 25 ? 1 : 2,
          },
        );
        if (selfHits.some((h) => isEndogenousSelfMatch(cand.assembledSeq.length, h))) {
          dropSet.add(cand);
          endogenousRefFiltered++;
        }
      }
      if (dropSet.size > 0) {
        keptCandidates = rawCandidates.filter((c) => !dropSet.has(c));
      }
    }
  }

  const annotated = annotateCrossStrainSpecificity(keptCandidates, targetIndices.length);
  reportBarcodeProgress('바코드 분석 완료', 100, annotated.length);

  return {
    scope: payload.scope,
    scannedLociCount: intervals.length,
    scannedBpTotal,
    elapsedMs: Math.round(performance.now() - t0),
    strainLabels: targetIndices.map((i) => strainSources[i].label),
    candidates: annotated,
    filteredStats: {
      repeatLociFiltered: telemetry.repeatLociFiltered,
      tandemOrArtifactFiltered: telemetry.tandemOrArtifactFiltered,
      endogenousRefFiltered,
    },
  };
}

self.onmessage = async (ev: MessageEvent<WorkerRequest>) => {
  const msg = ev.data;
  try {
    switch (msg.action) {
      case 'load_saccer3': {
        const url = msg.payload?.sacCer3Url ?? '/genomes/sacCer3.fa.gz';
        if (!cachedSacCer3Ref) {
          cachedSacCer3Ref = await MemoryFastaSource.fromGzipUrl(
            'sacCer3 (S. cerevisiae S288C)',
            url,
          );
        }
        genomeId = 'sacCer3';
        refSource = cachedSacCer3Ref;
        // 데모 모드에서 넘어온 경우 데모 트랙은 비우고 실제 로드된 로컬 BAM만 유지
        if (isDemoMode) {
          strainSources = [];
        }
        annotationIndex = getBuiltInAnnotation();
        annotationOrigin = 'builtin_s288c';
        isDemoMode = false;

        // chrV URA3 부근 기본 표시 (전체 17개 염색체 모두 드롭다운/all에서 선택 가능)
        replyOk(msg.reqId, await makeWorkspaceSummary('chrV', 115_600, 117_600));
        break;
      }

      case 'load_demo': {
        const demo = createYeastDemoDataset();
        genomeId = 'demo';
        refSource = demo.reference;
        strainSources = [demo.strain1, demo.strain2];
        annotationIndex = demo.annotation;
        annotationOrigin = 'demo';
        isDemoMode = true;

        replyOk(
          msg.reqId,
          await makeWorkspaceSummary(demo.defaultChrom, demo.defaultStart, demo.defaultEnd),
        );
        break;
      }

      case 'load_files': {
        const { fastaFile, faiFile, gffFile, strains: strainInputs } = msg.payload;
        genomeId = 'custom';
        refSource = await IndexedFastaSource.fromFiles(fastaFile, faiFile);
        strainSources = strainInputs.map(
          (st) => new GmodBamSource(st.id, st.label, st.bamFile, st.indexFile),
        );

        const contigs = refSource.getContigs();
        if (gffFile) {
          annotationIndex = new AnnotationIndex(await parseGff3File(gffFile));
          annotationOrigin = 'user_gff';
        } else if (looksLikeYeastGenome(contigs)) {
          annotationIndex = getBuiltInAnnotation();
          annotationOrigin = 'builtin_s288c';
        } else {
          annotationIndex = null;
          annotationOrigin = 'none';
        }

        isDemoMode = false;
        const firstChrom = contigs[0]?.name ?? '';
        const firstLen = contigs[0]?.length ?? 1000;
        replyOk(
          msg.reqId,
          await makeWorkspaceSummary(firstChrom, 0, Math.min(firstLen, 2_500)),
        );
        break;
      }

      case 'add_strains': {
        if (!refSource) throw new Error('먼저 Reference 지놈을 불러와주세요.');
        const newSources = msg.payload.strains.map(
          (st) => new GmodBamSource(st.id, st.label, st.bamFile, st.indexFile),
        );
        if (isDemoMode) {
          // 데모 트랙 대신 실제 사용자 트랙으로 전환 (Reference는 sacCer3로 자동 전환)
          if (cachedSacCer3Ref) {
            refSource = cachedSacCer3Ref;
            genomeId = 'sacCer3';
            annotationIndex = getBuiltInAnnotation();
            annotationOrigin = 'builtin_s288c';
          }
          strainSources = newSources;
          isDemoMode = false;
        } else {
          strainSources = [...strainSources, ...newSources];
        }
        const curChrom =
          refSource.getContigs().find((c) => c.name === msg.payload.currentChrom)?.name ??
          refSource.getContigs()[0]?.name ??
          '';
        replyOk(msg.reqId, await makeWorkspaceSummary(curChrom, 0, 2_500));
        break;
      }

      case 'remove_strain': {
        strainSources = strainSources.filter((s) => s.id !== msg.payload.id);
        const curChrom = msg.payload.currentChrom || refSource?.getContigs()[0]?.name || '';
        replyOk(msg.reqId, await makeWorkspaceSummary(curChrom, 0, 2_500));
        break;
      }

      case 'query_all_chroms': {
        replyOk(msg.reqId, await handleQueryAllChroms(msg.payload?.contigOverrides));
        break;
      }

      case 'query_window': {
        replyOk(msg.reqId, await handleQueryWindow(msg.payload));
        break;
      }

      case 'extract_region': {
        replyOk(msg.reqId, await handleExtractRegion(msg.payload));
        break;
      }

      case 'scan_chromosome': {
        replyOk(msg.reqId, await handleScanChromosome(msg.payload));
        break;
      }

      case 'cancel_scan': {
        activeScanToken++;
        replyOk(msg.reqId, true);
        break;
      }

      case 'search_genes': {
        const rawList = annotationIndex
          ? annotationIndex.searchGenes(msg.payload.query, msg.payload.limit ?? 25)
          : [];
        const mappedList = rawList.map((f) => ({
          ...f,
          chrom: refSource?.resolveContigName(f.chrom) ?? f.chrom,
        }));
        replyOk(msg.reqId, mappedList);
        break;
      }

      case 'run_blast': {
        replyOk(msg.reqId, await handleRunBlast(msg.payload));
        break;
      }

      case 'cancel_blast': {
        activeBlastToken++;
        replyOk(msg.reqId, true);
        break;
      }

      case 'scan_barcodes': {
        replyOk(msg.reqId, await handleScanBarcodes(msg.payload));
        break;
      }

      case 'cancel_barcode_scan': {
        activeBarcodeToken++;
        replyOk(msg.reqId, true);
        break;
      }
    }
  } catch (err) {
    replyErr(msg.reqId, err);
  }
};
