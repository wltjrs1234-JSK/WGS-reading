import { GmodBamSource, type AlignmentSource } from '../bio/bamSource';
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
import { AnnotationIndex, parseGff3File } from '../bio/gff3Parser';
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

let genomeId = 'sacCer3';
let refSource: ReferenceSequenceSource | null = null;
let cachedSacCer3Ref: MemoryFastaSource | null = null;
let cachedBuiltInAnnotation: AnnotationIndex | null = null;
let strainSources: AlignmentSource[] = [];
let annotationIndex: AnnotationIndex | null = null;
let annotationOrigin: AnnotationOrigin = 'none';
let isDemoMode = false;
let activeScanToken = 0;

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
    }
  } catch (err) {
    replyErr(msg.reqId, err);
  }
};
