import { forEachCigarOp } from './cigar';
import type {
  AlignedRead,
  ContigInfo,
  GffFeature,
  SvRawReadSignal,
} from './types';
import { isMitochondrialChrom, normalizeYeastChromName } from './yeastAliases';

export type CnvPloidyPreference = 'auto' | 1 | 2 | 4;
export type CnvScanMode = 'fast_sample' | 'deep_sv';

export interface ProbeBinMeasurement {
  chrom: string;
  start: number; // bin genomic start
  end: number;   // bin genomic end
  probeStart: number;
  probeEnd: number;
  rawDepth: number;
  sampledBp: number;
  hetSites: number;
  bafs: number[]; // non-ref allele frequencies at heterozygous sites
}

export interface NormalizedCnvBin {
  chrom: string;
  start: number;
  end: number;
  rawDepth: number;
  copyRatio: number;
  estCopies: number;
  log2Ratio: number;
  hetSites: number;
  sampledBp: number;
  hetPerKb: number;
  meanBaf: number | null;
}

export type ChromCopyStatus =
  | 'normal'
  | 'gain'
  | 'loss'
  | 'mitochondrial_low'
  | 'mitochondrial_high'
  | 'unmapped';

export interface ChromCopySummary {
  chrom: string;
  matchedChrom: string | null;
  length: number;
  medianDepth: number;
  meanDepth: number;
  copyRatio: number;
  estCopies: number;
  status: ChromCopyStatus;
  hetPerKb: number;
  bafMedian: number | null;
  note?: string;
}

export interface CnvSegmentEvent {
  id: string;
  strainId: string;
  strainLabel: string;
  chrom: string;
  start: number;
  end: number;
  lengthBp: number;
  meanCopyRatio: number;
  estCopies: number;
  direction: 'amplification' | 'deletion' | 'homozygous_deletion';
  overlappingGenes: string[];
}

export interface LohSegmentEvent {
  id: string;
  strainId: string;
  strainLabel: string;
  chrom: string;
  start: number;
  end: number;
  lengthBp: number;
  binCount: number;
  meanHetPerKb: number;
  genomeMedianHetPerKb: number;
  meanCopyRatio: number;
  lohType: 'cn_neutral_loh' | 'deletion_loh';
  overlappingGenes: string[];
}

export interface SvBreakpointEvent {
  id: string;
  strainId: string;
  strainLabel: string;
  svType: 'TRANSLOCATION' | 'LARGE_DELETION' | 'INVERSION_OR_DUP';
  chrom1: string;
  pos1: number; // 0-based representative breakpoint
  chrom2: string;
  pos2: number; // 0-based representative breakpoint
  spanBp: number | null; // null if inter-chromosomal
  supportPairs: number;
  splitReads: number;
  evidence: 'mate_pair' | 'split_read' | 'both';
  genes1: string[];
  genes2: string[];
  specificity: 'strain_specific' | 'shared' | 'single_strain';
  sharedWithStrains: string[];
}

export interface StrainCnvSvProfile {
  strainId: string;
  strainLabel: string;
  requestedPloidy: CnvPloidyPreference;
  effectivePloidy: 1 | 2 | 4;
  detectedPloidy: 1 | 2 | 4;
  ploidyConfidence: 'high' | 'medium' | 'low';
  ploidyReason: string;
  genomeMedianDepth: number;
  genomeMedianHetPerKb: number;
  mitoToNuclearRatio: number;
  totalHetSitesSampled: number;
  bafHistogram: number[]; // 20 bins covering [0.0..1.0) in 0.05 steps
  chromSummaries: ChromCopySummary[];
  bins: NormalizedCnvBin[];
  cnvSegments: CnvSegmentEvent[];
  lohSegments: LohSegmentEvent[];
  svEvents: SvBreakpointEvent[];
}

export interface CnvSvAnalysisResult {
  scanMode: CnvScanMode;
  ploidyPreference: CnvPloidyPreference;
  elapsedMs: number;
  contigs: ContigInfo[];
  strains: StrainCnvSvProfile[];
}

function medianOf(nums: number[]): number {
  if (nums.length === 0) return 0;
  const sorted = [...nums].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 1
    ? sorted[mid]
    : (sorted[mid - 1] + sorted[mid]) / 2;
}

function isRdnaRepeatWindow(chrom: string, start: number, end: number): boolean {
  return normalizeYeastChromName(chrom) === 'chrXII' && end >= 450_000 && start <= 492_000;
}

/**
 * 400~800 bp 프로브 구간에서 실제 매핑 Depth와 이형접합 B-Allele Frequency(BAF)를 정량 측정합니다.
 */
export function measureProbeWindow(
  chrom: string,
  binStart: number,
  binEnd: number,
  probeStart: number,
  probeEnd: number,
  refSeq: string,
  reads: AlignedRead[],
  minMapQ = 10,
): ProbeBinMeasurement {
  const len = Math.max(1, probeEnd - probeStart);
  const aArr = new Int32Array(len);
  const cArr = new Int32Array(len);
  const gArr = new Int32Array(len);
  const tArr = new Int32Array(len);
  const dArr = new Int32Array(len);

  for (const r of reads) {
    if ((r.flags & 0x704) !== 0 || r.mapq < minMapQ) continue;
    const seq = r.seq ? r.seq.toUpperCase() : '';
    let refPos = r.start;
    let queryPos = 0;

    forEachCigarOp(r.cigar, (op, opLen) => {
      switch (op) {
        case 'M':
        case '=':
        case 'X': {
          const ovS = refPos > probeStart ? refPos : probeStart;
          const ovE = refPos + opLen < probeEnd ? refPos + opLen : probeEnd;
          for (let rp = ovS; rp < ovE; rp++) {
            const qp = queryPos + (rp - refPos);
            const raw = qp < seq.length ? seq[qp] : 'N';
            const b = raw === '=' ? (refSeq[rp - probeStart] ?? 'N').toUpperCase() : raw;
            const idx = rp - probeStart;
            if (b === 'A') aArr[idx]++;
            else if (b === 'C') cArr[idx]++;
            else if (b === 'G') gArr[idx]++;
            else if (b === 'T') tArr[idx]++;
            dArr[idx]++;
          }
          refPos += opLen;
          queryPos += opLen;
          break;
        }
        case 'D':
        case 'N': {
          if (op === 'D') {
            const ovS = refPos > probeStart ? refPos : probeStart;
            const ovE = refPos + opLen < probeEnd ? refPos + opLen : probeEnd;
            for (let rp = ovS; rp < ovE; rp++) {
              dArr[rp - probeStart]++;
            }
          }
          refPos += opLen;
          break;
        }
        case 'I':
        case 'S':
          queryPos += opLen;
          break;
      }
    });
  }

  let depthSum = 0;
  let hetSites = 0;
  const bafs: number[] = [];

  for (let i = 0; i < len; i++) {
    const dp = dArr[i];
    depthSum += dp;
    if (dp < 10) continue;

    const acgtSum = aArr[i] + cArr[i] + gArr[i] + tArr[i];
    if (acgtSum < 10) continue;

    const counts = [
      { base: 'A', n: aArr[i] },
      { base: 'C', n: cArr[i] },
      { base: 'G', n: gArr[i] },
      { base: 'T', n: tArr[i] },
    ].sort((x, y) => y.n - x.n);

    // 2번째 다수 염기가 >=3 reads & >=15% 비율이면 진정한 이형접합 SNP로 판정
    if (counts[1].n >= 3 && counts[1].n / acgtSum >= 0.15) {
      hetSites++;
      const refB = (refSeq[i] ?? 'N').toUpperCase();
      const nonRefSum =
        refB === 'A'
          ? acgtSum - aArr[i]
          : refB === 'C'
            ? acgtSum - cArr[i]
            : refB === 'G'
              ? acgtSum - gArr[i]
              : refB === 'T'
                ? acgtSum - tArr[i]
                : counts[1].n;
      const baf = Math.min(0.99, Math.max(0.01, nonRefSum / acgtSum));
      bafs.push(baf);
    }
  }

  return {
    chrom,
    start: binStart,
    end: binEnd,
    probeStart,
    probeEnd,
    rawDepth: depthSum / len,
    sampledBp: len,
    hetSites,
    bafs,
  };
}

/**
 * 전체 핵 게놈의 이형접합 밀도 및 BAF 분포(0.5 vs 0.25/0.75)를 분석해 1n/2n/4n 배수성을 자동 추정합니다.
 */
export function estimatePloidyFromBafs(
  bafs: number[],
  totalSampledBp: number,
  meanDepth: number,
): { ploidy: 1 | 2 | 4; confidence: 'high' | 'medium' | 'low'; reason: string } {
  if (meanDepth < 8 || totalSampledBp < 10_000) {
    return {
      ploidy: 1,
      confidence: 'low',
      reason: `평균 커버리지(${meanDepth.toFixed(1)}x)가 낮아 기본 반수체(1n) 기준으로 계산했습니다.`,
    };
  }

  const hetPerKb = (bafs.length / Math.max(1, totalSampledBp)) * 1_000;
  if (bafs.length < 12 || hetPerKb < 0.05) {
    return {
      ploidy: 1,
      confidence: bafs.length < 4 ? 'high' : 'medium',
      reason: `이형접합 사이트가 거의 없음 (${bafs.length}개, ${hetPerKb.toFixed(2)}/kb) → 반수체(1n) 추정`,
    };
  }

  let halfCount = 0;
  let quarterCount = 0;
  for (const b of bafs) {
    if (b >= 0.42 && b <= 0.58) halfCount++;
    else if ((b >= 0.18 && b <= 0.34) || (b >= 0.66 && b <= 0.82)) quarterCount++;
  }

  const halfRatio = halfCount / bafs.length;
  const quarterRatio = quarterCount / bafs.length;

  // 사배체(4n): ABBB / AAAB 대립유전자 조합으로 0.25·0.75 피크가 0.50 피크와 비슷하거나 더 우세
  if (bafs.length >= 25 && quarterRatio >= 0.32 && quarterRatio >= halfRatio * 0.9) {
    return {
      ploidy: 4,
      confidence: quarterRatio >= 0.4 ? 'high' : 'medium',
      reason: `BAF 0.25·0.75 사배체 대립유전자 비율 우세 (${Math.round(quarterRatio * 100)}% vs 0.50 ${Math.round(halfRatio * 100)}%, ${hetPerKb.toFixed(2)} het/kb) → 사배체(4n) 추정`,
    };
  }

  return {
    ploidy: 2,
    confidence: halfRatio >= 0.35 ? 'high' : 'medium',
    reason: `이형접합 BAF 0.50 중심 피크 (${Math.round(halfRatio * 100)}%, ${hetPerKb.toFixed(2)} het/kb) → 이배체(2n) 추정`,
  };
}

/**
 * 측정된 빈(bin) 배열로부터 복제수 정규화, 염색체 이수성, CNV 구간, LOH 구간을 모두 계산합니다.
 */
export function analyzeCnvAndLohForStrain(
  strainId: string,
  strainLabel: string,
  contigs: ContigInfo[],
  matchedChromByRef: Record<string, string | null>,
  measurements: ProbeBinMeasurement[],
  ploidyPreference: CnvPloidyPreference,
  lookupGenes: (chrom: string, start: number, end: number) => GffFeature[],
): Omit<StrainCnvSvProfile, 'svEvents'> {
  const nuclearPositive = measurements.filter(
    (m) =>
      !isMitochondrialChrom(m.chrom) &&
      !isRdnaRepeatWindow(m.chrom, m.start, m.end) &&
      m.rawDepth > 0,
  );

  const genomeMedianDepth = Math.max(
    0,
    medianOf(nuclearPositive.map((m) => m.rawDepth)),
  );

  const nuclearBafs: number[] = [];
  let nuclearSampledBp = 0;
  for (const m of measurements) {
    if (isMitochondrialChrom(m.chrom) || isRdnaRepeatWindow(m.chrom, m.start, m.end)) continue;
    nuclearSampledBp += m.sampledBp;
    for (const b of m.bafs) nuclearBafs.push(b);
  }

  const ploidyEst = estimatePloidyFromBafs(
    nuclearBafs,
    nuclearSampledBp,
    genomeMedianDepth,
  );
  const effectivePloidy: 1 | 2 | 4 =
    ploidyPreference === 'auto' ? ploidyEst.ploidy : ploidyPreference;

  const bafHistogram = new Array<number>(20).fill(0);
  for (const b of nuclearBafs) {
    const idx = Math.max(0, Math.min(19, Math.floor(b * 20)));
    bafHistogram[idx]++;
  }

  const safeBaseline = genomeMedianDepth > 0 ? genomeMedianDepth : 1;

  const bins: NormalizedCnvBin[] = measurements.map((m) => {
    const copyRatio = genomeMedianDepth > 0 ? m.rawDepth / safeBaseline : 0;
    const estCopies = copyRatio * effectivePloidy;
    const log2Ratio =
      genomeMedianDepth > 0 ? Math.log2(Math.max(0.05, copyRatio)) : 0;
    const hetPerKb = (m.hetSites / Math.max(1, m.sampledBp)) * 1_000;
    const meanBaf =
      m.bafs.length > 0
        ? m.bafs.reduce((acc, v) => acc + v, 0) / m.bafs.length
        : null;
    return {
      chrom: m.chrom,
      start: m.start,
      end: m.end,
      rawDepth: Math.round(m.rawDepth * 10) / 10,
      copyRatio: Math.round(copyRatio * 100) / 100,
      estCopies: Math.round(estCopies * 100) / 100,
      log2Ratio: Math.round(log2Ratio * 100) / 100,
      hetSites: m.hetSites,
      sampledBp: m.sampledBp,
      hetPerKb: Math.round(hetPerKb * 100) / 100,
      meanBaf,
    };
  });

  // 염색체별 이수성 요약
  let mitoToNuclearRatio = 0;
  const chromSummaries: ChromCopySummary[] = contigs.map((c) => {
    const matchedChrom = matchedChromByRef[c.name] ?? null;
    const chrBins = bins.filter((b) => b.chrom === c.name);
    if (!matchedChrom || chrBins.length === 0) {
      return {
        chrom: c.name,
        matchedChrom,
        length: c.length,
        medianDepth: 0,
        meanDepth: 0,
        copyRatio: 0,
        estCopies: 0,
        status: 'unmapped',
        hetPerKb: 0,
        bafMedian: null,
      };
    }

    const cleanBins = isMitochondrialChrom(c.name)
      ? chrBins
      : chrBins.filter((b) => !isRdnaRepeatWindow(b.chrom, b.start, b.end));
    const useBins = cleanBins.length > 0 ? cleanBins : chrBins;

    const medDepth = medianOf(useBins.map((b) => b.rawDepth));
    const avgDepth =
      useBins.reduce((acc, b) => acc + b.rawDepth, 0) / useBins.length;
    const copyRatio = genomeMedianDepth > 0 ? medDepth / safeBaseline : 0;
    const estCopies = copyRatio * effectivePloidy;

    const chrMeas = measurements.filter(
      (m) =>
        m.chrom === c.name &&
        (isMitochondrialChrom(c.name) || !isRdnaRepeatWindow(m.chrom, m.start, m.end)),
    );
    const chrBafs: number[] = [];
    let chrHet = 0;
    let chrBp = 0;
    for (const m of chrMeas) {
      chrHet += m.hetSites;
      chrBp += m.sampledBp;
      for (const b of m.bafs) chrBafs.push(b);
    }
    const hetPerKb = Math.round(((chrHet / Math.max(1, chrBp)) * 1_000) * 100) / 100;
    const bafMedian = chrBafs.length > 0 ? Math.round(medianOf(chrBafs) * 100) / 100 : null;

    let status: ChromCopyStatus = 'normal';
    let note: string | undefined;

    if (isMitochondrialChrom(c.name)) {
      mitoToNuclearRatio = Math.round(copyRatio * 10) / 10;
      if (genomeMedianDepth >= 4 && copyRatio < 0.35) {
        status = 'mitochondrial_low';
        note = 'Petite / ρ0(mtDNA 소실·극저카피) 의심';
      } else if (genomeMedianDepth >= 4 && copyRatio > 80) {
        status = 'mitochondrial_high';
        note = 'mtDNA 고복제수';
      }
    } else if (genomeMedianDepth >= 4) {
      const delta = estCopies - effectivePloidy;
      if (delta >= 0.48) {
        status = 'gain';
        note = `이수성 증가 (${effectivePloidy}n → 약 ${Math.round(estCopies)} copy)`;
      } else if (delta <= -0.48) {
        status = 'loss';
        note = `이수성 감소 (${effectivePloidy}n → 약 ${Math.round(estCopies)} copy)`;
      }
    }

    return {
      chrom: c.name,
      matchedChrom,
      length: c.length,
      medianDepth: Math.round(medDepth * 10) / 10,
      meanDepth: Math.round(avgDepth * 10) / 10,
      copyRatio: Math.round(copyRatio * 100) / 100,
      estCopies: Math.round(estCopies * 10) / 10,
      status,
      hetPerKb,
      bafMedian,
      note,
    };
  });

  // 구간 CNV 세그먼트 (연속 >= 2개 빈, 최소 >= 15 kb 또는 단일 큰 결실/증폭 블록)
  const cnvSegments: CnvSegmentEvent[] = [];
  if (genomeMedianDepth >= 4) {
    for (const c of contigs) {
      if (isMitochondrialChrom(c.name)) continue;
      const cBins = bins.filter(
        (b) => b.chrom === c.name && !isRdnaRepeatWindow(b.chrom, b.start, b.end),
      );
      let curDir: CnvSegmentEvent['direction'] | null = null;
      let run: NormalizedCnvBin[] = [];

      const flushCnv = () => {
        if (!curDir || run.length < 2) {
          run = [];
          curDir = null;
          return;
        }
        const segStart = run[0].start;
        const segEnd = run[run.length - 1].end;
        if (segEnd - segStart < 10_000) {
          run = [];
          curDir = null;
          return;
        }
        const meanRatio = run.reduce((a, b) => a + b.copyRatio, 0) / run.length;
        const estCopies = meanRatio * effectivePloidy;
        const feats = lookupGenes(c.name, segStart, segEnd)
          .filter((f) => f.type === 'gene')
          .slice(0, 5)
          .map((f) => f.name);
        cnvSegments.push({
          id: `${strainId}_cnv_${c.name}_${segStart}_${segEnd}`,
          strainId,
          strainLabel,
          chrom: c.name,
          start: segStart,
          end: segEnd,
          lengthBp: segEnd - segStart,
          meanCopyRatio: Math.round(meanRatio * 100) / 100,
          estCopies: Math.round(estCopies * 10) / 10,
          direction: curDir,
          overlappingGenes: feats,
        });
        run = [];
        curDir = null;
      };

      for (const b of cBins) {
        const delta = b.estCopies - effectivePloidy;
        let dir: CnvSegmentEvent['direction'] | null = null;
        if (b.copyRatio <= 0.12) dir = 'homozygous_deletion';
        else if (delta <= -0.62) dir = 'deletion';
        else if (delta >= 0.65) dir = 'amplification';

        if (dir && dir === curDir) {
          run.push(b);
        } else {
          flushCnv();
          if (dir) {
            curDir = dir;
            run = [b];
          }
        }
      }
      flushCnv();
    }
  }

  // LOH 구간 탐색 (2n / 4n이면서 게놈 전체에 충분한 이형접합 배경이 있는 균주)
  const positiveHetBins = bins
    .filter(
      (b) =>
        !isMitochondrialChrom(b.chrom) &&
        !isRdnaRepeatWindow(b.chrom, b.start, b.end) &&
        b.hetPerKb > 0,
    )
    .map((b) => b.hetPerKb);
  const genomeMedianHetPerKb =
    Math.round(medianOf(positiveHetBins) * 100) / 100;

  const lohSegments: LohSegmentEvent[] = [];
  if (effectivePloidy >= 2 && genomeMedianHetPerKb >= 0.2 && nuclearBafs.length >= 25) {
    const lowHetCut = Math.max(0.02, genomeMedianHetPerKb * 0.15);
    for (const c of contigs) {
      if (isMitochondrialChrom(c.name)) continue;
      const cBins = bins.filter(
        (b) => b.chrom === c.name && !isRdnaRepeatWindow(b.chrom, b.start, b.end),
      );
      let run: NormalizedCnvBin[] = [];

      const flushLoh = () => {
        if (run.length < 3) {
          run = [];
          return;
        }
        const segStart = run[0].start;
        const segEnd = run[run.length - 1].end;
        const meanHet = run.reduce((a, b) => a + b.hetPerKb, 0) / run.length;
        const meanRatio = run.reduce((a, b) => a + b.copyRatio, 0) / run.length;
        const feats = lookupGenes(c.name, segStart, segEnd)
          .filter((f) => f.type === 'gene')
          .slice(0, 5)
          .map((f) => f.name);
        lohSegments.push({
          id: `${strainId}_loh_${c.name}_${segStart}_${segEnd}`,
          strainId,
          strainLabel,
          chrom: c.name,
          start: segStart,
          end: segEnd,
          lengthBp: segEnd - segStart,
          binCount: run.length,
          meanHetPerKb: Math.round(meanHet * 100) / 100,
          genomeMedianHetPerKb,
          meanCopyRatio: Math.round(meanRatio * 100) / 100,
          lohType: meanRatio < 0.72 ? 'deletion_loh' : 'cn_neutral_loh',
          overlappingGenes: feats,
        });
        run = [];
      };

      for (const b of cBins) {
        // 충분한 뎁스가 있는데도 이형접합이 완전히 사라진 연속 구간
        if (b.rawDepth >= Math.max(5, safeBaseline * 0.25) && b.hetPerKb <= lowHetCut) {
          run.push(b);
        } else {
          flushLoh();
        }
      }
      flushLoh();
    }
  }

  return {
    strainId,
    strainLabel,
    requestedPloidy: ploidyPreference,
    effectivePloidy,
    detectedPloidy: ploidyEst.ploidy,
    ploidyConfidence: ploidyEst.confidence,
    ploidyReason: ploidyEst.reason,
    genomeMedianDepth: Math.round(genomeMedianDepth * 10) / 10,
    genomeMedianHetPerKb,
    mitoToNuclearRatio,
    totalHetSitesSampled: nuclearBafs.length,
    bafHistogram,
    chromSummaries,
    bins,
    cnvSegments,
    lohSegments,
  };
}

interface RawSvEdge {
  readName: string;
  svType: SvBreakpointEvent['svType'];
  c1: string;
  p1: number;
  c2: string;
  p2: number;
  isSplit: boolean;
}

function canonicalEdge(
  readName: string,
  svType: SvBreakpointEvent['svType'],
  chrA: string,
  posA: number,
  chrB: string,
  posB: number,
  isSplit: boolean,
): RawSvEdge {
  if (chrA.localeCompare(chrB) < 0 || (chrA === chrB && posA <= posB)) {
    return { readName, svType, c1: chrA, p1: posA, c2: chrB, p2: posB, isSplit };
  }
  return { readName, svType, c1: chrB, p1: posB, c2: chrA, p2: posA, isSplit };
}

/**
 * Discordant Mate Pair 및 SA Split-Read 신호를 클러스터링하여 전좌·대결실·역위 이벤트를 추출합니다.
 */
export function clusterSvSignals(
  strainId: string,
  strainLabel: string,
  rawSignals: SvRawReadSignal[],
  resolveCanonicalChrom: (bamChrom: string) => string,
  lookupGenes: (chrom: string, start: number, end: number) => GffFeature[],
  options: { minSupport?: number; excludeRdnaSelf?: boolean } = {},
): SvBreakpointEvent[] {
  const minSupport = options.minSupport ?? 4;
  const excludeRdna = options.excludeRdnaSelf ?? true;
  const edges: RawSvEdge[] = [];

  for (const sig of rawSignals) {
    if (sig.mapq < 15) continue;
    const cA = resolveCanonicalChrom(sig.chrom);
    const pA = sig.strand === '+' ? sig.end : sig.start;

    // 1. Chimeric Split-Read (SA:Z)
    if (sig.saTag) {
      const entries = sig.saTag.split(';').filter(Boolean);
      for (const ent of entries) {
        const parts = ent.split(',');
        if (parts.length < 5) continue;
        const saChr = resolveCanonicalChrom(parts[0]);
        const saPos = Math.max(0, Number(parts[1]) - 1);
        const saMq = Number(parts[4]);
        if (!Number.isFinite(saPos) || !Number.isFinite(saMq) || saMq < 15) continue;
        if (saChr === cA && Math.abs(saPos - pA) < 2_000) continue;

        const svType: SvBreakpointEvent['svType'] =
          saChr !== cA
            ? 'TRANSLOCATION'
            : parts[2] !== sig.strand
              ? 'INVERSION_OR_DUP'
              : 'LARGE_DELETION';
        edges.push(canonicalEdge(sig.name, svType, cA, pA, saChr, saPos, true));
      }
    }

    // 2. Discordant Mate Pair
    if (sig.isPaired && sig.mateChrom) {
      const cB = resolveCanonicalChrom(sig.mateChrom);
      if (cB !== cA) {
        edges.push(canonicalEdge(sig.name, 'TRANSLOCATION', cA, sig.start, cB, sig.matePos, false));
      } else {
        const dist = Math.abs(sig.matePos - sig.start);
        if (sig.strand === sig.mateStrand && dist >= 1_000) {
          edges.push(
            canonicalEdge(sig.name, 'INVERSION_OR_DUP', cA, sig.start, cB, sig.matePos, false),
          );
        } else if (dist >= 2_500) {
          edges.push(
            canonicalEdge(sig.name, 'LARGE_DELETION', cA, sig.start, cB, sig.matePos, false),
          );
        }
      }
    }
  }

  interface Cluster {
    svType: SvBreakpointEvent['svType'];
    c1: string;
    c2: string;
    p1Sum: number;
    p2Sum: number;
    p1Avg: number;
    p2Avg: number;
    edgeCount: number;
    pairNames: Set<string>;
    splitNames: Set<string>;
  }

  const clusters: Cluster[] = [];
  const tol = 750;

  for (const e of edges) {
    if (
      excludeRdna &&
      isRdnaRepeatWindow(e.c1, e.p1, e.p1 + 1) &&
      isRdnaRepeatWindow(e.c2, e.p2, e.p2 + 1)
    ) {
      continue;
    }

    const match = clusters.find(
      (cl) =>
        cl.svType === e.svType &&
        cl.c1 === e.c1 &&
        cl.c2 === e.c2 &&
        Math.abs(cl.p1Avg - e.p1) <= tol &&
        Math.abs(cl.p2Avg - e.p2) <= tol,
    );
    if (match) {
      match.edgeCount++;
      match.p1Sum += e.p1;
      match.p2Sum += e.p2;
      match.p1Avg = Math.round(match.p1Sum / match.edgeCount);
      match.p2Avg = Math.round(match.p2Sum / match.edgeCount);
      if (e.isSplit) match.splitNames.add(e.readName);
      else match.pairNames.add(e.readName);
    } else {
      clusters.push({
        svType: e.svType,
        c1: e.c1,
        c2: e.c2,
        p1Sum: e.p1,
        p2Sum: e.p2,
        p1Avg: e.p1,
        p2Avg: e.p2,
        edgeCount: 1,
        pairNames: e.isSplit ? new Set() : new Set([e.readName]),
        splitNames: e.isSplit ? new Set([e.readName]) : new Set(),
      });
    }
  }

  const events: SvBreakpointEvent[] = [];
  for (const cl of clusters) {
    const uniqueReads = new Set([...cl.pairNames, ...cl.splitNames]).size;
    if (uniqueReads < minSupport) continue;

    const genes1 = lookupGenes(cl.c1, Math.max(0, cl.p1Avg - 500), cl.p1Avg + 500)
      .filter((f) => f.type === 'gene')
      .slice(0, 3)
      .map((f) => f.name);
    const genes2 = lookupGenes(cl.c2, Math.max(0, cl.p2Avg - 500), cl.p2Avg + 500)
      .filter((f) => f.type === 'gene')
      .slice(0, 3)
      .map((f) => f.name);

    const hasPairs = cl.pairNames.size > 0;
    const hasSplits = cl.splitNames.size > 0;

    events.push({
      id: `${strainId}_sv_${cl.svType}_${cl.c1}_${cl.p1Avg}_${cl.c2}_${cl.p2Avg}`,
      strainId,
      strainLabel,
      svType: cl.svType,
      chrom1: cl.c1,
      pos1: cl.p1Avg,
      chrom2: cl.c2,
      pos2: cl.p2Avg,
      spanBp: cl.c1 === cl.c2 ? Math.abs(cl.p2Avg - cl.p1Avg) : null,
      supportPairs: cl.pairNames.size,
      splitReads: cl.splitNames.size,
      evidence: hasPairs && hasSplits ? 'both' : hasSplits ? 'split_read' : 'mate_pair',
      genes1,
      genes2,
      specificity: 'single_strain',
      sharedWithStrains: [],
    });
  }

  return events.sort(
    (a, b) => b.supportPairs + b.splitReads - (a.supportPairs + a.splitReads),
  );
}

/**
 * 다수 균주 비교 시 SV 이벤트가 특정 균주에만 존재하는지 표시합니다.
 */
export function annotateCrossStrainSvSpecificity(
  profiles: StrainCnvSvProfile[],
): StrainCnvSvProfile[] {
  if (profiles.length <= 1) return profiles;
  const allSvs = profiles.flatMap((p) => p.svEvents);
  return profiles.map((p) => ({
    ...p,
    svEvents: p.svEvents.map((ev) => {
      const peers = new Set<string>();
      for (const other of allSvs) {
        if (other.strainId === ev.strainId) continue;
        if (
          other.svType === ev.svType &&
          other.chrom1 === ev.chrom1 &&
          other.chrom2 === ev.chrom2 &&
          Math.abs(other.pos1 - ev.pos1) <= 1_200 &&
          Math.abs(other.pos2 - ev.pos2) <= 1_200
        ) {
          peers.add(other.strainLabel);
        }
      }
      return {
        ...ev,
        specificity: peers.size === 0 ? 'strain_specific' : 'shared',
        sharedWithStrains: [...peers],
      };
    }),
  }));
}

export function exportCnvSvToTsv(result: CnvSvAnalysisResult): string {
  const lines: string[] = [
    '#Section\tStrain\tEffective_Ploidy\tDetected_Ploidy\tChrom1\tStart1_or_Pos1\tChrom2\tEnd1_or_Pos2\tType_or_Status\tValue_Copies_or_Pairs\tDetail_Genes',
  ];

  for (const st of result.strains) {
    for (const c of st.chromSummaries) {
      lines.push(
        [
          'CHROM_COPY',
          st.strainLabel,
          `${st.effectivePloidy}n`,
          `${st.detectedPloidy}n`,
          c.chrom,
          1,
          c.chrom,
          c.length,
          c.status,
          `${c.estCopies}x (ratio=${c.copyRatio})`,
          c.note ?? `het/kb=${c.hetPerKb}`,
        ].join('\t'),
      );
    }
    for (const seg of st.cnvSegments) {
      lines.push(
        [
          'CNV_SEGMENT',
          st.strainLabel,
          `${st.effectivePloidy}n`,
          `${st.detectedPloidy}n`,
          seg.chrom,
          seg.start + 1,
          seg.chrom,
          seg.end,
          seg.direction,
          `${seg.estCopies}x (${seg.lengthBp}bp)`,
          seg.overlappingGenes.join(';') || '-',
        ].join('\t'),
      );
    }
    for (const loh of st.lohSegments) {
      lines.push(
        [
          'LOH_SEGMENT',
          st.strainLabel,
          `${st.effectivePloidy}n`,
          `${st.detectedPloidy}n`,
          loh.chrom,
          loh.start + 1,
          loh.chrom,
          loh.end,
          loh.lohType,
          `${loh.meanHetPerKb}/kb (${loh.lengthBp}bp)`,
          loh.overlappingGenes.join(';') || '-',
        ].join('\t'),
      );
    }
    for (const sv of st.svEvents) {
      lines.push(
        [
          'SV_BREAKPOINT',
          st.strainLabel,
          `${st.effectivePloidy}n`,
          `${st.detectedPloidy}n`,
          sv.chrom1,
          sv.pos1 + 1,
          sv.chrom2,
          sv.pos2 + 1,
          sv.svType,
          `pairs=${sv.supportPairs},split=${sv.splitReads}`,
          `G1:${sv.genes1.join(';') || '-'}|G2:${sv.genes2.join(';') || '-'}`,
        ].join('\t'),
      );
    }
  }

  return '\uFEFF' + lines.join('\n') + '\n';
}

export function exportCnvSvToBed(result: CnvSvAnalysisResult): string {
  const lines = [
    '#chrom\tchromStart\tchromEnd\tname\tscore\tstrand',
  ];
  for (const st of result.strains) {
    for (const seg of st.cnvSegments) {
      lines.push(
        `${seg.chrom}\t${seg.start}\t${seg.end}\t${st.strainLabel}_${seg.direction}_${seg.estCopies}x\t${Math.min(1000, Math.round(seg.estCopies * 100))}\t.`,
      );
    }
    for (const loh of st.lohSegments) {
      lines.push(
        `${loh.chrom}\t${loh.start}\t${loh.end}\t${st.strainLabel}_${loh.lohType}\t500\t.`,
      );
    }
    for (const sv of st.svEvents) {
      lines.push(
        `${sv.chrom1}\t${Math.max(0, sv.pos1 - 250)}\t${sv.pos1 + 250}\t${st.strainLabel}_${sv.svType}_to_${sv.chrom2}:${sv.pos2 + 1}\t${Math.min(1000, (sv.supportPairs + sv.splitReads) * 50)}\t.`,
      );
    }
  }
  return lines.join('\n') + '\n';
}

