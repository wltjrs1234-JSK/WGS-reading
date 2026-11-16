import { forEachCigarOp } from './cigar';
import type {
  AlignedColumnMeta,
  AlignedRead,
  BaseCounts,
  ConsensusOptions,
  ConsensusPositionCall,
  DiffCategory,
  ThreeWayRegionAnalysis,
  VariantRecord,
} from './types';

const IUPAC_MAP: Record<string, string> = {
  A: 'A',
  C: 'C',
  G: 'G',
  T: 'T',
  AG: 'R',
  CT: 'Y',
  CG: 'S',
  AT: 'W',
  GT: 'K',
  AC: 'M',
  CGT: 'B',
  AGT: 'D',
  ACT: 'H',
  ACG: 'V',
  ACGT: 'N',
};

function allelesToIupac(alleles: string[]): string {
  const key = [...new Set(alleles)].sort().join('');
  return IUPAC_MAP[key] ?? 'N';
}

function emptyCounts(): BaseCounts {
  return { A: 0, C: 0, G: 0, T: 0, N: 0, del: 0 };
}

/**
 * 단일 균주의 특정 구간 [start, end)에 대한 Pileup 및 CIGAR 기반 Consensus 계산
 */
export function computeConsensusForRegion(
  refSeq: string,
  reads: AlignedRead[],
  start: number,
  end: number,
  options: ConsensusOptions,
): ConsensusPositionCall[] {
  const len = Math.max(0, end - start);
  const countsArr: BaseCounts[] = Array.from({ length: len }, emptyCounts);
  const insertionMaps: Array<Map<string, number>> = Array.from({ length: len }, () => new Map());

  for (const read of reads) {
    if (read.mapq < options.minMapQ) continue;
    // unmapped(0x4), secondary(0x100), QC-fail(0x200), duplicate(0x400) 필터링
    if ((read.flags & 0x704) !== 0) continue;
    if (read.end <= start || read.start >= end) continue;

    // per-base 객체/배열 할당 없이 직접 CIGAR를 순회해 고심도 GC 스톨을 제거합니다.
    const seq = read.seq && read.seq !== '*' ? read.seq.toUpperCase() : '';
    const qual = read.qual;
    const minBaseQ = options.minBaseQ;
    let refPos = read.start;
    let queryPos = 0;

    forEachCigarOp(read.cigar, (op, cigarLen) => {
      switch (op) {
        case 'M':
        case '=':
        case 'X': {
          const segStart = refPos;
          const overlapStart = segStart > start ? segStart : start;
          const overlapEnd = segStart + cigarLen < end ? segStart + cigarLen : end;
          for (let rp = overlapStart; rp < overlapEnd; rp++) {
            const qp = queryPos + (rp - segStart);
            if (qual && minBaseQ > 0 && qp < qual.length && qual[qp] < 255 && qual[qp] < minBaseQ) {
              continue;
            }
            const rawBase = qp < seq.length ? seq[qp] : 'N';
            const base =
              rawBase === '=' ? (refSeq[rp - start] ?? 'N').toUpperCase() : rawBase;
            const c = countsArr[rp - start];
            if (base === 'A') c.A++;
            else if (base === 'C') c.C++;
            else if (base === 'G') c.G++;
            else if (base === 'T') c.T++;
            else c.N++;
          }
          refPos += cigarLen;
          queryPos += cigarLen;
          break;
        }
        case 'I': {
          const anchorIdx = refPos - 1 - start;
          if (anchorIdx >= 0 && anchorIdx < len) {
            const insSeq =
              queryPos + cigarLen <= seq.length
                ? seq.slice(queryPos, queryPos + cigarLen)
                : 'N'.repeat(cigarLen);
            const m = insertionMaps[anchorIdx];
            m.set(insSeq, (m.get(insSeq) ?? 0) + 1);
          }
          queryPos += cigarLen;
          break;
        }
        case 'D': {
          const overlapStart = refPos > start ? refPos : start;
          const overlapEnd = refPos + cigarLen < end ? refPos + cigarLen : end;
          for (let rp = overlapStart; rp < overlapEnd; rp++) {
            countsArr[rp - start].del++;
          }
          refPos += cigarLen;
          break;
        }
        case 'N': {
          refPos += cigarLen;
          break;
        }
        case 'S': {
          queryPos += cigarLen;
          break;
        }
        case 'H':
        case 'P':
          break;
      }
    });
  }

  const results: ConsensusPositionCall[] = new Array(len);

  for (let i = 0; i < len; i++) {
    const pos = start + i;
    const refBase = (refSeq[i] ?? 'N').toUpperCase();
    const c = countsArr[i];
    const depth = c.A + c.C + c.G + c.T + c.N + c.del;

    let call = 'N';
    let isLowDepth = false;
    let isHeterozygous = false;

    if (depth < options.minDepth || depth === 0) {
      call = 'N';
      isLowDepth = true;
    } else if (c.del / depth >= options.indelThreshold) {
      call = '-';
    } else {
      const nucEntries: Array<[string, number]> = [
        ['A', c.A],
        ['C', c.C],
        ['G', c.G],
        ['T', c.T],
      ];
      nucEntries.sort((a, b) => b[1] - a[1]);

      if (nucEntries[0][1] === 0) {
        call = 'N';
      } else if (options.ploidyMode === 'diploid') {
        const passing = nucEntries
          .filter(([, cnt]) => cnt / depth >= options.hetThreshold)
          .map(([base]) => base);
        if (passing.length > 1) {
          call = allelesToIupac(passing);
          isHeterozygous = true;
        } else {
          call = passing[0] ?? nucEntries[0][0];
        }
      } else {
        call = nucEntries[0][0];
      }
    }

    // 삽입(Insertion) 합의 판정
    let insertion = '';
    if (!isLowDepth && insertionMaps[i].size > 0) {
      let bestIns = '';
      let bestCnt = 0;
      for (const [insSeq, cnt] of insertionMaps[i].entries()) {
        if (cnt > bestCnt) {
          bestCnt = cnt;
          bestIns = insSeq;
        }
      }
      if (bestCnt / depth >= options.indelThreshold) {
        insertion = bestIns;
      }
    }

    const isMismatch = !isLowDepth && (call !== refBase || insertion.length > 0);

    results[i] = {
      pos,
      refBase,
      call,
      insertion,
      depth,
      isLowDepth,
      isMismatch,
      isHeterozygous,
      counts: c,
    };
  }

  return results;
}

export function classifyPositionDiff(
  refBase: string,
  s1Call: ConsensusPositionCall | undefined,
  s2Call: ConsensusPositionCall | undefined,
  hasStrain2 = true,
): DiffCategory {
  const s1Low = !s1Call || s1Call.isLowDepth;
  const ref = refBase.toUpperCase();

  if (!hasStrain2) {
    if (s1Low) return 'gap_s1';
    const s1Same = s1Call.call === ref && s1Call.insertion === '';
    return s1Same ? 'match' : 'strain1_only';
  }

  const s2Low = !s2Call || s2Call.isLowDepth;

  if (s1Low && s2Low) return 'gap_both';
  if (s1Low) return 'gap_s1';
  if (s2Low) return 'gap_s2';

  const b1 = s1Call.call;
  const b2 = s2Call.call;

  const s1SameAsRef = b1 === ref && s1Call.insertion === '';
  const s2SameAsRef = b2 === ref && s2Call.insertion === '';

  if (s1SameAsRef && s2SameAsRef) return 'match';
  if (!s1SameAsRef && s2SameAsRef) return 'strain1_only';
  if (s1SameAsRef && !s2SameAsRef) return 'strain2_only';
  if (b1 === b2 && s1Call.insertion === s2Call.insertion) return 'shared';
  return 'divergent';
}

export interface RegionMappingStats {
  totalRefBp: number;
  coveredBp: number;
  zeroDepthBp: number;
  lowDepthBp: number;
  delBp: number;
  mismatchBp: number;
  meanDepth: number;
  minDepth: number;
  coveragePct: number;
}

/**
 * 선택 구간 전체에 대해 각 균주의 매핑 완전성(커버리지·평균 깊이·미매핑 N·결실 등)을 집계합니다.
 */
export function computeRegionMappingStats(
  calls: ConsensusPositionCall[],
  refSeq = '',
): RegionMappingStats {
  const totalRefBp = calls.length;
  if (totalRefBp === 0) {
    return {
      totalRefBp: 0,
      coveredBp: 0,
      zeroDepthBp: 0,
      lowDepthBp: 0,
      delBp: 0,
      mismatchBp: 0,
      meanDepth: 0,
      minDepth: 0,
      coveragePct: 0,
    };
  }

  let coveredBp = 0;
  let zeroDepthBp = 0;
  let lowDepthBp = 0;
  let delBp = 0;
  let mismatchBp = 0;
  let depthSum = 0;
  let minDepth = Number.POSITIVE_INFINITY;

  for (let i = 0; i < totalRefBp; i++) {
    const c = calls[i];
    depthSum += c.depth;
    if (c.depth < minDepth) minDepth = c.depth;
    if (c.depth === 0) zeroDepthBp++;
    if (c.isLowDepth) {
      lowDepthBp++;
    } else {
      coveredBp++;
    }
    if (c.call === '-') delBp++;
    if (c.isMismatch || (refSeq[i] && !c.isLowDepth && c.call !== refSeq[i].toUpperCase())) {
      mismatchBp++;
    }
  }

  return {
    totalRefBp,
    coveredBp,
    zeroDepthBp,
    lowDepthBp,
    delBp,
    mismatchBp,
    meanDepth: Math.round((depthSum / totalRefBp) * 10) / 10,
    minDepth: Number.isFinite(minDepth) ? minDepth : 0,
    coveragePct: Math.round((coveredBp / totalRefBp) * 1000) / 10,
  };
}

/**
 * 단일 균주의 Call 배열을 지정된 Exon/구간에 대해 ungapped·aligned 서열로 직렬화합니다.
 */
export function serializeCallsForIntervals(
  calls: ConsensusPositionCall[],
  baseStart: number,
  intervals: Array<{ start: number; end: number }>,
): { ungappedSeq: string; alignedSeq: string; subsetCalls: ConsensusPositionCall[] } {
  const ungappedChars: string[] = [];
  const alignedChars: string[] = [];
  const subsetCalls: ConsensusPositionCall[] = [];

  for (const iv of intervals) {
    for (let pos = iv.start; pos < iv.end; pos++) {
      const c = calls[pos - baseStart];
      if (!c) continue;
      subsetCalls.push(c);
      alignedChars.push(c.call);
      if (c.call !== '-') ungappedChars.push(c.call);
      if (c.insertion) {
        alignedChars.push(c.insertion);
        ungappedChars.push(c.insertion);
      }
    }
  }

  return {
    ungappedSeq: ungappedChars.join(''),
    alignedSeq: alignedChars.join(''),
    subsetCalls,
  };
}

function assembleIntervalsAnalysis(
  chrom: string,
  start: number,
  end: number,
  refSeq: string,
  s1Calls: ConsensusPositionCall[],
  s2Calls: ConsensusPositionCall[] | null,
  intervals: Array<{ start: number; end: number }>,
  variants: VariantRecord[],
): ThreeWayRegionAnalysis {
  const hasStrain2 = s2Calls !== null;
  const alignedRefChars: string[] = [];
  const alignedS1Chars: string[] = [];
  const alignedS2Chars: string[] = [];
  const ungappedRefChars: string[] = [];
  const ungappedS1Chars: string[] = [];
  const ungappedS2Chars: string[] = [];
  const refSubChars: string[] = [];
  const s1Sub: ConsensusPositionCall[] = [];
  const s2Sub: ConsensusPositionCall[] = [];
  const columns: AlignedColumnMeta[] = [];

  for (const iv of intervals) {
    for (let pos = iv.start; pos < iv.end; pos++) {
      const i = pos - start;
      const c1 = s1Calls[i];
      if (!c1) continue;
      const c2 = hasStrain2 ? s2Calls[i] : undefined;
      const refBase = (refSeq[i] ?? 'N').toUpperCase();
      const cat = classifyPositionDiff(refBase, c1, c2, hasStrain2);

      refSubChars.push(refBase);
      s1Sub.push(c1);
      if (c2) s2Sub.push(c2);

      alignedRefChars.push(refBase);
      alignedS1Chars.push(c1.call);
      if (c2) alignedS2Chars.push(c2.call);
      columns.push({ refPos: pos, isInsertion: false, category: cat });

      ungappedRefChars.push(refBase);
      if (c1.call !== '-') ungappedS1Chars.push(c1.call);
      if (c2 && c2.call !== '-') ungappedS2Chars.push(c2.call);

      const maxInsLen = Math.max(c1.insertion.length, c2 ? c2.insertion.length : 0);
      if (maxInsLen > 0) {
        for (let k = 0; k < maxInsLen; k++) {
          const i1 = c1.insertion[k] ?? '-';
          const i2 = c2 ? c2.insertion[k] ?? '-' : '-';
          let insCat: DiffCategory = 'match';
          if (!hasStrain2) {
            insCat = i1 !== '-' ? 'strain1_only' : 'match';
          } else {
            if (i1 !== '-' && i2 === '-') insCat = 'strain1_only';
            else if (i1 === '-' && i2 !== '-') insCat = 'strain2_only';
            else if (i1 === i2) insCat = 'shared';
            else insCat = 'divergent';
          }

          alignedRefChars.push('-');
          alignedS1Chars.push(i1);
          if (hasStrain2) alignedS2Chars.push(i2);
          columns.push({ refPos: pos, isInsertion: true, category: insCat });

          if (i1 !== '-') ungappedS1Chars.push(i1);
          if (hasStrain2 && i2 !== '-') ungappedS2Chars.push(i2);
        }
      }
    }
  }

  return {
    chrom,
    start,
    end,
    hasStrain2,
    refSeq: refSubChars.join(''),
    s1Calls: s1Sub,
    s2Calls: hasStrain2 ? s2Sub : [],
    ungappedRef: ungappedRefChars.join(''),
    ungappedS1: ungappedS1Chars.join(''),
    ungappedS2: hasStrain2 ? ungappedS2Chars.join('') : '',
    alignedRef: alignedRefChars.join(''),
    alignedS1: alignedS1Chars.join(''),
    alignedS2: hasStrain2 ? alignedS2Chars.join('') : '',
    columns,
    variants,
  };
}

/**
 * 미리 계산된 ConsensusPositionCall 배열로부터 ThreeWayRegionAnalysis(정렬·언갭·변이)를 조립합니다.
 */
export function buildAnalysisFromCalls(
  chrom: string,
  start: number,
  end: number,
  refSeq: string,
  s1Calls: ConsensusPositionCall[],
  s2Calls: ConsensusPositionCall[] | null,
  compactForWindow = false,
): ThreeWayRegionAnalysis {
  const hasStrain2 = s2Calls !== null;
  const variants = extractVariantRecords(
    chrom,
    start,
    refSeq,
    s1Calls,
    hasStrain2 ? s2Calls : null,
  );

  if (compactForWindow) {
    return {
      chrom,
      start,
      end,
      hasStrain2,
      refSeq: '',
      s1Calls,
      s2Calls: s2Calls ?? [],
      ungappedRef: '',
      ungappedS1: '',
      ungappedS2: '',
      alignedRef: '',
      alignedS1: '',
      alignedS2: '',
      columns: [],
      variants,
    };
  }

  return assembleIntervalsAnalysis(
    chrom,
    start,
    end,
    refSeq,
    s1Calls,
    s2Calls,
    [{ start, end }],
    variants,
  );
}

/**
 * 스플라이싱된 Exon 구간 목록(0-based half-open)만을 이어 붙여 ThreeWayRegionAnalysis 형태로 추출합니다.
 */
export function sliceSplicedAnalysis(
  analysis: ThreeWayRegionAnalysis,
  exons: Array<{ start: number; end: number }>,
): ThreeWayRegionAnalysis {
  const sorted = exons
    .map((ex) => ({
      start: Math.max(analysis.start, ex.start),
      end: Math.min(analysis.end, ex.end),
    }))
    .filter((ex) => ex.end > ex.start)
    .sort((a, b) => a.start - b.start);

  if (sorted.length === 0) return analysis;

  const exonVariants = analysis.variants.filter((v) =>
    sorted.some((iv) => v.end > iv.start && v.start < iv.end),
  );

  const spliced = assembleIntervalsAnalysis(
    analysis.chrom,
    analysis.start,
    analysis.end,
    analysis.refSeq,
    analysis.s1Calls,
    analysis.hasStrain2 ? analysis.s2Calls : null,
    sorted,
    exonVariants,
  );
  spliced.start = sorted[0].start;
  spliced.end = sorted[sorted.length - 1].end;
  return spliced;
}

/**
 * Reference · 균주1 (· 선택적 균주2) 정렬 비교 및 추출용 서열 생성
 * - s2Reads가 null이면 Reference vs 균주1 단일 비교(2-Way) 모드로 동작합니다.
 * - alignedRef / alignedS1 (/ alignedS2): 다중서열정렬처럼 동일 길이('-' 패딩 포함)
 * - ungappedRef / ungappedS1 (/ ungappedS2): 각 균주의 실제 염기 서열('-' 제거)
 */
export function buildThreeWayAnalysis(
  chrom: string,
  start: number,
  end: number,
  refSeq: string,
  s1Reads: AlignedRead[],
  s2Reads: AlignedRead[] | null,
  options: ConsensusOptions,
  compactForWindow = false,
): ThreeWayRegionAnalysis {
  const hasStrain2 = s2Reads !== null;
  const s1Calls = computeConsensusForRegion(refSeq, s1Reads, start, end, options);
  const s2Calls = hasStrain2
    ? computeConsensusForRegion(refSeq, s2Reads, start, end, options)
    : null;

  return buildAnalysisFromCalls(
    chrom,
    start,
    end,
    refSeq,
    s1Calls,
    s2Calls,
    compactForWindow,
  );
}

/**
 * 연속된 동일 카테고리 변이/결실 이벤트를 하나의 VariantRecord로 병합
 */
export function extractVariantRecords(
  chrom: string,
  start: number,
  refSeq: string,
  s1Calls: ConsensusPositionCall[],
  s2Calls: ConsensusPositionCall[] | null,
): VariantRecord[] {
  const hasStrain2 = s2Calls !== null;
  const variants: VariantRecord[] = [];
  let i = 0;

  while (i < s1Calls.length) {
    const c1 = s1Calls[i];
    const c2 = hasStrain2 ? s2Calls[i] : undefined;
    const refBase = (refSeq[i] ?? 'N').toUpperCase();
    const cat = classifyPositionDiff(refBase, c1, c2, hasStrain2);

    if (cat === 'match') {
      i++;
      continue;
    }

    const isGap = cat === 'gap_s1' || cat === 'gap_s2' || cat === 'gap_both';
    const isDel = c1.call === '-' || (c2 ? c2.call === '-' : false);
    let j = i + 1;

    if (isGap || isDel) {
      while (j < s1Calls.length) {
        const nc1 = s1Calls[j];
        const nc2 = hasStrain2 ? s2Calls[j] : undefined;
        const nRef = (refSeq[j] ?? 'N').toUpperCase();
        const nCat = classifyPositionDiff(nRef, nc1, nc2, hasStrain2);
        if (
          nCat === cat &&
          (isGap || nc1.call === '-' || (nc2 ? nc2.call === '-' : false))
        ) {
          j++;
        } else {
          break;
        }
      }
    }

    const refAllele = refSeq.slice(i, j).toUpperCase();
    const s1Allele = s1Calls
      .slice(i, j)
      .map((c) => c.call + c.insertion)
      .join('');
    const s2Allele = hasStrain2
      ? s2Calls
          .slice(i, j)
          .map((c) => c.call + c.insertion)
          .join('')
      : '';

    let type: VariantRecord['type'] = 'SNP';
    if (isGap) type = 'GAP';
    else if (s1Allele.includes('-') || s2Allele.includes('-')) type = 'DEL';
    else if (s1Calls[i].insertion || (c2 && c2.insertion)) type = 'INS';
    else if (j - i > 1) type = 'COMPLEX';

    const avgS1Depth = Math.round(
      s1Calls.slice(i, j).reduce((acc, c) => acc + c.depth, 0) / (j - i),
    );
    const avgS2Depth = hasStrain2
      ? Math.round(s2Calls.slice(i, j).reduce((acc, c) => acc + c.depth, 0) / (j - i))
      : 0;

    variants.push({
      chrom,
      start: start + i,
      end: start + j,
      type,
      category: cat,
      refAllele,
      s1Allele,
      s2Allele,
      s1Depth: avgS1Depth,
      s2Depth: avgS2Depth,
    });

    i = j;
  }

  return variants;
}

