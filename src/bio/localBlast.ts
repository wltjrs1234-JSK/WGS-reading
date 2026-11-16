import { reverseComplement } from './geneticCode';
import type { AlignedRead, ConsensusPositionCall } from './types';

export interface ParsedBlastQuery {
  id: string;
  seq: string;
}

export interface BlastContigSeq {
  name: string;
  seq: string; // uppercase genomic sequence
}

export interface BlastReferenceIndex {
  contigs: BlastContigSeq[];
  contigNames: string[];
  contigStarts: number[]; // concatenated global start offset per contig
  totalLength: number;
  kmerSize: number;
  counts: Int32Array;
  offsets: Int32Array;
  positions: Int32Array;
}

export interface BlastSearchOptions {
  maxHitsPerQuery?: number;
  minIdentityPct?: number;
  maxEValue?: number;
  shortPrimerMaxMismatches?: number; // applies automatically when queryLen <= 55
}

export interface BlastRawHit {
  queryId: string;
  queryLen: number;
  chrom: string;
  refStart: number; // 0-based inclusive on + genomic ref
  refEnd: number; // 0-based exclusive on + genomic ref
  qStart: number; // 0-based inclusive on user's original query
  qEnd: number; // 0-based exclusive on user's original query
  strand: '+' | '-';
  score: number;
  bitScore: number;
  eValue: number;
  identities: number;
  alignLen: number;
  mismatches: number;
  gaps: number;
  identityPct: number;
  queryCoveragePct: number;
  /** Chromosomal '+' order alignment strings (with '-' gaps), ready to compare against consensus calls and reverseComplement if strand === '-' */
  alignedQueryPlus: string;
  alignedRefPlus: string;
}

export interface BlastStrainComparison {
  id: string;
  label: string;
  matchedChrom: string | null;
  identityVsQueryPct: number;
  identitiesVsQuery: number;
  comparedCols: number;
  coveragePct: number;
  meanDepth: number;
  lowDepthBp: number;
  delBp: number;
  snpCountVsRef: number;
  softClipReads: number;
  insertionReads: number;
  /** Strain consensus aligned column-for-column with alignedQueryDisplay/alignedRefDisplay */
  alignedStrainDisplay: string;
}

export interface BlastForeignReadHit {
  strainId: string;
  strainLabel: string;
  queryId: string;
  chrom: string;
  anchorPos: number; // 0-based ref coordinate near clip/insertion/alignment
  clipSide: 'left_clip' | 'right_clip' | 'insertion' | 'unmapped_region';
  supportingReads: number;
  identityPct: number;
  matchLen: number;
  strand: '+' | '-';
  sampleReadName: string;
  sampleCigar: string;
  alignedQueryDisplay: string;
  alignedReadDisplay: string;
}

const KMER_K = 11;
const KMER_TABLE_SIZE = 1 << (KMER_K * 2); // 4^11 = 4,194,304 (~16MB Int32Array)
const KMER_MASK = KMER_TABLE_SIZE - 1;
const MAX_SEED_FREQ = 1200;

function baseTo2Bit(code: number): number {
  switch (code) {
    case 65: // A
    case 97: // a
      return 0;
    case 67: // C
    case 99: // c
      return 1;
    case 71: // G
    case 103: // g
      return 2;
    case 84: // T
    case 116: // t
    case 85: // U
    case 117: // u
      return 3;
    default:
      return -1;
  }
}

export function parseBlastQueries(rawInput: string, maxTotalLen = 60_000): ParsedBlastQuery[] {
  const lines = rawInput.split(/\r?\n/);
  const queries: ParsedBlastQuery[] = [];
  let curId: string | null = null;
  let chunks: string[] = [];

  const flush = () => {
    const seq = chunks
      .join('')
      .replace(/[\s0-9]/g, '')
      .toUpperCase()
      .replace(/U/g, 'T')
      .replace(/[^ACGTNRYSWKMBDHV]/g, '');
    if (seq.length > 0) {
      queries.push({
        id: curId ?? `Query_${queries.length + 1}`,
        seq,
      });
    }
    chunks = [];
  };

  for (const rawLine of lines) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      flush();
      curId = line.slice(1).trim().split(/\s+/)[0] || `Query_${queries.length + 1}`;
    } else {
      chunks.push(line);
    }
  }
  flush();

  if (queries.length === 0) {
    throw new Error('유효한 DNA 서열을 입력해주세요 (최소 10 bp 이상 권장).');
  }

  let totalLen = 0;
  for (const q of queries) {
    if (q.seq.length < 10) {
      throw new Error(`'${q.id}' 서열 길이가 너무 짧습니다 (${q.seq.length} bp). 최소 10 bp 이상 입력해주세요.`);
    }
    totalLen += q.seq.length;
  }
  if (totalLen > maxTotalLen) {
    throw new Error(
      `로컬 BLAST 1회 총 쿼리 길이는 ${maxTotalLen.toLocaleString()} bp 이하로 제한됩니다 (현재 ${totalLen.toLocaleString()} bp).`,
    );
  }
  return queries;
}

/**
 * Reference 전체 염색체에 대해 zero-GC CSR 11-mer 인덱스를 구축합니다.
 * (~12.16 Mb 효모 지놈 기준 ~64 MB 연속 TypedArray, 1회 구축 후 재사용)
 */
export function buildBlastReferenceIndex(contigs: BlastContigSeq[]): BlastReferenceIndex {
  const contigNames = contigs.map((c) => c.name);
  const contigStarts = new Array<number>(contigs.length);
  let totalLength = 0;
  for (let i = 0; i < contigs.length; i++) {
    contigStarts[i] = totalLength;
    totalLength += contigs[i].seq.length;
  }

  const counts = new Int32Array(KMER_TABLE_SIZE);

  // Pass 1: 각 11-mer 출현 빈도 집계
  for (let cIdx = 0; cIdx < contigs.length; cIdx++) {
    const seq = contigs[cIdx].seq;
    let code = 0;
    let validRun = 0;
    for (let i = 0; i < seq.length; i++) {
      const bit = baseTo2Bit(seq.charCodeAt(i));
      if (bit < 0) {
        code = 0;
        validRun = 0;
        continue;
      }
      code = ((code << 2) | bit) & KMER_MASK;
      validRun++;
      if (validRun >= KMER_K) {
        counts[code]++;
      }
    }
  }

  const offsets = new Int32Array(KMER_TABLE_SIZE + 1);
  let running = 0;
  for (let k = 0; k < KMER_TABLE_SIZE; k++) {
    offsets[k] = running;
    running += counts[k];
  }
  offsets[KMER_TABLE_SIZE] = running;

  const positions = new Int32Array(running);
  const writeCursor = new Int32Array(offsets.subarray(0, KMER_TABLE_SIZE));

  // Pass 2: 전역 시작 좌표 채우기
  for (let cIdx = 0; cIdx < contigs.length; cIdx++) {
    const seq = contigs[cIdx].seq;
    const baseGlobal = contigStarts[cIdx];
    let code = 0;
    let validRun = 0;
    for (let i = 0; i < seq.length; i++) {
      const bit = baseTo2Bit(seq.charCodeAt(i));
      if (bit < 0) {
        code = 0;
        validRun = 0;
        continue;
      }
      code = ((code << 2) | bit) & KMER_MASK;
      validRun++;
      if (validRun >= KMER_K) {
        const startPos = i - KMER_K + 1;
        const slot = writeCursor[code]++;
        positions[slot] = baseGlobal + startPos;
      }
    }
  }

  return {
    contigs,
    contigNames,
    contigStarts,
    totalLength,
    kmerSize: KMER_K,
    counts,
    offsets,
    positions,
  };
}

function locateGlobalContig(
  index: BlastReferenceIndex,
  globalPos: number,
): { contigIdx: number; localPos: number } {
  const starts = index.contigStarts;
  let lo = 0;
  let hi = starts.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (starts[mid] <= globalPos) lo = mid + 1;
    else hi = mid - 1;
  }
  const idx = Math.max(0, hi);
  return { contigIdx: idx, localPos: globalPos - starts[idx] };
}

interface SwAlignResult {
  score: number;
  qStart: number;
  qEnd: number;
  tStart: number;
  tEnd: number;
  alignedQuery: string;
  alignedTarget: string;
  identities: number;
  mismatches: number;
  gaps: number;
}

const MATCH_SCORE = 2;
const MISMATCH_PENALTY = -3;
const GAP_OPEN = -5;
const GAP_EXTEND = -2;

/**
 * 국소 Smith-Waterman (Linear-gap SW 또는 Banded SW) 정렬
 * 긴 쿼리(> 420 bp)에서는 주요 대각선(anchorDiag) 주변 Band 내에서만 계산해 수 ms 이내 완료
 */
function smithWatermanAlign(
  query: string,
  target: string,
  anchorDiag?: number,
): SwAlignResult | null {
  const m = query.length;
  const n = target.length;
  if (m === 0 || n === 0) return null;

  const useBanded = anchorDiag !== undefined && m > 420 && n > 420;
  const bandHalf = useBanded ? Math.min(180, Math.max(48, Math.floor(m * 0.1))) : Math.max(m, n) + 2;
  const diagCenter = anchorDiag ?? 0;

  // Affine-lite / linear gap-score SW with trace pointers: 1=diag, 2=up(gap in target), 3=left(gap in query)
  // Using affine approx: first gap costs (GAP_OPEN+GAP_EXTEND) by checking trace direction
  const stride = useBanded ? 2 * bandHalf + 1 : n + 1;
  const rows = m + 1;

  // Compact row storage when full matrix fits comfortably (<= 3.5M cells), else banded!
  const H = new Int32Array(rows * stride);
  const T = new Uint8Array(rows * stride); // 0=stop, 1=diag, 2=delInTarget(up), 3=insInTarget(left)

  let bestScore = 0;
  let bestI = 0;
  let bestCol = 0;

  if (!useBanded) {
    for (let i = 1; i <= m; i++) {
      const qCh = query.charCodeAt(i - 1);
      const rowOff = i * stride;
      const prevOff = (i - 1) * stride;
      for (let j = 1; j <= n; j++) {
        const tCh = target.charCodeAt(j - 1);
        const sub = qCh === tCh && qCh !== 78 ? MATCH_SCORE : MISMATCH_PENALTY;
        const diagVal = H[prevOff + j - 1] + sub;

        const upTrace = T[prevOff + j];
        const upVal = H[prevOff + j] + (upTrace === 2 ? GAP_EXTEND : GAP_OPEN);

        const leftTrace = T[rowOff + j - 1];
        const leftVal = H[rowOff + j - 1] + (leftTrace === 3 ? GAP_EXTEND : GAP_OPEN);

        let val = 0;
        let tr = 0;
        if (diagVal > val) {
          val = diagVal;
          tr = 1;
        }
        if (upVal > val) {
          val = upVal;
          tr = 2;
        }
        if (leftVal > val) {
          val = leftVal;
          tr = 3;
        }
        H[rowOff + j] = val;
        T[rowOff + j] = tr;
        if (val > bestScore) {
          bestScore = val;
          bestI = i;
          bestCol = j;
        }
      }
    }
  } else {
    for (let i = 1; i <= m; i++) {
      const qCh = query.charCodeAt(i - 1);
      const jMin = Math.max(1, i + diagCenter - bandHalf);
      const jMax = Math.min(n, i + diagCenter + bandHalf);
      const rowOff = i * stride;
      const prevOff = (i - 1) * stride;
      const bandBase = i + diagCenter - bandHalf;
      const prevBandBase = i - 1 + diagCenter - bandHalf;

      for (let j = jMin; j <= jMax; j++) {
        const col = j - bandBase; // [0 .. 2*bandHalf]

        // diag (i-1, j-1) -> prevCol = (j-1) - prevBandBase == col
        const diagScore = H[prevOff + col];
        const tCh = target.charCodeAt(j - 1);
        const sub = qCh === tCh && qCh !== 78 ? MATCH_SCORE : MISMATCH_PENALTY;
        const diagVal = diagScore + sub;

        // up (i-1, j) -> prevCol = j - prevBandBase == col + 1
        let upVal = 0;
        const upPrevCol = j - prevBandBase;
        if (upPrevCol >= 0 && upPrevCol < stride) {
          const upTr = T[prevOff + upPrevCol];
          upVal = H[prevOff + upPrevCol] + (upTr === 2 ? GAP_EXTEND : GAP_OPEN);
        }

        // left (i, j-1) -> curCol = col - 1
        let leftVal = 0;
        const leftCol = col - 1;
        if (leftCol >= 0) {
          const leftTr = T[rowOff + leftCol];
          leftVal = H[rowOff + leftCol] + (leftTr === 3 ? GAP_EXTEND : GAP_OPEN);
        }

        let val = 0;
        let tr = 0;
        if (diagVal > val) {
          val = diagVal;
          tr = 1;
        }
        if (upVal > val) {
          val = upVal;
          tr = 2;
        }
        if (leftVal > val) {
          val = leftVal;
          tr = 3;
        }
        H[rowOff + col] = val;
        T[rowOff + col] = tr;
        if (val > bestScore) {
          bestScore = val;
          bestI = i;
          bestCol = j; // store real j in bestCol
        }
      }
    }
  }

  if (bestScore <= 0 || bestI === 0 || bestCol === 0) return null;

  // Backtrack
  const qChars: string[] = [];
  const tChars: string[] = [];
  let ci = bestI;
  let cj = bestCol;
  let identities = 0;
  let mismatches = 0;
  let gaps = 0;

  while (ci > 0 && cj > 0) {
    const colIdx = useBanded ? cj - (ci + diagCenter - bandHalf) : cj;
    if (colIdx < 0 || colIdx >= stride) break;
    const tr = T[ci * stride + colIdx];
    if (tr === 0) break;

    if (tr === 1) {
      const qc = query[ci - 1];
      const tc = target[cj - 1];
      qChars.push(qc);
      tChars.push(tc);
      if (qc === tc && qc !== 'N') identities++;
      else mismatches++;
      ci--;
      cj--;
    } else if (tr === 2) {
      qChars.push(query[ci - 1]);
      tChars.push('-');
      gaps++;
      ci--;
    } else {
      qChars.push('-');
      tChars.push(target[cj - 1]);
      gaps++;
      cj--;
    }
  }

  qChars.reverse();
  tChars.reverse();

  return {
    score: bestScore,
    qStart: ci,
    qEnd: bestI,
    tStart: cj,
    tEnd: bestCol,
    alignedQuery: qChars.join(''),
    alignedTarget: tChars.join(''),
    identities,
    mismatches,
    gaps,
  };
}

function computeBlastStats(
  rawScore: number,
  queryLen: number,
  genomeLen: number,
): { bitScore: number; eValue: number } {
  // Standardized BLASTN lambda / K approximation
  const lambda = 0.625;
  const K = 0.41;
  const bitScore = Math.max(0, (lambda * rawScore - Math.log(K)) / Math.LN2);
  const effSearchSpace = Math.max(100, queryLen) * Math.max(1000, genomeLen * 2);
  const rawE = K * effSearchSpace * Math.exp(-lambda * rawScore);
  const eValue = Number.isFinite(rawE) ? Number(rawE.toExponential(2)) : 0;
  return {
    bitScore: Math.round(bitScore * 10) / 10,
    eValue,
  };
}

function collectSeedCandidatesForStrand(
  index: BlastReferenceIndex,
  orientedQuery: string,
): Array<{ contigIdx: number; refCenterStart: number; diag: number; weight: number }> {
  const qLen = orientedQuery.length;
  if (qLen < index.kmerSize) return [];

  const step = Math.max(1, Math.min(8, Math.floor(qLen / 80)));
  // Map from "contigIdx:diagBucket" -> { contigIdx, sumDiag, count, minRef, maxRef }
  const diagBuckets = new Map<
    string,
    { contigIdx: number; diagSum: number; count: number; minRef: number; maxRef: number }
  >();

  let highFreqSkipped = 0;

  for (let qi = 0; qi <= qLen - index.kmerSize; qi += step) {
    let code = 0;
    let valid = true;
    for (let k = 0; k < index.kmerSize; k++) {
      const bit = baseTo2Bit(orientedQuery.charCodeAt(qi + k));
      if (bit < 0) {
        valid = false;
        break;
      }
      code = (code << 2) | bit;
    }
    if (!valid) continue;

    const startOff = index.offsets[code];
    const endOff = index.offsets[code + 1];
    const freq = endOff - startOff;
    if (freq === 0) continue;
    if (freq > MAX_SEED_FREQ) {
      highFreqSkipped++;
      continue;
    }

    for (let p = startOff; p < endOff; p++) {
      const gPos = index.positions[p];
      const { contigIdx, localPos } = locateGlobalContig(index, gPos);
      const diag = localPos - qi;
      const bucketKey = `${contigIdx}:${Math.floor(diag / 24)}`;
      let entry = diagBuckets.get(bucketKey);
      if (!entry) {
        entry = {
          contigIdx,
          diagSum: diag,
          count: 1,
          minRef: localPos,
          maxRef: localPos + index.kmerSize,
        };
        diagBuckets.set(bucketKey, entry);
      } else {
        entry.diagSum += diag;
        entry.count++;
        if (localPos < entry.minRef) entry.minRef = localPos;
        if (localPos + index.kmerSize > entry.maxRef) entry.maxRef = localPos + index.kmerSize;
      }
    }
  }

  // Fallback for ultra-repetitive low-complexity queries if all seeds skipped
  if (diagBuckets.size === 0 && highFreqSkipped > 0) {
    for (let qi = 0; qi <= qLen - index.kmerSize; qi += Math.max(step, 5)) {
      let code = 0;
      let valid = true;
      for (let k = 0; k < index.kmerSize; k++) {
        const bit = baseTo2Bit(orientedQuery.charCodeAt(qi + k));
        if (bit < 0) {
          valid = false;
          break;
        }
        code = (code << 2) | bit;
      }
      if (!valid) continue;
      const startOff = index.offsets[code];
      const endOff = Math.min(index.offsets[code + 1], startOff + 60);
      for (let p = startOff; p < endOff; p++) {
        const gPos = index.positions[p];
        const { contigIdx, localPos } = locateGlobalContig(index, gPos);
        const diag = localPos - qi;
        const bucketKey = `${contigIdx}:${Math.floor(diag / 24)}`;
        const existing = diagBuckets.get(bucketKey);
        if (!existing) {
          diagBuckets.set(bucketKey, {
            contigIdx,
            diagSum: diag,
            count: 1,
            minRef: localPos,
            maxRef: localPos + index.kmerSize,
          });
        } else {
          existing.count++;
        }
      }
    }
  }

  // Merge adjacent buckets in same locus to avoid duplicate alignments of one region
  const sorted = [...diagBuckets.values()].sort((a, b) => b.count - a.count).slice(0, 60);
  const merged: Array<{ contigIdx: number; refCenterStart: number; diag: number; weight: number }> =
    [];

  for (const b of sorted) {
    const avgDiag = Math.round(b.diagSum / b.count);
    const estStart = avgDiag;
    const dup = merged.find(
      (m) => m.contigIdx === b.contigIdx && Math.abs(m.refCenterStart - estStart) < Math.max(40, qLen * 0.6),
    );
    if (!dup) {
      merged.push({
        contigIdx: b.contigIdx,
        refCenterStart: estStart,
        diag: avgDiag,
        weight: b.count,
      });
    }
  }
  return merged;
}

/**
 * 짧은 프라이머/sgRNA(12~55 bp)용 Pigeonhole 멀티-시드 불일치 허용 스캔
 * 최대 k개 mismatch에서도 반드시 1개 이상의 분할 시드가 정확히 일치하는 원리를 사용합니다.
 */
function findShortPrimerCandidates(
  contigs: BlastContigSeq[],
  orientedQuery: string,
  maxMismatches: number,
): Array<{ contigIdx: number; windowStart: number; windowEnd: number; diag: number }> {
  const qLen = orientedQuery.length;
  const effectiveM = Math.min(maxMismatches, Math.max(0, Math.floor((qLen - 10) / 5)));
  const numSeeds = effectiveM + 1;
  const seedLen = Math.max(5, Math.floor(qLen / numSeeds));

  const seenLocs = new Set<string>();
  const candidates: Array<{ contigIdx: number; windowStart: number; windowEnd: number; diag: number }> =
    [];

  for (let sIdx = 0; sIdx < numSeeds; sIdx++) {
    const qOff = Math.min(qLen - seedLen, sIdx * seedLen);
    const seed = orientedQuery.slice(qOff, qOff + seedLen);
    if (seed.includes('N')) continue;

    for (let cIdx = 0; cIdx < contigs.length; cIdx++) {
      const seq = contigs[cIdx].seq;
      let searchFrom = 0;
      let hitSafety = 0;
      while (searchFrom <= seq.length - seedLen && hitSafety < 2500) {
        const hitPos = seq.indexOf(seed, searchFrom);
        if (hitPos === -1) break;
        searchFrom = hitPos + 1;
        hitSafety++;

        const expectedStart = hitPos - qOff;
        if (expectedStart < -4 || expectedStart + qLen > seq.length + 4) continue;

        // Fast mismatch pre-filter allowing small edge overhang
        const winStart = Math.max(0, expectedStart - 6);
        const winEnd = Math.min(seq.length, expectedStart + qLen + 6);
        let approxMis = 0;
        const overlapS = Math.max(0, expectedStart);
        const overlapE = Math.min(seq.length, expectedStart + qLen);
        for (let rp = overlapS; rp < overlapE; rp++) {
          if (seq.charCodeAt(rp) !== orientedQuery.charCodeAt(rp - expectedStart)) {
            approxMis++;
            if (approxMis > effectiveM + 2) break;
          }
        }
        // Only full SW if close mismatch count OR near gap
        if (approxMis <= effectiveM + 3 || qLen >= 18) {
          if (approxMis <= effectiveM) {
            const key = `${cIdx}:${Math.floor(expectedStart / 8)}`;
            if (!seenLocs.has(key)) {
              seenLocs.add(key);
              candidates.push({
                contigIdx: cIdx,
                windowStart: winStart,
                windowEnd: winEnd,
                diag: expectedStart - winStart,
              });
            }
          }
        }
      }
    }
  }

  return candidates.slice(0, 80);
}

/**
 * 단일 Query를 Reference 전체(+/- 양가닥)에 대해 정렬 검색합니다.
 */
export function searchBlastReference(
  index: BlastReferenceIndex,
  query: ParsedBlastQuery,
  options: BlastSearchOptions = {},
): BlastRawHit[] {
  const maxHits = options.maxHitsPerQuery ?? 30;
  const minIdentityPct = options.minIdentityPct ?? 75;
  const maxEValue = options.maxEValue ?? 10;
  const maxPrimerMis = options.shortPrimerMaxMismatches ?? 2;

  const strands: Array<{ strand: '+' | '-'; orientedQ: string }> = [
    { strand: '+', orientedQ: query.seq },
    { strand: '-', orientedQ: reverseComplement(query.seq) },
  ];

  const allHits: BlastRawHit[] = [];
  const dedupSeen = new Set<string>();

  for (const { strand, orientedQ } of strands) {
    const qLen = orientedQ.length;
    const windows: Array<{
      contigIdx: number;
      winStart: number;
      winEnd: number;
      diagInWin: number;
    }> = [];

    if (qLen >= index.kmerSize) {
      const seeds = collectSeedCandidatesForStrand(index, orientedQ);
      const pad = Math.min(180, Math.max(30, Math.round(qLen * 0.25)));
      for (const sd of seeds) {
        const cLen = index.contigs[sd.contigIdx].seq.length;
        const winStart = Math.max(0, sd.refCenterStart - pad);
        const winEnd = Math.min(cLen, sd.refCenterStart + qLen + pad);
        if (winEnd > winStart) {
          windows.push({
            contigIdx: sd.contigIdx,
            winStart,
            winEnd,
            diagInWin: sd.diag - winStart,
          });
        }
      }
    }

    // 짧은 프라이머/sgRNA (10~55 bp)에서 11-mer 시드 히트가 없거나 적으면 Pigeonhole 불일치 스캔 보충
    if (qLen <= 55 && maxPrimerMis > 0) {
      const primerCands = findShortPrimerCandidates(index.contigs, orientedQ, maxPrimerMis);
      for (const pc of primerCands) {
        windows.push({
          contigIdx: pc.contigIdx,
          winStart: pc.windowStart,
          winEnd: pc.windowEnd,
          diagInWin: pc.diag,
        });
      }
    }

    for (const w of windows) {
      const contig = index.contigs[w.contigIdx];
      const targetSlice = contig.seq.slice(w.winStart, w.winEnd);
      const sw = smithWatermanAlign(orientedQ, targetSlice, w.diagInWin);
      if (!sw) continue;

      const alignLen = sw.alignedQuery.length;
      if (alignLen < Math.min(qLen, 10)) continue;

      const refStart = w.winStart + sw.tStart;
      const refEnd = w.winStart + sw.tEnd;
      const identityPct = Math.round((sw.identities / Math.max(1, alignLen)) * 1000) / 10;
      const alignedQuerySpan = sw.qEnd - sw.qStart;
      const queryCoveragePct = Math.round((alignedQuerySpan / qLen) * 1000) / 10;

      if (identityPct < minIdentityPct) continue;
      // Ignore tiny spurious hits on long queries (< 20 bp and < 30% query coverage)
      if (qLen > 30 && alignedQuerySpan < 18 && queryCoveragePct < 30) continue;

      const { bitScore, eValue } = computeBlastStats(sw.score, qLen, index.totalLength);
      if (eValue > maxEValue && alignedQuerySpan < qLen) continue;

      const dedupKey = `${strand}:${contig.name}:${refStart}:${refEnd}`;
      if (dedupSeen.has(dedupKey)) continue;
      dedupSeen.add(dedupKey);
      const qStartOrig = strand === '+' ? sw.qStart : qLen - sw.qEnd;
      const qEndOrig = strand === '+' ? sw.qEnd : qLen - sw.qStart;

      allHits.push({
        queryId: query.id,
        queryLen: qLen,
        chrom: contig.name,
        refStart,
        refEnd,
        qStart: qStartOrig,
        qEnd: qEndOrig,
        strand,
        score: sw.score,
        bitScore,
        eValue,
        identities: sw.identities,
        alignLen,
        mismatches: sw.mismatches,
        gaps: sw.gaps,
        identityPct,
        queryCoveragePct,
        alignedQueryPlus: sw.alignedQuery,
        alignedRefPlus: sw.alignedTarget,
      });
    }
  }

  allHits.sort((a, b) => b.score - a.score || a.eValue - b.eValue || b.identityPct - a.identityPct);
  return allHits.slice(0, maxHits);
}

/**
 * Reference Hit 구간 [hit.refStart, hit.refEnd)에 대해 특정 BAM 균주의 실제 Consensus·Depth·Soft-clip을 대조합니다.
 */
export function evaluateHitAgainstStrainCalls(
  hit: BlastRawHit,
  strainId: string,
  strainLabel: string,
  matchedChrom: string | null,
  calls: ConsensusPositionCall[],
  reads: AlignedRead[],
): BlastStrainComparison {
  if (!matchedChrom || calls.length === 0) {
    const emptyAligned = '-'.repeat(hit.alignedQueryPlus.length);
    return {
      id: strainId,
      label: strainLabel,
      matchedChrom,
      identityVsQueryPct: 0,
      identitiesVsQuery: 0,
      comparedCols: hit.alignedQueryPlus.length,
      coveragePct: 0,
      meanDepth: 0,
      lowDepthBp: hit.refEnd - hit.refStart,
      delBp: 0,
      snpCountVsRef: 0,
      softClipReads: 0,
      insertionReads: 0,
      alignedStrainDisplay: emptyAligned,
    };
  }

  let coveredBp = 0;
  let lowDepthBp = 0;
  let delBp = 0;
  let snpCountVsRef = 0;
  let depthSum = 0;

  for (const c of calls) {
    depthSum += c.depth;
    if (c.isLowDepth) lowDepthBp++;
    else coveredBp++;
    if (c.call === '-') delBp++;
    if (c.isMismatch) snpCountVsRef++;
  }

  let softClipReads = 0;
  let insertionReads = 0;
  for (const r of reads) {
    if ((r.flags & 0x704) !== 0) continue;
    if (r.cigar.includes('S')) softClipReads++;
    if (r.cigar.includes('I')) insertionReads++;
  }

  // alignedRefPlus 각 컬럼과 1:1 매칭되는 균주 Consensus 서열 생성
  const strainCharsPlus: string[] = [];
  let refCursor = 0;
  for (let i = 0; i < hit.alignedRefPlus.length; i++) {
    const refCh = hit.alignedRefPlus[i];
    if (refCh === '-') {
      // Query 삽입 위치: 직전 레퍼런스 위치에 균주 Insertion이 있다면 대응 염기 비교
      const prevIdx = Math.max(0, refCursor - 1);
      const ins = calls[prevIdx]?.insertion ?? '';
      // 몇 번째 연속 '-'인지 추적
      let insOffset = 0;
      for (let j = i - 1; j >= 0 && hit.alignedRefPlus[j] === '-'; j--) {
        insOffset++;
      }
      strainCharsPlus.push(ins[insOffset] ?? '-');
    } else {
      const c = calls[refCursor];
      strainCharsPlus.push(c ? c.call : 'N');
      refCursor++;
    }
  }

  const strainPlusStr = strainCharsPlus.join('');
  let identitiesVsQuery = 0;
  const totalCols = hit.alignedQueryPlus.length;
  for (let i = 0; i < totalCols; i++) {
    const qc = hit.alignedQueryPlus[i];
    const sc = strainPlusStr[i];
    if (qc !== '-' && sc !== '-' && sc !== 'N' && qc === sc) {
      identitiesVsQuery++;
    }
  }

  const alignedStrainDisplay =
    hit.strand === '-' ? reverseComplement(strainPlusStr) : strainPlusStr;

  const totalCalls = calls.length;
  return {
    id: strainId,
    label: strainLabel,
    matchedChrom,
    identityVsQueryPct: Math.round((identitiesVsQuery / Math.max(1, totalCols)) * 1000) / 10,
    identitiesVsQuery,
    comparedCols: totalCols,
    coveragePct: Math.round((coveredBp / Math.max(1, totalCalls)) * 1000) / 10,
    meanDepth: Math.round((depthSum / Math.max(1, totalCalls)) * 10) / 10,
    lowDepthBp,
    delBp,
    snpCountVsRef,
    softClipReads,
    insertionReads,
    alignedStrainDisplay,
  };
}

/**
 * Reference에 없는 외래 서열(KanMX, GFP, 벡터 카세트 등)을 BAM 리드의 Soft-clip / Insertion 구간에서 직접 탐색합니다.
 */
export function findForeignMatchesInReads(
  query: ParsedBlastQuery,
  strainId: string,
  strainLabel: string,
  chrom: string,
  reads: AlignedRead[],
  minMatchLen = 15,
): BlastForeignReadHit[] {
  const qLen = query.seq.length;
  const seedLen = Math.min(11, Math.max(8, Math.floor(qLen * 0.7)));
  if (qLen < seedLen) return [];

  const qForward = query.seq;
  const qRev = reverseComplement(query.seq);
  const forwardSeeds: string[] = [];
  const revSeeds: string[] = [];
  const step = Math.max(1, Math.min(6, Math.floor(seedLen / 2)));

  for (let i = 0; i <= qForward.length - seedLen; i += step) {
    forwardSeeds.push(qForward.slice(i, i + seedLen));
    revSeeds.push(qRev.slice(i, i + seedLen));
  }

  interface RawReadMatch {
    anchorPos: number;
    clipSide: BlastForeignReadHit['clipSide'];
    strand: '+' | '-';
    identityPct: number;
    matchLen: number;
    readName: string;
    cigar: string;
    alignedQuery: string;
    alignedRead: string;
  }

  const rawMatches: RawReadMatch[] = [];

  for (const r of reads) {
    if ((r.flags & 0x704) !== 0) continue;
    const hasInterestingCigar = r.cigar.includes('S') || r.cigar.includes('I') || r.cigar === '*';
    if (!hasInterestingCigar) continue;
    const seqUpper = r.seq?.toUpperCase() ?? '';
    if (seqUpper.length < seedLen) continue;

    let matchedStrand: '+' | '-' | null = null;
    for (let i = 0; i < forwardSeeds.length; i++) {
      if (seqUpper.includes(forwardSeeds[i])) {
        matchedStrand = '+';
        break;
      }
      if (seqUpper.includes(revSeeds[i])) {
        matchedStrand = '-';
        break;
      }
    }
    if (!matchedStrand) continue;

    const orientedQ = matchedStrand === '+' ? qForward : qRev;
    const sw = smithWatermanAlign(orientedQ, seqUpper);
    if (!sw) continue;

    const matchLen = sw.alignedQuery.length;
    const idPct = Math.round((sw.identities / Math.max(1, matchLen)) * 1000) / 10;
    if (matchLen < Math.min(minMatchLen, qLen) || idPct < 80) continue;

    let clipSide: BlastForeignReadHit['clipSide'] = 'insertion';
    let anchorPos = r.start;
    if (/^\d+S/.test(r.cigar)) {
      clipSide = 'left_clip';
      anchorPos = r.start;
    } else if (/\d+S$/.test(r.cigar)) {
      clipSide = 'right_clip';
      anchorPos = Math.max(r.start, r.end - 1);
    } else if (r.cigar === '*') {
      clipSide = 'unmapped_region';
    }

    rawMatches.push({
      anchorPos,
      clipSide,
      strand: matchedStrand,
      identityPct: idPct,
      matchLen,
      readName: r.name,
      cigar: r.cigar,
      alignedQuery: matchedStrand === '-' ? reverseComplement(sw.alignedQuery) : sw.alignedQuery,
      alignedRead: matchedStrand === '-' ? reverseComplement(sw.alignedTarget) : sw.alignedTarget,
    });
  }

  if (rawMatches.length === 0) return [];

  // 50 bp 이내의 동일 앵커를 하나의 삽입 후보 위치로 클러스터링
  rawMatches.sort((a, b) => a.anchorPos - b.anchorPos);
  const clusters: BlastForeignReadHit[] = [];

  for (const rm of rawMatches) {
    const existing = clusters.find(
      (c) => Math.abs(c.anchorPos - rm.anchorPos) <= 50 && c.clipSide === rm.clipSide,
    );
    if (existing) {
      existing.supportingReads++;
      if (rm.matchLen > existing.matchLen || rm.identityPct > existing.identityPct) {
        existing.matchLen = rm.matchLen;
        existing.identityPct = rm.identityPct;
        existing.sampleReadName = rm.readName;
        existing.sampleCigar = rm.cigar;
        existing.alignedQueryDisplay = rm.alignedQuery;
        existing.alignedReadDisplay = rm.alignedRead;
      }
    } else {
      clusters.push({
        strainId,
        strainLabel,
        queryId: query.id,
        chrom,
        anchorPos: rm.anchorPos,
        clipSide: rm.clipSide,
        supportingReads: 1,
        identityPct: rm.identityPct,
        matchLen: rm.matchLen,
        strand: rm.strand,
        sampleReadName: rm.readName,
        sampleCigar: rm.cigar,
        alignedQueryDisplay: rm.alignedQuery,
        alignedReadDisplay: rm.alignedRead,
      });
    }
  }

  clusters.sort((a, b) => b.supportingReads - a.supportingReads || b.matchLen - a.matchLen);
  return clusters.slice(0, 20);
}
