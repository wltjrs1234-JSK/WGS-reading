import { forEachCigarOp } from './cigar';
import { reverseComplement } from './geneticCode';
import type { AlignedRead, GffFeature } from './types';

export interface BarcodeScanOptions {
  minSupportReads?: number; // default 5
  minSupportRatio?: number; // default 0.12
  minInsertLen?: number; // default 8
  minClipLen?: number; // default 10
  minConsensusPct?: number; // default 72
  minMapQ?: number; // default 10 (유일 매핑 리드 우선, 멀티매핑 노이즈 배제)
  excludeRepeatLoci?: boolean; // default true (rDNA/Ty/텔로미어 반복 영역 제외)
  excludeEndogenousRefMatches?: boolean; // default true (sacCer3 내인성 서열 일치 제외)
  collectSubClipsForPairing?: boolean; // default true (양방향 클립 병합 시 minClipLen 달성 가능한 짧은 반쪽 클립 허용)
  customMotifSeq?: string; // optional extra sequence provided by user
}

export interface BarcodeBuildTelemetry {
  repeatLociFiltered: number;
  tandemOrArtifactFiltered: number;
}

export interface BarcodeRefContext {
  start: number;
  seq: string;
}

export interface BarcodeMotifDef {
  id: string;
  name: string;
  category: 'barseq_flank' | 'recombinase_scar' | 'cassette_marker' | 'cloning_site' | 'custom';
  seq: string;
  description: string;
}

export interface BarcodeMotifHit {
  motifId: string;
  name: string;
  category: BarcodeMotifDef['category'];
  strand: '+' | '-';
  matchedBases: number;
  isExact: boolean;
}

export interface RawBarcodeSignal {
  chrom: string;
  anchorPos: number; // 0-based genomic breakpoint coordinate
  signalType: 'insertion' | 'left_clip' | 'right_clip';
  seq: string; // insertion or clipped tail oriented away from genomic anchor for right/ins, or genomic left->right for left_clip
  readName: string;
}

export interface BarcodeCandidate {
  id: string;
  strainId: string;
  strainLabel: string;
  chrom: string;
  anchorPos: number; // 0-based
  signalType: 'insertion' | 'paired_insertion' | 'left_clip' | 'right_clip';
  classification: 'known_barcode_or_scar' | 'novel_insertion_barcode' | 'junction_breakpoint';
  specificity: 'strain_specific' | 'shared' | 'single_strain';
  assembledSeq: string;
  seqLength: number;
  extractedBarseq20mer?: string;
  supportReads: number;
  localDepth: number;
  supportRatioPct: number;
  consensusAgreementPct: number;
  motifHits: BarcodeMotifHit[];
  overlappingGenes: string[];
  sharedWithStrains: string[];
  sampleReadNames: string[];
}

export const BUILT_IN_BARCODE_MOTIFS: BarcodeMotifDef[] = [
  {
    id: 'yko_u1',
    name: 'Bar-seq UPTAG U1',
    category: 'barseq_flank',
    seq: 'GATGTCCACGAGGTCTCT',
    description: 'Yeast Deletion / Bar-seq UPTAG 5\' 유니버설 프라이머 플랭크',
  },
  {
    id: 'yko_u2_mx4',
    name: 'Bar-seq U2 / pFA6a-MX4 F',
    category: 'barseq_flank',
    seq: 'CGTACGCTGCAGGTCGAC',
    description: 'UPTAG 3\' 플랭크 겸 KanMX/NatMX/HphMX 카세트 5\' 링커',
  },
  {
    id: 'yko_d1',
    name: 'Bar-seq DNTAG D1',
    category: 'barseq_flank',
    seq: 'CGGTGTCGGTCTCGTAG',
    description: 'Yeast Deletion / Bar-seq DNTAG 5\' 유니버설 프라이머 플랭크',
  },
  {
    id: 'yko_d2_mxt',
    name: 'Bar-seq D2 / MX-Term R',
    category: 'barseq_flank',
    seq: 'ATCGATGAATTCGAGCTCG',
    description: 'DNTAG 3\' 플랭크 겸 pFA6a MX 터미네이터 링커',
  },
  {
    id: 'loxp',
    name: 'Cre loxP Scar',
    category: 'recombinase_scar',
    seq: 'ATAACTTCGTATAGCATACATTATACGAAGTTAT',
    description: 'Cre-loxP 마커 제거 후 남는 34 bp 재조합 스카',
  },
  {
    id: 'lox2272',
    name: 'Cre lox2272',
    category: 'recombinase_scar',
    seq: 'ATAACTTCGTATAGGATACTTTATACGAAGTTAT',
    description: '이형 lox2272 재조합 스카 서열',
  },
  {
    id: 'lox71_arm',
    name: 'lox71 / lox66 mutant arm',
    category: 'recombinase_scar',
    seq: 'TACCGTTCGTATAGCATACATTATACGAAGTTAT',
    description: '변이 lox71/lox66 스카 서열',
  },
  {
    id: 'frt',
    name: 'Flp FRT Scar',
    category: 'recombinase_scar',
    seq: 'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC',
    description: 'Flp-FRT 마커 제거 후 남는 34 bp 재조합 스카',
  },
  {
    id: 'iscei',
    name: 'I-SceI Site',
    category: 'cloning_site',
    seq: 'TAGGGATAACAGGGTAAT',
    description: 'I-SceI 메가뉴클레아제 인식 서열 (18 bp)',
  },
  {
    id: 'tef_prom_5p',
    name: 'A. gossypii TEF Promoter (MX)',
    category: 'cassette_marker',
    seq: 'GACATGGAGGCCCAGAATACCC',
    description: 'KanMX/NatMX/HphMX 카세트 내부 TEF 프로모터 시그니처',
  },
  {
    id: 'tef_term_3p',
    name: 'A. gossypii TEF Terminator (MX)',
    category: 'cassette_marker',
    seq: 'AGTATAGCGACCAGCATTCACATACG',
    description: 'KanMX/NatMX/HphMX 카세트 내부 TEF 터미네이터 시그니처',
  },
  {
    id: 't7_prom',
    name: 'T7 Promoter / Vector',
    category: 'cloning_site',
    seq: 'TAATACGACTCACTATAGGG',
    description: '플라스미드/합성 컨스트럭트 T7 프로모터 스카',
  },
  {
    id: 'm13_fwd',
    name: 'M13 Forward Vector',
    category: 'cloning_site',
    seq: 'GTAAAACGACGGCCAGT',
    description: '표준 클로닝 벡터 M13 Forward 링커',
  },
];

// Illumina / Nextera 시퀀싱 어댑터 오염 판별용 시그니처
const SEQUENCING_ADAPTER_SIGNATURES = [
  'AGATCGGAAGAGC', // TruSeq Universal / Read1 / Read2 Adapter core
  'CTGTCTCTTATACACATCT', // Nextera Transposase
  'TCGTCGGCAGCGTC', // Nextera i5
  'GTCTCGTGGGCTCGG', // Nextera i7
];

export const HOTSPOT_LOCUS_NAMES = [
  'HO',
  'YDL227C',
  'URA3',
  'YEL021W',
  'HIS3',
  'YOR202W',
  'LEU2',
  'YCL018W',
  'TRP1',
  'YDR007W',
  'CAN1',
  'YEL063C',
  'LYS2',
  'YBR115C',
  'MET15',
  'YLR303W',
  'ADE2',
  'YOR128C',
  'YPR015C', // YPRCΔ15 safe harbor 인접
];

const REPEAT_FEATURE_TYPES = new Set([
  'rRNA_gene',
  'transposable_element_gene',
  'telomere',
  'long_terminal_repeat',
]);

/**
 * rDNA 다중복제수 배열(chrXII RDN1/ASP3 구간), rRNA, Ty 트랜스포존, 텔로미어 등 내인성 반복 유전자좌를 판별합니다.
 */
export function isRepetitiveYeastRegion(
  chrom: string,
  pos: number,
  features: GffFeature[],
): boolean {
  if (chrom === 'chrXII' && pos >= 450_000 && pos <= 492_000) {
    return true;
  }
  for (const f of features) {
    if (f.chrom !== chrom) continue;
    if (
      REPEAT_FEATURE_TYPES.has(f.type) ||
      f.name.startsWith('RDN') ||
      f.name.startsWith('ETS') ||
      f.name.startsWith('ITS')
    ) {
      if (pos >= f.start - 1_500 && pos <= f.end + 1_500) {
        return true;
      }
    }
  }
  return false;
}

/**
 * 1~24 bp 주기의 직렬 연쇄 반복(VNTR/Minisatellite/SSR)이 지배적인 서열인지 탐지합니다.
 */
export function hasDominantTandemRepeat(seq: string): boolean {
  const upper = seq.toUpperCase();
  const len = upper.length;
  if (len < 12) return false;

  const maxPeriod = len < 16 ? Math.min(4, Math.floor(len / 2)) : Math.min(24, Math.floor(len / 2));
  for (let p = 1; p <= maxPeriod; p++) {
    const comparisons = len - p;
    if (comparisons < 8) continue;

    let matches = 0;
    let curRun = 0;
    let maxRun = 0;
    for (let i = 0; i < comparisons; i++) {
      if (upper.charCodeAt(i) === upper.charCodeAt(i + p)) {
        matches++;
        curRun++;
        if (curRun > maxRun) maxRun = curRun;
      } else {
        curRun = 0;
      }
    }

    const covered = maxRun + p;
    if (p <= 4 && covered >= 12 && covered / len >= 0.65) {
      return true;
    }
    if (covered >= Math.max(18, p * 2.2) && covered / len >= 0.55) {
      return true;
    }
    if (comparisons >= 12 && matches / comparisons >= 0.72) {
      return true;
    }
  }
  return false;
}

function matchesRefWithUpToOneMismatch(
  sub: string,
  refWindow: string,
  edge: 'exact_start' | 'exact_end',
): boolean {
  if (refWindow.includes(sub)) return true;
  if (sub.length < 15) return false;
  const seedLen = 11;
  const firstSeed = sub.slice(0, seedLen);
  let idx = refWindow.indexOf(firstSeed);
  while (idx !== -1) {
    if (idx + sub.length <= refWindow.length) {
      let mm = 0;
      for (let k = seedLen; k < sub.length; k++) {
        if (sub.charCodeAt(k) !== refWindow.charCodeAt(idx + k)) {
          mm++;
          if (mm > 1 || (edge === 'exact_end' && k >= sub.length - 4)) break;
        }
      }
      if (mm <= 1 && (edge !== 'exact_end' || sub.slice(-4) === refWindow.slice(idx + sub.length - 4, idx + sub.length))) {
        return true;
      }
    }
    idx = refWindow.indexOf(firstSeed, idx + 1);
  }
  const tailSeed = sub.slice(-seedLen);
  idx = refWindow.indexOf(tailSeed);
  while (idx !== -1) {
    const startIdx = idx - (sub.length - seedLen);
    if (startIdx >= 0) {
      let mm = 0;
      for (let k = 0; k < sub.length - seedLen; k++) {
        if (sub.charCodeAt(k) !== refWindow.charCodeAt(startIdx + k)) {
          mm++;
          if (mm > 1 || (edge === 'exact_start' && k < 4)) break;
        }
      }
      if (mm <= 1 && (edge !== 'exact_start' || sub.slice(0, 4) === refWindow.slice(startIdx, startIdx + 4))) {
        return true;
      }
    }
    idx = refWindow.indexOf(tailSeed, idx + 1);
  }
  return false;
}

/**
 * Soft-clip 합의 서열의 외곽 끝에 붙은 레퍼런스 플랭킹 염기(>=11 bp 일치)를 정밀 트리밍하여
 * 순수 비-레퍼런스 삽입 서열만 남깁니다.
 */
export function trimReferenceFlanksFromClip(
  seq: string,
  signalType: 'insertion' | 'left_clip' | 'right_clip',
  anchorPos: number,
  refWindowStart: number,
  refWindowSeq: string,
): string {
  if (signalType === 'insertion' || !refWindowSeq) return seq;
  const minRefMatch = 11;
  if (seq.length < minRefMatch) return seq;
  const relAnchor = anchorPos - refWindowStart;

  if (signalType === 'right_clip') {
    const wStart = Math.max(0, relAnchor - 20);
    const wEnd = Math.min(refWindowSeq.length, relAnchor + Math.max(180, seq.length + 60));
    if (wEnd - wStart < minRefMatch) return seq;
    const refAfter = refWindowSeq.slice(wStart, wEnd);

    for (let cut = 0; cut <= seq.length - minRefMatch; cut++) {
      const suffix = seq.slice(cut);
      if (matchesRefWithUpToOneMismatch(suffix, refAfter, 'exact_start')) {
        return seq.slice(0, cut);
      }
    }
    return seq;
  }

  const wStart = Math.max(0, relAnchor - Math.max(180, seq.length + 60));
  const wEnd = Math.min(refWindowSeq.length, relAnchor + 20);
  if (wEnd - wStart < minRefMatch) return seq;
  const refBefore = refWindowSeq.slice(wStart, wEnd);

  for (let cut = seq.length; cut >= minRefMatch; cut--) {
    const prefix = seq.slice(0, cut);
    if (matchesRefWithUpToOneMismatch(prefix, refBefore, 'exact_end')) {
      return seq.slice(cut);
    }
  }
  return seq;
}

/**
 * 염기 다양성(Shannon Entropy), 어댑터, 호모폴리머 및 1~24 bp VNTR 연쇄 반복 검사를 통해 생물학적/기술적 노이즈를 걸러냅니다.
 */
export function isLikelyArtifactSequence(seq: string): boolean {
  const upper = seq.toUpperCase();
  if (upper.length < 8) return true;

  // Illumina 어댑터 일치 여부
  for (const adp of SEQUENCING_ADAPTER_SIGNATURES) {
    if (upper.includes(adp) || reverseComplement(upper).includes(adp)) {
      return true;
    }
    // 클립 길이가 짧아 어댑터 앞 10 bp만 걸린 경우 체크
    if (upper.length >= 10 && (adp.startsWith(upper.slice(0, 10)) || adp.startsWith(upper.slice(-10)))) {
      return true;
    }
  }

  // 호모폴리머(동일 염기 9회 이상) 또는 N 다수
  const nCount = (upper.match(/N/g) ?? []).length;
  if (nCount / upper.length > 0.2) return true;
  if (/(A{9,}|C{9,}|G{9,}|T{9,})/.test(upper)) return true;

  // 디뉴클레오타이드 단순 반복 (ATx5, CAx5 등)
  if (/(AT){5,}|(TA){5,}|(CA){5,}|(AC){5,}|(GT){5,}|(TG){5,}|(GA){5,}|(AG){5,}|(CT){5,}|(TC){5,}/.test(upper)) {
    return true;
  }

  // 1~24 bp VNTR / 미니새틀라이트 직렬 반복
  if (hasDominantTandemRepeat(upper)) {
    return true;
  }

  // Shannon entropy (< 1.32 bits이면 극단적 편향 반복서열)
  const counts = { A: 0, C: 0, G: 0, T: 0 };
  let valid = 0;
  for (let i = 0; i < upper.length; i++) {
    const ch = upper[i];
    if (ch === 'A' || ch === 'C' || ch === 'G' || ch === 'T') {
      counts[ch]++;
      valid++;
    }
  }
  if (valid < 8) return true;
  let entropy = 0;
  for (const c of Object.values(counts)) {
    if (c === 0) continue;
    const p = c / valid;
    entropy -= p * Math.log2(p);
  }
  return entropy < 1.32;
}

function longestCommonSubstringLen(a: string, b: string): number {
  const k = Math.min(10, a.length, b.length);
  if (k < 8) return 0;
  let best = 0;
  for (let i = 0; i <= a.length - k; i++) {
    const sub = a.slice(i, i + k);
    let pos = b.indexOf(sub);
    while (pos !== -1) {
      let l = k;
      while (i + l < a.length && pos + l < b.length && a[i + l] === b[pos + l]) {
        l++;
      }
      if (l > best) best = l;
      pos = b.indexOf(sub, pos + 1);
    }
  }
  return best;
}

export function matchBarcodeMotifs(
  seq: string,
  extraMotif?: string,
): BarcodeMotifHit[] {
  const upper = seq.toUpperCase();
  const rc = reverseComplement(upper);
  const defs = [...BUILT_IN_BARCODE_MOTIFS];

  if (extraMotif) {
    const cleaned = extraMotif.replace(/[^ACGT]/gi, '').toUpperCase();
    if (cleaned.length >= 10) {
      defs.unshift({
        id: 'user_custom',
        name: '사용자 지정 의심 서열',
        category: 'custom',
        seq: cleaned,
        description: '사용자가 직접 입력한 의심 바코드/링커 서열',
      });
    }
  }

  const hits: BarcodeMotifHit[] = [];
  for (const def of defs) {
    const motif = def.seq;
    const fExact = upper.includes(motif);
    const rExact = rc.includes(motif);

    if (fExact || rExact) {
      hits.push({
        motifId: def.id,
        name: def.name,
        category: def.category,
        strand: fExact ? '+' : '-',
        matchedBases: motif.length,
        isExact: true,
      });
      continue;
    }

    const fPartial = longestCommonSubstringLen(upper, motif);
    const rPartial = longestCommonSubstringLen(rc, motif);
    const bestPartial = Math.max(fPartial, rPartial);
    if (bestPartial >= Math.min(12, Math.floor(motif.length * 0.7))) {
      hits.push({
        motifId: def.id,
        name: def.name,
        category: def.category,
        strand: fPartial >= rPartial ? '+' : '-',
        matchedBases: bestPartial,
        isExact: false,
      });
    }
  }

  return hits;
}

/**
 * U1..U2 또는 D1..D2 플랭크 사이에 끼어있는 Bar-seq 20 bp 분자 바코드를 자동 추출합니다.
 */
export function extractBarseqCore(seq: string): string | undefined {
  const checkOriented = (s: string): string | undefined => {
    const pairs: Array<[string, string]> = [
      ['GATGTCCACGAGGTCTCT', 'CGTACGCTGCAGGTCGAC'], // U1 .. U2
      ['CGGTGTCGGTCTCGTAG', 'ATCGATGAATTCGAGCTCG'], // D1 .. D2
    ];
    for (const [leftF, rightF] of pairs) {
      // 부분 경계(최소 8 bp 끝/시작) 포함 탐지
      const li = s.indexOf(leftF.slice(-10));
      const ri = s.indexOf(rightF.slice(0, 10));
      if (li !== -1 && ri !== -1 && ri > li + 10) {
        const innerStart = li + 10;
        const inner = s.slice(innerStart, ri);
        if (inner.length >= 15 && inner.length <= 26) {
          return inner;
        }
      }
    }
    return undefined;
  };
  return checkOriented(seq.toUpperCase()) ?? checkOriented(reverseComplement(seq.toUpperCase()));
}

/**
 * 단일 구간의 AlignedRead 목록에서 Insertion(I >= minInsertLen) 및 Soft-clip(S >= minClipLen) 신호를 추출합니다.
 */
export function collectBarcodeSignalsFromReads(
  chrom: string,
  reads: AlignedRead[],
  rangeStart: number,
  rangeEnd: number,
  options: BarcodeScanOptions = {},
): RawBarcodeSignal[] {
  const minIns = options.minInsertLen ?? 8;
  const reqMinClip = options.minClipLen ?? 10;
  const minClip =
    options.collectSubClipsForPairing !== false ? Math.min(reqMinClip, 11) : reqMinClip;
  const minMapQ = options.minMapQ ?? 10;
  const signals: RawBarcodeSignal[] = [];

  for (const r of reads) {
    if ((r.flags & 0x704) !== 0 || r.mapq < minMapQ) continue;
    const seq = r.seq && r.seq !== '*' ? r.seq.toUpperCase() : '';
    if (!seq) continue;

    // 1차 빠른 필터: 긴 S 또는 I 연산 가능성이 없는 순수 match 리드는 즉시 건너뜀
    if (!r.cigar.includes('S') && !r.cigar.includes('I')) continue;

    let refPos = r.start;
    let queryPos = 0;
    let opIndex = 0;
    let sawFirstH = false;

    forEachCigarOp(r.cigar, (op, len) => {
      switch (op) {
        case 'H': {
          if (opIndex === 0) sawFirstH = true;
          break;
        }
        case 'P':
          break;
        case 'M':
        case '=':
        case 'X': {
          refPos += len;
          queryPos += len;
          break;
        }
        case 'D':
        case 'N': {
          refPos += len;
          break;
        }
        case 'I': {
          if (len >= minIns) {
            const anchorPos = refPos - 1;
            if (anchorPos >= rangeStart && anchorPos < rangeEnd && queryPos + len <= seq.length) {
              signals.push({
                chrom,
                anchorPos,
                signalType: 'insertion',
                seq: seq.slice(queryPos, queryPos + len),
                readName: r.name,
              });
            }
          }
          queryPos += len;
          break;
        }
        case 'S': {
          const isLeading = opIndex === 0 || (opIndex === 1 && sawFirstH);
          if (len >= minClip && queryPos + len <= seq.length) {
            const anchorPos = isLeading ? r.start : refPos;
            if (anchorPos >= rangeStart && anchorPos < rangeEnd) {
              signals.push({
                chrom,
                anchorPos,
                signalType: isLeading ? 'left_clip' : 'right_clip',
                seq: seq.slice(queryPos, queryPos + len),
                readName: r.name,
              });
            }
          }
          queryPos += len;
          break;
        }
      }
      opIndex++;
    });
  }

  return signals;
}

/**
 * 앵커 경계부터 바깥쪽으로 염기 투표를 수행해 다수 리드 합의 클립/삽입 서열을 복원합니다.
 * - left_clip: 리드의 3' 쪽(breakpoint 인접 끝)부터 왼쪽으로 맞추어야 서로 다른 클립 길이도 정확히 정렬됨!
 * - right_clip / insertion: 리드의 5' 쪽(breakpoint 인접 끝)부터 오른쪽으로 정렬.
 */
function assembleConsensusFromBucket(
  seqs: string[],
  alignFromRightEnd: boolean,
): { consensus: string; agreementPct: number } {
  if (seqs.length === 0) return { consensus: '', agreementPct: 0 };

  const lengths = seqs.map((s) => s.length).sort((a, b) => a - b);
  // 최소 30% 이상 리드가 도달하는 길이까지만 안정적으로 합의 연장
  const maxStableLen = lengths[Math.max(0, Math.floor(lengths.length * 0.7) - 1)] || lengths[lengths.length - 1];
  const minActiveSupport = Math.max(2, Math.ceil(seqs.length * 0.3));

  const bases: string[] = [];
  let sumMajorityRatio = 0;

  for (let offset = 0; offset < maxStableLen; offset++) {
    let a = 0;
    let c = 0;
    let g = 0;
    let t = 0;
    let total = 0;

    for (const s of seqs) {
      if (offset >= s.length) continue;
      const ch = alignFromRightEnd ? s[s.length - 1 - offset] : s[offset];
      if (ch === 'A') a++;
      else if (ch === 'C') c++;
      else if (ch === 'G') g++;
      else if (ch === 'T') t++;
      else continue;
      total++;
    }

    if (total < minActiveSupport && bases.length >= 10) break;
    if (total === 0) break;

    const maxCnt = Math.max(a, c, g, t);
    const bestBase = maxCnt === a ? 'A' : maxCnt === c ? 'C' : maxCnt === g ? 'G' : 'T';
    bases.push(bestBase);
    sumMajorityRatio += maxCnt / total;
  }

  if (alignFromRightEnd) bases.reverse();
  const agreementPct = bases.length > 0 ? Math.round((sumMajorityRatio / bases.length) * 1000) / 10 : 0;
  return { consensus: bases.join(''), agreementPct };
}

function countReadsCoveringAnchor(reads: AlignedRead[], pos: number): number {
  let depth = 0;
  for (const r of reads) {
    if ((r.flags & 0x704) !== 0) continue;
    if (r.start <= pos && r.end > pos) depth++;
  }
  return depth;
}

/**
 * 좌측/우측 clip 또는 insertion 합의에서 반대측 clip과 10 bp 이상 겹치는 쌍을 찾아
 * 단일 완전 바코드(paired_insertion)로 병합 시도합니다.
 */
function tryMergeOverlappingRightLeft(rightSeq: string, leftSeq: string): string | null {
  const minOv = 10;
  const maxOv = Math.min(rightSeq.length, leftSeq.length);
  for (let ov = maxOv; ov >= minOv; ov--) {
    if (rightSeq.slice(-ov) === leftSeq.slice(0, ov)) {
      return rightSeq + leftSeq.slice(ov);
    }
  }
  return null;
}

export function buildBarcodeCandidatesForRegion(
  strainId: string,
  strainLabel: string,
  signals: RawBarcodeSignal[],
  readsForDepth: AlignedRead[],
  features: GffFeature[],
  options: BarcodeScanOptions = {},
  localRefContext?: BarcodeRefContext,
  telemetry?: BarcodeBuildTelemetry,
): BarcodeCandidate[] {
  const minSupport = options.minSupportReads ?? 5;
  const minRatio = options.minSupportRatio ?? 0.12;
  const minAgreement = options.minConsensusPct ?? 72;
  const reqMinIns = options.minInsertLen ?? 8;
  const reqMinClip = options.minClipLen ?? 10;
  const pairHalfMinClip =
    options.collectSubClipsForPairing !== false ? Math.min(reqMinClip, 11) : reqMinClip;

  // (chrom, signalType) 별 좌표 ±3 bp 내 클러스터링
  const sorted = [...signals].sort(
    (a, b) => a.chrom.localeCompare(b.chrom) || a.anchorPos - b.anchorPos,
  );

  interface Bucket {
    chrom: string;
    anchorPos: number;
    anchorSum: number;
    signalType: RawBarcodeSignal['signalType'];
    seqs: string[];
    readNames: Set<string>;
  }

  const buckets: Bucket[] = [];
  for (const sig of sorted) {
    const match = buckets.find(
      (b) =>
        b.chrom === sig.chrom &&
        b.signalType === sig.signalType &&
        Math.abs(b.anchorPos - sig.anchorPos) <= 3,
    );
    if (match) {
      match.seqs.push(sig.seq);
      match.readNames.add(sig.readName);
      match.anchorSum += sig.anchorPos;
      match.anchorPos = Math.round(match.anchorSum / match.seqs.length);
    } else {
      buckets.push({
        chrom: sig.chrom,
        anchorPos: sig.anchorPos,
        anchorSum: sig.anchorPos,
        signalType: sig.signalType,
        seqs: [sig.seq],
        readNames: new Set([sig.readName]),
      });
    }
  }

  interface AssembledCluster {
    chrom: string;
    anchorPos: number;
    signalType: RawBarcodeSignal['signalType'];
    consensus: string;
    agreementPct: number;
    supportReads: number;
    localDepth: number;
    readNames: string[];
    consumed: boolean;
  }

  const assembled: AssembledCluster[] = [];

  for (const b of buckets) {
    if (b.readNames.size < minSupport) continue;

    const { consensus: rawConsensus, agreementPct } = assembleConsensusFromBucket(
      b.seqs,
      b.signalType === 'left_clip',
    );
    if (agreementPct < minAgreement) continue;

    const trimmedConsensus = localRefContext
      ? trimReferenceFlanksFromClip(
          rawConsensus,
          b.signalType,
          b.anchorPos,
          localRefContext.start,
          localRefContext.seq,
        )
      : rawConsensus;
    const consensus = trimmedConsensus;

    const isInsertion = b.signalType === 'insertion';
    if (consensus.length < (isInsertion ? reqMinIns : pairHalfMinClip)) {
      continue;
    }
    const meetsFullLen = isInsertion || consensus.length >= reqMinClip;

    const preMotifHits = matchBarcodeMotifs(consensus, options.customMotifSeq);
    const preBarseq20 = extractBarseqCore(consensus);
    const hasKnownMotif = preMotifHits.length > 0 || Boolean(preBarseq20);

    if (!hasKnownMotif && isLikelyArtifactSequence(consensus)) {
      if (telemetry && meetsFullLen) telemetry.tandemOrArtifactFiltered++;
      continue;
    }

    if (
      options.excludeRepeatLoci !== false &&
      !hasKnownMotif &&
      isRepetitiveYeastRegion(b.chrom, b.anchorPos, features)
    ) {
      if (telemetry && meetsFullLen) telemetry.repeatLociFiltered++;
      continue;
    }

    const probePos =
      b.signalType === 'left_clip' ? b.anchorPos : Math.max(0, b.anchorPos - 1);
    const localDepth = Math.max(b.readNames.size, countReadsCoveringAnchor(readsForDepth, probePos));
    if (b.readNames.size / localDepth < minRatio) continue;

    assembled.push({
      chrom: b.chrom,
      anchorPos: b.anchorPos,
      signalType: b.signalType,
      consensus,
      agreementPct,
      supportReads: b.readNames.size,
      localDepth,
      readNames: [...b.readNames].slice(0, 6),
      consumed: false,
    });
  }

  // 동일 유전자좌(±20 bp)에 이미 완전한 리드 내부 Insertion이 존재하면 중복 Soft-clip 신호 제거
  const insertions = assembled.filter((c) => c.signalType === 'insertion');
  for (const cl of assembled) {
    if (cl.signalType === 'insertion') continue;
    const redundantIns = insertions.find(
      (ins) =>
        ins.chrom === cl.chrom &&
        Math.abs(ins.anchorPos - cl.anchorPos) <= 20 &&
        (Math.abs(ins.anchorPos - cl.anchorPos) <= 6 ||
          longestCommonSubstringLen(cl.consensus, ins.consensus) >= 10 ||
          longestCommonSubstringLen(cl.consensus, reverseComplement(ins.consensus)) >= 10),
    );
    if (redundantIns) {
      cl.consumed = true;
    }
  }

  const candidates: BarcodeCandidate[] = [];

  // 인접한 right_clip + left_clip 쌍(0~45 bp 이내) 병합 시도
  for (let i = 0; i < assembled.length; i++) {
    const a = assembled[i];
    if (a.consumed) continue;

    let mergedSeq: string | null = null;
    let partner: AssembledCluster | null = null;

    if (a.signalType === 'right_clip' || a.signalType === 'left_clip') {
      for (let j = i + 1; j < assembled.length; j++) {
        const b = assembled[j];
        if (b.consumed || b.chrom !== a.chrom) continue;
        if (b.signalType === a.signalType || b.signalType === 'insertion') continue;
        if (Math.abs(b.anchorPos - a.anchorPos) <= 45) {
          const right = a.signalType === 'right_clip' ? a : b;
          const left = a.signalType === 'left_clip' ? a : b;
          const trySeq = tryMergeOverlappingRightLeft(right.consensus, left.consensus);
          if (trySeq) {
            mergedSeq = trySeq;
            partner = b;
            b.consumed = true;
            break;
          }
        }
      }
    }

    const finalSeq = mergedSeq ?? a.consensus;
    const finalMinLen = a.signalType === 'insertion' ? reqMinIns : reqMinClip;
    if (finalSeq.length < finalMinLen) continue;
    const motifHits = matchBarcodeMotifs(finalSeq, options.customMotifSeq);
    const barseq20 = extractBarseqCore(finalSeq);
    const hasFinalMotif = motifHits.length > 0 || Boolean(barseq20);
    if (partner && !hasFinalMotif && isLikelyArtifactSequence(finalSeq)) {
      if (telemetry) telemetry.tandemOrArtifactFiltered++;
      continue;
    }

    const signalType: BarcodeCandidate['signalType'] = partner
      ? 'paired_insertion'
      : a.signalType;

    const classification: BarcodeCandidate['classification'] =
      hasFinalMotif
        ? 'known_barcode_or_scar'
        : signalType === 'insertion' || signalType === 'paired_insertion'
          ? 'novel_insertion_barcode'
          : 'junction_breakpoint';

    const supportReads = partner
      ? Math.max(a.supportReads, partner.supportReads)
      : a.supportReads;
    const localDepth = partner ? Math.max(a.localDepth, partner.localDepth) : a.localDepth;
    const agreementPct = partner
      ? Math.round(((a.agreementPct + partner.agreementPct) / 2) * 10) / 10
      : a.agreementPct;

    const ovGenes = features
      .filter((f) => f.start <= a.anchorPos + 50 && f.end >= a.anchorPos - 50)
      .slice(0, 3)
      .map((f) => `${f.name} (${f.strand})`);

    candidates.push({
      id: `${strainId}_${a.chrom}_${a.anchorPos}_${signalType}`,
      strainId,
      strainLabel,
      chrom: a.chrom,
      anchorPos: a.anchorPos,
      signalType,
      classification,
      specificity: 'single_strain',
      assembledSeq: finalSeq,
      seqLength: finalSeq.length,
      extractedBarseq20mer: barseq20,
      supportReads,
      localDepth,
      supportRatioPct: Math.min(100, Math.round((supportReads / Math.max(1, localDepth)) * 1000) / 10),
      consensusAgreementPct: agreementPct,
      motifHits,
      overlappingGenes: ovGenes,
      sharedWithStrains: [],
      sampleReadNames: a.readNames,
    });
  }

  // 좌표가 ±20 bp 이내에서 겹치는 유사 클립/삽입 중복 제거 (더 강한 지지 리드/긴 서열 유지)
  const deduped: BarcodeCandidate[] = [];
  for (const c of candidates) {
    const dup = deduped.find(
      (d) =>
        d.chrom === c.chrom &&
        Math.abs(d.anchorPos - c.anchorPos) <= 20 &&
        (d.signalType === c.signalType ||
          d.signalType === 'paired_insertion' ||
          c.signalType === 'paired_insertion') &&
        (longestCommonSubstringLen(d.assembledSeq, c.assembledSeq) >= 12 ||
          longestCommonSubstringLen(d.assembledSeq, reverseComplement(c.assembledSeq)) >= 12),
    );
    if (!dup) {
      deduped.push(c);
    } else if (
      c.supportReads > dup.supportReads ||
      (c.supportReads === dup.supportReads && c.seqLength > dup.seqLength)
    ) {
      Object.assign(dup, c);
    }
  }

  return deduped;
}

/**
 * 다수 균주(예: WT vs Mutant) 결과를 교차 비교해 균주 특이 바코드 후보를 우선 표시합니다.
 */
export function annotateCrossStrainSpecificity(
  allCandidates: BarcodeCandidate[],
  totalStrainsCount: number,
): BarcodeCandidate[] {
  if (totalStrainsCount <= 1) {
    return sortBarcodeCandidates(
      allCandidates.map((c) => ({ ...c, specificity: 'single_strain' as const })),
    );
  }

  const annotated = allCandidates.map((cand) => {
    const peers = new Set<string>();
    for (const other of allCandidates) {
      if (other.strainId === cand.strainId) continue;
      if (other.chrom !== cand.chrom) continue;
      if (Math.abs(other.anchorPos - cand.anchorPos) <= 25) {
        peers.add(other.strainLabel);
      }
    }
    return {
      ...cand,
      specificity: peers.size === 0 ? ('strain_specific' as const) : ('shared' as const),
      sharedWithStrains: [...peers],
    };
  });

  return sortBarcodeCandidates(annotated);
}

function sortBarcodeCandidates(list: BarcodeCandidate[]): BarcodeCandidate[] {
  const specificityRank = (s: BarcodeCandidate['specificity']) =>
    s === 'strain_specific' ? 3 : s === 'single_strain' ? 2 : 1;
  const classRank = (c: BarcodeCandidate['classification']) =>
    c === 'known_barcode_or_scar' ? 3 : c === 'novel_insertion_barcode' ? 2 : 1;

  return [...list].sort(
    (a, b) =>
      specificityRank(b.specificity) - specificityRank(a.specificity) ||
      classRank(b.classification) - classRank(a.classification) ||
      b.supportReads - a.supportReads ||
      b.seqLength - a.seqLength,
  );
}
