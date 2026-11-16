/**
 * DNA 역상보 변환 및 S. cerevisiae 핵(Table 1) / 미토콘드리아(Table 3) 유전부호 번역 유틸리티
 */

const COMPLEMENT_MAP: Record<string, string> = {
  A: 'T',
  T: 'A',
  U: 'A',
  C: 'G',
  G: 'C',
  R: 'Y',
  Y: 'R',
  S: 'S',
  W: 'W',
  K: 'M',
  M: 'K',
  B: 'V',
  V: 'B',
  D: 'H',
  H: 'D',
  N: 'N',
  '-': '-',
  a: 't',
  t: 'a',
  u: 'a',
  c: 'g',
  g: 'c',
  r: 'y',
  y: 'r',
  s: 's',
  w: 'w',
  k: 'm',
  m: 'k',
  b: 'v',
  v: 'b',
  d: 'h',
  h: 'd',
  n: 'n',
};

/**
 * 5'->3' 역상보 서열 반환 (IUPAC 모호 염기 및 갭 '-' 보존)
 */
export function reverseComplement(seq: string): string {
  const out = new Array<string>(seq.length);
  for (let i = 0; i < seq.length; i++) {
    const ch = seq[seq.length - 1 - i];
    out[i] = COMPLEMENT_MAP[ch] ?? ch;
  }
  return out.join('');
}

// NCBI Translation Table 1: Standard Code
const TABLE_1: Record<string, string> = {
  TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L',
  TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
  TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*',
  TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
  CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
  CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
  CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M',
  ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
  AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
  GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
  GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G',
};

// NCBI Translation Table 3: Yeast Mitochondrial Code
// 차이점: ATA->M, CTT/CTC/CTA/CTG->T, TGA->W, CGA/CGC->미할당(X)
const TABLE_3: Record<string, string> = {
  ...TABLE_1,
  ATA: 'M',
  CTT: 'T',
  CTC: 'T',
  CTA: 'T',
  CTG: 'T',
  TGA: 'W',
  CGA: 'X',
  CGC: 'X',
};

export type TranslationTableId = 1 | 3;

export function translateDnaWithCodons(
  dnaSeq: string,
  tableId: TranslationTableId = 1,
  frame: 0 | 1 | 2 = 0,
): { aa: string; codons: string[] } {
  // 갭('-') 제거 후 순수 염기서열 기준 번역
  const clean = dnaSeq.replace(/-/g, '').toUpperCase().replace(/U/g, 'T');
  const table = tableId === 3 ? TABLE_3 : TABLE_1;
  const aa: string[] = [];
  const codons: string[] = [];

  for (let i = frame; i + 2 < clean.length; i += 3) {
    const codon = clean.slice(i, i + 3);
    codons.push(codon);
    aa.push(table[codon] ?? 'X');
  }
  return { aa: aa.join(''), codons };
}

export function translateDna(
  dnaSeq: string,
  tableId: TranslationTableId = 1,
  frame: 0 | 1 | 2 = 0,
): string {
  return translateDnaWithCodons(dnaSeq, tableId, frame).aa;
}

export type ProteinEffectStatus =
  | 'Identical'
  | 'Synonymous'
  | 'Missense'
  | 'Nonsense'
  | 'Frameshift';

export interface AaChangeDetail {
  aaPos1: number; // 1-based residue index
  refAA: string;
  strainAA: string;
  refCodon: string;
  strainCodon: string;
  label: string;
  isAmbiguousGap: boolean;
  isStopGain: boolean;
  isStopLoss: boolean;
}

export interface StrainProteinEvaluation {
  aa: string;
  codons: string[];
  status: ProteinEffectStatus;
  changes: string[];
  details: AaChangeDetail[];
  synonymousPositions: number[]; // 1-based AA positions where codon changed silently
  realChangeCount: number;
  ambiguousXCount: number;
}

export interface ProteinComparison {
  refAA: string;
  s1AA: string;
  s2AA: string;
  refCodons: string[];
  s1Codons: string[];
  s2Codons: string[];
  s1Status: ProteinEffectStatus;
  s2Status: ProteinEffectStatus;
  s1Changes: string[];
  s2Changes: string[];
  s1Details: AaChangeDetail[];
  s2Details: AaChangeDetail[];
  s1SynonymousPositions: number[];
  s2SynonymousPositions: number[];
  s1RealChangeCount: number;
  s2RealChangeCount: number;
  s1AmbiguousXCount: number;
  s2AmbiguousXCount: number;
}

const PURE_ACGT_CODON = /^[ACGT]{3}$/;

export function evaluateStrainProtein(
  refDna: string,
  strainDna: string,
  tableId: TranslationTableId = 1,
  frame: 0 | 1 | 2 = 0,
  precomputedRef?: { aa: string; codons: string[]; clean: string },
): StrainProteinEvaluation {
  const refClean =
    precomputedRef?.clean ?? refDna.replace(/-/g, '').toUpperCase().replace(/U/g, 'T');
  const refTranslated = precomputedRef ?? translateDnaWithCodons(refClean, tableId, frame);

  const strainClean = strainDna.replace(/-/g, '').toUpperCase().replace(/U/g, 'T');
  const strainTranslated = translateDnaWithCodons(strainClean, tableId, frame);

  const refAA = refTranslated.aa;
  const strainAA = strainTranslated.aa;
  const changes: string[] = [];
  const details: AaChangeDetail[] = [];
  const synonymousPositions: number[] = [];
  const minLen = Math.min(refAA.length, strainAA.length);
  let hasStopGain = false;
  let realChangeCount = 0;
  let ambiguousXCount = 0;

  for (let i = 0; i < minLen; i++) {
    const rA = refAA[i];
    const sA = strainAA[i];
    const rCodon = refTranslated.codons[i] ?? '';
    const sCodon = strainTranslated.codons[i] ?? '';

    if (rA !== sA) {
      const label = `${rA}${i + 1}${sA}`;
      const isAmbiguousGap = sA === 'X' && !PURE_ACGT_CODON.test(sCodon);
      const isStopGain = sA === '*' && rA !== '*';
      const isStopLoss = rA === '*' && sA !== '*';

      changes.push(label);
      details.push({
        aaPos1: i + 1,
        refAA: rA,
        strainAA: sA,
        refCodon: rCodon,
        strainCodon: sCodon,
        label,
        isAmbiguousGap,
        isStopGain,
        isStopLoss,
      });

      if (isAmbiguousGap) {
        ambiguousXCount++;
      } else {
        realChangeCount++;
        if (isStopGain) hasStopGain = true;
      }
    } else if (rCodon && sCodon && rCodon !== sCodon) {
      synonymousPositions.push(i + 1);
    }
  }

  let status: ProteinEffectStatus = 'Identical';
  if (Math.abs(strainClean.length - refClean.length) % 3 !== 0) {
    status = 'Frameshift';
  } else if (hasStopGain) {
    status = 'Nonsense';
  } else if (changes.length > 0 || refAA.length !== strainAA.length) {
    status = 'Missense';
  } else if (refClean !== strainClean) {
    status = 'Synonymous';
  }

  return {
    aa: strainAA,
    codons: strainTranslated.codons,
    status,
    changes,
    details,
    synonymousPositions,
    realChangeCount,
    ambiguousXCount,
  };
}

export function compareExtractedProteins(
  refDna: string,
  s1Dna: string,
  s2Dna: string | null,
  tableId: TranslationTableId,
  frame: 0 | 1 | 2 = 0,
): ProteinComparison {
  const refClean = refDna.replace(/-/g, '').toUpperCase().replace(/U/g, 'T');
  const refTranslated = translateDnaWithCodons(refClean, tableId, frame);
  const preRef = { ...refTranslated, clean: refClean };

  const s1Eval = evaluateStrainProtein(refDna, s1Dna, tableId, frame, preRef);
  const s2Eval: StrainProteinEvaluation =
    s2Dna !== null
      ? evaluateStrainProtein(refDna, s2Dna, tableId, frame, preRef)
      : {
          aa: '',
          codons: [],
          status: 'Identical',
          changes: [],
          details: [],
          synonymousPositions: [],
          realChangeCount: 0,
          ambiguousXCount: 0,
        };

  return {
    refAA: refTranslated.aa,
    s1AA: s1Eval.aa,
    s2AA: s2Eval.aa,
    refCodons: refTranslated.codons,
    s1Codons: s1Eval.codons,
    s2Codons: s2Eval.codons,
    s1Status: s1Eval.status,
    s2Status: s2Eval.status,
    s1Changes: s1Eval.changes,
    s2Changes: s2Eval.changes,
    s1Details: s1Eval.details,
    s2Details: s2Eval.details,
    s1SynonymousPositions: s1Eval.synonymousPositions,
    s2SynonymousPositions: s2Eval.synonymousPositions,
    s1RealChangeCount: s1Eval.realChangeCount,
    s2RealChangeCount: s2Eval.realChangeCount,
    s1AmbiguousXCount: s1Eval.ambiguousXCount,
    s2AmbiguousXCount: s2Eval.ambiguousXCount,
  };
}

