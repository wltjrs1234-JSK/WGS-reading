import {
  expandIupacBase,
  reverseComplement,
  translateCodon,
  type TranslationTableId,
} from './geneticCode';
import type { AnnotationIndex } from './gff3Parser';
import type {
  ContigInfo,
  GffFeature,
  VariantAnnotationMeta,
  VariantConsequence,
  VariantImpact,
  VariantRecord,
} from './types';
import { isMitochondrialChrom } from './yeastAliases';

export interface ExonInterval {
  start: number; // 0-based inclusive
  end: number; // 0-based exclusive
}

export interface GeneVariantSummary {
  geneKey: string;
  geneName: string;
  systematicName?: string;
  chrom: string;
  start: number;
  end: number;
  strand: '+' | '-' | '.';
  featureType: string;
  highestImpact: VariantImpact;
  variantCount: number;
  highCount: number;
  moderateCount: number;
  lowCount: number;
  modifierCount: number;
  topExamples: string[]; // e.g. ["p.T511I", "p.Q142*"]
  categories: string[];
}

export const IMPACT_ORDER: Record<VariantImpact, number> = {
  HIGH: 4,
  MODERATE: 3,
  LOW: 2,
  MODIFIER: 1,
};

/**
 * GFF Feature에서 5'->3' 코딩 순서(Exon 배열)를 반환합니다.
 */
export function getOrderedCodingExons(feat: GffFeature): ExonInterval[] {
  const cds = feat.subFeatures
    .filter((sf) => sf.type === 'CDS' || sf.type === 'exon')
    .map((sf) => ({ start: sf.start, end: sf.end }));
  const list = cds.length > 0 ? cds : [{ start: feat.start, end: feat.end }];
  list.sort((a, b) => a.start - b.start);
  if (feat.strand === '-') {
    list.reverse();
  }
  return list;
}

function getTargetAllelesToEvaluate(v: VariantRecord): Array<{ allele: string; isS2: boolean }> {
  if (v.category === 'strain2_only') {
    return v.s2Allele ? [{ allele: v.s2Allele, isS2: true }] : [];
  }
  if (v.category === 'divergent') {
    const out: Array<{ allele: string; isS2: boolean }> = [];
    if (v.s1Allele) out.push({ allele: v.s1Allele, isS2: false });
    if (v.s2Allele) out.push({ allele: v.s2Allele, isS2: true });
    return out;
  }
  return v.s1Allele ? [{ allele: v.s1Allele, isS2: false }] : [];
}

/**
 * 5'->3' 정렬된 엑손 목록에서 특정 0-based 게놈 좌표가 CDS 내 몇 번째 염기(0-based)인지 계산합니다.
 */
function locatePositionInExons(
  genomicPos: number,
  orderedExons: ExonInterval[],
  strand: '+' | '-' | '.',
): { cdsOffset: number; exonNumber: number; totalCdsLen: number } | null {
  let cumulative = 0;
  let foundOffset: number | null = null;
  let foundExonNum = 0;

  for (let i = 0; i < orderedExons.length; i++) {
    const ex = orderedExons[i];
    const exLen = Math.max(0, ex.end - ex.start);
    if (foundOffset === null && genomicPos >= ex.start && genomicPos < ex.end) {
      foundOffset =
        strand === '-'
          ? cumulative + (ex.end - 1 - genomicPos)
          : cumulative + (genomicPos - ex.start);
      foundExonNum = i + 1;
    }
    cumulative += exLen;
  }

  if (foundOffset === null) return null;
  return {
    cdsOffset: foundOffset,
    exonNumber: foundExonNum,
    totalCdsLen: cumulative,
  };
}

/**
 * 정렬된 엑손 목록에서 특정 CDS 오프셋(0-based)에 해당하는 게놈 좌표(0-based)를 반환합니다.
 */
function genomicPosForCdsOffset(
  cdsOffset: number,
  orderedExons: ExonInterval[],
  strand: '+' | '-' | '.',
): number | null {
  if (cdsOffset < 0) return null;
  let rem = cdsOffset;
  for (const ex of orderedExons) {
    const exLen = Math.max(0, ex.end - ex.start);
    if (rem < exLen) {
      return strand === '-' ? ex.end - 1 - rem : ex.start + rem;
    }
    rem -= exLen;
  }
  return null;
}

function isNearExonBoundary(pos: number, exons: ExonInterval[], maxDist = 2): boolean {
  for (const ex of exons) {
    if (pos >= ex.start - maxDist && pos < ex.start) return true;
    if (pos >= ex.end && pos < ex.end + maxDist) return true;
  }
  return false;
}

function classifySnpInCodon(
  refCodon: string,
  codonPos: number,
  codingAltBase: string,
  aaPos1: number,
  cdsPos1: number,
  refCodingBase: string,
  tableId: TranslationTableId,
): {
  consequence: VariantConsequence;
  impact: VariantImpact;
  hgvsP: string;
  hgvsC: string;
  codonChange: string;
  isHeterozygous: boolean;
} {
  const refAA = translateCodon(refCodon, tableId);
  const candidates = expandIupacBase(codingAltBase);
  const isHet = candidates.length > 1;
  const nonRefCandidates = candidates.filter((b) => b !== refCodingBase && b !== 'N');
  const effectiveAlts = nonRefCandidates.length > 0 ? nonRefCandidates : candidates;

  let bestImpact: VariantImpact = 'LOW';
  let bestConsequence: VariantConsequence = 'synonymous';
  let bestAltCodon = refCodon;
  const aaOutcomes: string[] = [];

  for (const altB of effectiveAlts) {
    if (altB === 'N') continue;
    const altCodon =
      refCodon.slice(0, codonPos) + altB + refCodon.slice(codonPos + 1);
    const altAA = translateCodon(altCodon, tableId);
    aaOutcomes.push(altAA);

    let cons: VariantConsequence = 'synonymous';
    let imp: VariantImpact = 'LOW';

    if (refAA === 'M' && aaPos1 === 1 && altAA !== 'M') {
      cons = 'start_lost';
      imp = 'HIGH';
    } else if (refAA !== '*' && altAA === '*') {
      cons = 'nonsense';
      imp = 'HIGH';
    } else if (refAA === '*' && altAA !== '*') {
      cons = 'stop_lost';
      imp = 'HIGH';
    } else if (altAA !== refAA) {
      cons = 'missense';
      imp = 'MODERATE';
    }

    if (IMPACT_ORDER[imp] > IMPACT_ORDER[bestImpact]) {
      bestImpact = imp;
      bestConsequence = cons;
      bestAltCodon = altCodon;
    } else if (bestAltCodon === refCodon) {
      bestAltCodon = altCodon;
    }
  }

  const primaryAltAA =
    aaOutcomes.find((a) => a === '*') ??
    aaOutcomes.find((a) => a !== refAA) ??
    aaOutcomes[0] ??
    refAA;
  const hetSuffix = isHet && bestConsequence !== 'synonymous' ? `/${refAA}` : '';
  const hgvsP =
    bestConsequence === 'synonymous'
      ? `p.${refAA}${aaPos1}=`
      : bestConsequence === 'stop_lost'
        ? `p.*${aaPos1}${primaryAltAA}ext${hetSuffix}`
        : `p.${refAA}${aaPos1}${primaryAltAA}${hetSuffix}`;

  const hgvsC = `c.${cdsPos1}${refCodingBase}>${codingAltBase}`;
  return {
    consequence: bestConsequence,
    impact: bestImpact,
    hgvsP,
    hgvsC,
    codonChange: `${refCodon}>${bestAltCodon}`,
    isHeterozygous: isHet,
  };
}

function featurePriorityRank(feat: GffFeature, vStart: number, vEnd: number): number {
  const exons = getOrderedCodingExons(feat);
  const inExon = exons.some((ex) => vStart < ex.end && vEnd > ex.start);
  if (feat.type === 'gene' && inExon) return 100;
  if (feat.type === 'gene') return 80;
  if (
    feat.type === 'tRNA_gene' ||
    feat.type === 'rRNA_gene' ||
    feat.type === 'ncRNA_gene' ||
    feat.type === 'snRNA_gene' ||
    feat.type === 'snoRNA_gene'
  ) {
    return 60;
  }
  if (feat.type === 'pseudogene' || feat.type === 'transposable_element_gene') return 40;
  if (feat.type === 'dubious_gene') return inExon ? 30 : 20;
  return 10;
}

function computeCdsOverlapBp(vStart: number, vEnd: number, exons: ExonInterval[]): number {
  let total = 0;
  for (const ex of exons) {
    const ov = Math.min(vEnd, ex.end) - Math.max(vStart, ex.start);
    if (ov > 0) total += ov;
  }
  return total;
}

/**
 * 단일 VariantRecord에 대해 유전자·코돈·아미노산 효과를 정밀 계산합니다.
 */
export function annotateSingleVariant(
  v: VariantRecord,
  nearbyFeatures: GffFeature[],
  getRefBaseAt: (pos: number) => string,
  tableId: TranslationTableId = isMitochondrialChrom(v.chrom) ? 3 : 1,
): VariantAnnotationMeta {
  const vcfAnchorBase =
    (v.type === 'DEL' || v.type === 'GAP') && v.start > 0
      ? getRefBaseAt(v.start - 1).toUpperCase() || 'N'
      : undefined;

  // 1. 유전자 본체와 직접 겹치는 Feature 탐색
  const overlapping = nearbyFeatures
    .filter((f) => f.chrom === v.chrom && v.start < f.end && v.end > f.start)
    .sort((a, b) => featurePriorityRank(b, v.start, v.end) - featurePriorityRank(a, v.start, v.end));

  const primary = overlapping[0];
  if (!primary) {
    // 상류 500 bp 프로모터 검사
    const promoters = nearbyFeatures
      .filter((f) => {
        if (f.chrom !== v.chrom || f.type === 'dubious_gene') return false;
        if (f.strand === '-') {
          return v.start >= f.end && v.start < f.end + 500;
        }
        return v.end <= f.start && v.end > f.start - 500;
      })
      .sort((a, b) => {
        const dA = a.strand === '-' ? v.start - a.end : a.start - v.end;
        const dB = b.strand === '-' ? v.start - b.end : b.start - v.end;
        return dA - dB;
      });

    const prom = promoters[0];
    if (prom) {
      const dist = prom.strand === '-' ? v.start - prom.end + 1 : prom.start - v.start;
      return {
        geneId: prom.id,
        geneName: prom.name,
        systematicName: prom.systematicName,
        strand: prom.strand,
        featureType: prom.type,
        consequence: 'promoter_500bp',
        impact: 'MODIFIER',
        hgvsC: `c.-${Math.max(1, dist)}`,
        vcfAnchorBase,
      };
    }

    return {
      geneName: 'Intergenic',
      featureType: 'intergenic',
      consequence: 'intergenic',
      impact: 'MODIFIER',
      vcfAnchorBase,
    };
  }

  // 2. 비코딩 RNA 유전자 등 non-CDS 특성 처리
  const isCodingFeature =
    primary.type === 'gene' ||
    primary.type === 'dubious_gene' ||
    primary.type === 'transposable_element_gene' ||
    primary.subFeatures.some((sf) => sf.type === 'CDS');

  if (!isCodingFeature) {
    return {
      geneId: primary.id,
      geneName: primary.name,
      systematicName: primary.systematicName,
      strand: primary.strand,
      featureType: primary.type,
      consequence: 'ncrna_or_trna',
      impact: 'MODIFIER',
      vcfAnchorBase,
    };
  }

  const exons = getOrderedCodingExons(primary);
  const cdsOvBp = computeCdsOverlapBp(v.start, v.end, exons);

  // 3. GAP (커버리지 결실/미커버)
  if (v.type === 'GAP') {
    if (cdsOvBp === 0) {
      const nearSplice = isNearExonBoundary(v.start, exons) || isNearExonBoundary(v.end - 1, exons);
      return {
        geneId: primary.id,
        geneName: primary.name,
        systematicName: primary.systematicName,
        strand: primary.strand,
        featureType: primary.type,
        consequence: nearSplice ? 'splice_site' : 'intron',
        impact: nearSplice ? 'HIGH' : 'MODIFIER',
        vcfAnchorBase,
      };
    }
    const isMajorLoss = cdsOvBp >= 10;
    return {
      geneId: primary.id,
      geneName: primary.name,
      systematicName: primary.systematicName,
      strand: primary.strand,
      featureType: primary.type,
      consequence: isMajorLoss ? 'gene_loss' : 'low_depth_gap',
      impact: isMajorLoss ? 'HIGH' : 'MODIFIER',
      hgvsP: isMajorLoss ? `결실/미커버 (${cdsOvBp} bp)` : `저심도 (${cdsOvBp} bp)`,
      vcfAnchorBase,
    };
  }

  // 4. CDS 엑손 내부가 아닐 때 -> Splice site vs Intron
  const probePos =
    v.type === 'INS' && primary.strand === '-' && locatePositionInExons(v.start, exons, primary.strand) === null
      ? v.start + 1
      : v.start;
  const exonHit = locatePositionInExons(probePos, exons, primary.strand);

  if (!exonHit && cdsOvBp === 0) {
    const nearSplice = isNearExonBoundary(v.start, exons) || isNearExonBoundary(v.end - 1, exons);
    return {
      geneId: primary.id,
      geneName: primary.name,
      systematicName: primary.systematicName,
      strand: primary.strand,
      featureType: primary.type,
      consequence: nearSplice ? 'splice_site' : 'intron',
      impact: nearSplice ? 'HIGH' : 'MODIFIER',
      hgvsP: nearSplice ? 'splice±2bp' : undefined,
      vcfAnchorBase,
    };
  }

  // 5. 코딩 서열 내 변이 평가
  const targets = getTargetAllelesToEvaluate(v);

  // 5-A. 단일 염기 SNP 정밀 코돈 번역
  if (v.type === 'SNP' && v.end - v.start === 1 && exonHit && targets.length > 0) {
    const codonIdx = Math.floor(exonHit.cdsOffset / 3);
    const codonStartOffset = codonIdx * 3;
    const codonPos = exonHit.cdsOffset % 3;
    const g0 = genomicPosForCdsOffset(codonStartOffset, exons, primary.strand);
    const g1 = genomicPosForCdsOffset(codonStartOffset + 1, exons, primary.strand);
    const g2 = genomicPosForCdsOffset(codonStartOffset + 2, exons, primary.strand);

    if (g0 !== null && g1 !== null && g2 !== null) {
      const orient = (b: string) => (primary.strand === '-' ? reverseComplement(b) : b.toUpperCase());
      const refCodon = orient(getRefBaseAt(g0)) + orient(getRefBaseAt(g1)) + orient(getRefBaseAt(g2));
      const refCodingBase = orient(v.refAllele[0] ?? getRefBaseAt(v.start));
      const aaPos1 = codonIdx + 1;
      const cdsPos1 = exonHit.cdsOffset + 1;

      if (/^[ACGT]{3}$/.test(refCodon)) {
        let bestMeta: VariantAnnotationMeta | null = null;
        for (const t of targets) {
          const altGenomicBase = t.allele.replace(/-/g, '')[0] ?? 'N';
          const codingAlt = orient(altGenomicBase);
          const res = classifySnpInCodon(
            refCodon,
            codonPos,
            codingAlt,
            aaPos1,
            cdsPos1,
            refCodingBase,
            tableId,
          );
          if (!bestMeta || IMPACT_ORDER[res.impact] > IMPACT_ORDER[bestMeta.impact]) {
            bestMeta = {
              geneId: primary.id,
              geneName: primary.name,
              systematicName: primary.systematicName,
              strand: primary.strand,
              featureType: primary.type,
              consequence: res.consequence,
              impact: res.impact,
              hgvsP: res.hgvsP,
              hgvsC: res.hgvsC,
              codonChange: res.codonChange,
              exonNumber: exonHit.exonNumber,
              totalExons: exons.length,
              isHeterozygous: res.isHeterozygous,
              vcfAnchorBase,
            };
          }
        }
        if (bestMeta) return bestMeta;
      }
    }
  }

  // 5-B. Insertion / Deletion / Complex indels
  const refOffsetForCds =
    exonHit ??
    locatePositionInExons(
      Math.max(v.start, exons[0]?.start ?? v.start),
      exons,
      primary.strand,
    );
  const cdsPos1 = refOffsetForCds ? refOffsetForCds.cdsOffset + 1 : 1;
  const aaPos1 = refOffsetForCds ? Math.floor(refOffsetForCds.cdsOffset / 3) + 1 : 1;
  const refG0 = refOffsetForCds
    ? genomicPosForCdsOffset(Math.floor(refOffsetForCds.cdsOffset / 3) * 3, exons, primary.strand)
    : null;
  const refG1 = refOffsetForCds
    ? genomicPosForCdsOffset(Math.floor(refOffsetForCds.cdsOffset / 3) * 3 + 1, exons, primary.strand)
    : null;
  const refG2 = refOffsetForCds
    ? genomicPosForCdsOffset(Math.floor(refOffsetForCds.cdsOffset / 3) * 3 + 2, exons, primary.strand)
    : null;
  let anchorAA = '';
  if (refG0 !== null && refG1 !== null && refG2 !== null) {
    const orient = (b: string) => (primary.strand === '-' ? reverseComplement(b) : b.toUpperCase());
    anchorAA = translateCodon(
      orient(getRefBaseAt(refG0)) + orient(getRefBaseAt(refG1)) + orient(getRefBaseAt(refG2)),
      tableId,
    );
  }

  let bestIndelImpact: VariantImpact = 'MODERATE';
  let bestIndelCons: VariantConsequence = 'inframe_indel';
  let bestHgvsP = `p.${anchorAA || 'X'}${aaPos1}indel`;
  let bestHgvsC = `c.${cdsPos1}`;

  for (const t of targets) {
    const ungappedAltLen = t.allele.replace(/-/g, '').length;
    const refLen = v.refAllele.length;
    const netDelta = ungappedAltLen - refLen;

    if (v.type === 'INS' || netDelta > 0) {
      const insLen = Math.max(1, netDelta);
      const isFs = insLen % 3 !== 0;
      const imp: VariantImpact = isFs ? 'HIGH' : 'MODERATE';
      const cons: VariantConsequence = isFs ? 'frameshift' : 'inframe_indel';
      if (IMPACT_ORDER[imp] >= IMPACT_ORDER[bestIndelImpact]) {
        bestIndelImpact = imp;
        bestIndelCons = cons;
        bestHgvsP = isFs
          ? `p.${anchorAA || 'X'}${aaPos1}fs (+${insLen}bp)`
          : `p.${anchorAA || 'X'}${aaPos1}ins${insLen / 3}aa`;
        bestHgvsC = `c.${cdsPos1}ins+${insLen}`;
      }
    } else if (v.type === 'DEL' || netDelta < 0) {
      const delLen = Math.max(1, -netDelta);
      const effectiveDelCds = exonHit && cdsOvBp > 0 ? cdsOvBp : delLen;
      const isFs = effectiveDelCds % 3 !== 0;
      const imp: VariantImpact = isFs ? 'HIGH' : 'MODERATE';
      const cons: VariantConsequence = isFs ? 'frameshift' : 'inframe_indel';
      if (IMPACT_ORDER[imp] >= IMPACT_ORDER[bestIndelImpact]) {
        bestIndelImpact = imp;
        bestIndelCons = cons;
        bestHgvsP = isFs
          ? `p.${anchorAA || 'X'}${aaPos1}fs (-${effectiveDelCds}bp)`
          : `p.${anchorAA || 'X'}${aaPos1}del${effectiveDelCds / 3}aa`;
        bestHgvsC = `c.${cdsPos1}del${effectiveDelCds}`;
      }
    } else {
      bestIndelImpact = 'MODERATE';
      bestIndelCons = 'missense';
      bestHgvsP = `p.${anchorAA || 'X'}${aaPos1}delins`;
      bestHgvsC = `c.${cdsPos1}delins`;
    }
  }

  return {
    geneId: primary.id,
    geneName: primary.name,
    systematicName: primary.systematicName,
    strand: primary.strand,
    featureType: primary.type,
    consequence: bestIndelCons,
    impact: bestIndelImpact,
    hgvsP: bestHgvsP,
    hgvsC: bestHgvsC,
    exonNumber: refOffsetForCds?.exonNumber,
    totalExons: exons.length,
    vcfAnchorBase,
  };
}

/**
 * 변이 목록 전체에 대해 주석을 계산하고 annotation이 부착된 새 배열을 반환합니다.
 * getRefBaseAt는 구간 외부 코돈 경계(±2 bp) 또는 VCF 앵커 염기(-1 bp) 조회에도 사용됩니다.
 */
export function annotateVariants(
  variants: VariantRecord[],
  annotationIndex: AnnotationIndex | null,
  getRefBaseAt: (chrom: string, pos: number) => string,
): VariantRecord[] {
  if (variants.length === 0) return variants;
  return variants.map((v) => {
    const features = annotationIndex
      ? annotationIndex.getFeaturesInRange(v.chrom, Math.max(0, v.start - 550), v.end + 550)
      : [];
    const annotation = annotateSingleVariant(
      v,
      features,
      (pos) => getRefBaseAt(v.chrom, pos),
      isMitochondrialChrom(v.chrom) ? 3 : 1,
    );
    return { ...v, annotation };
  });
}

/**
 * 주석된 변이를 유전자 단위로 집계해 영향도가 큰 유전자부터 정렬합니다.
 */
export function summarizeVariantsByGene(variants: VariantRecord[]): GeneVariantSummary[] {
  const map = new Map<string, GeneVariantSummary>();
  for (const v of variants) {
    const ann = v.annotation;
    if (!ann || ann.featureType === 'intergenic') continue;
    const key = ann.geneId || `${v.chrom}:${ann.geneName}`;
    let entry = map.get(key);
    if (!entry) {
      entry = {
        geneKey: key,
        geneName: ann.geneName,
        systematicName: ann.systematicName,
        chrom: v.chrom,
        start: v.start,
        end: v.end,
        strand: ann.strand ?? '.',
        featureType: ann.featureType,
        highestImpact: ann.impact,
        variantCount: 0,
        highCount: 0,
        moderateCount: 0,
        lowCount: 0,
        modifierCount: 0,
        topExamples: [],
        categories: [],
      };
      map.set(key, entry);
    }

    entry.start = Math.min(entry.start, v.start);
    entry.end = Math.max(entry.end, v.end);
    entry.variantCount++;
    if (ann.impact === 'HIGH') entry.highCount++;
    else if (ann.impact === 'MODERATE') entry.moderateCount++;
    else if (ann.impact === 'LOW') entry.lowCount++;
    else entry.modifierCount++;

    if (IMPACT_ORDER[ann.impact] > IMPACT_ORDER[entry.highestImpact]) {
      entry.highestImpact = ann.impact;
    }

    const exStr = ann.hgvsP || ann.hgvsC || `${v.type}@${v.start + 1}`;
    if (entry.topExamples.length < 4 && !entry.topExamples.includes(exStr)) {
      entry.topExamples.push(exStr);
    }
    if (!entry.categories.includes(v.category)) {
      entry.categories.push(v.category);
    }
  }

  return [...map.values()].sort(
    (a, b) =>
      IMPACT_ORDER[b.highestImpact] - IMPACT_ORDER[a.highestImpact] ||
      b.highCount - a.highCount ||
      b.moderateCount - a.moderateCount ||
      b.variantCount - a.variantCount ||
      a.geneName.localeCompare(b.geneName),
  );
}

function resolveVcfRefAltGt(
  v: VariantRecord,
  rawAllele: string,
): { pos1: number; ref: string; alt: string; gt: string } | null {
  if (!rawAllele || v.type === 'GAP') return null;

  if (v.type === 'INS') {
    const cleanAlt = rawAllele.toUpperCase().replace(/[^ACGTN]/g, '');
    const cleanRef = (v.refAllele[0] ?? 'N').toUpperCase();
    if (!cleanAlt || cleanAlt === cleanRef) return null;
    return { pos1: v.start + 1, ref: cleanRef, alt: cleanAlt, gt: '1/1' };
  }

  if (v.type === 'DEL' || rawAllele.includes('-')) {
    const anchor = v.annotation?.vcfAnchorBase ?? 'N';
    const altRemaining = rawAllele.replace(/-/g, '');
    if (v.start > 0) {
      return {
        pos1: v.start, // 1-based anchor position = (v.start - 1) + 1
        ref: anchor + v.refAllele.toUpperCase(),
        alt: anchor + altRemaining.toUpperCase(),
        gt: '1/1',
      };
    }
    return {
      pos1: 1,
      ref: v.refAllele.toUpperCase(),
      alt: altRemaining.toUpperCase() || '.',
      gt: '1/1',
    };
  }

  // SNP or multi-base substitution
  const refUpper = v.refAllele.toUpperCase();
  const altUpper = rawAllele.toUpperCase();
  if (refUpper.length === 1 && altUpper.length === 1) {
    const expanded = expandIupacBase(altUpper).filter((b) => b !== refUpper && b !== 'N');
    if (expanded.length === 0) return null;
    const isHet = expandIupacBase(altUpper).length > 1;
    return {
      pos1: v.start + 1,
      ref: refUpper,
      alt: expanded[0],
      gt: isHet ? '0/1' : '1/1',
    };
  }
  return { pos1: v.start + 1, ref: refUpper, alt: altUpper, gt: '1/1' };
}

/**
 * 표준 VCF v4.2 포맷 문자열을 생성합니다.
 */
export function exportVariantsToVcf(
  variants: VariantRecord[],
  refContigs: ContigInfo[],
  s1Label: string,
  s2Label: string | null,
): string {
  const lines: string[] = [
    '##fileformat=VCFv4.2',
    '##source=YeastWGS_Studio',
    ...refContigs.map((c) => `##contig=<ID=${c.name},length=${c.length}>`),
    '##INFO=<ID=VTYPE,Number=1,Type=String,Description="Variant type">',
    '##INFO=<ID=CAT,Number=1,Type=String,Description="Comparison category">',
    '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotation: Gene|Consequence|Impact|HGVS_c|HGVS_p">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
  ];

  const sanitizeSample = (s: string) => s.replace(/\s+/g, '_');
  const sampleCols = s2Label
    ? [sanitizeSample(s1Label), sanitizeSample(s2Label)]
    : [sanitizeSample(s1Label)];
  lines.push(`#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sampleCols.join('\t')}`);

  for (const v of variants) {
    if (v.type === 'GAP') continue;
    const primaryAllele =
      v.category === 'strain2_only' ? v.s2Allele : v.s1Allele || v.s2Allele;
    const rec = resolveVcfRefAltGt(v, primaryAllele);
    if (!rec || rec.alt === '.') continue;

    // 2균주 비교 시 ALTs 병합 처리
    let altCol = rec.alt;
    let s1Gt = v.category === 'strain2_only' ? '0/0' : rec.gt;
    let s2Gt = '0/0';
    if (s2Label) {
      if (v.category === 'shared') {
        s2Gt = rec.gt;
      } else if (v.category === 'strain2_only') {
        s2Gt = rec.gt;
      } else if (v.category === 'divergent') {
        const rec2 = resolveVcfRefAltGt(v, v.s2Allele);
        if (rec2 && rec2.alt !== '.' && rec2.alt !== rec.alt) {
          altCol = `${rec.alt},${rec2.alt}`;
          s2Gt = rec2.gt === '0/1' ? '0/2' : '2/2';
        } else if (rec2) {
          s2Gt = rec2.gt;
        }
      }
    }

    const ann = v.annotation;
    const annField = ann
      ? `${ann.geneName}|${ann.consequence}|${ann.impact}|${ann.hgvsC ?? ''}|${ann.hgvsP ?? ''}`
      : 'Intergenic|intergenic|MODIFIER||';
    const info = `VTYPE=${v.type};CAT=${v.category};ANN=${annField.replace(/[\s;]/g, '_')}`;
    const sampleData = s2Label
      ? `${s1Gt}:${v.s1Depth}\t${s2Gt}:${v.s2Depth}`
      : `${s1Gt}:${v.s1Depth}`;
    lines.push(`${v.chrom}\t${rec.pos1}\t.\t${rec.ref}\t${altCol}\t.\tPASS\t${info}\tGT:DP\t${sampleData}`);
  }

  return lines.join('\n') + '\n';
}

/**
 * Excel 호환 UTF-8 BOM TSV 문자열을 생성합니다.
 */
export function exportVariantsToTsv(
  variants: VariantRecord[],
  s1Label: string,
  s2Label: string | null,
): string {
  const header = [
    'Chrom',
    'Start_1based',
    'End_1based',
    'Type',
    'Category',
    'Gene',
    'Systematic_ID',
    'Strand',
    'Impact',
    'Consequence',
    'HGVS_c',
    'HGVS_p',
    'Codon_Change',
    'Ref',
    s1Label,
    `${s1Label}_Depth`,
    ...(s2Label ? [s2Label, `${s2Label}_Depth`] : []),
  ].join('\t');

  const rows = variants.map((v) => {
    const a = v.annotation;
    return [
      v.chrom,
      v.start + 1,
      v.end,
      v.type,
      v.category,
      a?.geneName ?? '-',
      a?.systematicName ?? '-',
      a?.strand ?? '-',
      a?.impact ?? 'MODIFIER',
      a?.consequence ?? '-',
      a?.hgvsC ?? '-',
      a?.hgvsP ?? '-',
      a?.codonChange ?? '-',
      v.refAllele,
      v.s1Allele,
      v.s1Depth,
      ...(s2Label ? [v.s2Allele, v.s2Depth] : []),
    ].join('\t');
  });

  return '\uFEFF' + [header, ...rows].join('\n') + '\n';
}

