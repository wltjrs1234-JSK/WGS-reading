import { normalizeYeastChromName, resolveMatchingRefName } from './yeastAliases';
import type { GffFeature } from './types';

function decodeGffAttr(value: string): string {
  try {
    return decodeURIComponent(value);
  } catch {
    return value;
  }
}

function parseAttributes(attrCol: string): Record<string, string> {
  const attrs: Record<string, string> = {};
  for (const part of attrCol.split(';')) {
    const trimmed = part.trim();
    if (!trimmed) continue;
    const eqIdx = trimmed.indexOf('=');
    if (eqIdx === -1) continue;
    const key = trimmed.slice(0, eqIdx).trim();
    const val = decodeGffAttr(trimmed.slice(eqIdx + 1).trim());
    attrs[key] = val;
  }
  return attrs;
}

const GENE_PARENT_TYPES = new Set([
  'gene',
  'pseudogene',
  'ncRNA_gene',
  'tRNA_gene',
  'rRNA_gene',
  'snoRNA_gene',
  'snRNA_gene',
  'transposable_element_gene',
  'ARS',
  'centromere',
  'telomere',
]);

const SUBFEATURE_TYPES = new Set(['CDS', 'exon', 'noncoding_exon', 'five_prime_UTR_intron', 'intron']);

export async function parseGff3File(file: File): Promise<GffFeature[]> {
  let text: string;
  if (file.name.toLowerCase().endsWith('.gz')) {
    const stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
    text = await new Response(stream).text();
  } else {
    text = await file.text();
  }
  return parseGff3Text(text);
}

export function parseGff3Text(text: string): GffFeature[] {
  const featuresById = new Map<string, GffFeature>();
  const mrnaToGeneId = new Map<string, string>();
  const topFeatures: GffFeature[] = [];
  const pendingSubs: Array<{
    parentIds: string[];
    type: string;
    start: number;
    end: number;
    isDubious: boolean;
  }> = [];

  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line === '##FASTA') break;
    if (line.startsWith('#')) continue;

    const cols = line.split('\t');
    if (cols.length < 9) continue;

    const chrom = cols[0];
    const featureType = cols[2];
    const start1 = Number(cols[3]);
    const end1 = Number(cols[4]);
    const strand = (cols[6] === '+' || cols[6] === '-' ? cols[6] : '.') as '+' | '-' | '.';
    const attrs = parseAttributes(cols[8]);

    if (Number.isNaN(start1) || Number.isNaN(end1)) continue;
    const start = Math.max(0, start1 - 1); // 0-based
    const end = Math.max(start, end1);

    const id = attrs.ID || `${chrom}:${start}-${end}:${featureType}`;

    if (
      featureType === 'mRNA' ||
      featureType === 'transcript' ||
      featureType === 'tRNA' ||
      featureType === 'rRNA' ||
      featureType === 'snoRNA' ||
      featureType === 'snRNA' ||
      featureType === 'ncRNA'
    ) {
      if (attrs.Parent) {
        mrnaToGeneId.set(id, attrs.Parent.split(',')[0]);
      }
      continue;
    }

    if (GENE_PARENT_TYPES.has(featureType)) {
      if (featureType === 'ARS' && id.endsWith('_ARS')) continue;
      const systematic = attrs.Name || attrs.ID || id;
      const standard = attrs.gene || systematic;
      const aliases = (attrs.Alias ?? '')
        .split(',')
        .map((s) => s.trim())
        .filter((s) => s && s !== standard && s !== systematic && s.length <= 24);
      const baseNote = attrs.Note || attrs.display || '';
      const note =
        aliases.length > 0 ? `[${aliases.slice(0, 3).join(', ')}] ${baseNote}`.trim() : baseNote || undefined;

      const feat: GffFeature = {
        id,
        name: standard,
        systematicName: systematic !== standard ? systematic : undefined,
        chrom,
        type: featureType,
        start,
        end,
        strand,
        note,
        subFeatures: [],
      };
      featuresById.set(id, feat);
      topFeatures.push(feat);
    } else if (SUBFEATURE_TYPES.has(featureType) && attrs.Parent) {
      pendingSubs.push({
        parentIds: attrs.Parent.split(',').map((s) => s.trim()),
        type: featureType,
        start,
        end,
        isDubious: attrs.orf_classification === 'Dubious',
      });
    }
  }

  for (const sub of pendingSubs) {
    for (const pid of sub.parentIds) {
      const directParent = featuresById.get(pid);
      const geneParent = directParent ?? featuresById.get(mrnaToGeneId.get(pid) ?? '');
      if (geneParent) {
        if (sub.isDubious && geneParent.type === 'gene') {
          geneParent.type = 'dubious_gene';
        }
        const exists = geneParent.subFeatures.some(
          (f) => f.type === sub.type && f.start === sub.start && f.end === sub.end,
        );
        if (!exists) {
          geneParent.subFeatures.push({
            type: sub.type,
            start: sub.start,
            end: sub.end,
          });
        }
      }
    }
  }

  // SGD GFF3의 5'/3' UTR 포함 gene 전사체 범위 대신 실제 코딩 CDS 범위로 보정 (IGV RefSeq 일치)
  for (const feat of topFeatures) {
    feat.subFeatures.sort((a, b) => a.start - b.start);
    const cdsSubs = feat.subFeatures.filter((s) => s.type === 'CDS');
    if (cdsSubs.length > 0) {
      feat.start = Math.min(...cdsSubs.map((c) => c.start));
      feat.end = Math.max(...cdsSubs.map((c) => c.end));
      feat.subFeatures = cdsSubs;
    }
  }

  topFeatures.sort((a, b) => a.start - b.start);
  return topFeatures;
}

export class AnnotationIndex {
  private features: GffFeature[];
  private chromList: string[];

  constructor(features: GffFeature[]) {
    this.features = features;
    this.chromList = [...new Set(features.map((f) => f.chrom))];
  }

  get totalCount(): number {
    return this.features.length;
  }

  getFeaturesInRange(targetChrom: string, start: number, end: number): GffFeature[] {
    const matchedChrom = resolveMatchingRefName(targetChrom, this.chromList);
    if (!matchedChrom) {
      const targetNorm = normalizeYeastChromName(targetChrom);
      return this.features.filter(
        (f) =>
          normalizeYeastChromName(f.chrom) === targetNorm && f.end > start && f.start < end,
      );
    }
    return this.features.filter(
      (f) => f.chrom === matchedChrom && f.end > start && f.start < end,
    );
  }

  searchGenes(query: string, limit = 25): GffFeature[] {
    const q = query.trim().toLowerCase();
    if (!q) return [];
    const exact: GffFeature[] = [];
    const prefix: GffFeature[] = [];
    const substring: GffFeature[] = [];

    for (const f of this.features) {
      const nameL = f.name.toLowerCase();
      const sysL = f.systematicName?.toLowerCase() ?? '';
      const noteL = f.note?.toLowerCase() ?? '';

      if (nameL === q || sysL === q) {
        exact.push(f);
      } else if (nameL.startsWith(q) || sysL.startsWith(q)) {
        prefix.push(f);
      } else if (nameL.includes(q) || sysL.includes(q) || noteL.includes(q)) {
        substring.push(f);
      }
    }

    return [...exact, ...prefix, ...substring].slice(0, limit);
  }
}

