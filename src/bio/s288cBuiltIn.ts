import { AnnotationIndex } from './gff3Parser';
import rawGenes from './s288cGenes.json';
import type { GffFeature } from './types';

type CompactGeneTuple = [
  string,          // chrom (e.g. 'chrI', 'chrM')
  number,          // 0-based start
  number,          // 0-based end
  '+' | '-' | '.', // strand
  string,          // standard gene name (or systematic if none)
  string,          // systematicName ('' if same as name)
  string,          // type ('gene', 'tRNA_gene', ...)
  string,          // short note
  [number, number][]?, // optional exon/CDS intervals
];

export function createBuiltInS288CAnnotationIndex(): AnnotationIndex {
  const tuples = rawGenes as unknown as CompactGeneTuple[];
  const features: GffFeature[] = tuples.map((t, idx) => {
    const [chrom, start, end, strand, name, systematic, type, note, subs] = t;
    return {
      id: systematic || `${name}_${idx}`,
      name,
      systematicName: systematic || undefined,
      chrom,
      type,
      start,
      end,
      strand,
      note: note || undefined,
      subFeatures: subs
        ? subs.map(([s, e]) => ({ type: 'CDS', start: s, end: e }))
        : [],
    };
  });
  return new AnnotationIndex(features);
}

