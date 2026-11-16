export type CigarOperator = 'M' | 'I' | 'D' | 'N' | 'S' | 'H' | 'P' | '=' | 'X';

export interface CigarOp {
  op: CigarOperator;
  len: number;
}

export interface AlignedRead {
  id: string;
  name: string;
  start: number; // 0-based inclusive
  end: number;   // 0-based exclusive
  strand: '+' | '-';
  seq: string;
  cigar: string;
  mapq: number;
  qual?: Uint8Array | number[];
  flags: number;
}

export interface SoftClipInfo {
  left: number;
  right: number;
}

export interface BaseCounts {
  A: number;
  C: number;
  G: number;
  T: number;
  N: number;
  del: number;
}

export interface PositionPileup {
  pos: number; // 0-based reference position
  depth: number;
  counts: BaseCounts;
  insertions: Record<string, number>;
}

export type PloidyMode = 'haploid' | 'diploid';

export interface ConsensusOptions {
  minDepth: number;
  minMapQ: number;
  minBaseQ: number;
  ploidyMode: PloidyMode;
  hetThreshold: number;   // e.g. 0.25 for diploid IUPAC call
  indelThreshold: number; // e.g. 0.50 for calling deletion or insertion
}

export const DEFAULT_CONSENSUS_OPTIONS: ConsensusOptions = {
  minDepth: 5,
  minMapQ: 0,
  minBaseQ: 0,
  ploidyMode: 'haploid',
  hetThreshold: 0.25,
  indelThreshold: 0.5,
};

export interface ConsensusPositionCall {
  pos: number;
  refBase: string;
  call: string;          // 'A'|'C'|'G'|'T'|'-'|'N'|IUPAC
  insertion: string;     // consensus insertion immediately after pos, or ''
  depth: number;
  isLowDepth: boolean;
  isMismatch: boolean;
  isHeterozygous: boolean;
  counts: BaseCounts;
}

export type DiffCategory =
  | 'match'
  | 'strain1_only'
  | 'strain2_only'
  | 'shared'
  | 'divergent'
  | 'gap_s1'
  | 'gap_s2'
  | 'gap_both';

export interface VariantRecord {
  chrom: string;
  start: number; // 0-based
  end: number;   // 0-based exclusive
  type: 'SNP' | 'INS' | 'DEL' | 'GAP' | 'COMPLEX';
  category: DiffCategory;
  refAllele: string;
  s1Allele: string;
  s2Allele: string;
  s1Depth: number;
  s2Depth: number;
}

export interface AlignedColumnMeta {
  refPos: number; // 0-based coordinate (for insertion column, anchor refPos)
  isInsertion: boolean;
  category: DiffCategory;
}

export interface ThreeWayRegionAnalysis {
  chrom: string;
  start: number;
  end: number;
  hasStrain2: boolean;
  refSeq: string;
  s1Calls: ConsensusPositionCall[];
  s2Calls: ConsensusPositionCall[];
  ungappedRef: string;
  ungappedS1: string;
  ungappedS2: string;
  alignedRef: string;
  alignedS1: string;
  alignedS2: string;
  columns: AlignedColumnMeta[];
  variants: VariantRecord[];
}

export interface FastaIndexEntry {
  name: string;
  length: number;
  offset: number;    // byte offset of first sequence character
  lineBases: number; // bases per line
  lineWidth: number; // bytes per line including newline
}

export interface GffSubFeature {
  type: string;
  start: number; // 0-based inclusive
  end: number;   // 0-based exclusive
}

export interface GffFeature {
  id: string;
  name: string;
  systematicName?: string;
  chrom: string;
  type: string;
  start: number; // 0-based inclusive
  end: number;   // 0-based exclusive
  strand: '+' | '-' | '.';
  note?: string;
  subFeatures: GffSubFeature[];
}

export interface ContigInfo {
  name: string;
  length: number;
}

