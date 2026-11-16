import type { RegionMappingStats } from '../bio/consensus';
import type { AlleleTick, PackedRead } from '../bio/readPacking';
import type {
  AlignedRead,
  ConsensusOptions,
  ContigInfo,
  GffFeature,
  ThreeWayRegionAnalysis,
} from '../bio/types';

export interface StrainFileInput {
  id: string;
  label: string;
  bamFile: File;
  indexFile: File;
}

export interface LoadLocalFilesPayload {
  fastaFile: File;
  faiFile?: File;
  gffFile?: File;
  strains: StrainFileInput[];
}

export interface TrackSummaryMeta {
  id: string;
  label: string;
  contigs: ContigInfo[];
  matchedChrom: string | null;
  mappings: Record<string, string | null>;
}

export type AnnotationOrigin = 'user_gff' | 'builtin_s288c' | 'demo' | 'none';

export interface WorkspaceLoadedSummary {
  genomeId: string; // 'sacCer3' | 'custom' | 'demo'
  referenceLabel: string;
  contigs: ContigInfo[];
  gffCount: number;
  annotationOrigin: AnnotationOrigin;
  strains: TrackSummaryMeta[];
  defaultChrom: string;
  defaultStart: number;
  defaultEnd: number;
  isDemo: boolean;
}

export interface StrainWindowData {
  id: string;
  label: string;
  matchedChrom: string | null;
  reads: AlignedRead[];
  packedReads?: PackedRead[];
  totalRows?: number;
  alleleTicks?: AlleleTick[];
  overviewBins: number[];
  meanDepth: number;
}

export interface WindowQueryResult {
  chrom: string;
  start: number;
  end: number;
  isOverviewOnly: boolean;
  refSeq: string;
  features: GffFeature[];
  strains: StrainWindowData[];
  analysis: ThreeWayRegionAnalysis | null;
}

export interface AllChromsStrainProfile {
  id: string;
  label: string;
  byChrom: Array<{
    chrom: string;
    matchedChrom: string | null;
    bins: number[];
    meanDepth: number;
  }>;
}

export interface AllChromsOverviewResult {
  contigs: ContigInfo[];
  strains: AllChromsStrainProfile[];
}

export interface ExtractionPayload {
  chrom: string;
  start: number;
  end: number;
  primaryStrainIndex: number;
  secondaryStrainIndex: number | null;
  options: ConsensusOptions;
  exons?: Array<{ start: number; end: number }>;
}

export interface ExtractedExtraStrain {
  label: string;
  ungappedSeq: string;
  alignedSeq: string;
  splicedUngappedSeq?: string;
  splicedAlignedSeq?: string;
  mappingStats?: RegionMappingStats;
  splicedMappingStats?: RegionMappingStats;
}

export interface ExtractionResult {
  analysis: ThreeWayRegionAnalysis;
  hasStrain1: boolean;
  strain1Label: string;
  strain2Label: string | null;
  extraStrains: ExtractedExtraStrain[];
}

export interface ScanProgressEvent {
  type: 'scan_progress';
  chrom: string;
  scannedBp: number;
  totalBp: number;
  variantsSoFar: number;
}

export interface LoadProgressEvent {
  type: 'load_progress';
  message: string;
}

export type WorkerRequest =
  | { reqId: number; action: 'load_saccer3'; payload?: { sacCer3Url: string } }
  | { reqId: number; action: 'load_demo' }
  | { reqId: number; action: 'load_files'; payload: LoadLocalFilesPayload }
  | { reqId: number; action: 'add_strains'; payload: { strains: StrainFileInput[]; currentChrom: string } }
  | { reqId: number; action: 'remove_strain'; payload: { id: string; currentChrom: string } }
  | { reqId: number; action: 'query_all_chroms'; payload?: { contigOverrides?: Record<string, Record<string, string>> } }
  | {
      reqId: number;
      action: 'query_window';
      payload: {
        chrom: string;
        start: number;
        end: number;
        viewSpan?: number;
        primaryStrainIndex: number;
        secondaryStrainIndex: number | null;
        options: ConsensusOptions;
        contigOverrides?: Record<string, Record<string, string>>;
      };
    }
  | {
      reqId: number;
      action: 'extract_region';
      payload: ExtractionPayload & {
        contigOverrides?: Record<string, Record<string, string>>;
      };
    }
  | {
      reqId: number;
      action: 'scan_chromosome';
      payload: {
        chrom: string;
        primaryStrainIndex: number;
        secondaryStrainIndex: number | null;
        options: ConsensusOptions;
        contigOverrides?: Record<string, Record<string, string>>;
      };
    }
  | { reqId: number; action: 'cancel_scan' }
  | { reqId: number; action: 'search_genes'; payload: { query: string; limit?: number } };

export type WorkerResponse =
  | { reqId: number; ok: true; data: unknown }
  | { reqId: number; ok: false; error: string }
  | ScanProgressEvent
  | LoadProgressEvent;
