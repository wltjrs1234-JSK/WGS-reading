import type { BarcodeCandidate, BarcodeScanOptions } from '../bio/barcodeScanner';
import type {
  CnvPloidyPreference,
  CnvScanMode,
} from '../bio/cnvSvAnalyzer';
import type { RegionMappingStats } from '../bio/consensus';
import type {
  BlastForeignReadHit,
  BlastRawHit,
  BlastStrainComparison,
} from '../bio/localBlast';
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

export type VariantScanScope = 'current_chrom' | 'all_chroms';

export interface ScanProgressEvent {
  type: 'scan_progress';
  chrom: string;
  scannedBp: number;
  totalBp: number;
  variantsSoFar: number;
  chromIndex?: number;
  totalChroms?: number;
}

export interface LoadProgressEvent {
  type: 'load_progress';
  message: string;
}

export interface BlastProgressEvent {
  type: 'blast_progress';
  message: string;
  percent: number;
}

export interface BlastRunPayload {
  rawInput: string;
  maxHitsPerQuery?: number;
  minIdentityPct?: number;
  shortPrimerMaxMismatches?: number;
  foreignBamScope?: 'off' | 'current_chrom' | 'all_chroms';
  currentChrom?: string;
  options: ConsensusOptions;
  contigOverrides?: Record<string, Record<string, string>>;
}

export interface BlastOverlappingGene {
  id: string;
  name: string;
  systematicName?: string;
  strand: '+' | '-' | '.';
  start: number;
  end: number;
  relation: 'exact_or_contains' | 'inside_gene' | 'partial_overlap';
}

export interface BlastHitWithStrains extends BlastRawHit {
  alignedQueryDisplay: string;
  alignedRefDisplay: string;
  overlappingGenes: BlastOverlappingGene[];
  strains: BlastStrainComparison[];
}

export interface BlastQueryResult {
  queryId: string;
  queryLength: number;
  hits: BlastHitWithStrains[];
  foreignHits: BlastForeignReadHit[];
}

export interface BlastSearchResult {
  referenceLabel: string;
  totalRefLength: number;
  elapsedMs: number;
  queries: BlastQueryResult[];
  strainLabels: string[];
}

export type BarcodeScanScope = 'hotspots' | 'current_chrom' | 'all_chroms';

export interface BarcodeScanPayload {
  scope: BarcodeScanScope;
  currentChrom?: string;
  targetStrainIndex?: number | null; // null = scan all loaded strains and compare
  options?: BarcodeScanOptions;
  contigOverrides?: Record<string, Record<string, string>>;
}

export interface BarcodeProgressEvent {
  type: 'barcode_progress';
  message: string;
  percent: number;
  foundSoFar: number;
}

export interface BarcodeScanFilteredStats {
  repeatLociFiltered: number;
  tandemOrArtifactFiltered: number;
  endogenousRefFiltered: number;
}

export interface BarcodeScanResult {
  scope: BarcodeScanScope;
  scannedLociCount: number;
  scannedBpTotal: number;
  elapsedMs: number;
  strainLabels: string[];
  candidates: BarcodeCandidate[];
  filteredStats?: BarcodeScanFilteredStats;
}

export interface CnvSvPayload {
  scanMode: CnvScanMode;
  ploidyPreference: CnvPloidyPreference;
  targetStrainIndex?: number | null;
  excludeRdnaSv?: boolean;
  minSvSupport?: number;
  contigOverrides?: Record<string, Record<string, string>>;
}

export interface CnvProgressEvent {
  type: 'cnv_progress';
  message: string;
  percent: number;
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
        scope?: VariantScanScope;
        chrom: string;
        primaryStrainIndex: number;
        secondaryStrainIndex: number | null;
        options: ConsensusOptions;
        contigOverrides?: Record<string, Record<string, string>>;
      };
    }
  | { reqId: number; action: 'cancel_scan' }
  | { reqId: number; action: 'cancel_window' }
  | { type: 'read_blob_reply'; relayId: number; ok: boolean; buffer?: ArrayBuffer; error?: string }
  | { reqId: number; action: 'search_genes'; payload: { query: string; limit?: number } }
  | { reqId: number; action: 'run_blast'; payload: BlastRunPayload }
  | { reqId: number; action: 'cancel_blast' }
  | { reqId: number; action: 'scan_barcodes'; payload: BarcodeScanPayload }
  | { reqId: number; action: 'cancel_barcode_scan' }
  | { reqId: number; action: 'analyze_cnv_sv'; payload: CnvSvPayload }
  | { reqId: number; action: 'cancel_cnv_sv' };

export interface BlobRelayRequestEvent {
  type: 'read_blob_slice';
  relayId: number;
  fileKey: string;
  start: number;
  end: number;
}

export type WorkerResponse =
  | { reqId: number; ok: true; data: unknown }
  | { reqId: number; ok: false; error: string; aborted?: boolean }
  | BlobRelayRequestEvent
  | ScanProgressEvent
  | LoadProgressEvent
  | BlastProgressEvent
  | BarcodeProgressEvent
  | CnvProgressEvent;
