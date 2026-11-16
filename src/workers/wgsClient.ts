import type { ConsensusOptions, GffFeature, VariantRecord } from '../bio/types';
import type {
  AllChromsOverviewResult,
  BarcodeProgressEvent,
  BarcodeScanPayload,
  BarcodeScanResult,
  BlastProgressEvent,
  BlastRunPayload,
  BlastSearchResult,
  ExtractionPayload,
  ExtractionResult,
  LoadLocalFilesPayload,
  LoadProgressEvent,
  ScanProgressEvent,
  StrainFileInput,
  WindowQueryResult,
  WorkspaceLoadedSummary,
  WorkerResponse,
} from './wgsProtocol';

interface PendingCall {
  resolve: (data: unknown) => void;
  reject: (err: Error) => void;
}

export class WgsWorkerClient {
  private worker: Worker;
  private nextReqId = 1;
  private pending = new Map<number, PendingCall>();
  private onScanProgress?: (ev: ScanProgressEvent) => void;
  private onLoadProgress?: (ev: LoadProgressEvent) => void;
  private onBlastProgress?: (ev: BlastProgressEvent) => void;
  private onBarcodeProgress?: (ev: BarcodeProgressEvent) => void;

  constructor(
    onScanProgress?: (ev: ScanProgressEvent) => void,
    onLoadProgress?: (ev: LoadProgressEvent) => void,
    onBlastProgress?: (ev: BlastProgressEvent) => void,
    onBarcodeProgress?: (ev: BarcodeProgressEvent) => void,
  ) {
    this.onScanProgress = onScanProgress;
    this.onLoadProgress = onLoadProgress;
    this.onBlastProgress = onBlastProgress;
    this.onBarcodeProgress = onBarcodeProgress;
    this.worker = new Worker(new URL('./wgsWorker.ts', import.meta.url), {
      type: 'module',
    });

    this.worker.onmessage = (ev: MessageEvent<WorkerResponse>) => {
      const data = ev.data;
      if ('type' in data) {
        if (data.type === 'scan_progress') {
          this.onScanProgress?.(data);
          return;
        }
        if (data.type === 'load_progress') {
          this.onLoadProgress?.(data);
          return;
        }
        if (data.type === 'blast_progress') {
          this.onBlastProgress?.(data);
          return;
        }
        if (data.type === 'barcode_progress') {
          this.onBarcodeProgress?.(data);
          return;
        }
      }
      const res = data as Exclude<
        WorkerResponse,
        ScanProgressEvent | LoadProgressEvent | BlastProgressEvent | BarcodeProgressEvent
      >;
      const call = this.pending.get(res.reqId);
      if (!call) return;
      this.pending.delete(res.reqId);
      if (res.ok) {
        call.resolve(res.data);
      } else {
        call.reject(new Error(res.error));
      }
    };

    this.worker.onerror = (ev) => {
      for (const call of this.pending.values()) {
        call.reject(new Error(ev.message || 'Worker runtime error'));
      }
      this.pending.clear();
    };
  }

  private call<T>(action: string, payload?: unknown): Promise<T> {
    const reqId = this.nextReqId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(reqId, {
        resolve: resolve as (d: unknown) => void,
        reject,
      });
      if (payload !== undefined) {
        this.worker.postMessage({ reqId, action, payload });
      } else {
        this.worker.postMessage({ reqId, action });
      }
    });
  }

  loadSacCer3(sacCer3Url: string): Promise<WorkspaceLoadedSummary> {
    return this.call<WorkspaceLoadedSummary>('load_saccer3', { sacCer3Url });
  }

  loadDemo(): Promise<WorkspaceLoadedSummary> {
    return this.call<WorkspaceLoadedSummary>('load_demo');
  }

  loadFiles(payload: LoadLocalFilesPayload): Promise<WorkspaceLoadedSummary> {
    return this.call<WorkspaceLoadedSummary>('load_files', payload);
  }

  addStrains(strains: StrainFileInput[], currentChrom: string): Promise<WorkspaceLoadedSummary> {
    return this.call<WorkspaceLoadedSummary>('add_strains', { strains, currentChrom });
  }

  removeStrain(id: string, currentChrom: string): Promise<WorkspaceLoadedSummary> {
    return this.call<WorkspaceLoadedSummary>('remove_strain', { id, currentChrom });
  }

  queryAllChroms(
    contigOverrides?: Record<string, Record<string, string>>,
  ): Promise<AllChromsOverviewResult> {
    return this.call<AllChromsOverviewResult>('query_all_chroms', { contigOverrides });
  }

  queryWindow(payload: {
    chrom: string;
    start: number;
    end: number;
    viewSpan?: number;
    primaryStrainIndex: number;
    secondaryStrainIndex: number | null;
    options: ConsensusOptions;
    contigOverrides?: Record<string, Record<string, string>>;
  }): Promise<WindowQueryResult> {
    return this.call<WindowQueryResult>('query_window', payload);
  }

  extractRegion(
    payload: ExtractionPayload & {
      contigOverrides?: Record<string, Record<string, string>>;
    },
  ): Promise<ExtractionResult> {
    return this.call<ExtractionResult>('extract_region', payload);
  }

  scanChromosome(payload: {
    chrom: string;
    primaryStrainIndex: number;
    secondaryStrainIndex: number | null;
    options: ConsensusOptions;
    contigOverrides?: Record<string, Record<string, string>>;
  }): Promise<VariantRecord[]> {
    return this.call<VariantRecord[]>('scan_chromosome', payload);
  }

  cancelScan(): Promise<boolean> {
    return this.call<boolean>('cancel_scan');
  }

  searchGenes(query: string, limit = 25): Promise<GffFeature[]> {
    return this.call<GffFeature[]>('search_genes', { query, limit });
  }

  runBlast(payload: BlastRunPayload): Promise<BlastSearchResult> {
    return this.call<BlastSearchResult>('run_blast', payload);
  }

  cancelBlast(): Promise<boolean> {
    return this.call<boolean>('cancel_blast');
  }

  scanBarcodes(payload: BarcodeScanPayload): Promise<BarcodeScanResult> {
    return this.call<BarcodeScanResult>('scan_barcodes', payload);
  }

  cancelBarcodeScan(): Promise<boolean> {
    return this.call<boolean>('cancel_barcode_scan');
  }

  dispose() {
    this.worker.terminate();
  }
}
