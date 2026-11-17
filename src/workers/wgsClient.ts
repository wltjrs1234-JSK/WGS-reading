import type { CnvSvAnalysisResult } from '../bio/cnvSvAnalyzer';
import type { ConsensusOptions, GffFeature, VariantRecord } from '../bio/types';
import type {
  AllChromsOverviewResult,
  BarcodeProgressEvent,
  BarcodeScanPayload,
  BarcodeScanResult,
  BlastProgressEvent,
  BlastRunPayload,
  BlastSearchResult,
  BlobRelayRequestEvent,
  CnvProgressEvent,
  CnvSvPayload,
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
  private worker!: Worker;
  private nextReqId = 1;
  private pending = new Map<number, PendingCall>();
  private blobRegistry = new Map<string, File>();
  private onScanProgress?: (ev: ScanProgressEvent) => void;
  private onLoadProgress?: (ev: LoadProgressEvent) => void;
  private onBlastProgress?: (ev: BlastProgressEvent) => void;
  private onBarcodeProgress?: (ev: BarcodeProgressEvent) => void;
  private onCnvProgress?: (ev: CnvProgressEvent) => void;

  constructor(
    onScanProgress?: (ev: ScanProgressEvent) => void,
    onLoadProgress?: (ev: LoadProgressEvent) => void,
    onBlastProgress?: (ev: BlastProgressEvent) => void,
    onBarcodeProgress?: (ev: BarcodeProgressEvent) => void,
    onCnvProgress?: (ev: CnvProgressEvent) => void,
  ) {
    this.onScanProgress = onScanProgress;
    this.onLoadProgress = onLoadProgress;
    this.onBlastProgress = onBlastProgress;
    this.onBarcodeProgress = onBarcodeProgress;
    this.onCnvProgress = onCnvProgress;
    this.spawnWorker();
  }

  private spawnWorker(): void {
    this.worker = new Worker(new URL('./wgsWorker.ts', import.meta.url), {
      type: 'module',
    });

    this.worker.onmessage = (ev: MessageEvent<WorkerResponse>) => {
      const data = ev.data;
      if ('type' in data) {
        if (data.type === 'read_blob_slice') {
          void this.handleBlobRelayRequest(data);
          return;
        }
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
        if (data.type === 'cnv_progress') {
          this.onCnvProgress?.(data);
          return;
        }
      }
      const res = data as Exclude<
        WorkerResponse,
        | BlobRelayRequestEvent
        | ScanProgressEvent
        | LoadProgressEvent
        | BlastProgressEvent
        | BarcodeProgressEvent
        | CnvProgressEvent
      >;
      const call = this.pending.get(res.reqId);
      if (!call) return;
      this.pending.delete(res.reqId);
      if (res.ok) {
        call.resolve(res.data);
      } else {
        const err = new Error(res.error);
        if (res.aborted) err.name = 'AbortError';
        call.reject(err);
      }
    };

    const failPending = (reason: string) => {
      for (const call of this.pending.values()) {
        call.reject(new Error(reason));
      }
      this.pending.clear();
    };

    this.worker.onerror = (ev) => {
      failPending(ev.message || '워커 런타임 오류가 발생했습니다.');
    };
    this.worker.onmessageerror = () => {
      failPending('워커 메시지 직렬화 오류가 발생했습니다.');
    };
  }

  private async handleBlobRelayRequest(req: BlobRelayRequestEvent): Promise<void> {
    const targetWorker = this.worker;
    const file = this.blobRegistry.get(req.fileKey);
    const postSafe = (payload: unknown, transfer?: Transferable[]) => {
      if (this.worker !== targetWorker) return;
      try {
        if (transfer) {
          targetWorker.postMessage(payload, transfer);
        } else {
          targetWorker.postMessage(payload);
        }
      } catch {
        // worker already terminated
      }
    };
    if (!file) {
      postSafe({
        type: 'read_blob_reply',
        relayId: req.relayId,
        ok: false,
        error: 'Main-thread file reference unavailable',
      });
      return;
    }
    try {
      const ab = await file.slice(req.start, req.end).arrayBuffer();
      postSafe({ type: 'read_blob_reply', relayId: req.relayId, ok: true, buffer: ab }, [ab]);
    } catch (err) {
      postSafe({
        type: 'read_blob_reply',
        relayId: req.relayId,
        ok: false,
        error: err instanceof Error ? err.message : String(err),
      });
    }
  }

  private registerStrains(strains: StrainFileInput[]): void {
    for (const st of strains) {
      this.blobRegistry.set(`${st.id}:bam`, st.bamFile);
      this.blobRegistry.set(`${st.id}:bai`, st.indexFile);
    }
  }

  /**
   * OS I/O 정체 또는 워커 비정상 상태를 사용자가 즉시 중단하고 새 워커로 복구합니다.
   */
  restartWorker(): void {
    try {
      this.worker.terminate();
    } catch {
      // ignore
    }
    for (const call of this.pending.values()) {
      const err = new Error('작업이 사용자에 의해 중단되었습니다.');
      err.name = 'AbortError';
      call.reject(err);
    }
    this.pending.clear();
    this.spawnWorker();
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
    this.blobRegistry.set('custom:fasta', payload.fastaFile);
    this.registerStrains(payload.strains);
    return this.call<WorkspaceLoadedSummary>('load_files', payload);
  }

  addStrains(strains: StrainFileInput[], currentChrom: string): Promise<WorkspaceLoadedSummary> {
    this.registerStrains(strains);
    return this.call<WorkspaceLoadedSummary>('add_strains', { strains, currentChrom });
  }

  cancelWindow(): Promise<boolean> {
    return this.call<boolean>('cancel_window');
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
    scope?: 'current_chrom' | 'all_chroms';
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

  analyzeCnvSv(payload: CnvSvPayload): Promise<CnvSvAnalysisResult> {
    return this.call<CnvSvAnalysisResult>('analyze_cnv_sv', payload);
  }

  cancelCnvSv(): Promise<boolean> {
    return this.call<boolean>('cancel_cnv_sv');
  }

  dispose() {
    this.worker.terminate();
  }
}
