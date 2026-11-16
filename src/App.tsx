import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import {
  AlignJustify,
  Crosshair,
  Dna,
  FolderOpen,
  Hand,
  Maximize2,
  PanelRight,
  Plus,
  QrCode,
  Scissors,
  Search,
  SlidersHorizontal,
  X,
  ZoomIn,
  ZoomOut,
} from 'lucide-react';
import {
  DEFAULT_CONSENSUS_OPTIONS,
  type ConsensusOptions,
  type GffFeature,
  type VariantRecord,
} from './bio/types';
import { AllChromsOverviewCanvas } from './components/AllChromsOverviewCanvas';
import { BarcodeScannerModal } from './components/BarcodeScannerModal';
import { BlastModal } from './components/BlastModal';
import { ChromosomeIdeogram } from './components/ChromosomeIdeogram';
import { FileLoaderModal } from './components/FileLoaderModal';
import { GenomeCanvasViewer, type InteractionMode } from './components/GenomeCanvasViewer';
import { SequenceExtractorModal } from './components/SequenceExtractorModal';
import { VariantSidebar } from './components/VariantSidebar';
import { matchDroppedFiles } from './utils/fileMatcher';
import { WgsWorkerClient } from './workers/wgsClient';
import type {
  AllChromsOverviewResult,
  BarcodeProgressEvent,
  BarcodeScanPayload,
  BarcodeScanResult,
  BlastProgressEvent,
  BlastRunPayload,
  BlastSearchResult,
  ExtractionResult,
  LoadLocalFilesPayload,
  ScanProgressEvent,
  StrainFileInput,
  WindowQueryResult,
  WorkspaceLoadedSummary,
} from './workers/wgsProtocol';

function parseCoordinateInput(raw: string): { chrom?: string; start: number; end: number } | null {
  const clean = raw.trim().replace(/,/g, '');
  const m = clean.match(/^(?:([^:]+):)?(\d+)\s*[-–~]\s*(\d+)$/);
  if (!m) return null;
  const s1 = Number(m[2]);
  const e1 = Number(m[3]);
  if (Number.isNaN(s1) || Number.isNaN(e1)) return null;
  const start = Math.max(0, Math.min(s1, e1) - 1);
  const end = Math.max(start + 10, Math.max(s1, e1));
  return { chrom: m[1]?.trim(), start, end };
}

function computeBufferedWindowRange(targetStart: number, targetEnd: number, chromLen: number) {
  const safeLen = Math.max(1, chromLen);
  const vSpan = Math.max(1, Math.min(safeLen, Math.round(targetEnd - targetStart)));
  const vStart = Math.max(0, Math.min(safeLen - vSpan, Math.round(targetStart)));
  const vEnd = Math.min(safeLen, Math.max(vStart + 1, vStart + vSpan));
  const pad =
    vSpan <= 6_000
      ? Math.min(2_500, Math.max(600, Math.round(vSpan * 0.55)))
      : vSpan <= 30_000
        ? Math.min(2_000, Math.max(0, Math.floor((30_000 - vSpan) / 2)))
        : Math.round(vSpan * 0.2);
  const bufStart = Math.max(0, vStart - pad);
  const bufEnd = Math.min(safeLen, vEnd + pad);
  return { viewStart: vStart, viewEnd: vEnd, bufStart, bufEnd, viewSpan: vEnd - vStart };
}

export default function App() {
  const workerRef = useRef<WgsWorkerClient | null>(null);
  const latestQueryRef = useRef(0);
  const searchDropdownRef = useRef<HTMLDivElement>(null);
  const appRootRef = useRef<HTMLDivElement>(null);

  const [workspace, setWorkspace] = useState<WorkspaceLoadedSummary | null>(null);
  const [chrom, setChrom] = useState('chrV');
  const [start, setStart] = useState(115_600);
  const [end, setEnd] = useState(117_600);
  const [navVersion, setNavVersion] = useState(0);

  // IGV 통합 위치/유전자 검색창 상태
  const [locusInput, setLocusInput] = useState('');
  const [geneSuggestions, setGeneSuggestions] = useState<GffFeature[]>([]);
  const [showSearchDropdown, setShowSearchDropdown] = useState(false);
  const [highlightedGene, setHighlightedGene] = useState<GffFeature | null>(null);

  const [windowData, setWindowData] = useState<WindowQueryResult | null>(null);
  const [allChromsData, setAllChromsData] = useState<AllChromsOverviewResult | null>(null);

  const [primaryIdx, setPrimaryIdx] = useState(0);
  const [secondaryIdx, setSecondaryIdx] = useState<number | null>(null);
  const [options, setOptions] = useState<ConsensusOptions>(DEFAULT_CONSENSUS_OPTIONS);
  const [contigOverrides, setContigOverrides] = useState<Record<string, Record<string, string>>>({});

  const [interactionMode, setInteractionMode] = useState<InteractionMode>('pan');
  const [compactMode, setCompactMode] = useState(false);
  const [showAllBases, setShowAllBases] = useState(false);
  const [showCompareStrip, setShowCompareStrip] = useState(true);
  const [showCrosshairs, setShowCrosshairs] = useState(false);
  const [showCenterLine, setShowCenterLine] = useState(false);
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);

  const [isLoaderOpen, setIsLoaderOpen] = useState(false);
  const [loaderMode, setLoaderMode] = useState<'tracks' | 'genome'>('tracks');
  const [isLoadingFiles, setIsLoadingFiles] = useState(false);
  const [loadingMessage, setLoadingMessage] = useState<string | null>(null);
  const [loaderError, setLoaderError] = useState<string | null>(null);
  const [dropError, setDropError] = useState<string | null>(null);
  const [isWindowDragOver, setIsWindowDragOver] = useState(false);
  const [dropNotice, setDropNotice] = useState<string | null>(null);
  const initPromiseRef = useRef<Promise<void> | null>(null);

  const [activeSelection, setActiveSelection] = useState<{ start: number; end: number } | null>(null);
  const [isExtractorOpen, setIsExtractorOpen] = useState(false);
  const [isExtracting, setIsExtracting] = useState(false);
  const [extractError, setExtractError] = useState<string | null>(null);
  const [extractResult, setExtractResult] = useState<ExtractionResult | null>(null);
  const [extractFeature, setExtractFeature] = useState<GffFeature | null>(null);

  const [scannedVariants, setScannedVariants] = useState<VariantRecord[]>([]);
  const [isScanning, setIsScanning] = useState(false);
  const [scanProgress, setScanProgress] = useState<ScanProgressEvent | null>(null);

  const [isBlastOpen, setIsBlastOpen] = useState(false);
  const [isBlasting, setIsBlasting] = useState(false);
  const [blastProgress, setBlastProgress] = useState<BlastProgressEvent | null>(null);
  const [blastError, setBlastError] = useState<string | null>(null);
  const [blastResult, setBlastResult] = useState<BlastSearchResult | null>(null);
  const [blastInitialQuery, setBlastInitialQuery] = useState('');

  const [isBarcodeOpen, setIsBarcodeOpen] = useState(false);
  const [isScanningBarcodes, setIsScanningBarcodes] = useState(false);
  const [barcodeProgress, setBarcodeProgress] = useState<BarcodeProgressEvent | null>(null);
  const [barcodeError, setBarcodeError] = useState<string | null>(null);
  const [barcodeResult, setBarcodeResult] = useState<BarcodeScanResult | null>(null);

  const currentChromLength = useMemo(() => {
    if (!workspace) return 576874;
    return workspace.contigs.find((c) => c.name === chrom)?.length ?? 576874;
  }, [workspace, chrom]);

  const workspaceRef = useRef(workspace);
  workspaceRef.current = workspace;
  const windowDataRef = useRef(windowData);
  windowDataRef.current = windowData;
  const chromRef = useRef(chrom);
  chromRef.current = chrom;
  const strainParamsRef = useRef({ primaryIdx, secondaryIdx, options, contigOverrides });
  strainParamsRef.current = { primaryIdx, secondaryIdx, options, contigOverrides };

  const isPrefetchingRef = useRef(false);
  const pendingPrefetchRef = useRef<{ chrom: string; start: number; end: number } | null>(null);
  const prefetchTimerRef = useRef<number | null>(null);
  const skipNextEffectQueryRef = useRef(false);

  const pumpPrefetchQueue = useCallback(async () => {
    if (isPrefetchingRef.current) return;
    const req = pendingPrefetchRef.current;
    if (!req) return;
    pendingPrefetchRef.current = null;

    const client = workerRef.current;
    const ws = workspaceRef.current;
    if (!client || !ws || req.chrom === 'all') return;

    const cLen = ws.contigs.find((c) => c.name === req.chrom)?.length ?? Math.max(req.end, 576874);
    const buf = computeBufferedWindowRange(req.start, req.end, cLen);
    const ticket = latestQueryRef.current;
    const { primaryIdx: pIdx, secondaryIdx: sIdx, options: opts, contigOverrides: overrides } =
      strainParamsRef.current;

    isPrefetchingRef.current = true;
    try {
      const res = await client.queryWindow({
        chrom: req.chrom,
        start: buf.bufStart,
        end: buf.bufEnd,
        viewSpan: buf.viewSpan,
        primaryStrainIndex: pIdx,
        secondaryStrainIndex: sIdx,
        options: opts,
        contigOverrides: overrides,
      });
      if (ticket === latestQueryRef.current && chromRef.current === req.chrom) {
        setWindowData(res);
      }
    } catch (err) {
      console.error(err);
    } finally {
      isPrefetchingRef.current = false;
      if (pendingPrefetchRef.current) {
        void pumpPrefetchQueue();
      }
    }
  }, []);

  const handleViewerViewportChange = useCallback(
    (nextStart: number, nextEnd: number, urgent = false) => {
      const curChrom = chromRef.current;
      const ws = workspaceRef.current;
      if (!curChrom || curChrom === 'all' || !ws) return;

      const cLen = ws.contigs.find((c) => c.name === curChrom)?.length ?? Math.max(nextEnd, 576874);
      const s = Math.max(0, Math.round(nextStart));
      const e = Math.min(cLen, Math.max(s + 1, Math.round(nextEnd)));

      setStart(s);
      setEnd(e);
      setLocusInput(`${curChrom}:${(s + 1).toLocaleString()}-${e.toLocaleString()}`);

      const wd = windowDataRef.current;
      const vSpan = Math.max(1, e - s);
      const wantOverview = vSpan > 30_000;
      const margin = Math.min(vSpan * 0.2, 1200);
      const needsRefresh =
        !wd ||
        wd.chrom !== curChrom ||
        wd.isOverviewOnly !== wantOverview ||
        (s < wd.start + margin && wd.start > 0) ||
        (e > wd.end - margin && wd.end < cLen) ||
        (wantOverview && (vSpan < (wd.end - wd.start) * 0.5 || vSpan > (wd.end - wd.start) * 1.15));

      if (!needsRefresh) return;

      pendingPrefetchRef.current = { chrom: curChrom, start: s, end: e };
      if (prefetchTimerRef.current !== null) {
        window.clearTimeout(prefetchTimerRef.current);
        prefetchTimerRef.current = null;
      }

      if (urgent && (!wd || wd.isOverviewOnly === wantOverview)) {
        void pumpPrefetchQueue();
      } else {
        prefetchTimerRef.current = window.setTimeout(() => {
          prefetchTimerRef.current = null;
          void pumpPrefetchQueue();
        }, 55);
      }
    },
    [pumpPrefetchQueue],
  );

  const fetchWindow = useCallback(
    async (
      targetChrom: string,
      targetStart: number,
      targetEnd: number,
      pIdx = primaryIdx,
      sIdx = secondaryIdx,
      opts = options,
      overrides = contigOverrides,
    ) => {
      const client = workerRef.current;
      if (!client || !targetChrom) return;
      pendingPrefetchRef.current = null;
      if (prefetchTimerRef.current !== null) {
        window.clearTimeout(prefetchTimerRef.current);
        prefetchTimerRef.current = null;
      }
      const ticket = ++latestQueryRef.current;
      if (targetChrom === 'all') {
        try {
          const allRes = await client.queryAllChroms(overrides);
          if (ticket === latestQueryRef.current) {
            setChrom('all');
            setAllChromsData(allRes);
            setLocusInput('all');
          }
        } catch (err) {
          console.error(err);
        }
        return;
      }

      const cLen =
        workspaceRef.current?.contigs.find((c) => c.name === targetChrom)?.length ??
        Math.max(targetEnd, 576874);
      const buf = computeBufferedWindowRange(targetStart, targetEnd, cLen);

      try {
        const res = await client.queryWindow({
          chrom: targetChrom,
          start: buf.bufStart,
          end: buf.bufEnd,
          viewSpan: buf.viewSpan,
          primaryStrainIndex: pIdx,
          secondaryStrainIndex: sIdx,
          options: opts,
          contigOverrides: overrides,
        });
        if (ticket === latestQueryRef.current) {
          setWindowData(res);
          setChrom(res.chrom);
          setStart(buf.viewStart);
          setEnd(buf.viewEnd);
          setNavVersion((v) => v + 1);
          setLocusInput(
            `${res.chrom}:${(buf.viewStart + 1).toLocaleString()}-${buf.viewEnd.toLocaleString()}`,
          );
        }
      } catch (err) {
        console.error(err);
      }
    },
    [primaryIdx, secondaryIdx, options, contigOverrides],
  );

  // 시작 시 sacCer3 (전체 17개 염색체 + 7,728 S288C 유전자) 자동 로드
  useEffect(() => {
    const client = new WgsWorkerClient(
      (prog) => {
        setScanProgress(prog);
      },
      (loadEv) => {
        setLoadingMessage(loadEv.message);
      },
      (blastEv) => {
        setBlastProgress(blastEv);
      },
      (barEv) => {
        setBarcodeProgress(barEv);
      },
    );
    workerRef.current = client;

    const init = async () => {
      setIsLoadingFiles(true);
      setLoadingMessage('sacCer3 기본 지놈 불러오는 중...');
      try {
        const url = new URL('/genomes/sacCer3.fa.gz', window.location.origin).toString();
        const summary = await client.loadSacCer3(url);
        workspaceRef.current = summary;
        setWorkspace(summary);
        setPrimaryIdx(0);
        setSecondaryIdx(null);
        const defLen =
          summary.contigs.find((c) => c.name === summary.defaultChrom)?.length ??
          summary.defaultEnd;
        const buf = computeBufferedWindowRange(summary.defaultStart, summary.defaultEnd, defLen);
        const res = await client.queryWindow({
          chrom: summary.defaultChrom,
          start: buf.bufStart,
          end: buf.bufEnd,
          viewSpan: buf.viewSpan,
          primaryStrainIndex: 0,
          secondaryStrainIndex: null,
          options: DEFAULT_CONSENSUS_OPTIONS,
          contigOverrides: {},
        });
        setWindowData(res);
        setChrom(res.chrom);
        setStart(buf.viewStart);
        setEnd(buf.viewEnd);
        setNavVersion((v) => v + 1);
        setLocusInput(
          `${res.chrom}:${(buf.viewStart + 1).toLocaleString()}-${buf.viewEnd.toLocaleString()}`,
        );
      } catch (err) {
        const msg = err instanceof Error ? err.message : String(err);
        setLoaderError(msg);
        setDropError(msg);
        setIsLoaderOpen(true);
      } finally {
        setIsLoadingFiles(false);
        setLoadingMessage(null);
      }
    };
    const p = init();
    initPromiseRef.current = p;
    void p;

    return () => {
      if (prefetchTimerRef.current !== null) {
        window.clearTimeout(prefetchTimerRef.current);
      }
      client.dispose();
      workerRef.current = null;
    };
  }, []);

  // 브라우저 창 전체에서 Ctrl+휠 페이지 확대 차단
  useEffect(() => {
    const el = appRootRef.current;
    if (!el) return;
    const preventPageZoom = (e: WheelEvent) => {
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
      }
    };
    el.addEventListener('wheel', preventPageZoom, { passive: false });
    return () => el.removeEventListener('wheel', preventPageZoom);
  }, []);

  // 전역 Esc 키 처리
  useEffect(() => {
    const handleEsc = (e: KeyboardEvent) => {
      if (e.key !== 'Escape') return;
      if (showSearchDropdown) {
        setShowSearchDropdown(false);
        return;
      }
      if (isBlastOpen) {
        setIsBlastOpen(false);
        return;
      }
      if (isBarcodeOpen) {
        setIsBarcodeOpen(false);
        return;
      }
      if (isExtractorOpen) {
        setIsExtractorOpen(false);
        return;
      }
      if (isLoaderOpen) {
        setIsLoaderOpen(false);
        return;
      }
      if (activeSelection) {
        setActiveSelection(null);
        return;
      }
      if (interactionMode === 'select') {
        setInteractionMode('pan');
      }
    };
    window.addEventListener('keydown', handleEsc);
    return () => window.removeEventListener('keydown', handleEsc);
  }, [showSearchDropdown, isBlastOpen, isBarcodeOpen, isExtractorOpen, isLoaderOpen, activeSelection, interactionMode]);

  useEffect(() => {
    const onDocClick = (ev: MouseEvent) => {
      if (searchDropdownRef.current && !searchDropdownRef.current.contains(ev.target as Node)) {
        setShowSearchDropdown(false);
      }
    };
    document.addEventListener('mousedown', onDocClick);
    return () => document.removeEventListener('mousedown', onDocClick);
  }, []);

  useEffect(() => {
    if (skipNextEffectQueryRef.current) {
      skipNextEffectQueryRef.current = false;
      return;
    }
    if (workspace && chrom && chrom !== 'all') {
      void fetchWindow(chrom, start, end, primaryIdx, secondaryIdx, options, contigOverrides);
    } else if (workspace && chrom === 'all') {
      void fetchWindow('all', 0, 1);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [primaryIdx, secondaryIdx, options, contigOverrides]);

  const applyWorkspaceSummary = async (
    summary: WorkspaceLoadedSummary,
    keepCurrentLocus = false,
  ) => {
    const s2Idx = summary.strains.length > 1 ? 1 : null;
    if (primaryIdx !== 0 || secondaryIdx !== s2Idx) {
      skipNextEffectQueryRef.current = true;
    }
    workspaceRef.current = summary;
    setWorkspace(summary);
    setPrimaryIdx(0);
    setSecondaryIdx(s2Idx);
    setScannedVariants([]);
    setActiveSelection(null);

    if (keepCurrentLocus && chrom) {
      if (chrom === 'all') {
        await fetchWindow('all', 0, 1, 0, s2Idx, options, contigOverrides);
      } else {
        const matched = summary.contigs.find((c) => c.name === chrom);
        if (matched) {
          await fetchWindow(chrom, start, Math.min(matched.length, end), 0, s2Idx, options, contigOverrides);
          return;
        }
      }
    }
    setHighlightedGene(null);
    await fetchWindow(
      summary.defaultChrom,
      summary.defaultStart,
      summary.defaultEnd,
      0,
      s2Idx,
      options,
      contigOverrides,
    );
  };

  const handleSelectSacCer3 = async () => {
    const client = workerRef.current;
    if (!client) return;
    setIsLoadingFiles(true);
    setLoaderError(null);
    try {
      const url = new URL('/genomes/sacCer3.fa.gz', window.location.origin).toString();
      const summary = await client.loadSacCer3(url);
      setContigOverrides({});
      setIsLoaderOpen(false);
      await applyWorkspaceSummary(summary, false);
    } catch (err) {
      setLoaderError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsLoadingFiles(false);
    }
  };

  const handleLoadDemo = async () => {
    const client = workerRef.current;
    if (!client) return;
    setIsLoadingFiles(true);
    setLoaderError(null);
    try {
      const summary = await client.loadDemo();
      setContigOverrides({});
      setIsLoaderOpen(false);
      await applyWorkspaceSummary(summary, false);
    } catch (err) {
      setLoaderError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsLoadingFiles(false);
    }
  };

  const handleLoadCustomGenome = async (payload: LoadLocalFilesPayload): Promise<boolean> => {
    const client = workerRef.current;
    if (!client) return false;
    setIsLoadingFiles(true);
    setLoadingMessage('커스텀 Reference · BAM/BAI 불러오는 중...');
    setLoaderError(null);
    setDropError(null);
    try {
      const summary = await client.loadFiles(payload);
      setContigOverrides({});
      await applyWorkspaceSummary(summary, false);
      setIsLoaderOpen(false);
      return true;
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      setLoaderError(msg);
      setDropError(msg);
      return false;
    } finally {
      setIsLoadingFiles(false);
      setLoadingMessage(null);
    }
  };

  const handleAddTracksOnly = async (newStrains: StrainFileInput[]): Promise<boolean> => {
    const client = workerRef.current;
    if (!client || newStrains.length === 0) return false;
    setIsLoadingFiles(true);
    setLoadingMessage(`${newStrains.length}개 BAM/BAI 트랙 불러오는 중...`);
    setLoaderError(null);
    setDropError(null);
    try {
      if (!workspaceRef.current && initPromiseRef.current) {
        setLoadingMessage('sacCer3 기본 지놈 초기화 대기 중...');
        await initPromiseRef.current;
      }
      const summary = await client.addStrains(newStrains, chrom === 'all' ? 'chrV' : chrom);
      await applyWorkspaceSummary(summary, true);
      setIsLoaderOpen(false);

      const totalRef = summary.contigs.length;
      const refLenByChrom = new Map(summary.contigs.map((c) => [c.name, c.length]));
      const addedStrains = summary.strains.slice(-newStrains.length);
      const summaryMsgs: string[] = [];
      let hasWarning = false;

      for (const st of addedStrains) {
        const mappedEntries = Object.entries(st.mappings).filter(
          (entry): entry is [string, string] => Boolean(entry[1]),
        );
        const mappedCount = mappedEntries.length;
        const lenDiffs = mappedEntries.filter(([refChr, bamChr]) => {
          const rLen = refLenByChrom.get(refChr) ?? 0;
          const bLen = st.contigs.find((c) => c.name === bamChr)?.length ?? 0;
          return rLen > 0 && bLen > 0 && Math.abs(rLen - bLen) > 100;
        });

        if (mappedCount === 0) {
          hasWarning = true;
          summaryMsgs.push(
            `'${st.label}': ${totalRef}개 중 매핑된 염색체 0개 (BAM 헤더: ${st.contigs.map((c) => c.name).slice(0, 4).join(', ')}... — 우측 [매핑·설정]에서 직접 지정 필요)`,
          );
        } else if (lenDiffs.length > 0) {
          hasWarning = true;
          summaryMsgs.push(
            `'${st.label}' 추가됨 (${totalRef}개 중 ${mappedCount}개 매핑, 주의: BAM 염색체 길이가 현재 Reference와 ${lenDiffs.length}개 다름)`,
          );
        } else if (mappedCount < totalRef) {
          summaryMsgs.push(
            `'${st.label}' 추가됨: ${totalRef}개 중 ${mappedCount}개 염색체 매핑 (BAM 헤더: ${st.contigs.length}개)`,
          );
        } else {
          summaryMsgs.push(`'${st.label}' 추가됨: ${totalRef}개 전체 염색체 매핑 완료`);
        }
      }

      if (summaryMsgs.length > 0) {
        setDropNotice(summaryMsgs.join(' | '));
        if (!hasWarning) {
          setTimeout(() => setDropNotice(null), 7000);
        }
      }
      return true;
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      setLoaderError(msg);
      setDropError(msg);
      return false;
    } finally {
      setIsLoadingFiles(false);
      setLoadingMessage(null);
    }
  };

  const handleRemoveStrain = async (id: string) => {
    const client = workerRef.current;
    if (!client) return;
    const summary = await client.removeStrain(id, chrom === 'all' ? 'chrV' : chrom);
    await applyWorkspaceSummary(summary, true);
  };

  // 메인 화면 아무 곳에나 BAM+BAI(또는 FASTA/GFF3) 드래그 앤 드롭 처리
  const handleGlobalFileDrop = async (e: React.DragEvent) => {
    e.preventDefault();
    setIsWindowDragOver(false);
    const files = Array.from(e.dataTransfer.files);
    if (files.length === 0) return;
    const matched = matchDroppedFiles(files);

    if (matched.fastaFile) {
      await handleLoadCustomGenome({
        fastaFile: matched.fastaFile,
        faiFile: matched.faiFile,
        gffFile: matched.gffFile,
        strains: matched.pairedStrains.map((s) => ({
          id: s.id,
          label: s.label,
          bamFile: s.bamFile,
          indexFile: s.indexFile,
        })),
      });
      return;
    }

    if (matched.pairedStrains.length > 0) {
      await handleAddTracksOnly(
        matched.pairedStrains.map((s) => ({
          id: s.id,
          label: s.label,
          bamFile: s.bamFile,
          indexFile: s.indexFile,
        })),
      );
      if (matched.unpairedBams.length > 0) {
        setDropNotice(
          `인덱스(.bai)가 없는 BAM(${matched.unpairedBams.map((f) => f.name).join(', ')})은 제외되었습니다.`,
        );
        setTimeout(() => setDropNotice(null), 5000);
      }
      return;
    }

    if (matched.unpairedBams.length > 0) {
      setDropNotice('.bam 파일과 동일한 이름의 .bai (또는 .csi) 인덱스 파일을 함께 드롭해주세요.');
      setTimeout(() => setDropNotice(null), 6000);
      return;
    }

    if (matched.gffFile || matched.faiFile) {
      setDropNotice(
        '사용자 GFF3 / FAI 파일은 상단 [Genome 변경] 메뉴에서 Reference FASTA와 함께 불러와주세요.',
      );
      setTimeout(() => setDropNotice(null), 6000);
      return;
    }

    if (matched.unpairedIndexes.length > 0) {
      setDropNotice('인덱스(.bai/.csi)만 드롭되었습니다. 짝이 되는 .bam 파일을 함께 드롭해주세요.');
      setTimeout(() => setDropNotice(null), 6000);
      return;
    }

    if (matched.unrecognized.length > 0) {
      setDropError(
        `지원하지 않는 파일 형식입니다 (${matched.unrecognized.map((f) => f.name).join(', ')}). FASTA(.fa/.fasta/.fa.gz), GFF3, BAM+BAI(.bam/.bai/.csi)를 사용해주세요.`,
      );
    }
  };

  const handleChromChange = (nextChrom: string) => {
    if (!workspace) return;
    setScannedVariants([]);
    setActiveSelection(null);
    setHighlightedGene(null);
    if (nextChrom === 'all') {
      void fetchWindow('all', 0, 1);
      return;
    }
    const len = workspace.contigs.find((c) => c.name === nextChrom)?.length ?? 2000;
    const initSpan = Math.min(len, 2_500);
    const initStart = nextChrom === 'chrV' ? Math.min(len - initSpan, 115_600) : 0;
    void fetchWindow(nextChrom, initStart, initStart + initSpan);
  };

  const handleZoom = (factor: number) => {
    if (chrom === 'all') return;
    const span = end - start;
    const center = (start + end) / 2;
    const newSpan = Math.max(30, Math.min(currentChromLength, Math.round(span * factor)));
    let ns = Math.round(center - newSpan / 2);
    ns = Math.max(0, Math.min(currentChromLength - newSpan, ns));
    void fetchWindow(chrom, ns, ns + newSpan);
  };

  const jumpToGene = async (gene: GffFeature) => {
    if (!workspace) return;
    const targetContig = workspace.contigs.find((c) => c.name === gene.chrom);
    const len = targetContig?.length ?? currentChromLength;
    const pad = Math.max(150, Math.round((gene.end - gene.start) * 0.25));
    const s = Math.max(0, gene.start - pad);
    const e = Math.min(len, gene.end + pad);
    if (gene.chrom !== chrom) setScannedVariants([]);
    setHighlightedGene(gene);
    setShowSearchDropdown(false);
    await fetchWindow(gene.chrom, s, e);
  };

  const handleLocusChange = async (val: string) => {
    setLocusInput(val);
    if (!workerRef.current || !workspace || workspace.gffCount === 0) return;
    const trimmed = val.trim();
    // 좌표 형식이면 드롭다운 닫기, 문자면 유전자 자동완성 조회
    if (!trimmed || parseCoordinateInput(trimmed) !== null || /^chr/i.test(trimmed) && trimmed.includes(':')) {
      setShowSearchDropdown(false);
      return;
    }
    const hits = await workerRef.current.searchGenes(trimmed, 12);
    setGeneSuggestions(hits);
    setShowSearchDropdown(hits.length > 0);
  };

  const handleLocusSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!workspace || !workerRef.current) return;
    const parsed = parseCoordinateInput(locusInput);
    if (parsed) {
      const targetChrom =
        parsed.chrom && workspace.contigs.some((c) => c.name === parsed.chrom)
          ? parsed.chrom
          : chrom === 'all'
            ? workspace.contigs[0]?.name ?? 'chrI'
            : chrom;
      const len = workspace.contigs.find((c) => c.name === targetChrom)?.length ?? currentChromLength;
      const s = Math.max(0, Math.min(len - 30, parsed.start));
      const en = Math.min(len, Math.max(s + 30, parsed.end));
      if (targetChrom !== chrom) setScannedVariants([]);
      setHighlightedGene(null);
      setShowSearchDropdown(false);
      await fetchWindow(targetChrom, s, en);
      return;
    }

    if (locusInput.trim().toLowerCase() === 'all') {
      handleChromChange('all');
      return;
    }

    if (workspace.gffCount > 0 && locusInput.trim()) {
      const hits =
        geneSuggestions.length > 0
          ? geneSuggestions
          : await workerRef.current.searchGenes(locusInput.trim(), 1);
      if (hits.length > 0) {
        await jumpToGene(hits[0]);
      }
    }
  };

  const runExtract = async (
    targetChrom: string,
    selStart: number,
    selEnd: number,
    nextOptions = options,
    featureOverride?: GffFeature | null,
  ) => {
    const client = workerRef.current;
    if (!client || targetChrom === 'all') return;

    let resolvedFeature: GffFeature | null =
      featureOverride !== undefined ? featureOverride : null;
    if (featureOverride === undefined) {
      if (
        highlightedGene &&
        highlightedGene.chrom === targetChrom &&
        highlightedGene.start === selStart &&
        highlightedGene.end === selEnd
      ) {
        resolvedFeature = highlightedGene;
      } else {
        resolvedFeature =
          windowDataRef.current?.features.find(
            (f) => f.chrom === targetChrom && f.start === selStart && f.end === selEnd,
          ) ?? null;
      }
    }

    const cdsExons = resolvedFeature
      ? resolvedFeature.subFeatures
          .filter((sf) => sf.type === 'CDS' || sf.type === 'exon')
          .map((sf) => ({ start: sf.start, end: sf.end }))
      : [];

    setExtractFeature(resolvedFeature);
    setActiveSelection({ start: selStart, end: selEnd });
    setIsExtractorOpen(true);
    setIsExtracting(true);
    setExtractError(null);
    try {
      const res = await client.extractRegion({
        chrom: targetChrom,
        start: selStart,
        end: selEnd,
        primaryStrainIndex: primaryIdx,
        secondaryStrainIndex: secondaryIdx,
        options: nextOptions,
        exons: cdsExons.length > 1 ? cdsExons : undefined,
        contigOverrides,
      });
      setExtractResult(res);
    } catch (err) {
      setExtractError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsExtracting(false);
    }
  };

  const handleCancelSelection = () => {
    setActiveSelection(null);
    setExtractFeature(null);
    setIsExtractorOpen(false);
  };

  const handleSelectFeature = (feat: GffFeature) => {
    setHighlightedGene(feat);
    void runExtract(feat.chrom, feat.start, feat.end, options, feat);
  };

  const handleStartScan = async () => {
    const client = workerRef.current;
    if (!client || !chrom || chrom === 'all') return;
    setIsScanning(true);
    setScanProgress(null);
    try {
      const vars = await client.scanChromosome({
        chrom,
        primaryStrainIndex: primaryIdx,
        secondaryStrainIndex: secondaryIdx,
        options,
        contigOverrides,
      });
      setScannedVariants(vars);
    } catch (err) {
      console.warn(err);
    } finally {
      setIsScanning(false);
    }
  };

  const handleCancelScan = async () => {
    await workerRef.current?.cancelScan();
    setIsScanning(false);
  };

  const handleRunBlast = async (payload: BlastRunPayload) => {
    const client = workerRef.current;
    if (!client) return;
    setIsBlasting(true);
    setBlastError(null);
    setBlastProgress({ type: 'blast_progress', message: '로컬 BLAST 준비 중...', percent: 2 });
    try {
      const res = await client.runBlast({
        ...payload,
        contigOverrides,
      });
      setBlastResult(res);
    } catch (err) {
      setBlastError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsBlasting(false);
    }
  };

  const handleCancelBlast = async () => {
    await workerRef.current?.cancelBlast();
    setIsBlasting(false);
  };

  const handleJumpToBlastHit = (hitChrom: string, hitStart: number, hitEnd: number) => {
    const targetLen =
      workspace?.contigs.find((c) => c.name === hitChrom)?.length ?? currentChromLength;
    const span = Math.max(1, hitEnd - hitStart);
    const pad = Math.max(80, Math.round(span * 0.4));
    const viewS = Math.max(0, hitStart - pad);
    const viewE = Math.min(targetLen, hitEnd + pad);
    if (hitChrom !== chrom) setScannedVariants([]);
    setHighlightedGene(null);
    setActiveSelection({ start: hitStart, end: hitEnd });
    setIsBlastOpen(false);
    void fetchWindow(hitChrom, viewS, viewE);
  };

  const handleStartBarcodeScan = async (payload: BarcodeScanPayload) => {
    const client = workerRef.current;
    if (!client) return;
    setIsScanningBarcodes(true);
    setBarcodeError(null);
    setBarcodeProgress({
      type: 'barcode_progress',
      message: '바코드 · 삽입 스카 탐색 준비 중...',
      percent: 2,
      foundSoFar: 0,
    });
    try {
      const res = await client.scanBarcodes({
        ...payload,
        contigOverrides,
      });
      setBarcodeResult(res);
    } catch (err) {
      setBarcodeError(err instanceof Error ? err.message : String(err));
    } finally {
      setIsScanningBarcodes(false);
    }
  };

  const handleCancelBarcodeScan = async () => {
    await workerRef.current?.cancelBarcodeScan();
    setIsScanningBarcodes(false);
  };

  const handleJumpToBarcodeBreakpoint = (bpChrom: string, anchorPos: number) => {
    const targetLen =
      workspace?.contigs.find((c) => c.name === bpChrom)?.length ?? currentChromLength;
    const viewS = Math.max(0, anchorPos - 150);
    const viewE = Math.min(targetLen, anchorPos + 150);
    if (bpChrom !== chrom) setScannedVariants([]);
    setHighlightedGene(null);
    setActiveSelection({ start: Math.max(0, anchorPos - 10), end: Math.min(targetLen, anchorPos + 10) });
    setIsBarcodeOpen(false);
    void fetchWindow(bpChrom, viewS, viewE);
  };

  const hasStrain2 =
    secondaryIdx !== null &&
    workspace !== null &&
    secondaryIdx >= 0 &&
    secondaryIdx < workspace.strains.length;

  return (
    <div
      ref={appRootRef}
      className="relative flex h-screen flex-col overflow-hidden bg-white text-slate-900"
      onDragOver={(e) => {
        if (e.dataTransfer.types.includes('Files')) {
          e.preventDefault();
          setIsWindowDragOver(true);
        }
      }}
      onDragLeave={(e) => {
        if (e.currentTarget === e.target) setIsWindowDragOver(false);
      }}
      onDrop={(e) => void handleGlobalFileDrop(e)}
    >
      {/* Row 1: IGV Classic Top Menu Bar */}
      <header className="flex flex-wrap items-center justify-between gap-2 border-b border-slate-300 bg-slate-100 px-3 py-1.5 text-xs">
        <div className="flex flex-wrap items-center gap-2">
          <div className="mr-1 flex items-center gap-1.5 font-extrabold tracking-tight text-slate-800">
            <Dna className="h-4 w-4 text-sky-700" />
            <span>YeastWGS</span>
          </div>

          {/* Genome Menu Pill */}
          <div className="flex items-center overflow-hidden rounded border border-slate-300 bg-white">
            <span className="border-r border-slate-200 bg-slate-50 px-2 py-1 font-semibold text-slate-600">
              Genome
            </span>
            <span className="px-2 py-1 font-mono font-bold text-slate-900">
              {workspace?.genomeId === 'sacCer3'
                ? 'sacCer3'
                : workspace?.isDemo
                  ? 'Demo'
                  : workspace?.referenceLabel.slice(0, 18) ?? 'Loading...'}
            </span>
            <button
              type="button"
              onClick={() => {
                setLoaderMode('genome');
                setIsLoaderOpen(true);
              }}
              className="border-l border-slate-200 bg-slate-50 px-2 py-1 text-[11px] font-medium text-sky-700 hover:bg-sky-50"
              title="다른 FASTA 레퍼런스로 변경 또는 sacCer3 복귀"
            >
              변경
            </button>
          </div>

          {/* Tracks (+ BAM/BAI) */}
          <button
            type="button"
            onClick={() => {
              setLoaderMode('tracks');
              setIsLoaderOpen(true);
            }}
            className="inline-flex items-center gap-1 rounded border border-sky-600 bg-sky-600 px-2.5 py-1 font-semibold text-white hover:bg-sky-700"
          >
            <Plus className="h-3.5 w-3.5" />
            BAM/BAI 트랙 추가
          </button>

          {/* Loaded Track Chips */}
          {workspace && workspace.strains.length > 0 && (
            <div className="flex flex-wrap items-center gap-1">
              {workspace.strains.map((s, idx) => (
                <span
                  key={s.id}
                  className="inline-flex items-center gap-1 rounded border border-slate-300 bg-white px-2 py-0.5 text-[11px] font-medium text-slate-800"
                >
                  <span
                    className={`h-2 w-2 rounded-full ${
                      idx === primaryIdx
                        ? 'bg-sky-600'
                        : idx === secondaryIdx
                          ? 'bg-emerald-600'
                          : 'bg-slate-400'
                    }`}
                  />
                  {s.label}
                  <button
                    type="button"
                    onClick={() => void handleRemoveStrain(s.id)}
                    className="ml-0.5 rounded text-slate-400 hover:text-rose-600"
                    title="트랙 제거"
                  >
                    <X className="h-3 w-3" />
                  </button>
                </span>
              ))}
            </div>
          )}
        </div>

        {/* Right IGV Tools */}
        <div className="flex flex-wrap items-center gap-1.5">
          {/* Pan vs Drag Extract */}
          <div className="flex overflow-hidden rounded border border-slate-300 bg-white">
            <button
              type="button"
              onClick={() => setInteractionMode('pan')}
              className={`inline-flex items-center gap-1 px-2 py-1 font-medium ${
                interactionMode === 'pan'
                  ? 'bg-slate-800 text-white'
                  : 'text-slate-700 hover:bg-slate-100'
              }`}
              title="마우스 드래그/휠 좌우 이동"
            >
              <Hand className="h-3 w-3" /> 이동
            </button>
            <button
              type="button"
              onClick={() => setInteractionMode('select')}
              className={`inline-flex items-center gap-1 border-l border-slate-200 px-2 py-1 font-medium ${
                interactionMode === 'select'
                  ? 'bg-sky-600 text-white'
                  : 'text-slate-700 hover:bg-slate-100'
              }`}
              title="드래그해 선택한 구간의 실제 균주 서열 추출 (Shift+드래그 가능)"
            >
              <Scissors className="h-3 w-3" /> 드래그 서열 추출
            </button>
          </div>

          {activeSelection && (
            <div className="flex items-center overflow-hidden rounded border border-emerald-600 bg-white">
              <button
                type="button"
                onClick={() => void runExtract(chrom, activeSelection.start, activeSelection.end)}
                className="bg-emerald-600 px-2 py-1 font-semibold text-white hover:bg-emerald-700"
              >
                추출 ({(activeSelection.end - activeSelection.start).toLocaleString()} bp)
              </button>
              <button
                type="button"
                onClick={handleCancelSelection}
                className="px-1.5 py-1 text-slate-600 hover:bg-rose-50 hover:text-rose-700"
                title="선택 해제 (Esc)"
              >
                <X className="h-3.5 w-3.5" />
              </button>
            </div>
          )}

          <button
            type="button"
            onClick={() => setShowCrosshairs((v) => !v)}
            className={`inline-flex items-center gap-1 rounded border px-2 py-1 font-medium ${
              showCrosshairs
                ? 'border-sky-500 bg-sky-100 text-sky-900'
                : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'
            }`}
          >
            <Crosshair className="h-3 w-3" /> Crosshairs
          </button>

          <button
            type="button"
            onClick={() => setShowCenterLine((v) => !v)}
            className={`rounded border px-2 py-1 font-medium ${
              showCenterLine
                ? 'border-sky-500 bg-sky-100 text-sky-900'
                : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'
            }`}
          >
            Center Line
          </button>

          <button
            type="button"
            onClick={() => setShowCompareStrip((v) => !v)}
            className={`inline-flex items-center gap-1 rounded border px-2 py-1 font-medium ${
              showCompareStrip
                ? 'border-sky-500 bg-sky-100 text-sky-900'
                : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'
            }`}
            title="Reference vs 균주 Consensus 비교 띠 표시/숨김"
          >
            <SlidersHorizontal className="h-3 w-3" /> {hasStrain2 ? '3-Way 비교 띠' : '2-Way 비교 띠'}
          </button>

          <button
            type="button"
            onClick={() => setCompactMode((v) => !v)}
            className={`inline-flex items-center gap-1 rounded border px-2 py-1 font-medium ${
              compactMode
                ? 'border-amber-500 bg-amber-100 text-amber-900'
                : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'
            }`}
            title="리드 접기/펼치기"
          >
            <AlignJustify className="h-3 w-3" /> 리드 접기
          </button>

          <button
            type="button"
            onClick={() => setShowAllBases((v) => !v)}
            className={`rounded border px-2 py-1 font-medium ${
              showAllBases
                ? 'border-sky-500 bg-sky-100 text-sky-900'
                : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'
            }`}
          >
            전체 염기
          </button>

          <button
            type="button"
            onClick={() => setIsBlastOpen(true)}
            className="inline-flex items-center gap-1 rounded border border-indigo-600 bg-indigo-600 px-2.5 py-1 font-semibold text-white hover:bg-indigo-700"
            title="입력한 DNA/프라이머 서열을 Reference와 로드된 BAM 균주에 동시 BLAST 매핑"
          >
            <Dna className="h-3.5 w-3.5" />
            로컬 BLAST
          </button>

          <button
            type="button"
            onClick={() => setIsBarcodeOpen(true)}
            className="inline-flex items-center gap-1 rounded border border-amber-600 bg-amber-600 px-2.5 py-1 font-semibold text-white hover:bg-amber-700"
            title="외부 삽입 바코드(UPTAG/DNTAG·loxP·FRT) 및 미지 워터마크 스카 자동 스캔"
          >
            <QrCode className="h-3.5 w-3.5" />
            바코드 스캔
          </button>

          <button
            type="button"
            onClick={() => setIsSidebarOpen((v) => !v)}
            className={`inline-flex items-center gap-1 rounded border px-2.5 py-1 font-semibold ${
              isSidebarOpen
                ? 'border-sky-600 bg-sky-600 text-white'
                : 'border-slate-300 bg-white text-slate-800 hover:bg-slate-50'
            }`}
          >
            <PanelRight className="h-3.5 w-3.5" />
            변이·매핑 패널
          </button>
        </div>
      </header>

      {/* Row 2: IGV Location & Unified Gene/Locus Search Bar */}
      {workspace && (
        <div className="flex flex-wrap items-center justify-between gap-2 border-b border-slate-300 bg-white px-3 py-1.5 text-xs">
          <div className="flex flex-wrap items-center gap-2">
            {/* Chromosome Selector starting with 'all' + chrI..chrXVI + chrM */}
            <select
              value={chrom}
              onChange={(e) => handleChromChange(e.target.value)}
              className="w-36 rounded border border-slate-300 bg-slate-50 px-2 py-1 font-sans text-xs font-semibold text-slate-900 focus:border-sky-500 focus:outline-none"
            >
              <option value="all">all (전체 염색체)</option>
              {workspace.contigs.map((c) => (
                <option key={c.name} value={c.name}>
                  {c.name}
                </option>
              ))}
            </select>

            {/* Unified Locus / Gene Search Input (IGV style) */}
            <div ref={searchDropdownRef} className="relative">
              <form onSubmit={handleLocusSubmit} className="flex items-center gap-1.5">
                <div className="relative">
                  <input
                    type="text"
                    value={locusInput}
                    onChange={(e) => void handleLocusChange(e.target.value)}
                    placeholder="좌표(chrV:115,600-117,600) 또는 유전자명(URA3, YER134C)"
                    className="w-72 rounded border border-slate-300 bg-white py-1 pl-2.5 pr-7 font-mono text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                  />
                  <button
                    type="submit"
                    className="absolute right-1.5 top-1 text-slate-500 hover:text-sky-700"
                    title="이동 / 유전자 검색"
                  >
                    <Search className="h-3.5 w-3.5" />
                  </button>
                </div>
              </form>

              {showSearchDropdown && (
                <div className="absolute left-0 top-full z-40 mt-1 max-h-72 w-84 overflow-y-auto rounded-xl border border-slate-300 bg-white p-1.5 shadow-xl">
                  {geneSuggestions.map((g) => (
                    <div
                      key={g.id}
                      className="flex items-center justify-between gap-2 rounded-lg px-2.5 py-2 hover:bg-sky-50"
                    >
                      <button
                        type="button"
                        onClick={() => void jumpToGene(g)}
                        className="min-w-0 flex-1 text-left"
                      >
                        <div className="flex items-center gap-1.5 text-xs font-bold text-slate-900">
                          <span className="text-blue-800">{g.name}</span>
                          {g.systematicName && (
                            <span className="font-mono text-[10px] font-normal text-slate-500">
                              ({g.systematicName})
                            </span>
                          )}
                        </div>
                        <div className="font-mono text-[10px] text-slate-600">
                          {g.chrom}:{(g.start + 1).toLocaleString()}–{g.end.toLocaleString()} ({g.strand})
                        </div>
                      </button>
                      <button
                        type="button"
                        onClick={() => {
                          setShowSearchDropdown(false);
                          void jumpToGene(g);
                          void runExtract(g.chrom, g.start, g.end, options, g);
                        }}
                        className="shrink-0 rounded bg-sky-100 px-2 py-1 text-[10px] font-semibold text-sky-800 hover:bg-sky-200"
                      >
                        추출
                      </button>
                    </div>
                  ))}
                </div>
              )}
            </div>

            {chrom !== 'all' && (
              <span className="font-mono text-xs text-slate-600">
                {(end - start).toLocaleString()} bp
              </span>
            )}
          </div>

          <div className="flex items-center gap-3">
            <span className="hidden text-[11px] text-slate-500 lg:inline">
              {interactionMode === 'select'
                ? '드래그로 서열 선택 중 (우클릭/Esc 취소)'
                : '마우스 휠: 좌우 이동 · Ctrl+휠: 확대/축소 · Shift+드래그: 서열 추출'}
            </span>

            {chrom !== 'all' && (
              <div className="flex items-center overflow-hidden rounded border border-slate-300 bg-white">
                <button
                  type="button"
                  onClick={() => handleZoom(0.5)}
                  className="p-1 text-slate-700 hover:bg-slate-100"
                  title="확대 (Ctrl+휠)"
                >
                  <ZoomIn className="h-3.5 w-3.5" />
                </button>
                <button
                  type="button"
                  onClick={() => handleZoom(2)}
                  className="border-l border-slate-200 p-1 text-slate-700 hover:bg-slate-100"
                  title="축소 (Ctrl+휠)"
                >
                  <ZoomOut className="h-3.5 w-3.5" />
                </button>
                <button
                  type="button"
                  onClick={() => void fetchWindow(chrom, 0, currentChromLength)}
                  className="border-l border-slate-200 p-1 text-slate-700 hover:bg-slate-100"
                  title="현재 염색체 전체"
                >
                  <Maximize2 className="h-3.5 w-3.5" />
                </button>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Row 3: IGV Chromosome Ideogram */}
      {workspace && (
        <ChromosomeIdeogram
          chrom={chrom}
          chromLength={currentChromLength}
          start={start}
          end={end}
          contigs={workspace.contigs}
          onJumpToPosition={(ns, ne) => void fetchWindow(chrom, ns, ne)}
        />
      )}

      {isLoadingFiles && (
        <div className="flex items-center gap-2 border-b border-amber-200 bg-amber-50 px-4 py-1.5 text-xs font-semibold text-amber-900">
          <span className="inline-block h-2.5 w-2.5 animate-pulse rounded-full bg-amber-500" />
          <span>{loadingMessage || '파일 불러오는 중...'}</span>
        </div>
      )}

      {dropError && (
        <div className="flex items-center justify-between gap-2 border-b border-rose-200 bg-rose-50 px-4 py-1.5 text-xs font-medium text-rose-900">
          <span>오류: {dropError}</span>
          <button
            type="button"
            onClick={() => setDropError(null)}
            className="rounded p-0.5 text-rose-700 hover:bg-rose-100"
          >
            <X className="h-3.5 w-3.5" />
          </button>
        </div>
      )}

      {dropNotice && (
        <div className="flex items-center justify-between gap-2 border-b border-sky-200 bg-sky-50 px-4 py-1.5 text-xs font-medium text-sky-900">
          <span>{dropNotice}</span>
          <div className="flex items-center gap-2">
            <button
              type="button"
              onClick={() => setIsSidebarOpen(true)}
              className="rounded border border-sky-300 bg-white px-2 py-0.5 text-[11px] font-semibold text-sky-800 hover:bg-sky-100"
            >
              매핑 패널 열기
            </button>
            <button
              type="button"
              onClick={() => setDropNotice(null)}
              className="rounded p-0.5 text-sky-700 hover:bg-sky-100"
            >
              <X className="h-3.5 w-3.5" />
            </button>
          </div>
        </div>
      )}

      {chrom !== 'all' &&
        windowData &&
        windowData.strains.some((s) => s.matchedChrom === null) && (
          <div className="flex flex-wrap items-center justify-between gap-2 border-b border-amber-200 bg-amber-50 px-4 py-1.5 text-xs text-amber-900">
            <span>
              현재 염색체(<span className="font-mono font-bold">{chrom}</span>)와 자동 매칭되지 않은 BAM 트랙:{' '}
              <span className="font-semibold">
                {windowData.strains
                  .filter((s) => s.matchedChrom === null)
                  .map((s) => s.label)
                  .join(', ')}
              </span>
            </span>
            <button
              type="button"
              onClick={() => setIsSidebarOpen(true)}
              className="rounded border border-amber-300 bg-white px-2 py-0.5 text-[11px] font-semibold text-amber-900 hover:bg-amber-100"
            >
              [변이·매핑 패널]에서 컨티그 직접 연결
            </button>
          </div>
        )}

      {/* Main Content Area */}
      <div className="flex flex-1 overflow-hidden">
        <main className="flex flex-1 flex-col overflow-hidden">
          {chrom === 'all' ? (
            <AllChromsOverviewCanvas
              data={allChromsData}
              onSelectChrom={(c) => handleChromChange(c)}
            />
          ) : (
            <GenomeCanvasViewer
              data={windowData}
              viewStart={start}
              viewEnd={end}
              navVersion={navVersion}
              primaryStrainIndex={primaryIdx}
              secondaryStrainIndex={secondaryIdx}
              interactionMode={interactionMode}
              compactMode={compactMode}
              showAllBases={showAllBases}
              showCompareStrip={showCompareStrip && (workspace?.strains.length ?? 0) > 0}
              showCrosshairs={showCrosshairs}
              showCenterLine={showCenterLine}
              activeSelection={activeSelection}
              highlightedGene={highlightedGene}
              scannedVariants={scannedVariants}
              chromLength={currentChromLength}
              onViewportChange={handleViewerViewportChange}
              onSelectRegion={(selS, selE) => void runExtract(chrom, selS, selE)}
              onCancelSelection={handleCancelSelection}
              onSelectFeature={handleSelectFeature}
            />
          )}
        </main>

        {workspace && isSidebarOpen && (
          <VariantSidebar
            chrom={chrom === 'all' ? workspace.contigs[0]?.name ?? 'chrI' : chrom}
            windowVariants={windowData?.analysis?.variants ?? []}
            scannedVariants={scannedVariants}
            isScanning={isScanning}
            scanProgress={scanProgress}
            hasGff={workspace.gffCount > 0}
            annotationOrigin={workspace.annotationOrigin}
            strains={workspace.strains}
            primaryStrainIndex={primaryIdx}
            secondaryStrainIndex={secondaryIdx}
            options={options}
            contigOverrides={contigOverrides}
            refContigs={workspace.contigs}
            onStartScan={() => void handleStartScan()}
            onCancelScan={() => void handleCancelScan()}
            onJumpToRange={(targetChr, s, e) => {
              if (targetChr !== chrom) setScannedVariants([]);
              void fetchWindow(targetChr, s, e);
            }}
            onJumpToChrom={handleChromChange}
            onExtractFeature={handleSelectFeature}
            onSearchGenes={async (q) => (await workerRef.current?.searchGenes(q, 25)) ?? []}
            onChangePrimaryIndex={setPrimaryIdx}
            onChangeSecondaryIndex={setSecondaryIdx}
            onChangeOptions={setOptions}
            onChangeContigOverride={(strainId, refChr, bamChr) => {
              setContigOverrides((prev) => ({
                ...prev,
                [strainId]: {
                  ...(prev[strainId] ?? {}),
                  [refChr]: bamChr,
                },
              }));
            }}
          />
        )}
      </div>

      {/* Global Drag-Over Overlay */}
      {isWindowDragOver && (
        <div className="pointer-events-none fixed inset-0 z-50 flex items-center justify-center bg-sky-600/15 backdrop-blur-xs">
          <div className="rounded-2xl border-2 border-dashed border-sky-600 bg-white px-8 py-6 text-center shadow-2xl">
            <FolderOpen className="mx-auto mb-2 h-8 w-8 text-sky-600" />
            <p className="text-sm font-bold text-slate-900">
              여기에 BAM + BAI (또는 FASTA) 파일을 놓으면 바로 불러옵니다
            </p>
          </div>
        </div>
      )}

      {/* Modals */}
      <FileLoaderModal
        isOpen={isLoaderOpen}
        mode={loaderMode}
        activeGenomeLabel={workspace?.referenceLabel ?? 'sacCer3'}
        isLoading={isLoadingFiles}
        loadingMsg={loadingMessage}
        errorMsg={loaderError}
        onClose={() => setIsLoaderOpen(false)}
        onAddTracksOnly={(s) => handleAddTracksOnly(s)}
        onLoadCustomGenome={(p) => handleLoadCustomGenome(p)}
        onSelectSacCer3={() => void handleSelectSacCer3()}
        onLoadDemo={() => void handleLoadDemo()}
      />

      <SequenceExtractorModal
        isOpen={isExtractorOpen}
        isLoading={isExtracting}
        errorMsg={extractError}
        result={extractResult}
        selectedFeature={extractFeature}
        chromLength={currentChromLength}
        options={options}
        onClose={() => setIsExtractorOpen(false)}
        onClearAndClose={handleCancelSelection}
        onSendToLocalBlast={(qFasta) => {
          setIsExtractorOpen(false);
          setBlastInitialQuery(qFasta);
          setIsBlastOpen(true);
        }}
        onRecompute={(c, s, e, nextOpts) => {
          setOptions(nextOpts);
          const keepFeat =
            extractFeature &&
            extractFeature.chrom === c &&
            s <= extractFeature.start &&
            e >= extractFeature.end
              ? extractFeature
              : null;
          void runExtract(c, s, e, nextOpts, keepFeat);
        }}
      />

      <BlastModal
        isOpen={isBlastOpen}
        isRunning={isBlasting}
        progress={blastProgress}
        errorMsg={blastError}
        result={blastResult}
        initialQuery={blastInitialQuery}
        referenceLabel={workspace?.referenceLabel ?? 'sacCer3'}
        strainLabels={workspace?.strains.map((s) => s.label) ?? []}
        currentChrom={chrom}
        options={options}
        onClose={() => setIsBlastOpen(false)}
        onRunBlast={(p) => void handleRunBlast(p)}
        onCancelBlast={() => void handleCancelBlast()}
        onJumpToHit={handleJumpToBlastHit}
      />

      <BarcodeScannerModal
        isOpen={isBarcodeOpen}
        isScanning={isScanningBarcodes}
        progress={barcodeProgress}
        errorMsg={barcodeError}
        result={barcodeResult}
        strainLabels={workspace?.strains.map((s) => s.label) ?? []}
        currentChrom={chrom}
        onClose={() => setIsBarcodeOpen(false)}
        onStartScan={(p) => void handleStartBarcodeScan(p)}
        onCancelScan={() => void handleCancelBarcodeScan()}
        onJumpToBreakpoint={handleJumpToBarcodeBreakpoint}
        onSendToBlast={(qFasta) => {
          setIsBarcodeOpen(false);
          setBlastInitialQuery(qFasta);
          setIsBlastOpen(true);
        }}
      />
    </div>
  );
}
