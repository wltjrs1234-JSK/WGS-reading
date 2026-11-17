import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import {
  expandReadAlignment,
  type ReadInsertionMarker,
  type ReadSegment,
} from '../bio/cigar';
import { packGeneFeatures, type PackedGene } from '../bio/geneLayout';
import {
  MAX_VISIBLE_PILEUP_ROWS,
  packReadsWithHeap,
  type AlleleTick,
  type PackedRead,
} from '../bio/readPacking';
import type {
  AlignedRead,
  BaseCounts,
  DiffCategory,
  GffFeature,
  VariantRecord,
} from '../bio/types';
import type { WindowQueryResult } from '../workers/wgsProtocol';
import {
  FONT_CHANGE_EVENT,
  formatCanvasFont,
  getActiveCanvasFontConfig,
} from '../utils/fontSettings';

export type InteractionMode = 'pan' | 'select';

interface PreparedBase {
  refPos: number;
  base: string;
  qual?: number;
}

interface PreparedRead {
  read: PackedRead;
  segments: ReadSegment[];
  insertions: ReadInsertionMarker[];
  mismatches: PreparedBase[];
  allBases: PreparedBase[];
}

interface StrainLayout {
  id: string;
  label: string;
  matchedChrom: string | null;
  meanDepth: number;
  overviewBins: number[];
  maxBinDepth: number;
  preparedReads: PreparedRead[];
  readsByRow: Map<number, PackedRead[]>;
  totalRows: number;
  covY: number;
  covH: number;
  pileupY: number;
  pileupH: number;
  alleleTicks: AlleleTick[] | null;
}

interface ViewerLayout {
  rulerY: number;
  refY: number;
  strains: StrainLayout[];
  emptyHintY: number;
  showCompareTrack: boolean;
  compareY: number;
  compareH: number;
  hasFeatures: boolean;
  packedGenes: { packed: PackedGene[]; totalRows: number };
  geneY: number;
  geneH: number;
  totalHeight: number;
}

interface GenomeCanvasViewerProps {
  data: WindowQueryResult | null;
  viewStart: number;
  viewEnd: number;
  navVersion: number;
  primaryStrainIndex: number;
  secondaryStrainIndex: number | null;
  interactionMode: InteractionMode;
  compactMode: boolean;
  showAllBases: boolean;
  showCompareStrip: boolean;
  showCrosshairs: boolean;
  showCenterLine: boolean;
  activeSelection: { start: number; end: number } | null;
  highlightedGene: GffFeature | null;
  scannedVariants: VariantRecord[];
  chromLength: number;
  onViewportChange: (newStart: number, newEnd: number, urgentBufferCheck?: boolean) => void;
  onSelectRegion: (start: number, end: number) => void;
  onCancelSelection: () => void;
  onSelectFeature: (feat: GffFeature) => void;
}

interface HoverInfo {
  x: number;
  y: number;
  pos: number;
  refBase?: string;
  s1Call?: string;
  s1Depth?: number;
  s2Call?: string;
  s2Depth?: number;
  hasStrain2?: boolean;
  read?: AlignedRead;
  feature?: GffFeature;
}

type DragState =
  | {
      mode: 'pan';
      startX: number;
      origStart: number;
      origEnd: number;
    }
  | {
      mode: 'select';
      anchorX: number;
      anchorPos: number;
      currentPos: number;
    };

const IGV_BASE_COLORS: Record<string, string> = {
  A: '#009e73', // green
  C: '#1d4ed8', // blue
  G: '#d97706', // amber/brown-orange
  T: '#dc2626', // red
  N: '#94a3b8',
  '-': '#0f172a',
};

const DIFF_BG_COLORS: Record<DiffCategory, string> = {
  match: 'transparent',
  strain1_only: '#e0f2fe',
  strain2_only: '#d1fae5',
  shared: '#ffedd5',
  divergent: '#ffe4e6',
  gap_s1: '#f1f5f9',
  gap_s2: '#f1f5f9',
  gap_both: '#e2e8f0',
};

const DIFF_MARKER_COLORS: Record<DiffCategory, string> = {
  match: 'transparent',
  strain1_only: '#0284c7',
  strain2_only: '#059669',
  shared: '#ea580c',
  divergent: '#e11d48',
  gap_s1: '#64748b',
  gap_s2: '#64748b',
  gap_both: '#475569',
};

const LEFT_LABEL_W = 112;
const RULER_H = 26;
const REF_TRACK_H = 18;
const COV_TRACK_H = 50;
const READ_ROW_H = 9;
const READ_ROW_GAP = 2;
const MAX_VISIBLE_ROWS = MAX_VISIBLE_PILEUP_ROWS;
const COMPARE_TWO_WAY_H = 38;
const COMPARE_THREE_WAY_H = 56;
const GENE_LANE_H = 34;
const GENE_MIN_H = 48;

function niceTickStep(span: number, targetTicks: number): number {
  const raw = span / Math.max(1, targetTicks);
  const mag = Math.pow(10, Math.floor(Math.log10(Math.max(1, raw))));
  const norm = raw / mag;
  if (norm < 1.5) return mag;
  if (norm < 3.5) return 2 * mag;
  if (norm < 7.5) return 5 * mag;
  return 10 * mag;
}

/**
 * IGV 스타일 커버리지 위 변이 비율 틱을 그리기 위한 포지션별 대립유전자 집계
 */
function computeAlleleTicks(
  refSeq: string,
  reads: AlignedRead[],
  start: number,
  end: number,
): Array<{ pos: number; depth: number; counts: BaseCounts }> {
  const span = Math.max(0, end - start);
  const countsArr: BaseCounts[] = Array.from({ length: span }, () => ({
    A: 0,
    C: 0,
    G: 0,
    T: 0,
    N: 0,
    del: 0,
  }));

  for (const r of reads) {
    if ((r.flags & 0x704) !== 0) continue;
    const exp = expandReadAlignment(r, start, end);
    for (const b of exp.bases) {
      if (b.qual !== undefined && b.qual < 255 && b.qual < 5) continue;
      const idx = b.refPos - start;
      if (idx < 0 || idx >= span) continue;
      const c = countsArr[idx];
      if (b.base === 'A') c.A++;
      else if (b.base === 'C') c.C++;
      else if (b.base === 'G') c.G++;
      else if (b.base === 'T') c.T++;
      else c.N++;
    }
  }

  const ticks: Array<{ pos: number; depth: number; counts: BaseCounts }> = [];
  for (let i = 0; i < span; i++) {
    const c = countsArr[i];
    const depth = c.A + c.C + c.G + c.T + c.N;
    if (depth < 2) continue;
    const refBase = (refSeq[i] ?? 'N').toUpperCase();
    const nonRef =
      depth -
      (refBase === 'A'
        ? c.A
        : refBase === 'C'
          ? c.C
          : refBase === 'G'
            ? c.G
            : refBase === 'T'
              ? c.T
              : 0);
    // IGV 기본 allele frequency 임계치(0.20) 일치
    if (nonRef >= 2 && nonRef / depth >= 0.2) {
      ticks.push({ pos: start + i, depth, counts: c });
    }
  }
  return ticks;
}

function checkNearBufferEdge(
  vStart: number,
  vEnd: number,
  data: WindowQueryResult,
  chromLength: number,
): boolean {
  const vSpan = Math.max(1, vEnd - vStart);
  const wantOverview = vSpan > 30_000;
  if (data.isOverviewOnly !== wantOverview) return true;
  if (vStart < data.start + vSpan * 0.2 && data.start > 0) return true;
  if (vEnd > data.end - vSpan * 0.2 && data.end < chromLength) return true;
  if (wantOverview) {
    const dataSpan = Math.max(1, data.end - data.start);
    if (vSpan < dataSpan * 0.5 || vSpan > dataSpan * 1.15) return true;
  }
  return false;
}

export const GenomeCanvasViewer: React.FC<GenomeCanvasViewerProps> = ({
  data,
  viewStart,
  viewEnd,
  navVersion,
  primaryStrainIndex,
  secondaryStrainIndex,
  interactionMode,
  compactMode,
  showAllBases,
  showCompareStrip,
  showCrosshairs,
  showCenterLine,
  activeSelection,
  highlightedGene,
  scannedVariants,
  chromLength,
  onViewportChange,
  onSelectRegion,
  onCancelSelection,
  onSelectFeature,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const measureCtxRef = useRef<CanvasRenderingContext2D | null>(null);
  const [size, setSize] = useState({ width: 960, height: 560 });

  // 60fps 클라이언트 뷰포트(부동소수점)와 목표 뷰포트(휠 관성 보간용)
  const liveViewRef = useRef<{ chrom: string; start: number; end: number }>({
    chrom: data?.chrom ?? '',
    start: viewStart,
    end: viewEnd,
  });
  const targetViewRef = useRef<{ chrom: string; start: number; end: number }>({
    chrom: data?.chrom ?? '',
    start: viewStart,
    end: viewEnd,
  });
  const lastSyncedNavRef = useRef<{ navVersion: number; chrom: string }>({
    navVersion: -1,
    chrom: '',
  });

  const dragRef = useRef<DragState | null>(null);
  const hoverRef = useRef<HoverInfo | null>(null);
  const [hoverInfo, setHoverInfo] = useState<HoverInfo | null>(null);

  const drawRafRef = useRef<number | null>(null);
  const glideRafRef = useRef<number | null>(null);
  const syncTimerRef = useRef<number | null>(null);
  const lastSyncAtRef = useRef<number>(0);
  const lastHoverCommitRef = useRef<number>(0);

  const getMeasureCtx = (): CanvasRenderingContext2D | null => {
    const fontCfg = getActiveCanvasFontConfig();
    if (!measureCtxRef.current && typeof document !== 'undefined') {
      const c = document.createElement('canvas');
      const mctx = c.getContext('2d');
      if (mctx) {
        measureCtxRef.current = mctx;
      }
    }
    if (measureCtxRef.current) {
      measureCtxRef.current.font = formatCanvasFont(
        fontCfg.uiFontFamily,
        10.5 + fontCfg.canvasDeltaPx,
        600,
      );
    }
    return measureCtxRef.current;
  };

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setSize({
          width: Math.max(360, Math.floor(entry.contentRect.width)),
          height: Math.max(320, Math.floor(entry.contentRect.height)),
        });
      }
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  const hasStrain2 =
    secondaryStrainIndex !== null &&
    data !== null &&
    secondaryStrainIndex >= 0 &&
    secondaryStrainIndex < data.strains.length;

  const compareTrackHeight = hasStrain2 ? COMPARE_THREE_WAY_H : COMPARE_TWO_WAY_H;

  // 유전자 다중 레인 배치 (IGV RefSeq처럼 Dubious ORF는 검색 하이라이트 시에만 표시, 좌우 스크롤 중 레인 안정 유지)
  const viewSpanRounded = Math.max(30, Math.round(viewEnd - viewStart));
  const packedGenes = useMemo(() => {
    if (!data || data.features.length === 0) {
      return { packed: [] as PackedGene[], totalRows: 0 };
    }
    const visibleFeatures = data.features.filter(
      (f) => f.type !== 'dubious_gene' || (highlightedGene !== null && f.id === highlightedGene.id),
    );
    if (visibleFeatures.length === 0) {
      return { packed: [] as PackedGene[], totalRows: 0 };
    }
    const trackW = Math.max(100, size.width - LEFT_LABEL_W);
    const pxPerBp = trackW / viewSpanRounded;
    const mctx = getMeasureCtx();
    return packGeneFeatures(
      visibleFeatures,
      pxPerBp,
      mctx ? (name) => mctx.measureText(name).width : undefined,
    );
  }, [data, size.width, viewSpanRounded, highlightedGene]);

  const layout = useMemo<ViewerLayout | null>(() => {
    if (!data) return null;
    let y = 0;
    const rulerY = y;
    y += RULER_H;

    const refY = y;
    y += REF_TRACK_H;

    const strains: StrainLayout[] = [];
    for (let idx = 0; idx < data.strains.length; idx++) {
      const s = data.strains[idx];
      const covY = y;
      y += COV_TRACK_H;

      const showPileup = !compactMode && !data.isOverviewOnly;
      const packedResult = showPileup
        ? s.packedReads
          ? { visiblePacked: s.packedReads, totalRows: s.totalRows ?? s.packedReads.length }
          : packReadsWithHeap(s.reads, MAX_VISIBLE_ROWS)
        : { visiblePacked: [] as PackedRead[], totalRows: s.totalRows ?? 0 };
      const visibleRows = Math.min(
        MAX_VISIBLE_ROWS,
        Math.max(packedResult.totalRows, showPileup ? 4 : 0),
      );
      const pileupY = y;
      const pileupH = showPileup ? visibleRows * (READ_ROW_H + READ_ROW_GAP) + 10 : 0;
      y += pileupH + 4;

      // 가시 행(최대 28줄)만 CIGAR 전개 결과를 1회 사전 계산하여 60fps 렌더 중 재파싱 제거
      const preparedReads: PreparedRead[] = [];
      const readsByRow = new Map<number, PackedRead[]>();
      if (showPileup) {
        const refSeq = data.refSeq;
        const dStart = data.start;
        for (const r of packedResult.visiblePacked) {
          if (r.row >= MAX_VISIBLE_ROWS) continue;
          const exp = expandReadAlignment(r, data.start, data.end);
          const mismatches: PreparedBase[] = [];
          const normalizedBases: PreparedBase[] = exp.bases.map((b) => {
            if (!refSeq) return b;
            const refIdx = b.refPos - dStart;
            const refBase = refIdx >= 0 && refIdx < refSeq.length ? refSeq[refIdx] : 'N';
            const resolvedBase = b.base === '=' ? refBase : b.base;
            const norm = resolvedBase === b.base ? b : { ...b, base: resolvedBase };
            if (resolvedBase !== refBase) {
              mismatches.push(norm);
            }
            return norm;
          });
          preparedReads.push({
            read: r,
            segments: exp.segments,
            insertions: exp.insertions,
            mismatches,
            allBases: normalizedBases,
          });
          let bucket = readsByRow.get(r.row);
          if (!bucket) {
            bucket = [];
            readsByRow.set(r.row, bucket);
          }
          bucket.push(r);
        }
      }

      const alleleTicks =
        !data.isOverviewOnly && s.alleleTicks
          ? s.alleleTicks
          : !data.isOverviewOnly && data.refSeq && s.reads.length > 0
            ? computeAlleleTicks(data.refSeq, s.reads, data.start, data.end)
            : null;

      const maxBinDepth = Math.max(5, Math.ceil(Math.max(...s.overviewBins, 1)));
      strains.push({
        id: s.id,
        label: s.label,
        matchedChrom: s.matchedChrom,
        meanDepth: s.meanDepth,
        overviewBins: s.overviewBins,
        maxBinDepth,
        preparedReads,
        readsByRow,
        totalRows: packedResult.totalRows,
        covY,
        covH: COV_TRACK_H,
        pileupY,
        pileupH,
        alleleTicks,
      });
    }

    const emptyHintY = y;
    if (data.strains.length === 0) {
      y += 96;
    }

    const hasCompare =
      showCompareStrip &&
      data.strains.length > 0 &&
      ((!data.isOverviewOnly && data.analysis !== null) ||
        (data.isOverviewOnly && scannedVariants.length > 0));
    const compareY = y;
    if (hasCompare) {
      y += compareTrackHeight;
    }

    const hasFeatures = packedGenes.packed.length > 0;
    const geneY = y;
    const geneH = hasFeatures ? Math.max(GENE_MIN_H, packedGenes.totalRows * GENE_LANE_H + 8) : 0;
    if (hasFeatures) {
      y += geneH;
    }

    return {
      rulerY,
      refY,
      strains,
      emptyHintY,
      showCompareTrack: hasCompare,
      compareY,
      compareH: compareTrackHeight,
      hasFeatures,
      packedGenes,
      geneY,
      geneH,
      totalHeight: Math.max(size.height, y + 18),
    };
  }, [
    data,
    compactMode,
    showCompareStrip,
    scannedVariants.length,
    size.height,
    compareTrackHeight,
    packedGenes,
  ]);

  // 최신 props/상태 ref (requestAnimationFrame 내 클로저 무효화 방지)
  const depsRef = useRef({
    data,
    layout,
    size,
    chromLength,
    hasStrain2,
    activeSelection,
    highlightedGene,
    scannedVariants,
    showAllBases,
    showCrosshairs,
    showCenterLine,
    onViewportChange,
    onSelectRegion,
  });
  depsRef.current = {
    data,
    layout,
    size,
    chromLength,
    hasStrain2,
    activeSelection,
    highlightedGene,
    scannedVariants,
    showAllBases,
    showCrosshairs,
    showCenterLine,
    onViewportChange,
    onSelectRegion,
  };

  const drawCanvas = useCallback(() => {
    const canvas = canvasRef.current;
    const {
      data: curData,
      layout: curLayout,
      size: curSize,
      hasStrain2: curHasS2,
      activeSelection: curActiveSel,
      highlightedGene: curHiGene,
      scannedVariants: curScanned,
      showAllBases: curShowAll,
      showCrosshairs: curCrosshairs,
      showCenterLine: curCenterLine,
    } = depsRef.current;

    if (!canvas || !curData || !curLayout) return;

    const dpr = window.devicePixelRatio || 1;
    const w = curSize.width;
    const h = curLayout.totalHeight;
    const targetPixW = Math.floor(w * dpr);
    const targetPixH = Math.floor(h * dpr);
    if (canvas.width !== targetPixW || canvas.height !== targetPixH) {
      canvas.width = targetPixW;
      canvas.height = targetPixH;
      canvas.style.width = `${w}px`;
      canvas.style.height = `${h}px`;
    }

    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const view =
      liveViewRef.current.chrom === curData.chrom
        ? liveViewRef.current
        : { chrom: curData.chrom, start: curData.start, end: curData.end };
    const effectiveStart = view.start;
    const effectiveEnd = view.end;

    const trackW = Math.max(100, w - LEFT_LABEL_W);
    const span = Math.max(1, effectiveEnd - effectiveStart);
    const pxPerBp = trackW / span;
    const toX = (pos: number) => LEFT_LABEL_W + (pos - effectiveStart) * pxPerBp;
    const fontCfg = getActiveCanvasFontConfig();
    const cUi = (px: number, weight: 'normal' | 'bold' | number = 'normal') =>
      formatCanvasFont(fontCfg.uiFontFamily, px + fontCfg.canvasDeltaPx, weight);
    const cMono = (px: number, weight: 'normal' | 'bold' | number = 'normal') =>
      formatCanvasFont(fontCfg.monoFontFamily, px + fontCfg.canvasDeltaPx, weight);
    const readBaseFontPx = Math.min(9.5, Math.max(8, 8.5 + fontCfg.canvasDeltaPx * 0.4));

    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, w, h);

    // Left Gutter
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, LEFT_LABEL_W, h);
    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    ctx.moveTo(LEFT_LABEL_W - 0.5, 0);
    ctx.lineTo(LEFT_LABEL_W - 0.5, h);
    ctx.stroke();

    // 1. Ruler Track (IGV style)
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(LEFT_LABEL_W, curLayout.rulerY, trackW, RULER_H);
    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    ctx.moveTo(0, RULER_H - 0.5);
    ctx.lineTo(w, RULER_H - 0.5);
    ctx.stroke();

    ctx.fillStyle = '#475569';
    ctx.font = cUi(10, 600);
    ctx.textBaseline = 'middle';
    ctx.fillText('bp', 8, curLayout.rulerY + 14);

    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_LABEL_W, 0, trackW, h);
    ctx.clip();

    const tickStep = niceTickStep(span, Math.floor(trackW / 110));
    const firstTick = Math.ceil(effectiveStart / tickStep) * tickStep;
    ctx.font = cMono(10, 500);
    ctx.textBaseline = 'top';

    for (let t = firstTick; t <= effectiveEnd; t += tickStep) {
      const x = Math.round(toX(t)) + 0.5;
      ctx.strokeStyle = '#cbd5e1';
      ctx.beginPath();
      ctx.moveTo(x, curLayout.rulerY + 18);
      ctx.lineTo(x, curLayout.rulerY + RULER_H);
      ctx.stroke();

      ctx.fillStyle = '#334155';
      ctx.fillText(`${(t + 1).toLocaleString()} bp`, x + 4, curLayout.rulerY + 4);
    }

    // 2. Reference Track (IGV crisp nucleotide color bar)
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(LEFT_LABEL_W, curLayout.refY, trackW, REF_TRACK_H);
    if (!curData.isOverviewOnly && curData.refSeq) {
      const visStart = Math.max(curData.start, Math.floor(effectiveStart));
      const visEnd = Math.min(curData.end, Math.ceil(effectiveEnd));
      const baseStep = pxPerBp < 0.7 ? Math.max(1, Math.floor(1 / pxPerBp)) : 1;
      for (let pos = visStart; pos < visEnd; pos += baseStep) {
        const base = curData.refSeq[pos - curData.start] ?? 'N';
        const x0 = toX(pos);
        const bw = Math.max(1, pxPerBp * baseStep);
        ctx.fillStyle = IGV_BASE_COLORS[base] ?? '#94a3b8';
        if (pxPerBp >= 8) {
          ctx.fillRect(x0 + 0.3, curLayout.refY + 2, Math.max(1, bw - 0.6), REF_TRACK_H - 4);
          ctx.fillStyle = '#ffffff';
          ctx.font = cMono(10, 'bold');
          ctx.textBaseline = 'middle';
          ctx.textAlign = 'center';
          ctx.fillText(base, x0 + bw / 2, curLayout.refY + REF_TRACK_H / 2 + 0.5);
          ctx.textAlign = 'left';
        } else {
          ctx.fillRect(x0, curLayout.refY + 4, Math.ceil(bw), REF_TRACK_H - 8);
        }
      }
    } else {
      ctx.fillStyle = '#64748b';
      ctx.font = cUi(10.5);
      ctx.textBaseline = 'middle';
      ctx.fillText(
        '조감 모드 (Ctrl+휠로 30 kb 이하 확대 시 염기 및 리드 표시)',
        LEFT_LABEL_W + 10,
        curLayout.refY + REF_TRACK_H / 2,
      );
    }

    // 3. Strain Tracks (IGV Coverage + Allele Ticks + Neutral Gray Reads)
    const dataSpan = Math.max(1, curData.end - curData.start);
    const visibleCovMax = new Map<string, number>();
    for (const s of curLayout.strains) {
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(LEFT_LABEL_W, s.covY, trackW, s.covH);

      const nBins = s.overviewBins.length;
      let covScaleMax = s.maxBinDepth;
      if (nBins > 0) {
        const binBp = dataSpan / nBins;
        const bStart = Math.max(0, Math.floor((effectiveStart - curData.start) / binBp));
        const bEnd = Math.min(nBins - 1, Math.ceil((effectiveEnd - curData.start) / binBp));

        if (bEnd >= bStart) {
          let localMax = 0;
          for (let b = bStart; b <= bEnd; b++) {
            if (s.overviewBins[b] > localMax) localMax = s.overviewBins[b];
          }
          covScaleMax = Math.max(5, Math.ceil(localMax));
          visibleCovMax.set(s.id, covScaleMax);

          ctx.fillStyle = '#9ca3af';
          ctx.beginPath();
          let started = false;
          let lastX = LEFT_LABEL_W;
          for (let b = bStart; b <= bEnd; b++) {
            const bx0 = toX(curData.start + b * binBp);
            const bx1 = toX(curData.start + (b + 1) * binBp);
            const ratio = Math.min(1, s.overviewBins[b] / covScaleMax);
            const topY = s.covY + s.covH - ratio * (s.covH - 6);
            if (!started) {
              ctx.moveTo(bx0, s.covY + s.covH);
              started = true;
            }
            ctx.lineTo(bx0, topY);
            ctx.lineTo(bx1, topY);
            lastX = bx1;
          }
          if (started) {
            ctx.lineTo(lastX, s.covY + s.covH);
            ctx.closePath();
            ctx.fill();
          }
        }

        if (s.alleleTicks) {
          for (const tk of s.alleleTicks) {
            if (tk.pos + 1 < effectiveStart || tk.pos > effectiveEnd) continue;
            const tx = toX(tk.pos);
            const tw = Math.max(1.2, pxPerBp);
            const totalH = Math.min(1, tk.depth / covScaleMax) * (s.covH - 6);
            let yBottom = s.covY + s.covH;
            const entries: Array<[string, number]> = [
              ['A', tk.counts.A],
              ['C', tk.counts.C],
              ['G', tk.counts.G],
              ['T', tk.counts.T],
            ];
            for (const [base, cnt] of entries) {
              if (cnt === 0) continue;
              const segH = (cnt / tk.depth) * totalH;
              ctx.fillStyle = IGV_BASE_COLORS[base];
              ctx.fillRect(tx, yBottom - segH, tw, segH);
              yBottom -= segH;
            }
          }
        }
      }

      // Track Title Pill
      ctx.font = cUi(11, 600);
      const pillText = s.matchedChrom ? s.label : `${s.label} (미매칭)`;
      const pillW = Math.min(220, ctx.measureText(pillText).width + 14);
      ctx.fillStyle = 'rgba(255, 255, 255, 0.92)';
      ctx.fillRect(LEFT_LABEL_W + 6, s.covY + 4, pillW, 18);
      ctx.strokeStyle = s.matchedChrom ? '#64748b' : '#dc2626';
      ctx.strokeRect(LEFT_LABEL_W + 6.5, s.covY + 4.5, pillW - 1, 17);
      ctx.fillStyle = s.matchedChrom ? '#0f172a' : '#b91c1c';
      ctx.textBaseline = 'middle';
      ctx.fillText(pillText, LEFT_LABEL_W + 13, s.covY + 13.5, pillW - 12);

      // Divider below coverage
      ctx.strokeStyle = '#e2e8f0';
      ctx.beginPath();
      ctx.moveTo(LEFT_LABEL_W, s.covY + s.covH + 0.5);
      ctx.lineTo(w, s.covY + s.covH + 0.5);
      ctx.stroke();

      // Read Pileup (사전 전개된 가시 행만 직접 렌더링 — 프레임당 CIGAR 파싱 0회)
      if (s.pileupH > 0) {
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(LEFT_LABEL_W, s.pileupY, trackW, s.pileupH);

        for (const pr of s.preparedReads) {
          const r = pr.read;
          if (r.end <= effectiveStart || r.start >= effectiveEnd) continue;

          const ry = s.pileupY + 3 + r.row * (READ_ROW_H + READ_ROW_GAP);

          for (const seg of pr.segments) {
            if (seg.end <= effectiveStart || seg.start >= effectiveEnd) continue;
            const sx0 = Math.max(LEFT_LABEL_W, toX(seg.start));
            const sx1 = Math.min(w, toX(seg.end));
            const sw = Math.max(1, sx1 - sx0);

            if (seg.kind === 'match') {
              if (r.mapq === 0) {
                ctx.fillStyle = '#e2e8f0';
                ctx.fillRect(sx0, ry, sw, READ_ROW_H);
                ctx.strokeStyle = '#94a3b8';
                ctx.lineWidth = 0.7;
                ctx.strokeRect(sx0 + 0.5, ry + 0.5, Math.max(1, sw - 1), READ_ROW_H - 1);
                ctx.lineWidth = 1;
              } else {
                ctx.fillStyle = '#b4b8bf';
                ctx.fillRect(sx0, ry, sw, READ_ROW_H);
              }
            } else if (seg.kind === 'del') {
              ctx.strokeStyle = '#0f172a';
              ctx.lineWidth = 1.2;
              ctx.beginPath();
              ctx.moveTo(sx0, ry + READ_ROW_H / 2);
              ctx.lineTo(sx0 + sw, ry + READ_ROW_H / 2);
              ctx.stroke();
              ctx.lineWidth = 1;
            } else if (seg.kind === 'skip') {
              ctx.setLineDash([2, 2]);
              ctx.strokeStyle = '#64748b';
              ctx.beginPath();
              ctx.moveTo(sx0, ry + READ_ROW_H / 2);
              ctx.lineTo(sx0 + sw, ry + READ_ROW_H / 2);
              ctx.stroke();
              ctx.setLineDash([]);
            }
          }

          if (curData.refSeq) {
            const drawFullBases = curShowAll || pxPerBp >= 7;
            const basesToDraw = drawFullBases ? pr.allBases : pr.mismatches;
            for (const b of basesToDraw) {
              if (b.refPos + 1 < effectiveStart || b.refPos > effectiveEnd) continue;
              const refIdx = b.refPos - curData.start;
              const refBase =
                refIdx >= 0 && refIdx < curData.refSeq.length ? curData.refSeq[refIdx] : 'N';
              const isMismatch = !drawFullBases || b.base !== refBase;

              const bx = toX(b.refPos);
              const bw = Math.max(1.4, pxPerBp);
              const qualAlpha =
                b.qual !== undefined && b.qual < 255
                  ? b.qual < 5
                    ? 0.15
                    : b.qual >= 20
                      ? 1
                      : 0.15 + ((b.qual - 5) / 15) * 0.85
                  : 1;
              if (isMismatch) {
                ctx.globalAlpha = qualAlpha;
                ctx.fillStyle = IGV_BASE_COLORS[b.base] ?? '#d97706';
                ctx.fillRect(bx, ry, bw, READ_ROW_H);
                ctx.globalAlpha = 1;
              }
              if (pxPerBp >= 7) {
                ctx.globalAlpha = isMismatch ? qualAlpha : 1;
                ctx.fillStyle = isMismatch ? '#ffffff' : '#1e293b';
                ctx.font = formatCanvasFont(
                  fontCfg.monoFontFamily,
                  readBaseFontPx,
                  isMismatch ? 'bold' : 'normal',
                );
                ctx.textAlign = 'center';
                ctx.textBaseline = 'middle';
                ctx.fillText(b.base, bx + bw / 2, ry + READ_ROW_H / 2 + 0.5);
                ctx.textAlign = 'left';
                ctx.globalAlpha = 1;
              }
            }
          }

          for (const ins of pr.insertions) {
            if (ins.anchorPos + 1 < effectiveStart || ins.anchorPos > effectiveEnd) continue;
            const ix = toX(ins.anchorPos + 1);
            ctx.fillStyle = '#9333ea';
            ctx.fillRect(ix - 1, ry, 2, READ_ROW_H);
            ctx.fillRect(ix - 2.5, ry, 5, 1.5);
            ctx.fillRect(ix - 2.5, ry + READ_ROW_H - 1.5, 5, 1.5);
          }
        }

        if (s.totalRows > MAX_VISIBLE_ROWS) {
          ctx.fillStyle = 'rgba(255, 255, 255, 0.92)';
          ctx.fillRect(w - 138, s.pileupY + s.pileupH - 19, 132, 17);
          ctx.strokeStyle = '#cbd5e1';
          ctx.strokeRect(w - 138, s.pileupY + s.pileupH - 19, 132, 17);
          ctx.fillStyle = '#475569';
          ctx.font = cUi(10, 600);
          ctx.fillText(
            `+${s.totalRows - MAX_VISIBLE_ROWS}개 리드 더 있음`,
            w - 132,
            s.pileupY + s.pileupH - 8,
          );
        }

        ctx.strokeStyle = '#cbd5e1';
        ctx.beginPath();
        ctx.moveTo(0, s.pileupY + s.pileupH + 0.5);
        ctx.lineTo(w, s.pileupY + s.pileupH + 0.5);
        ctx.stroke();
      }
    }

    // Empty Track Guidance when no BAM loaded yet
    if (curData.strains.length === 0) {
      ctx.fillStyle = '#f8fafc';
      ctx.fillRect(LEFT_LABEL_W + 16, curLayout.emptyHintY + 12, trackW - 32, 72);
      ctx.strokeStyle = '#cbd5e1';
      ctx.setLineDash([5, 4]);
      ctx.strokeRect(LEFT_LABEL_W + 16, curLayout.emptyHintY + 12, trackW - 32, 72);
      ctx.setLineDash([]);
      ctx.fillStyle = '#334155';
      ctx.font = cUi(12, 600);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(
        'sacCer3 (S. cerevisiae 17개 전체 염색체)가 준비되었습니다. 상단 [+ BAM/BAI 추가] 또는 화면에 BAM+BAI 파일을 드롭하세요.',
        LEFT_LABEL_W + trackW / 2,
        curLayout.emptyHintY + 48,
      );
      ctx.textAlign = 'left';
    }

    // 4. Optional Compact Comparison Strip (2-Way / 3-Way)
    if (curLayout.showCompareTrack) {
      const cy = curLayout.compareY;
      ctx.fillStyle = '#f8fafc';
      ctx.fillRect(LEFT_LABEL_W, cy, trackW, curLayout.compareH);

      if (!curData.isOverviewOnly && curData.analysis) {
        const { s1Calls, s2Calls, hasStrain2: analysisHasS2 } = curData.analysis;
        const rowH = 16;
        const refRowY = cy + 3;
        const s1RowY = cy + 20;
        const s2RowY = cy + 37;

        const iStart = Math.max(0, Math.floor(effectiveStart) - curData.start);
        const iEnd = Math.min(s1Calls.length - 1, Math.ceil(effectiveEnd) - curData.start);

        for (let i = iStart; i <= iEnd; i++) {
          const pos = curData.start + i;
          const x0 = toX(pos);
          const bw = Math.max(1, pxPerBp);
          const c1 = s1Calls[i];
          const c2 = analysisHasS2 ? s2Calls[i] : undefined;
          const refBase = c1?.refBase ?? 'N';

          let cat: DiffCategory = 'match';
          const s1Low = !c1 || c1.isLowDepth;
          if (!analysisHasS2) {
            if (s1Low) cat = 'gap_s1';
            else cat = c1.call === refBase && c1.insertion === '' ? 'match' : 'strain1_only';
          } else {
            const s2Low = !c2 || c2.isLowDepth;
            if (s1Low && s2Low) cat = 'gap_both';
            else if (s1Low) cat = 'gap_s1';
            else if (s2Low) cat = 'gap_s2';
            else {
              const s1Same = c1.call === refBase && c1.insertion === '';
              const s2Same = c2.call === refBase && c2.insertion === '';
              if (!s1Same && s2Same) cat = 'strain1_only';
              else if (s1Same && !s2Same) cat = 'strain2_only';
              else if (!s1Same && !s2Same) {
                cat =
                  c1.call === c2.call && c1.insertion === c2.insertion ? 'shared' : 'divergent';
              }
            }
          }

          if (cat !== 'match') {
            ctx.fillStyle = DIFF_BG_COLORS[cat];
            ctx.fillRect(x0, cy, bw, curLayout.compareH);
          }

          if (pxPerBp >= 7) {
            ctx.font = cMono(9.5, 'bold');
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';

            ctx.fillStyle = '#475569';
            ctx.fillText(refBase, x0 + bw / 2, refRowY + rowH / 2);

            if (c1) {
              ctx.fillStyle = c1.isLowDepth
                ? '#94a3b8'
                : c1.isMismatch
                  ? '#0369a1'
                  : '#1e293b';
              ctx.fillText(c1.call, x0 + bw / 2, s1RowY + rowH / 2);
              if (c1.insertion) {
                ctx.fillStyle = '#9333ea';
                ctx.fillRect(x0 + bw - 2, s1RowY + 2, 2, rowH - 4);
              }
            }
            if (analysisHasS2 && c2) {
              ctx.fillStyle = c2.isLowDepth
                ? '#94a3b8'
                : c2.isMismatch
                  ? '#047857'
                  : '#1e293b';
              ctx.fillText(c2.call, x0 + bw / 2, s2RowY + rowH / 2);
              if (c2.insertion) {
                ctx.fillStyle = '#9333ea';
                ctx.fillRect(x0 + bw - 2, s2RowY + 2, 2, rowH - 4);
              }
            }
            ctx.textAlign = 'left';
          } else if (cat !== 'match') {
            ctx.fillStyle = DIFF_MARKER_COLORS[cat];
            ctx.fillRect(x0, s1RowY + 2, Math.max(2, bw), analysisHasS2 ? rowH * 2 - 4 : rowH - 4);
          }
        }
      } else {
        for (const v of curScanned) {
          if (v.end < effectiveStart || v.start > effectiveEnd) continue;
          const vx = toX(v.start);
          const vw = Math.max(2, (v.end - v.start) * pxPerBp);
          ctx.fillStyle = DIFF_MARKER_COLORS[v.category] ?? '#d97706';
          ctx.fillRect(vx, cy + 6, vw, curLayout.compareH - 12);
        }
      }
    }

    // 5. IGV-Style RefSeq / SGD Genes Track (다중 레인 배치로 중첩/인접 유전자 충돌 방지)
    if (curLayout.hasFeatures) {
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(LEFT_LABEL_W, curLayout.geneY, trackW, curLayout.geneH);

      for (const pg of curLayout.packedGenes.packed) {
        const feat = pg.feature;
        if (feat.end <= effectiveStart || feat.start >= effectiveEnd) continue;

        const rawX0 = toX(feat.start);
        const rawX1 = toX(feat.end);
        const fx0 = Math.max(LEFT_LABEL_W, rawX0);
        const fx1 = Math.min(w, rawX1);
        const fw = Math.max(2, fx1 - fx0);

        const laneTop = curLayout.geneY + 4 + pg.row * GENE_LANE_H;
        const cy = laneTop + 10;

        const isHighlighted =
          curHiGene !== null &&
          curHiGene.name === feat.name &&
          curHiGene.start === feat.start;

        if (isHighlighted) {
          ctx.fillStyle = 'rgba(250, 204, 21, 0.25)';
          ctx.fillRect(fx0 - 2, laneTop, fw + 4, GENE_LANE_H - 2);
        }

        const geneBlue = isHighlighted ? '#d97706' : '#00008b';
        ctx.strokeStyle = geneBlue;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(fx0, cy);
        ctx.lineTo(fx0 + fw, cy);
        ctx.stroke();
        ctx.lineWidth = 1;

        const blocks =
          feat.subFeatures.length > 0
            ? feat.subFeatures
            : [{ type: feat.type, start: feat.start, end: feat.end }];

        for (const blk of blocks) {
          if (blk.end <= effectiveStart || blk.start >= effectiveEnd) continue;
          const bx0 = Math.max(LEFT_LABEL_W, toX(blk.start));
          const bx1 = Math.min(w, toX(blk.end));
          const bw = Math.max(2, bx1 - bx0);

          // 인접 블록 간 1px 흰색 경계선으로 시각적 분리
          ctx.fillStyle = geneBlue;
          ctx.fillRect(bx0, cy - 5.5, bw, 11);
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1;
          ctx.strokeRect(bx0 - 0.5, cy - 6, bw + 1, 12);

          if (bw >= 14 && (feat.strand === '+' || feat.strand === '-')) {
            ctx.strokeStyle = '#ffffff';
            ctx.lineWidth = 1.2;
            const step = 18;
            for (let cx = bx0 + 8; cx < bx0 + bw - 6; cx += step) {
              ctx.beginPath();
              if (feat.strand === '+') {
                ctx.moveTo(cx - 2.5, cy - 3);
                ctx.lineTo(cx + 2, cy);
                ctx.lineTo(cx - 2.5, cy + 3);
              } else {
                ctx.moveTo(cx + 2.5, cy - 3);
                ctx.lineTo(cx - 2, cy);
                ctx.lineTo(cx + 2.5, cy + 3);
              }
              ctx.stroke();
            }
            ctx.lineWidth = 1;
          }
        }

        if (pg.drawLabel) {
          const trueMidX = (rawX0 + rawX1) / 2;
          const visMidX =
            rawX0 < LEFT_LABEL_W || rawX1 > w
              ? Math.max(fx0 + 28, Math.min(fx1 - 28, (fx0 + fx1) / 2))
              : trueMidX;
          ctx.fillStyle = isHighlighted ? '#92400e' : '#00008b';
          ctx.font = cUi(10.5, 600);
          ctx.textAlign = 'center';
          ctx.textBaseline = 'top';
          ctx.fillText(feat.name, visMidX, laneTop + 18, Math.max(60, fw + 30));
          ctx.textAlign = 'left';
        }
      }
    }

    // Center Line Guide
    if (curCenterLine) {
      const centerPos = Math.floor((effectiveStart + effectiveEnd) / 2);
      const cx0 = toX(centerPos);
      const cx1 = toX(centerPos + 1);
      ctx.strokeStyle = 'rgba(15, 23, 42, 0.35)';
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(cx0, RULER_H);
      ctx.lineTo(cx0, h);
      if (cx1 - cx0 >= 3) {
        ctx.moveTo(cx1, RULER_H);
        ctx.lineTo(cx1, h);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Crosshairs
    const curHover = hoverRef.current;
    const curDrag = dragRef.current;
    if (curCrosshairs && curHover && !curDrag) {
      ctx.strokeStyle = 'rgba(2, 132, 199, 0.45)';
      ctx.setLineDash([2, 2]);
      ctx.beginPath();
      ctx.moveTo(curHover.x + 0.5, 0);
      ctx.lineTo(curHover.x + 0.5, h);
      ctx.moveTo(LEFT_LABEL_W, curHover.y + 0.5);
      ctx.lineTo(w, curHover.y + 0.5);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Highlighted gene dashed vertical marker
    if (curHiGene && curHiGene.end > effectiveStart && curHiGene.start < effectiveEnd) {
      const gx0 = Math.max(LEFT_LABEL_W, toX(curHiGene.start));
      const gx1 = Math.min(w, toX(curHiGene.end));
      if (gx1 > gx0) {
        ctx.strokeStyle = 'rgba(217, 119, 6, 0.55)';
        ctx.setLineDash([4, 3]);
        ctx.strokeRect(gx0, RULER_H, gx1 - gx0, h - RULER_H);
        ctx.setLineDash([]);
      }
    }

    // Active / Dragging Selection Overlay
    const sel =
      curDrag?.mode === 'select'
        ? {
            start: Math.min(curDrag.anchorPos, curDrag.currentPos),
            end: Math.max(curDrag.anchorPos, curDrag.currentPos),
          }
        : curActiveSel;

    if (sel && sel.end > sel.start) {
      const sx0 = Math.max(LEFT_LABEL_W, toX(sel.start));
      const sx1 = Math.min(w, toX(sel.end));
      if (sx1 > sx0) {
        ctx.fillStyle = 'rgba(14, 165, 233, 0.14)';
        ctx.fillRect(sx0, 0, sx1 - sx0, h);
        ctx.strokeStyle = '#0284c7';
        ctx.lineWidth = 1.5;
        ctx.strokeRect(sx0, 0.5, sx1 - sx0, h - 1);
        ctx.lineWidth = 1;

        const badgeText =
          curDrag?.mode === 'select'
            ? `${(sel.start + 1).toLocaleString()}–${sel.end.toLocaleString()} (${(
                sel.end - sel.start
              ).toLocaleString()} bp · Esc/우클릭 취소)`
            : `${(sel.start + 1).toLocaleString()}–${sel.end.toLocaleString()} (${(
                sel.end - sel.start
              ).toLocaleString()} bp)`;
        ctx.font = cMono(10, 600);
        const bw = ctx.measureText(badgeText).width + 16;
        const bx = Math.min(w - bw - 6, Math.max(LEFT_LABEL_W + 4, sx0));
        ctx.fillStyle = '#0284c7';
        ctx.fillRect(bx, 2, bw, 18);
        ctx.fillStyle = '#ffffff';
        ctx.textBaseline = 'middle';
        ctx.fillText(badgeText, bx + 8, 11);
      }
    }

    ctx.restore();

    // Left Gutter Labels & IGV Coverage Y-Axis Scale
    ctx.fillStyle = '#475569';
    ctx.font = cUi(10.5, 600);
    ctx.textBaseline = 'middle';
    ctx.fillText('Reference', 8, curLayout.refY + REF_TRACK_H / 2);

    for (const s of curLayout.strains) {
      const axisX = LEFT_LABEL_W - 14;
      ctx.strokeStyle = '#64748b';
      ctx.beginPath();
      ctx.moveTo(axisX - 4, s.covY + 5);
      ctx.lineTo(axisX, s.covY + 5);
      ctx.lineTo(axisX, s.covY + s.covH - 3);
      ctx.lineTo(axisX - 4, s.covY + s.covH - 3);
      ctx.stroke();

      ctx.fillStyle = '#475569';
      ctx.font = cMono(10);
      ctx.textAlign = 'right';
      ctx.textBaseline = 'top';
      ctx.fillText(String(visibleCovMax.get(s.id) ?? s.maxBinDepth), axisX - 7, s.covY + 3);
      ctx.textBaseline = 'bottom';
      ctx.fillText('0', axisX - 7, s.covY + s.covH - 1);
      ctx.textAlign = 'left';

      if (s.pileupH > 0) {
        ctx.fillStyle = '#64748b';
        ctx.font = cMono(10);
        ctx.textBaseline = 'middle';
        ctx.fillText(`${s.meanDepth.toFixed(1)}x`, 10, s.pileupY + 14);
      }
    }

    if (curLayout.showCompareTrack) {
      ctx.fillStyle = '#0f172a';
      ctx.font = cUi(10.5, 600);
      ctx.textBaseline = 'middle';
      ctx.fillText(
        curHasS2 ? '3-Way 비교' : '2-Way 비교',
        8,
        curLayout.compareY + curLayout.compareH / 2,
      );
    }

    if (curLayout.hasFeatures) {
      ctx.fillStyle = '#334155';
      ctx.font = cUi(10.5, 600);
      ctx.textBaseline = 'middle';
      ctx.fillText('RefSeq / SGD', 8, curLayout.geneY + 18);
    }
  }, []);

  const scheduleDraw = useCallback(() => {
    if (drawRafRef.current !== null) return;
    drawRafRef.current = window.requestAnimationFrame(() => {
      drawRafRef.current = null;
      drawCanvas();
    });
  }, [drawCanvas]);

  const stopGlideLoop = useCallback(() => {
    if (glideRafRef.current !== null) {
      window.cancelAnimationFrame(glideRafRef.current);
      glideRafRef.current = null;
    }
  }, []);

  const startGlideLoop = useCallback(() => {
    if (glideRafRef.current !== null) return;
    const tick = () => {
      const live = liveViewRef.current;
      const target = targetViewRef.current;
      if (live.chrom !== target.chrom) {
        liveViewRef.current = { ...target };
        glideRafRef.current = null;
        drawCanvas();
        return;
      }

      const dStart = target.start - live.start;
      const dEnd = target.end - live.end;
      const span = Math.max(1, target.end - target.start);
      const eps = Math.max(0.005, span * 0.0004);
      if (Math.abs(dStart) <= eps && Math.abs(dEnd) <= eps) {
        liveViewRef.current = { ...target };
        glideRafRef.current = null;
        drawCanvas();
        return;
      }

      const alpha = 0.34;
      liveViewRef.current = {
        chrom: live.chrom,
        start: live.start + dStart * alpha,
        end: live.end + dEnd * alpha,
      };
      drawCanvas();
      glideRafRef.current = window.requestAnimationFrame(tick);
    };
    glideRafRef.current = window.requestAnimationFrame(tick);
  }, [drawCanvas]);

  const emitViewport = useCallback(
    (startVal: number, endVal: number, force = false) => {
      const cLen = depsRef.current.chromLength;
      const s = Math.max(0, Math.round(startVal));
      const e = Math.min(cLen, Math.max(s + 1, Math.round(endVal)));
      const now = performance.now();
      if (force || now - lastSyncAtRef.current >= 55) {
        lastSyncAtRef.current = now;
        if (syncTimerRef.current !== null) {
          window.clearTimeout(syncTimerRef.current);
          syncTimerRef.current = null;
        }
        depsRef.current.onViewportChange(s, e, force);
      } else if (syncTimerRef.current === null) {
        syncTimerRef.current = window.setTimeout(() => {
          syncTimerRef.current = null;
          lastSyncAtRef.current = performance.now();
          const cur = targetViewRef.current;
          const cl = depsRef.current.chromLength;
          const cs = Math.max(0, Math.round(cur.start));
          const ce = Math.min(cl, Math.max(cs + 1, Math.round(cur.end)));
          depsRef.current.onViewportChange(cs, ce, true);
        }, 65);
      }
    },
    [],
  );

  // 외부 내비게이션(초기 로드, 염색체 변경, 유전자 검색 이동, 좌표 입력, Ideogram 클릭, 버튼 줌) 시 뷰포트 동기화
  useEffect(() => {
    if (!data) return;
    const isExternalNav =
      lastSyncedNavRef.current.navVersion !== navVersion ||
      lastSyncedNavRef.current.chrom !== data.chrom;
    if (isExternalNav) {
      lastSyncedNavRef.current = { navVersion, chrom: data.chrom };
      stopGlideLoop();
      liveViewRef.current = { chrom: data.chrom, start: viewStart, end: viewEnd };
      targetViewRef.current = { chrom: data.chrom, start: viewStart, end: viewEnd };
    }
    scheduleDraw();
  }, [data, viewStart, viewEnd, navVersion, stopGlideLoop, scheduleDraw]);

  // 기타 렌더 의존성 변경 시 즉시 캔버스 갱신
  useEffect(() => {
    scheduleDraw();
  }, [
    layout,
    size,
    activeSelection,
    highlightedGene,
    scannedVariants,
    showAllBases,
    showCrosshairs,
    showCenterLine,
    primaryStrainIndex,
    secondaryStrainIndex,
    scheduleDraw,
  ]);

  // 전역 Esc 키로 진행 중인 드래그/선택 즉시 취소, 창 바깥에서 마우스를 떼어도 드래그 정상 종료 및 타이머 정리
  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape' && dragRef.current) {
        dragRef.current = null;
        scheduleDraw();
      }
    };
    const onWindowMouseUp = () => {
      const drag = dragRef.current;
      if (!drag) return;
      dragRef.current = null;
      if (drag.mode === 'pan') {
        const cur = liveViewRef.current;
        emitViewport(cur.start, cur.end, true);
      } else {
        const s = Math.min(drag.anchorPos, drag.currentPos);
        const en = Math.max(drag.anchorPos, drag.currentPos);
        if (en - s >= 2) {
          depsRef.current.onSelectRegion?.(s, en);
        }
        scheduleDraw();
      }
    };
    window.addEventListener('keydown', onKeyDown);
    window.addEventListener('mouseup', onWindowMouseUp);
    window.addEventListener(FONT_CHANGE_EVENT, scheduleDraw);
    document.fonts?.addEventListener?.('loadingdone', scheduleDraw);
    return () => {
      window.removeEventListener('keydown', onKeyDown);
      window.removeEventListener('mouseup', onWindowMouseUp);
      window.removeEventListener(FONT_CHANGE_EVENT, scheduleDraw);
      document.fonts?.removeEventListener?.('loadingdone', scheduleDraw);
      stopGlideLoop();
      if (drawRafRef.current !== null) window.cancelAnimationFrame(drawRafRef.current);
      if (syncTimerRef.current !== null) window.clearTimeout(syncTimerRef.current);
    };
  }, [scheduleDraw, stopGlideLoop, emitViewport]);

  // Native non-passive Wheel 리스너: Ctrl+휠 브라우저 페이지 확대 차단 및 부드러운 휠 좌우 스크롤
  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;

    const onNativeWheel = (e: WheelEvent) => {
      const currentData = depsRef.current.data;
      if (!currentData) return;
      if (e.altKey) return;

      e.preventDefault();
      e.stopPropagation();

      const curTarget =
        targetViewRef.current.chrom === currentData.chrom
          ? targetViewRef.current
          : {
              chrom: currentData.chrom,
              start: liveViewRef.current.start,
              end: liveViewRef.current.end,
            };
      const cLen = depsRef.current.chromLength;
      const tSpan = Math.max(1, curTarget.end - curTarget.start);
      const trackW = Math.max(100, depsRef.current.size.width - LEFT_LABEL_W);

      const normY =
        e.deltaMode === 1 ? e.deltaY * 32 : e.deltaMode === 2 ? e.deltaY * 400 : e.deltaY;
      const normX =
        e.deltaMode === 1 ? e.deltaX * 32 : e.deltaMode === 2 ? e.deltaX * 400 : e.deltaX;

      if (e.ctrlKey || e.metaKey) {
        if (normY === 0) return;
        const rect = el.getBoundingClientRect();
        const relX = e.clientX - rect.left;
        const ratio = Math.max(0, Math.min(1, (relX - LEFT_LABEL_W) / trackW));

        const stepMag = Math.min(280, Math.abs(normY));
        const rawFactor = 1 + stepMag * 0.0022;
        const factor = normY > 0 ? rawFactor : 1 / rawFactor;
        const newSpan = Math.max(30, Math.min(cLen, tSpan * factor));

        const anchor = curTarget.start + ratio * tSpan;
        let newStart = anchor - ratio * newSpan;
        newStart = Math.max(0, Math.min(cLen - newSpan, newStart));
        const newEnd = newStart + newSpan;

        targetViewRef.current = { chrom: currentData.chrom, start: newStart, end: newEnd };
        hoverRef.current = null;
        setHoverInfo(null);
        startGlideLoop();
        emitViewport(newStart, newEnd, true);
        return;
      }

      const rawDelta = Math.abs(normX) > Math.abs(normY) ? normX : normY;
      if (rawDelta === 0) return;

      const shiftBp = (rawDelta / Math.max(200, trackW)) * tSpan * 0.9;
      const newStart = Math.max(0, Math.min(cLen - tSpan, curTarget.start + shiftBp));
      const newEnd = newStart + tSpan;
      if (Math.abs(newStart - curTarget.start) < 1e-6) return;

      targetViewRef.current = { chrom: currentData.chrom, start: newStart, end: newEnd };
      hoverRef.current = null;
      setHoverInfo(null);
      startGlideLoop();

      const nearEdge = checkNearBufferEdge(newStart, newEnd, currentData, cLen);
      emitViewport(newStart, newEnd, nearEdge);
    };

    el.addEventListener('wheel', onNativeWheel, { passive: false });
    return () => el.removeEventListener('wheel', onNativeWheel);
  }, [data !== null, startGlideLoop, emitViewport]);

  const xToPos = (clientX: number): number => {
    const canvas = canvasRef.current;
    if (!canvas || !data) return 0;
    const rect = canvas.getBoundingClientRect();
    const x = clientX - rect.left;
    const trackW = Math.max(100, size.width - LEFT_LABEL_W);
    const ratio = (x - LEFT_LABEL_W) / trackW;
    const v =
      liveViewRef.current.chrom === data.chrom
        ? liveViewRef.current
        : { chrom: data.chrom, start: data.start, end: data.end };
    const span = Math.max(1, v.end - v.start);
    return Math.max(0, Math.min(chromLength, Math.floor(v.start + ratio * span)));
  };

  const posToCanvasX = (pos: number): number => {
    const trackW = Math.max(100, size.width - LEFT_LABEL_W);
    const v =
      data && liveViewRef.current.chrom === data.chrom
        ? liveViewRef.current
        : { start: viewStart, end: viewEnd };
    const span = Math.max(1, v.end - v.start);
    return LEFT_LABEL_W + ((pos - v.start) / span) * trackW;
  };

  const handleMouseDown = (e: React.MouseEvent) => {
    if (!data || e.button !== 0) return;
    const rect = e.currentTarget.getBoundingClientRect();
    const relX = e.clientX - rect.left;
    if (relX < LEFT_LABEL_W) return;

    stopGlideLoop();
    const useSelect = interactionMode === 'select' || e.shiftKey;
    hoverRef.current = null;
    setHoverInfo(null);

    if (useSelect) {
      const pos = xToPos(e.clientX);
      dragRef.current = {
        mode: 'select',
        anchorX: e.clientX,
        anchorPos: pos,
        currentPos: pos,
      };
      scheduleDraw();
    } else {
      const v =
        liveViewRef.current.chrom === data.chrom
          ? liveViewRef.current
          : { chrom: data.chrom, start: data.start, end: data.end };
      dragRef.current = {
        mode: 'pan',
        startX: e.clientX,
        origStart: v.start,
        origEnd: v.end,
      };
    }
  };

  const hitTestGene = (clientX: number, clientY: number, containerEl: EventTarget & Element): GffFeature | undefined => {
    if (!layout || !layout.hasFeatures) return undefined;
    const rect = containerEl.getBoundingClientRect();
    const relX = clientX - rect.left;
    const relY = clientY - rect.top;
    if (relY < layout.geneY || relY >= layout.geneY + layout.geneH) return undefined;
    const row = Math.floor((relY - layout.geneY - 4) / GENE_LANE_H);
    return layout.packedGenes.packed.find((pg) => {
      if (pg.row !== row) return false;
      const gx0 = posToCanvasX(pg.feature.start);
      const gx1 = posToCanvasX(pg.feature.end);
      const mid = (gx0 + gx1) / 2;
      const hitLeft = Math.min(gx0, mid - pg.labelWidthPx / 2) - 4;
      const hitRight = Math.max(gx1, mid + pg.labelWidthPx / 2) + 4;
      return relX >= hitLeft && relX <= hitRight;
    })?.feature;
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (!data || !layout) return;
    const rect = e.currentTarget.getBoundingClientRect();
    const relX = e.clientX - rect.left;
    const relY = e.clientY - rect.top;

    const drag = dragRef.current;
    if (drag) {
      if (drag.mode === 'select') {
        drag.currentPos = xToPos(e.clientX);
        scheduleDraw();
      } else {
        const trackW = Math.max(100, size.width - LEFT_LABEL_W);
        const span = drag.origEnd - drag.origStart;
        const dx = e.clientX - drag.startX;
        const deltaBp = (-dx / trackW) * span;
        const newStart = Math.max(0, Math.min(chromLength - span, drag.origStart + deltaBp));
        const newEnd = newStart + span;
        liveViewRef.current = { chrom: data.chrom, start: newStart, end: newEnd };
        targetViewRef.current = { chrom: data.chrom, start: newStart, end: newEnd };
        scheduleDraw();
        const nearEdge = checkNearBufferEdge(newStart, newEnd, data, chromLength);
        emitViewport(newStart, newEnd, nearEdge);
      }
      return;
    }

    if (relX < LEFT_LABEL_W) {
      hoverRef.current = null;
      setHoverInfo(null);
      if (showCrosshairs) scheduleDraw();
      return;
    }

    const pos = xToPos(e.clientX);
    const feature = hitTestGene(e.clientX, e.clientY, e.currentTarget);

    let read: AlignedRead | undefined;
    for (const s of layout.strains) {
      if (s.pileupH > 0 && relY >= s.pileupY && relY <= s.pileupY + s.pileupH) {
        const row = Math.floor((relY - s.pileupY - 3) / (READ_ROW_H + READ_ROW_GAP));
        const rowReads = s.readsByRow.get(row);
        if (rowReads) {
          read = rowReads.find((r) => pos >= r.start && pos < r.end);
        }
        break;
      }
    }

    const callIdx = pos - data.start;
    const c1 = data.analysis?.s1Calls[callIdx];
    const c2 = data.analysis?.hasStrain2 ? data.analysis.s2Calls[callIdx] : undefined;
    const refBase =
      data.refSeq && callIdx >= 0 && callIdx < data.refSeq.length ? data.refSeq[callIdx] : undefined;

    const nextHover: HoverInfo = {
      x: relX,
      y: relY,
      pos,
      refBase,
      s1Call: c1?.call,
      s1Depth: c1?.depth,
      s2Call: c2?.call,
      s2Depth: c2?.depth,
      hasStrain2: data.analysis?.hasStrain2 ?? false,
      read,
      feature,
    };
    hoverRef.current = nextHover;
    const now = performance.now();
    if (now - lastHoverCommitRef.current >= 32) {
      lastHoverCommitRef.current = now;
      setHoverInfo(nextHover);
    }
    if (showCrosshairs) scheduleDraw();
  };

  const finishDrag = (e: React.MouseEvent) => {
    const drag = dragRef.current;
    if (!drag || !data) return;
    dragRef.current = null;

    if (drag.mode === 'select') {
      const s = Math.min(drag.anchorPos, drag.currentPos);
      const en = Math.max(drag.anchorPos, drag.currentPos);
      const dragPx = Math.abs(e.clientX - drag.anchorX);
      scheduleDraw();
      if (dragPx >= 5 && en - s >= 1) {
        onSelectRegion(s, en);
      } else if (dragPx < 5) {
        if (activeSelection) {
          onCancelSelection();
        } else {
          const hit = hitTestGene(e.clientX, e.clientY, e.currentTarget);
          if (hit) onSelectFeature(hit);
        }
      }
    } else {
      const moved = Math.abs(e.clientX - drag.startX) > 4;
      if (moved) {
        const cur = liveViewRef.current;
        emitViewport(cur.start, cur.end, true);
      } else {
        const hit = hitTestGene(e.clientX, e.clientY, e.currentTarget);
        if (hit) onSelectFeature(hit);
      }
    }
  };

  const handleContextMenu = (e: React.MouseEvent) => {
    if (dragRef.current) {
      e.preventDefault();
      dragRef.current = null;
      scheduleDraw();
    } else if (activeSelection) {
      e.preventDefault();
      onCancelSelection();
    }
  };

  if (!data) {
    return (
      <div
        ref={containerRef}
        className="flex flex-1 items-center justify-center bg-white text-sm text-slate-500"
      >
        sacCer3 지놈을 불러오는 중입니다...
      </div>
    );
  }

  return (
    <div
      ref={containerRef}
      className="relative flex-1 overflow-y-auto overflow-x-hidden bg-white select-none"
    >
      <canvas
        ref={canvasRef}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={finishDrag}
        onContextMenu={handleContextMenu}
        onMouseLeave={() => {
          if (!dragRef.current) {
            hoverRef.current = null;
            setHoverInfo(null);
            if (showCrosshairs) scheduleDraw();
          }
        }}
        className={`block ${
          interactionMode === 'select' ? 'cursor-crosshair' : 'cursor-grab active:cursor-grabbing'
        }`}
      />

      {hoverInfo && (
        <div
          className="pointer-events-none absolute z-20 rounded-lg border border-slate-300 bg-white/95 px-3 py-2 text-xs text-slate-800 shadow-lg backdrop-blur"
          style={{
            left: Math.min(size.width - 260, hoverInfo.x + 14),
            top: Math.max(8, hoverInfo.y + 14),
          }}
        >
          <div className="font-mono font-bold text-sky-700">
            {data.chrom}:{(hoverInfo.pos + 1).toLocaleString()}
          </div>
          {hoverInfo.refBase && (
            <div className="mt-0.5 flex gap-2 font-mono text-[11px]">
              <span className="text-slate-700">Ref:{hoverInfo.refBase}</span>
              {hoverInfo.s1Call !== undefined && (
                <span className="font-semibold text-sky-700">
                  균주1:{hoverInfo.s1Call}({hoverInfo.s1Depth ?? 0}x)
                </span>
              )}
              {hoverInfo.hasStrain2 && (
                <span className="font-semibold text-emerald-700">
                  균주2:{hoverInfo.s2Call ?? '-'}({hoverInfo.s2Depth ?? 0}x)
                </span>
              )}
            </div>
          )}
          {hoverInfo.feature && (
            <div className="mt-1 border-t border-slate-200 pt-1 text-[11px] text-blue-800">
              유전자: <b>{hoverInfo.feature.name}</b> ({hoverInfo.feature.strand},{' '}
              {(hoverInfo.feature.start + 1).toLocaleString()}–
              {hoverInfo.feature.end.toLocaleString()}) — 클릭해 서열 추출
            </div>
          )}
          {hoverInfo.read && (
            <div className="mt-1 border-t border-slate-200 pt-1 font-mono text-[11px] text-slate-600">
              <div>{hoverInfo.read.name}</div>
              <div>
                CIGAR: {hoverInfo.read.cigar} · MAPQ: {hoverInfo.read.mapq} · {hoverInfo.read.strand}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
};
