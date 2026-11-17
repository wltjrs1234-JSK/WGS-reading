import React, { useEffect, useMemo, useRef, useState } from 'react';
import type { CnvSvAnalysisResult } from '../bio/cnvSvAnalyzer';
import { isMitochondrialChrom, normalizeYeastChromName } from '../bio/yeastAliases';
import type { AllChromsOverviewResult } from '../workers/wgsProtocol';
import {
  FONT_CHANGE_EVENT,
  formatCanvasFont,
  getActiveCanvasFontConfig,
} from '../utils/fontSettings';

interface AllChromsOverviewCanvasProps {
  data: AllChromsOverviewResult | null;
  cnvData?: CnvSvAnalysisResult | null;
  onSelectChrom: (chrom: string) => void;
  onOpenCnvModal?: () => void;
}

const LEFT_LABEL_W = 144;
const HEADER_H = 34;
const STRAIN_ROW_H = 68;

export const AllChromsOverviewCanvas: React.FC<AllChromsOverviewCanvasProps> = ({
  data,
  cnvData,
  onSelectChrom,
  onOpenCnvModal,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [size, setSize] = useState({ width: 960, height: 500 });
  const [hoveredChrom, setHoveredChrom] = useState<string | null>(null);
  const [viewMode, setViewMode] = useState<'raw_overview' | 'copy_number'>('copy_number');
  const [fontTick, setFontTick] = useState(0);

  useEffect(() => {
    const bump = () => setFontTick((v) => v + 1);
    window.addEventListener(FONT_CHANGE_EVENT, bump);
    document.fonts?.addEventListener?.('loadingdone', bump);
    return () => {
      window.removeEventListener(FONT_CHANGE_EVENT, bump);
      document.fonts?.removeEventListener?.('loadingdone', bump);
    };
  }, []);

  const showCnvOverlay = Boolean(cnvData && cnvData.strains.length > 0 && viewMode === 'copy_number');

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setSize({
          width: Math.max(400, Math.floor(entry.contentRect.width)),
          height: Math.max(320, Math.floor(entry.contentRect.height)),
        });
      }
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const onWheel = (e: WheelEvent) => {
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
      }
    };
    el.addEventListener('wheel', onWheel, { passive: false });
    return () => el.removeEventListener('wheel', onWheel);
  }, []);

  const layout = useMemo(() => {
    if (!data) return null;
    const trackW = Math.max(200, size.width - LEFT_LABEL_W - 12);
    const totalLen = Math.max(
      1,
      data.contigs.reduce((acc, c) => acc + c.length, 0),
    );

    // 각 염색체에 최소 폭(26px)을 보장하면서 비례 배분
    const minW = 26;
    const rawWidths = data.contigs.map((c) => Math.max(minW, (c.length / totalLen) * trackW));
    const sumRaw = rawWidths.reduce((a, b) => a + b, 0);
    const scale = trackW / sumRaw;

    let xCursor = LEFT_LABEL_W;
    const cols = data.contigs.map((c, idx) => {
      const w = rawWidths[idx] * scale;
      const col = { chrom: c.name, length: c.length, x: xCursor, w };
      xCursor += w;
      return col;
    });

    const totalH = Math.max(
      size.height,
      HEADER_H + Math.max(1, data.strains.length) * STRAIN_ROW_H + 30,
    );

    return { cols, totalH };
  }, [data, size]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !data || !layout) return;

    const dpr = window.devicePixelRatio || 1;
    const w = size.width;
    const h = layout.totalH;
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
    canvas.style.width = `${w}px`;
    canvas.style.height = `${h}px`;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, w, h);

    // Header bar
    ctx.fillStyle = '#f1f5f9';
    ctx.fillRect(0, 0, w, HEADER_H);
    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    ctx.moveTo(0, HEADER_H + 0.5);
    ctx.lineTo(w, HEADER_H + 0.5);
    ctx.stroke();

    const fontCfg = getActiveCanvasFontConfig();
    const cUi = (px: number, weight: 'normal' | 'bold' | number = 'normal') =>
      formatCanvasFont(fontCfg.uiFontFamily, px + fontCfg.canvasDeltaPx, weight);
    const cMono = (px: number, weight: 'normal' | 'bold' | number = 'normal') =>
      formatCanvasFont(fontCfg.monoFontFamily, px + fontCfg.canvasDeltaPx, weight);

    ctx.fillStyle = '#0f172a';
    ctx.font = cUi(11, 700);
    ctx.textBaseline = 'middle';
    ctx.fillText('전체 지놈 (all)', 10, HEADER_H / 2);

    for (let i = 0; i < layout.cols.length; i++) {
      const col = layout.cols[i];
      const isHovered = col.chrom === hoveredChrom;
      if (isHovered) {
        ctx.fillStyle = 'rgba(14, 165, 233, 0.1)';
        ctx.fillRect(col.x, 0, col.w, h);
      }

      if (i % 2 === 1 && !isHovered) {
        ctx.fillStyle = '#f8fafc';
        ctx.fillRect(col.x, HEADER_H, col.w, h - HEADER_H);
      }

      ctx.strokeStyle = '#cbd5e1';
      ctx.beginPath();
      ctx.moveTo(col.x + 0.5, 0);
      ctx.lineTo(col.x + 0.5, h);
      ctx.stroke();

      const shortLabel = normalizeYeastChromName(col.chrom).replace(/^chr/, '');
      ctx.fillStyle = isHovered ? '#0369a1' : '#334155';
      ctx.font = cMono(11, 700);
      ctx.textAlign = 'center';
      ctx.fillText(shortLabel, col.x + col.w / 2, HEADER_H / 2, col.w - 4);
      ctx.textAlign = 'left';
    }

    if (data.strains.length === 0) {
      ctx.fillStyle = '#64748b';
      ctx.font = cUi(13);
      ctx.textAlign = 'center';
      ctx.fillText(
        '상단 [+ BAM/BAI 추가] 또는 화면에 BAM+BAI 파일을 드롭하면 17개 전체 염색체 커버리지가 표시됩니다.',
        LEFT_LABEL_W + (w - LEFT_LABEL_W) / 2,
        HEADER_H + 80,
      );
      ctx.textAlign = 'left';
      return;
    }

    for (let sIdx = 0; sIdx < data.strains.length; sIdx++) {
      const strain = data.strains[sIdx];
      const cnvProfile = showCnvOverlay
        ? cnvData?.strains.find((p) => p.strainId === strain.id || p.strainLabel === strain.label)
        : undefined;
      const rowY = HEADER_H + sIdx * STRAIN_ROW_H;

      // Left label
      ctx.fillStyle = '#f8fafc';
      ctx.fillRect(0, rowY, LEFT_LABEL_W, STRAIN_ROW_H);
      ctx.strokeStyle = '#cbd5e1';
      ctx.beginPath();
      ctx.moveTo(LEFT_LABEL_W - 0.5, rowY);
      ctx.lineTo(LEFT_LABEL_W - 0.5, rowY + STRAIN_ROW_H);
      ctx.moveTo(0, rowY + STRAIN_ROW_H + 0.5);
      ctx.lineTo(w, rowY + STRAIN_ROW_H + 0.5);
      ctx.stroke();

      const matchedCount = strain.byChrom.filter((b) => b.matchedChrom !== null).length;
      ctx.fillStyle = '#0f172a';
      ctx.font = cUi(11, 600);
      ctx.fillText(strain.label, 10, rowY + 22, LEFT_LABEL_W - 16);
      ctx.fillStyle = cnvProfile
        ? '#0f766e'
        : matchedCount === layout.cols.length
          ? '#059669'
          : '#d97706';
      ctx.font = cMono(10);
      ctx.fillText(
        cnvProfile
          ? `${cnvProfile.effectivePloidy}n · 중앙 ${cnvProfile.genomeMedianDepth}x`
          : `매핑: ${matchedCount}/${layout.cols.length} chr`,
        10,
        rowY + 42,
        LEFT_LABEL_W - 16,
      );

      // CNV baseline dashed horizontal line at basePloidy (mid-row)
      if (cnvProfile) {
        const baseY = rowY + STRAIN_ROW_H - 6 - 0.5 * (STRAIN_ROW_H - 18);
        ctx.save();
        ctx.setLineDash([3, 3]);
        ctx.strokeStyle = 'rgba(15, 118, 110, 0.35)';
        ctx.beginPath();
        ctx.moveTo(LEFT_LABEL_W, baseY);
        ctx.lineTo(w, baseY);
        ctx.stroke();
        ctx.restore();
      }

      let maxVal = 1;
      if (!cnvProfile) {
        for (const bc of strain.byChrom) {
          for (const v of bc.bins) {
            if (Number.isFinite(v) && v > maxVal) maxVal = v;
          }
        }
      }

      for (let cIdx = 0; cIdx < layout.cols.length; cIdx++) {
        const col = layout.cols[cIdx];
        const info = strain.byChrom[cIdx];
        if (!info) continue;

        if (!info.matchedChrom) {
          ctx.fillStyle = 'rgba(244, 63, 94, 0.12)';
          ctx.fillRect(col.x + 1, rowY + 6, Math.max(1, col.w - 2), STRAIN_ROW_H - 12);
          continue;
        }

        if (cnvProfile) {
          const chrBins = cnvProfile.bins.filter((b) => b.chrom === col.chrom);
          const chrSummary = cnvProfile.chromSummaries.find((cs) => cs.chrom === col.chrom);
          const isMito = isMitochondrialChrom(col.chrom);
          const maxCopies = cnvProfile.effectivePloidy * 2;
          const binW = col.w / Math.max(1, chrBins.length);

          for (let b = 0; b < chrBins.length; b++) {
            const val = chrBins[b];
            const ratioH = isMito
              ? Math.min(1, Math.log10(Math.max(1, val.copyRatio)) / 2.2)
              : Math.min(1, Math.max(0.05, val.estCopies / maxCopies));
            const barH = ratioH * (STRAIN_ROW_H - 18);
            const delta = val.estCopies - cnvProfile.effectivePloidy;
            ctx.fillStyle = isMito
              ? '#0d9488'
              : delta >= 0.48
                ? '#e11d48'
                : delta <= -0.48
                  ? '#4f46e5'
                  : '#38bdf8';
            ctx.fillRect(
              col.x + b * binW,
              rowY + STRAIN_ROW_H - 6 - barH,
              Math.max(1.2, binW + 0.3),
              barH,
            );
          }

          if (chrSummary && (chrSummary.status === 'gain' || chrSummary.status === 'loss')) {
            ctx.fillStyle = chrSummary.status === 'gain' ? '#be123c' : '#3730a3';
            ctx.font = cMono(9, 700);
            ctx.textAlign = 'center';
            ctx.fillText(
              `${chrSummary.estCopies}n`,
              col.x + col.w / 2,
              rowY + 12,
              col.w - 2,
            );
            ctx.textAlign = 'left';
          }
          continue;
        }

        const bins = info.bins;
        const binW = col.w / Math.max(1, bins.length);
        ctx.fillStyle = '#9ca3af';
        for (let b = 0; b < bins.length; b++) {
          const val = Number.isFinite(bins[b]) ? bins[b] : 0;
          const barH = (val / maxVal) * (STRAIN_ROW_H - 16);
          ctx.fillRect(
            col.x + b * binW,
            rowY + STRAIN_ROW_H - 6 - barH,
            Math.max(1, binW + 0.3),
            barH,
          );
        }
      }
    }
  }, [data, cnvData, showCnvOverlay, layout, size, hoveredChrom, fontTick]);

  const colAtX = (clientX: number): string | null => {
    const canvas = canvasRef.current;
    if (!canvas || !layout) return null;
    const rect = canvas.getBoundingClientRect();
    const x = clientX - rect.left;
    const found = layout.cols.find((c) => x >= c.x && x < c.x + c.w);
    return found?.chrom ?? null;
  };

  if (!data) {
    return (
      <div className="flex flex-1 items-center justify-center bg-white text-xs text-slate-500">
        전체 지놈 커버리지를 불러오는 중입니다...
      </div>
    );
  }

  return (
    <div className="flex flex-1 flex-col overflow-hidden bg-white">
      {data.strains.length > 0 && (
        <div className="flex flex-wrap items-center justify-between gap-2 border-b border-slate-200 bg-slate-50 px-3 py-1.5 text-xs">
          <div className="flex items-center gap-2">
            <span className="font-semibold text-slate-700">전체 염색체 표시 방식:</span>
            <div className="inline-flex rounded border border-slate-300 bg-white p-0.5 text-[11px]">
              <button
                type="button"
                onClick={() => setViewMode('raw_overview')}
                className={`rounded px-2 py-0.5 font-medium ${
                  !showCnvOverlay ? 'bg-slate-800 text-white' : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                상대 밀도 (BAI 인덱스)
              </button>
              <button
                type="button"
                onClick={() => {
                  if (!cnvData && onOpenCnvModal) onOpenCnvModal();
                  else setViewMode('copy_number');
                }}
                className={`rounded px-2 py-0.5 font-medium ${
                  showCnvOverlay ? 'bg-teal-600 text-white' : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                정규화 복제수 (CNV · 이수성)
              </button>
            </div>
            {showCnvOverlay && (
              <span className="text-[11px] text-slate-600">
                (<span className="font-semibold text-rose-600">빨강</span>=증가 ·{' '}
                <span className="font-semibold text-indigo-600">파랑</span>=결실 · 점선=기준 배수)
              </span>
            )}
          </div>
          {onOpenCnvModal && (
            <button
              type="button"
              onClick={onOpenCnvModal}
              className="rounded border border-teal-300 bg-teal-50 px-2.5 py-0.5 text-[11px] font-semibold text-teal-800 hover:bg-teal-100"
            >
              {cnvData ? 'CNV · LOH · 전좌 상세 보기' : '정밀 CNV · LOH · 전좌 스캔 열기'}
            </button>
          )}
        </div>
      )}
      <div ref={containerRef} className="relative flex-1 overflow-y-auto overflow-x-hidden bg-white">
        <canvas
          ref={canvasRef}
          onMouseMove={(e) => setHoveredChrom(colAtX(e.clientX))}
          onMouseLeave={() => setHoveredChrom(null)}
          onClick={(e) => {
            const c = colAtX(e.clientX);
            if (c) onSelectChrom(c);
          }}
          className="block cursor-pointer"
        />
      </div>
    </div>
  );
};

