import React, { useEffect, useMemo, useRef, useState } from 'react';
import { normalizeYeastChromName } from '../bio/yeastAliases';
import type { AllChromsOverviewResult } from '../workers/wgsProtocol';

interface AllChromsOverviewCanvasProps {
  data: AllChromsOverviewResult | null;
  onSelectChrom: (chrom: string) => void;
}

const LEFT_LABEL_W = 144;
const HEADER_H = 34;
const STRAIN_ROW_H = 68;

export const AllChromsOverviewCanvas: React.FC<AllChromsOverviewCanvasProps> = ({
  data,
  onSelectChrom,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [size, setSize] = useState({ width: 960, height: 500 });
  const [hoveredChrom, setHoveredChrom] = useState<string | null>(null);

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

    ctx.fillStyle = '#0f172a';
    ctx.font = '700 11px ui-sans-serif, system-ui, sans-serif';
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
      ctx.font = '700 11px ui-monospace, monospace';
      ctx.textAlign = 'center';
      ctx.fillText(shortLabel, col.x + col.w / 2, HEADER_H / 2, col.w - 4);
      ctx.textAlign = 'left';
    }

    if (data.strains.length === 0) {
      ctx.fillStyle = '#64748b';
      ctx.font = '13px ui-sans-serif, system-ui, sans-serif';
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
      ctx.font = '600 11px ui-sans-serif, system-ui, sans-serif';
      ctx.fillText(strain.label, 10, rowY + 22, LEFT_LABEL_W - 16);
      ctx.fillStyle = matchedCount === layout.cols.length ? '#059669' : '#d97706';
      ctx.font = '10px ui-monospace, monospace';
      ctx.fillText(
        `매핑: ${matchedCount}/${layout.cols.length} chr`,
        10,
        rowY + 42,
        LEFT_LABEL_W - 16,
      );

      let maxVal = 1;
      for (const bc of strain.byChrom) {
        for (const v of bc.bins) {
          if (Number.isFinite(v) && v > maxVal) maxVal = v;
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
  }, [data, layout, size, hoveredChrom]);

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
  );
};

