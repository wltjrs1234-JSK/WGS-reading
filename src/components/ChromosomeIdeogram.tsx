import React, { useRef, useState } from 'react';
import type { ContigInfo } from '../bio/types';

interface ChromosomeIdeogramProps {
  chrom: string;
  chromLength: number;
  start: number;
  end: number;
  contigs: ContigInfo[];
  onJumpToPosition: (newStart: number, newEnd: number) => void;
}

export const ChromosomeIdeogram: React.FC<ChromosomeIdeogramProps> = ({
  chrom,
  chromLength,
  start,
  end,
  contigs,
  onJumpToPosition,
}) => {
  const barRef = useRef<HTMLDivElement>(null);
  const [dragging, setDragging] = useState(false);

  if (chrom === 'all') {
    return (
      <div className="flex items-center justify-between border-b border-slate-300 bg-slate-50 px-3 py-1.5 text-[11px] text-slate-700">
        <span>
          전체 지놈 보기 (<b>{contigs.length}</b>개 염색체) — 아래 각 염색체 영역을 클릭하면 해당 염색체로 확대 이동합니다.
        </span>
      </div>
    );
  }

  const span = Math.max(1, end - start);
  const leftPct = Math.max(0, Math.min(100, (start / Math.max(1, chromLength)) * 100));
  const widthPct = Math.max(0.6, Math.min(100 - leftPct, (span / Math.max(1, chromLength)) * 100));

  const moveFromClientX = (clientX: number) => {
    const el = barRef.current;
    if (!el || chromLength <= 0) return;
    const rect = el.getBoundingClientRect();
    const ratio = Math.max(0, Math.min(1, (clientX - rect.left) / Math.max(1, rect.width)));
    const center = Math.round(ratio * chromLength);
    const ns = Math.max(0, Math.min(chromLength - span, Math.round(center - span / 2)));
    onJumpToPosition(ns, ns + span);
  };

  return (
    <div className="flex items-center gap-3 border-b border-slate-300 bg-white px-3 py-1.5">
      <div className="w-24 shrink-0 truncate font-mono text-[10px] font-semibold text-slate-600">
        {chrom} ({(chromLength / 1000).toFixed(1)} kb)
      </div>

      <div
        ref={barRef}
        onMouseDown={(e) => {
          setDragging(true);
          moveFromClientX(e.clientX);
        }}
        onMouseMove={(e) => {
          if (dragging) moveFromClientX(e.clientX);
        }}
        onMouseUp={() => setDragging(false)}
        onMouseLeave={() => setDragging(false)}
        className="relative h-3.5 flex-1 cursor-pointer overflow-hidden rounded-full border border-slate-300 bg-gradient-to-b from-slate-100 to-slate-200"
        title="클릭하거나 드래그해 염색체 내 원하는 위치로 즉시 이동"
      >
        {/* Subtle chromosome band ticks */}
        {[20, 40, 60, 80].map((pct) => (
          <div
            key={pct}
            className="pointer-events-none absolute top-0 bottom-0 w-px bg-slate-300/70"
            style={{ left: `${pct}%` }}
          />
        ))}

        {/* IGV Red Viewport Box */}
        <div
          className="pointer-events-none absolute top-0 bottom-0 border-2 border-red-600 bg-red-500/20"
          style={{
            left: `${leftPct}%`,
            width: `${widthPct}%`,
          }}
        />
      </div>
    </div>
  );
};

