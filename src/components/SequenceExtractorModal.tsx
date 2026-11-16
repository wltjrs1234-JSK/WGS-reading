import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  AlertTriangle,
  ArrowLeftRight,
  Check,
  CheckCircle2,
  ChevronLeft,
  ChevronRight,
  Copy,
  Download,
  ExternalLink,
  FlaskConical,
  Scissors,
  X,
} from 'lucide-react';
import {
  computeRegionMappingStats,
  sliceSplicedAnalysis,
  type RegionMappingStats,
} from '../bio/consensus';
import {
  evaluateStrainProtein,
  reverseComplement,
  translateDnaWithCodons,
  type AaChangeDetail,
  type StrainProteinEvaluation,
  type TranslationTableId,
} from '../bio/geneticCode';
import { isMitochondrialChrom } from '../bio/yeastAliases';
import type { ConsensusOptions, ConsensusPositionCall, GffFeature } from '../bio/types';
import type { ExtractionResult } from '../workers/wgsProtocol';

interface SequenceExtractorModalProps {
  isOpen: boolean;
  isLoading: boolean;
  errorMsg: string | null;
  result: ExtractionResult | null;
  selectedFeature?: GffFeature | null;
  chromLength: number;
  options: ConsensusOptions;
  onClose: () => void;
  onClearAndClose: () => void;
  onSendToLocalBlast?: (queryFasta: string) => void;
  onRecompute: (
    chrom: string,
    start: number,
    end: number,
    nextOptions: ConsensusOptions,
  ) => void;
}

interface AlignmentRowSpec {
  key: string;
  label: string;
  seq: string;
  labelClass: string;
  classifyChar: (ch: string, refCh: string, colIdx0: number) => {
    className: string;
    title?: string;
  };
}

interface CoalescedRun {
  text: string;
  className: string;
  title?: string;
}

const STRAIN_CARD_STYLES = [
  {
    accent: 'border-sky-200 bg-sky-50/50',
    badge: 'text-sky-800',
    labelClass: 'text-sky-700',
    diffClass: 'bg-sky-200 font-bold text-sky-950',
  },
  {
    accent: 'border-emerald-200 bg-emerald-50/50',
    badge: 'text-emerald-800',
    labelClass: 'text-emerald-700',
    diffClass: 'bg-emerald-200 font-bold text-emerald-950',
  },
  {
    accent: 'border-violet-200 bg-violet-50/50',
    badge: 'text-violet-800',
    labelClass: 'text-violet-700',
    diffClass: 'bg-violet-200 font-bold text-violet-950',
  },
  {
    accent: 'border-indigo-200 bg-indigo-50/50',
    badge: 'text-indigo-800',
    labelClass: 'text-indigo-700',
    diffClass: 'bg-indigo-200 font-bold text-indigo-950',
  },
];

function formatFastaRecord(header: string, seq: string, lineLen = 80): string {
  const lines: string[] = [`>${header}`];
  for (let i = 0; i < seq.length; i += lineLen) {
    lines.push(seq.slice(i, i + lineLen));
  }
  return lines.join('\n');
}

function downloadTextFile(filename: string, content: string, mime = 'text/plain;charset=utf-8') {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function coalesceRuns(
  seq: string,
  refSlice: string,
  startCol0: number,
  classify: AlignmentRowSpec['classifyChar'],
): CoalescedRun[] {
  if (seq.length === 0) return [];
  const runs: CoalescedRun[] = [];
  let curStart = 0;
  let curMeta = classify(seq[0], refSlice[0] ?? '', startCol0);

  for (let i = 1; i < seq.length; i++) {
    const meta = classify(seq[i], refSlice[i] ?? '', startCol0 + i);
    if (meta.className === curMeta.className && meta.title === curMeta.title) {
      continue;
    }
    runs.push({
      text: seq.slice(curStart, i),
      className: curMeta.className,
      title: curMeta.title,
    });
    curStart = i;
    curMeta = meta;
  }
  runs.push({
    text: seq.slice(curStart),
    className: curMeta.className,
    title: curMeta.title,
  });
  return runs;
}

const SequenceBlockAlignment: React.FC<{
  title: string;
  unitLabel: string;
  refLabel: string;
  refSeq: string;
  strainRows: AlignmentRowSpec[];
  blockSize?: number;
  pageSizeUnits: number;
  diffOnly: boolean;
  onToggleDiffOnly: () => void;
  page: number;
  onPageChange: (p: number) => void;
  focusedCol0: number | null;
  formatMatchSymbol?: (colIdx0: number, refCh: string, strainChars: string[]) => string;
  formatBlockCoord?: (startCol0: number, endCol0: number) => string;
  emptyMessage: string;
}> = ({
  refLabel,
  refSeq,
  strainRows,
  blockSize = 60,
  pageSizeUnits,
  diffOnly,
  onToggleDiffOnly,
  page,
  onPageChange,
  focusedCol0,
  formatMatchSymbol,
  formatBlockCoord,
  emptyMessage,
}) => {
  const scrollContainerRef = useRef<HTMLDivElement | null>(null);
  const blockRefs = useRef<Map<number, HTMLDivElement>>(new Map());

  const totalLen = useMemo(() => {
    let maxL = refSeq.length;
    for (const r of strainRows) {
      if (r.seq.length > maxL) maxL = r.seq.length;
    }
    return maxL;
  }, [refSeq, strainRows]);

  const allBlocks = useMemo(() => {
    const blocks: Array<{ start: number; end: number; hasDiff: boolean }> = [];
    for (let s = 0; s < totalLen; s += blockSize) {
      const e = Math.min(totalLen, s + blockSize);
      let hasDiff = false;
      for (let i = s; i < e; i++) {
        const rc = refSeq[i] ?? '';
        for (const row of strainRows) {
          const sc = row.seq[i];
          if (sc !== undefined && sc !== rc) {
            hasDiff = true;
            break;
          }
        }
        if (hasDiff) break;
      }
      blocks.push({ start: s, end: e, hasDiff });
    }
    return blocks;
  }, [totalLen, blockSize, refSeq, strainRows]);

  const visibleBlocks = useMemo(
    () => (diffOnly ? allBlocks.filter((b) => b.hasDiff) : allBlocks),
    [allBlocks, diffOnly],
  );

  const blocksPerPage = Math.max(1, Math.floor(pageSizeUnits / blockSize));
  const totalPages = Math.max(1, Math.ceil(visibleBlocks.length / blocksPerPage));
  const safePage = Math.min(Math.max(0, page), totalPages - 1);
  const pageBlocks = visibleBlocks.slice(safePage * blocksPerPage, (safePage + 1) * blocksPerPage);

  useEffect(() => {
    if (focusedCol0 === null) return;
    const blockStart = Math.floor(focusedCol0 / blockSize) * blockSize;
    const el = blockRefs.current.get(blockStart);
    if (el) {
      el.scrollIntoView({ behavior: 'smooth', block: 'center' });
    }
  }, [focusedCol0, safePage, diffOnly, blockSize]);

  if (totalLen === 0) {
    return (
      <div className="rounded-lg border border-slate-200 bg-white p-3 text-center text-xs text-slate-500">
        {emptyMessage}
      </div>
    );
  }

  return (
    <div className="space-y-2">
      {/* Toolbar & Pagination */}
      <div className="flex flex-wrap items-center justify-between gap-2 text-[11px]">
        <div className="flex items-center gap-2">
          <button
            type="button"
            onClick={onToggleDiffOnly}
            className={`rounded-lg border px-2.5 py-1 font-semibold transition ${
              diffOnly
                ? 'border-sky-500 bg-sky-600 text-white'
                : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-100'
            }`}
          >
            {diffOnly
              ? `차이 포함 블록만 보는 중 (${visibleBlocks.length}/${allBlocks.length} 블록)`
              : `전체 서열 (${totalLen.toLocaleString()}개) · 차이 블록 필터`}
          </button>
        </div>

        {totalPages > 1 && (
          <div className="flex items-center gap-1.5">
            <button
              type="button"
              disabled={safePage === 0}
              onClick={() => onPageChange(safePage - 1)}
              className="inline-flex items-center rounded border border-slate-300 bg-white px-1.5 py-0.5 text-slate-700 hover:bg-slate-100 disabled:opacity-40"
            >
              <ChevronLeft className="h-3.5 w-3.5" />
            </button>
            <select
              value={safePage}
              onChange={(e) => onPageChange(Number(e.target.value))}
              className="rounded border border-slate-300 bg-white px-2 py-0.5 font-mono text-[11px] text-slate-800"
            >
              {Array.from({ length: totalPages }, (_, pIdx) => {
                const first = visibleBlocks[pIdx * blocksPerPage];
                const last =
                  visibleBlocks[Math.min(visibleBlocks.length - 1, (pIdx + 1) * blocksPerPage - 1)];
                return (
                  <option key={pIdx} value={pIdx}>
                    페이지 {pIdx + 1}/{totalPages} ({(first.start + 1).toLocaleString()}–
                    {last.end.toLocaleString()})
                  </option>
                );
              })}
            </select>
            <button
              type="button"
              disabled={safePage >= totalPages - 1}
              onClick={() => onPageChange(safePage + 1)}
              className="inline-flex items-center rounded border border-slate-300 bg-white px-1.5 py-0.5 text-slate-700 hover:bg-slate-100 disabled:opacity-40"
            >
              <ChevronRight className="h-3.5 w-3.5" />
            </button>
          </div>
        )}
      </div>

      {/* Block List */}
      <div
        ref={scrollContainerRef}
        className="max-h-72 space-y-2.5 overflow-x-auto overflow-y-auto rounded-lg border border-slate-200 bg-white p-3 font-mono text-[11px] leading-[18px]"
      >
        {pageBlocks.length === 0 ? (
          <div className="py-4 text-center text-xs text-slate-500">
            일치하지 않는(차이/변이) 블록이 없습니다. 전체 서열 보기 모드로 전환해 전체를 확인할 수 있습니다.
          </div>
        ) : (
          pageBlocks.map((blk) => {
            const refSlice = refSeq.slice(blk.start, blk.end);
            const isFocusedBlock =
              focusedCol0 !== null && focusedCol0 >= blk.start && focusedCol0 < blk.end;
            const matchChars =
              formatMatchSymbol && strainRows.length > 0
                ? refSlice
                    .split('')
                    .map((rc, offset) =>
                      formatMatchSymbol(
                        blk.start + offset,
                        rc,
                        strainRows.map((r) => r.seq[blk.start + offset] ?? ''),
                      ),
                    )
                    .join('')
                : null;

            return (
              <div
                key={blk.start}
                ref={(el) => {
                  if (el) blockRefs.current.set(blk.start, el);
                  else blockRefs.current.delete(blk.start);
                }}
                className={`rounded px-1.5 py-1 ${
                  isFocusedBlock
                    ? 'bg-amber-50/70 ring-1 ring-amber-400'
                    : blk.hasDiff
                      ? 'bg-slate-50/70'
                      : ''
                }`}
              >
                <div className="mb-0.5 flex items-center justify-between text-[10px] text-slate-400">
                  <span>
                    Pos {(blk.start + 1).toLocaleString()}–{blk.end.toLocaleString()}
                    {formatBlockCoord ? ` · ${formatBlockCoord(blk.start, blk.end)}` : ''}
                  </span>
                  {blk.hasDiff && (
                    <span className="rounded bg-amber-100 px-1.5 font-sans font-semibold text-amber-800">
                      변이/차이 구간
                    </span>
                  )}
                </div>

                {/* Reference Row */}
                <div className="whitespace-nowrap">
                  <span className="mr-2 inline-block w-28 truncate align-bottom font-semibold text-slate-600">
                    {refLabel}
                  </span>
                  {focusedCol0 !== null && isFocusedBlock ? (
                    <>
                      <span className="text-slate-800">
                        {refSlice.slice(0, focusedCol0 - blk.start)}
                      </span>
                      <span className="bg-amber-300 font-bold text-slate-950">
                        {refSlice[focusedCol0 - blk.start]}
                      </span>
                      <span className="text-slate-800">
                        {refSlice.slice(focusedCol0 - blk.start + 1)}
                      </span>
                    </>
                  ) : (
                    <span className="text-slate-800">{refSlice}</span>
                  )}
                </div>

                {/* Match Indicator Row */}
                {matchChars !== null && (
                  <div className="whitespace-nowrap text-slate-400 select-none">
                    <span className="mr-2 inline-block w-28 truncate align-bottom text-[10px]">
                      일치
                    </span>
                    <span>{matchChars}</span>
                  </div>
                )}

                {/* Strain Rows */}
                {strainRows.map((row) => {
                  const sSlice = row.seq.slice(blk.start, blk.end);
                  const runs = coalesceRuns(sSlice, refSlice, blk.start, row.classifyChar);
                  return (
                    <div key={row.key} className="whitespace-nowrap">
                      <span
                        className={`mr-2 inline-block w-28 truncate align-bottom font-semibold ${row.labelClass}`}
                        title={row.label}
                      >
                        {row.label}
                      </span>
                      {runs.map((run, rIdx) => (
                        <span key={rIdx} className={run.className} title={run.title}>
                          {run.text}
                        </span>
                      ))}
                    </div>
                  );
                })}
              </div>
            );
          })
        )}
      </div>
    </div>
  );
};

export const SequenceExtractorModal: React.FC<SequenceExtractorModalProps> = ({
  isOpen,
  isLoading,
  errorMsg,
  result,
  selectedFeature,
  chromLength,
  options,
  onClose,
  onClearAndClose,
  onSendToLocalBlast,
  onRecompute,
}) => {
  const [seqFormatMode, setSeqFormatMode] = useState<'ungapped' | 'aligned'>('ungapped');
  const [isRevComp, setIsRevComp] = useState(false);
  const [useSplicedCds, setUseSplicedCds] = useState(false);
  const [frame, setFrame] = useState<0 | 1 | 2>(0);
  const [manualTableId, setManualTableId] = useState<TranslationTableId | null>(null);
  const [copiedKey, setCopiedKey] = useState<string | null>(null);
  const [expandedCards, setExpandedCards] = useState(false);

  const [editStart1, setEditStart1] = useState<number | null>(null);
  const [editEnd1, setEditEnd1] = useState<number | null>(null);

  const [dnaPage, setDnaPage] = useState(0);
  const [dnaDiffOnly, setDnaDiffOnly] = useState(false);
  const [focusedDnaCol, setFocusedDnaCol] = useState<number | null>(null);

  const [aaPage, setAaPage] = useState(0);
  const [aaDiffOnly, setAaDiffOnly] = useState(false);
  const [focusedAaPos, setFocusedAaPos] = useState<number | null>(null);

  const rawAnalysis = result?.analysis ?? null;

  const featureExons = useMemo(() => {
    if (!selectedFeature || !rawAnalysis) return [];
    return selectedFeature.subFeatures
      .filter((sf) => sf.type === 'CDS' || sf.type === 'exon')
      .map((sf) => ({
        start: Math.max(rawAnalysis.start, sf.start),
        end: Math.min(rawAnalysis.end, sf.end),
      }))
      .filter((ex) => ex.end > ex.start)
      .sort((a, b) => a.start - b.start);
  }, [selectedFeature, rawAnalysis]);

  const isExactFeatureBounds = Boolean(
    selectedFeature &&
      rawAnalysis &&
      rawAnalysis.chrom === selectedFeature.chrom &&
      rawAnalysis.start === selectedFeature.start &&
      rawAnalysis.end === selectedFeature.end,
  );
  const hasSplicedExons = featureExons.length > 1;

  const selectionResetKey = rawAnalysis
    ? `${rawAnalysis.chrom}:${rawAnalysis.start}-${rawAnalysis.end}:${selectedFeature?.id ?? ''}`
    : '';

  useEffect(() => {
    if (!rawAnalysis) return;
    setIsRevComp(selectedFeature?.strand === '-');
    setUseSplicedCds(hasSplicedExons && isExactFeatureBounds);
    setFrame(0);
    setManualTableId(null);
    setEditStart1(null);
    setEditEnd1(null);
    setDnaPage(0);
    setAaPage(0);
    setFocusedDnaCol(null);
    setFocusedAaPos(null);
    setDnaDiffOnly(false);
    setAaDiffOnly(false);
  }, [selectionResetKey]); // eslint-disable-line react-hooks/exhaustive-deps

  const activeAnalysis = useMemo(() => {
    if (!rawAnalysis) return null;
    if (useSplicedCds && hasSplicedExons) {
      return sliceSplicedAnalysis(rawAnalysis, featureExons);
    }
    return rawAnalysis;
  }, [rawAnalysis, useSplicedCds, hasSplicedExons, featureExons]);

  const hasStrain1 = Boolean(result?.hasStrain1 ?? (result && result.strain1Label !== 'Reference Only'));
  const hasStrain2 = Boolean(rawAnalysis?.hasStrain2 && result?.strain2Label);

  const autoTableId: TranslationTableId = rawAnalysis
    ? isMitochondrialChrom(rawAnalysis.chrom)
      ? 3
      : 1
    : 1;
  const activeTableId = manualTableId ?? autoTableId;

  const extractedModel = useMemo(() => {
    if (!result || !activeAnalysis || !rawAnalysis) return null;

    const orient = (s: string) => (isRevComp ? reverseComplement(s) : s);

    // 1. 카드 및 FASTA 추출용 서열 (선택된 Real/Aligned & 정방향/역상보 반영)
    const refCardSeq = orient(
      seqFormatMode === 'ungapped' ? activeAnalysis.ungappedRef : activeAnalysis.alignedRef,
    );
    const alignedRefOriented = orient(activeAnalysis.alignedRef);
    const ungappedRefOriented = orient(activeAnalysis.ungappedRef);

    // 2. Column metadata (strand 방향 반영)
    const orderedCols = isRevComp
      ? [...activeAnalysis.columns].reverse()
      : activeAnalysis.columns;

    // Reference 번역 사전 계산
    const refTrans = translateDnaWithCodons(ungappedRefOriented, activeTableId, frame);
    const refClean = ungappedRefOriented.replace(/-/g, '').toUpperCase().replace(/U/g, 'T');
    const preRef = { ...refTrans, clean: refClean };

    interface StrainEntry {
      key: string;
      title: string;
      cardSeq: string;
      ungappedOriented: string;
      alignedOriented: string;
      mappingStats: RegionMappingStats | null;
      lowDepthFlagsByAlignedCol?: boolean[];
      proteinEval: StrainProteinEvaluation;
      style: (typeof STRAIN_CARD_STYLES)[number];
    }

    const strains: StrainEntry[] = [];

    const buildLowDepthFlags = (calls: ConsensusPositionCall[]): boolean[] => {
      // aligned 컬럼 기준으로 해당 컬럼의 refPos call이 lowDepth인지 매핑
      const callByPos = new Map<number, ConsensusPositionCall>();
      for (const c of calls) callByPos.set(c.pos, c);
      return orderedCols.map((col) => Boolean(callByPos.get(col.refPos)?.isLowDepth));
    };

    if (hasStrain1) {
      const s1CardSeq = orient(
        seqFormatMode === 'ungapped' ? activeAnalysis.ungappedS1 : activeAnalysis.alignedS1,
      );
      const s1Ungapped = orient(activeAnalysis.ungappedS1);
      const s1Aligned = orient(activeAnalysis.alignedS1);
      strains.push({
        key: 's1',
        title: result.strain1Label,
        cardSeq: s1CardSeq,
        ungappedOriented: s1Ungapped,
        alignedOriented: s1Aligned,
        mappingStats: computeRegionMappingStats(activeAnalysis.s1Calls, activeAnalysis.refSeq),
        lowDepthFlagsByAlignedCol: buildLowDepthFlags(activeAnalysis.s1Calls),
        proteinEval: evaluateStrainProtein(
          ungappedRefOriented,
          s1Ungapped,
          activeTableId,
          frame,
          preRef,
        ),
        style: STRAIN_CARD_STYLES[0],
      });
    }

    if (hasStrain2 && result.strain2Label) {
      const s2CardSeq = orient(
        seqFormatMode === 'ungapped' ? activeAnalysis.ungappedS2 : activeAnalysis.alignedS2,
      );
      const s2Ungapped = orient(activeAnalysis.ungappedS2);
      const s2Aligned = orient(activeAnalysis.alignedS2);
      strains.push({
        key: 's2',
        title: result.strain2Label,
        cardSeq: s2CardSeq,
        ungappedOriented: s2Ungapped,
        alignedOriented: s2Aligned,
        mappingStats: computeRegionMappingStats(activeAnalysis.s2Calls, activeAnalysis.refSeq),
        lowDepthFlagsByAlignedCol: buildLowDepthFlags(activeAnalysis.s2Calls),
        proteinEval: evaluateStrainProtein(
          ungappedRefOriented,
          s2Ungapped,
          activeTableId,
          frame,
          preRef,
        ),
        style: STRAIN_CARD_STYLES[1],
      });
    }

    result.extraStrains.forEach((ex, idx) => {
      const useSp = useSplicedCds && hasSplicedExons && ex.splicedUngappedSeq !== undefined;
      const rawU = useSp ? ex.splicedUngappedSeq! : ex.ungappedSeq;
      const rawA = useSp ? ex.splicedAlignedSeq! : ex.alignedSeq;
      const uSeq = orient(rawU);
      const aSeq = orient(rawA);
      const cardSeq = seqFormatMode === 'ungapped' ? uSeq : aSeq;
      const mStats =
        useSp && ex.splicedMappingStats ? ex.splicedMappingStats : (ex.mappingStats ?? null);

      strains.push({
        key: `ex_${idx}`,
        title: ex.label,
        cardSeq,
        ungappedOriented: uSeq,
        alignedOriented: aSeq,
        mappingStats: mStats,
        proteinEval: evaluateStrainProtein(
          ungappedRefOriented,
          uSeq,
          activeTableId,
          frame,
          preRef,
        ),
        style: STRAIN_CARD_STYLES[(2 + idx) % STRAIN_CARD_STYLES.length],
      });
    });

    // DNA 변이 컬럼 빠른 이동 목록
    const dnaDiffMarkers: Array<{ colIdx0: number; refPos1: number; summary: string; isGap: boolean }> =
      [];
    for (let i = 0; i < alignedRefOriented.length; i++) {
      const rc = alignedRefOriented[i];
      let hasDiff = false;
      let allGap = false;
      const parts: string[] = [];
      for (const s of strains) {
        const sc = s.alignedOriented[i];
        if (sc !== undefined && sc !== rc) {
          hasDiff = true;
          if (sc === 'N' && s.lowDepthFlagsByAlignedCol?.[i]) {
            allGap = true;
          }
          parts.push(`${s.title.slice(0, 10)}:${sc}`);
        }
      }
      if (hasDiff) {
        const colMeta = orderedCols[i];
        dnaDiffMarkers.push({
          colIdx0: i,
          refPos1: (colMeta?.refPos ?? rawAnalysis.start) + 1,
          summary: `${rc}→${parts[0] ?? '?'}`,
          isGap: allGap,
        });
      }
    }

    const strandTag = isRevComp ? 'revcomp' : 'sense';
    const cdsTag = useSplicedCds && hasSplicedExons ? '_CDS' : '';
    const modeTag = seqFormatMode === 'ungapped' ? 'real_dna' : 'aligned';
    const coordTag = `${rawAnalysis.chrom}:${rawAnalysis.start + 1}-${rawAnalysis.end}${cdsTag}_${strandTag}_${modeTag}`;
    const proteinTag = `${rawAnalysis.chrom}:${rawAnalysis.start + 1}-${rawAnalysis.end}${cdsTag}_${strandTag}_frame+${frame + 1}_table${activeTableId}`;

    const fastaEntries = [formatFastaRecord(`Reference|${coordTag}`, refCardSeq)];
    const proteinFastaEntries = [formatFastaRecord(`Reference|${proteinTag}`, refTrans.aa)];
    for (const s of strains) {
      fastaEntries.push(formatFastaRecord(`${s.title.replace(/\s+/g, '_')}|${coordTag}`, s.cardSeq));
      proteinFastaEntries.push(
        formatFastaRecord(`${s.title.replace(/\s+/g, '_')}|${proteinTag}`, s.proteinEval.aa),
      );
    }

    const refNBp = activeAnalysis.refSeq.split('').filter((b) => b.toUpperCase() === 'N').length;

    return {
      refCardSeq,
      alignedRefOriented,
      ungappedRefOriented,
      orderedCols,
      refAA: refTrans.aa,
      refCodons: refTrans.codons,
      refNBp,
      strains,
      dnaDiffMarkers,
      coordTag,
      proteinTag,
      multiFasta: fastaEntries.join('\n'),
      proteinMultiFasta: proteinFastaEntries.join('\n'),
    };
  }, [
    result,
    activeAnalysis,
    rawAnalysis,
    seqFormatMode,
    isRevComp,
    useSplicedCds,
    hasSplicedExons,
    hasStrain1,
    hasStrain2,
    activeTableId,
    frame,
  ]);

  if (!isOpen) return null;

  const copyText = async (key: string, text: string) => {
    try {
      await navigator.clipboard.writeText(text);
      setCopiedKey(key);
      setTimeout(() => setCopiedKey(null), 1800);
    } catch {
      const ta = document.createElement('textarea');
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      document.body.removeChild(ta);
      setCopiedKey(key);
      setTimeout(() => setCopiedKey(null), 1800);
    }
  };

  const currentStart1 = editStart1 ?? (rawAnalysis ? rawAnalysis.start + 1 : 1);
  const currentEnd1 = editEnd1 ?? (rawAnalysis ? rawAnalysis.end : 1);

  const applyRange = (nextStart0: number, nextEnd0: number, nextOpts = options) => {
    if (!rawAnalysis) return;
    const s = Math.max(0, Math.min(chromLength - 1, nextStart0));
    const e = Math.min(chromLength, Math.max(s + 1, nextEnd0));
    setEditStart1(null);
    setEditEnd1(null);
    onRecompute(rawAnalysis.chrom, s, e, nextOpts);
  };

  const jumpToDnaCol = (col0: number) => {
    setDnaDiffOnly(false);
    const pageIdx = Math.floor(col0 / 6000);
    setDnaPage(pageIdx);
    setFocusedDnaCol(col0);
  };

  const jumpToAaPos = (aaPos1: number, alsoJumpDna = false) => {
    const col0 = Math.max(0, aaPos1 - 1);
    setAaDiffOnly(false);
    setAaPage(Math.floor(col0 / 1800));
    setFocusedAaPos(col0);
    if (alsoJumpDna) {
      const dnaCol0 = frame + col0 * 3;
      jumpToDnaCol(dnaCol0);
    }
  };

  const DNA_PAGE_SIZE = 6000;
  const AA_PAGE_SIZE = 1800;

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/45 p-4 backdrop-blur-xs"
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div className="flex max-h-[93vh] w-full max-w-6xl flex-col overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-2xl">
        {/* Header */}
        <div className="flex items-center justify-between gap-3 border-b border-slate-200 bg-slate-50 px-6 py-3.5">
          <div className="flex items-center gap-3">
            <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-xl bg-sky-100 text-sky-700">
              <Scissors className="h-5 w-5" />
            </div>
            <div>
              <div className="flex flex-wrap items-center gap-2">
                <h2 className="text-base font-bold text-slate-900">
                  선택 구간 균주별 합의 서열(Consensus) 및 단백질 전체 비교
                </h2>
                {selectedFeature && (
                  <span className="inline-flex items-center gap-1 rounded-md border border-indigo-200 bg-indigo-50 px-2 py-0.5 text-xs font-bold text-indigo-800">
                    {selectedFeature.name}
                    {selectedFeature.systematicName &&
                    selectedFeature.systematicName !== selectedFeature.name
                      ? ` (${selectedFeature.systematicName})`
                      : ''}{' '}
                    · {selectedFeature.strand === '-' ? '(-) 역방향' : '(+) 정방향'}
                    {hasSplicedExons ? ` · Exon ${featureExons.length}개` : ''}
                  </span>
                )}
              </div>
              {rawAnalysis && (
                <p className="font-mono text-xs text-slate-600">
                  {rawAnalysis.chrom}:{(rawAnalysis.start + 1).toLocaleString()}–
                  {rawAnalysis.end.toLocaleString()} (게놈{' '}
                  {(rawAnalysis.end - rawAnalysis.start).toLocaleString()} bp
                  {useSplicedCds && activeAnalysis
                    ? ` · CDS 스플라이싱 ${activeAnalysis.refSeq.length.toLocaleString()} bp`
                    : ''}
                  {' · '}변이 {rawAnalysis.variants.length}건 ·{' '}
                  {extractedModel?.strains.length
                    ? `Ref + 균주 ${extractedModel.strains.length}개 비교`
                    : 'Reference 단독'})
                </p>
              )}
            </div>
          </div>
          <div className="flex items-center gap-2">
            <button
              type="button"
              onClick={onClearAndClose}
              className="rounded-lg border border-slate-300 bg-white px-3 py-1.5 text-xs font-medium text-slate-700 hover:bg-slate-100"
            >
              선택 취소 후 닫기
            </button>
            <button
              onClick={onClose}
              className="rounded-lg p-1.5 text-slate-500 hover:bg-slate-200 hover:text-slate-900"
              title="창 닫기 (Esc)"
            >
              <X className="h-5 w-5" />
            </button>
          </div>
        </div>

        {/* Body */}
        <div className="flex-1 space-y-4 overflow-y-auto p-5">
          {errorMsg && (
            <div className="rounded-xl border border-rose-300 bg-rose-50 p-3 text-xs text-rose-800">
              {errorMsg}
            </div>
          )}

          {isLoading && (
            <div className="rounded-xl border border-sky-200 bg-sky-50 p-4 text-center text-xs font-medium text-sky-800">
              선택한 구간의 전체 리드 Pileup 및 CIGAR 기반 Consensus 서열을 계산하고 있습니다...
            </div>
          )}

          {rawAnalysis && activeAnalysis && extractedModel && (
            <>
              {/* Control Bar */}
              <div className="flex flex-wrap items-center justify-between gap-2.5 rounded-xl border border-slate-200 bg-slate-50 p-3">
                <div className="flex flex-wrap items-center gap-2 text-xs">
                  <span className="font-medium text-slate-700">좌표(1-based):</span>
                  <input
                    type="number"
                    min={1}
                    max={chromLength}
                    value={currentStart1}
                    onChange={(e) => setEditStart1(Number(e.target.value))}
                    className="w-24 rounded-lg border border-slate-300 bg-white px-2 py-1 font-mono text-slate-900 focus:border-sky-500 focus:outline-none"
                  />
                  <span className="text-slate-400">–</span>
                  <input
                    type="number"
                    min={1}
                    max={chromLength}
                    value={currentEnd1}
                    onChange={(e) => setEditEnd1(Number(e.target.value))}
                    className="w-24 rounded-lg border border-slate-300 bg-white px-2 py-1 font-mono text-slate-900 focus:border-sky-500 focus:outline-none"
                  />
                  <button
                    type="button"
                    onClick={() => applyRange(currentStart1 - 1, currentEnd1)}
                    className="rounded-lg bg-sky-100 px-2.5 py-1 font-semibold text-sky-800 hover:bg-sky-200"
                  >
                    적용
                  </button>
                  <button
                    type="button"
                    onClick={() => applyRange(rawAnalysis.start - 50, rawAnalysis.end + 50)}
                    className="rounded-lg border border-slate-300 bg-white px-2 py-1 font-medium text-slate-700 hover:bg-slate-100"
                    title="프라이머 설계용 앞뒤 50bp 확장"
                  >
                    ±50bp
                  </button>
                  <button
                    type="button"
                    onClick={() => applyRange(rawAnalysis.start - 200, rawAnalysis.end + 200)}
                    className="rounded-lg border border-slate-300 bg-white px-2 py-1 font-medium text-slate-700 hover:bg-slate-100"
                  >
                    ±200bp
                  </button>
                  {selectedFeature && !isExactFeatureBounds && (
                    <button
                      type="button"
                      onClick={() => applyRange(selectedFeature.start, selectedFeature.end)}
                      className="rounded-lg border border-indigo-300 bg-indigo-50 px-2 py-1 font-semibold text-indigo-800 hover:bg-indigo-100"
                    >
                      {selectedFeature.name} CDS 경계 복원
                    </button>
                  )}
                </div>

                <div className="flex flex-wrap items-center gap-2 text-xs">
                  <div className="flex overflow-hidden rounded-lg border border-slate-300 bg-white">
                    <button
                      type="button"
                      onClick={() => setSeqFormatMode('ungapped')}
                      className={`px-2.5 py-1 font-semibold ${
                        seqFormatMode === 'ungapped'
                          ? 'bg-sky-600 text-white'
                          : 'text-slate-700 hover:bg-slate-100'
                      }`}
                    >
                      실제 서열 (Ungapped)
                    </button>
                    <button
                      type="button"
                      onClick={() => setSeqFormatMode('aligned')}
                      className={`border-l border-slate-200 px-2.5 py-1 font-semibold ${
                        seqFormatMode === 'aligned'
                          ? 'bg-sky-600 text-white'
                          : 'text-slate-700 hover:bg-slate-100'
                      }`}
                    >
                      비교 정렬 (Aligned)
                    </button>
                  </div>

                  {hasSplicedExons && (
                    <button
                      type="button"
                      onClick={() => setUseSplicedCds((v) => !v)}
                      className={`rounded-lg border px-2.5 py-1 font-semibold ${
                        useSplicedCds
                          ? 'border-indigo-400 bg-indigo-100 text-indigo-900'
                          : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-100'
                      }`}
                      title="Intron이 포함된 유전자에서 Exon(CDS)만 연결해 DNA/아미노산 번역 정확도를 유지합니다."
                    >
                      {useSplicedCds
                        ? `CDS 스플라이싱 (${featureExons.length} Exons)`
                        : '전체 게놈 (Intron 포함)'}
                    </button>
                  )}

                  <button
                    type="button"
                    onClick={() => setIsRevComp((v) => !v)}
                    className={`inline-flex items-center gap-1 rounded-lg border px-2.5 py-1 font-semibold ${
                      isRevComp
                        ? 'border-violet-400 bg-violet-100 text-violet-900'
                        : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-100'
                    }`}
                  >
                    <ArrowLeftRight className="h-3.5 w-3.5" />
                    {isRevComp ? "역상보 5'→3' (-)" : "정방향 5'→3' (+)"}
                  </button>

                  <button
                    type="button"
                    onClick={() =>
                      applyRange(rawAnalysis.start, rawAnalysis.end, {
                        ...options,
                        ploidyMode: options.ploidyMode === 'haploid' ? 'diploid' : 'haploid',
                      })
                    }
                    className={`rounded-lg border px-2.5 py-1 font-semibold ${
                      options.ploidyMode === 'diploid'
                        ? 'border-amber-400 bg-amber-100 text-amber-900'
                        : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-100'
                    }`}
                  >
                    {options.ploidyMode === 'haploid' ? 'Haploid' : 'Diploid (IUPAC)'}
                  </button>
                </div>
              </div>

              {/* Full-region Mapping Verification Banner */}
              {extractedModel.strains.length > 0 && (() => {
                const gaps = extractedModel.strains.filter(
                  (s) => s.mappingStats && s.mappingStats.lowDepthBp > 0,
                );
                const totalBp = activeAnalysis.refSeq.length;
                if (gaps.length === 0) {
                  return (
                    <div className="flex flex-wrap items-center justify-between gap-2 rounded-xl border border-emerald-200 bg-emerald-50/70 px-3.5 py-2 text-xs text-emerald-900">
                      <span className="inline-flex items-center gap-1.5 font-semibold">
                        <CheckCircle2 className="h-4 w-4 text-emerald-600" />
                        선택 구간 전체 매핑 완료: Reference와 로드된 균주{' '}
                        {extractedModel.strains.length}개 모두 {totalBp.toLocaleString()} bp 전
                        구간이 기준 Depth({options.minDepth}x+) 이상으로 완전히 매핑되었습니다.
                      </span>
                      <span className="font-mono text-[11px] text-emerald-800">
                        {extractedModel.strains
                          .map((s) =>
                            s.mappingStats
                              ? `${s.title}: 평균 ${s.mappingStats.meanDepth}x`
                              : s.title,
                          )
                          .join(' · ')}
                      </span>
                    </div>
                  );
                }
                return (
                  <div className="flex flex-wrap items-center justify-between gap-2 rounded-xl border border-amber-300 bg-amber-50 px-3.5 py-2 text-xs text-amber-900">
                    <span className="inline-flex items-center gap-1.5 font-semibold">
                      <AlertTriangle className="h-4 w-4 text-amber-600" />
                      일부 구간 저심도/미매핑(N) 발견 (Depth &lt; {options.minDepth}x):
                    </span>
                    <span className="font-mono text-[11px]">
                      {gaps
                        .map((g) => {
                          const st = g.mappingStats!;
                          return `${g.title} ${st.coveragePct}% (N ${st.lowDepthBp} bp · 평균 ${st.meanDepth}x)`;
                        })
                        .join(' / ')}
                    </span>
                  </div>
                );
              })()}

              {/* Strain Cards + Full Scrollable Sequences */}
              <div className="space-y-2">
                <div className="flex items-center justify-between text-[11px] text-slate-600">
                  <span className="font-semibold">
                    Reference 및 매핑 균주 서열 ({1 + extractedModel.strains.length}개 트랙)
                  </span>
                  <button
                    type="button"
                    onClick={() => setExpandedCards((v) => !v)}
                    className="font-semibold text-sky-700 hover:underline"
                  >
                    {expandedCards ? '카드 서열 접기' : '카드 내 전체 서열 펼치기'}
                  </button>
                </div>

                <div
                  className={`grid gap-3 ${
                    extractedModel.strains.length === 0
                      ? 'md:grid-cols-1'
                      : extractedModel.strains.length === 1
                        ? 'md:grid-cols-2'
                        : 'md:grid-cols-3'
                  }`}
                >
                  {[
                    {
                      key: 'ref',
                      title: 'Reference',
                      seq: extractedModel.refCardSeq,
                      accent: 'border-emerald-200 bg-emerald-50/50',
                      badge: 'text-emerald-800',
                      statsText:
                        extractedModel.refNBp > 0
                          ? `Reference N ${extractedModel.refNBp} bp`
                          : 'Reference 100%',
                      isComplete: extractedModel.refNBp === 0,
                    },
                    ...extractedModel.strains.map((s) => ({
                      key: s.key,
                      title: s.title,
                      seq: s.cardSeq,
                      accent: s.style.accent,
                      badge: s.style.badge,
                      statsText: s.mappingStats
                        ? `매핑 ${s.mappingStats.coveragePct}% (${s.mappingStats.coveredBp.toLocaleString()}/${s.mappingStats.totalRefBp.toLocaleString()} bp) · ${s.mappingStats.meanDepth}x`
                        : '매핑 통계 없음',
                      isComplete: s.mappingStats ? s.mappingStats.lowDepthBp === 0 : true,
                    })),
                  ].map((item) => {
                    const previewSeq =
                      expandedCards && item.seq.length <= 30_000
                        ? item.seq
                        : item.seq.slice(0, 120);
                    return (
                      <div
                        key={item.key}
                        className={`flex flex-col justify-between rounded-xl border p-3 ${item.accent}`}
                      >
                        <div className="flex items-center justify-between gap-2">
                          <span className={`truncate text-xs font-bold ${item.badge}`} title={item.title}>
                            {item.title}
                          </span>
                          <span className="shrink-0 font-mono text-xs font-semibold text-slate-700">
                            {item.seq.replace(/-/g, '').length.toLocaleString()} bp
                          </span>
                        </div>
                        <div
                          className={`mt-1 text-[10px] font-medium ${
                            item.isComplete ? 'text-emerald-700' : 'text-amber-700'
                          }`}
                        >
                          {item.statsText}
                        </div>
                        <div
                          className={`mt-2 rounded border border-slate-200 bg-white px-2 py-1.5 font-mono text-[11px] text-slate-800 ${
                            expandedCards
                              ? 'max-h-32 overflow-y-auto break-all leading-4'
                              : 'truncate'
                          }`}
                        >
                          {previewSeq}
                          {(!expandedCards && item.seq.length > 120) ||
                          (expandedCards && item.seq.length > 30_000)
                            ? '...'
                            : ''}
                        </div>
                        <div className="mt-2.5 flex gap-2">
                          <button
                            type="button"
                            onClick={() => copyText(`${item.key}_raw`, item.seq)}
                            className="inline-flex flex-1 items-center justify-center gap-1 rounded-lg border border-slate-300 bg-white px-2 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                          >
                            {copiedKey === `${item.key}_raw` ? (
                              <Check className="h-3 w-3 text-emerald-600" />
                            ) : (
                              <Copy className="h-3 w-3" />
                            )}
                            순수 서열 복사
                          </button>
                          <button
                            type="button"
                            onClick={() =>
                              copyText(
                                `${item.key}_fa`,
                                formatFastaRecord(
                                  `${item.title.replace(/\s+/g, '_')}|${extractedModel.coordTag}`,
                                  item.seq,
                                ),
                              )
                            }
                            className="inline-flex flex-1 items-center justify-center gap-1 rounded-lg border border-slate-300 bg-white px-2 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                          >
                            {copiedKey === `${item.key}_fa` ? (
                              <Check className="h-3 w-3 text-emerald-600" />
                            ) : (
                              <Copy className="h-3 w-3" />
                            )}
                            FASTA 복사
                          </button>
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>

              {/* Full Nucleotide Alignment Comparison */}
              <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
                <div className="mb-2 flex flex-wrap items-center justify-between gap-2">
                  <div>
                    <span className="text-xs font-bold text-slate-800">
                      염기 단위 전체 정렬 비교 (
                      {extractedModel.alignedRefOriented.length.toLocaleString()} 컬럼 ·{' '}
                      {isRevComp ? "역상보 5'→3' (-)" : "정방향 5'→3' (+)"})
                    </span>
                    <span className="ml-2 text-[11px] text-slate-500">
                      차이/변이 {extractedModel.dnaDiffMarkers.length}컬럼
                    </span>
                  </div>
                  <div className="flex flex-wrap items-center gap-2 text-[11px] text-slate-600">
                    {extractedModel.strains.map((s) => (
                      <span key={s.key} className={`font-semibold ${s.style.labelClass}`}>
                        ■ {s.title.slice(0, 12)} 변이
                      </span>
                    ))}
                    <span className="font-semibold text-rose-700">■ 결실(-)</span>
                    <span className="font-semibold text-amber-700">■ 저심도(N)</span>
                  </div>
                </div>

                {extractedModel.dnaDiffMarkers.length > 0 && (
                  <div className="mb-2 flex flex-wrap items-center gap-1">
                    <span className="mr-1 text-[10px] font-semibold text-slate-500">
                      변이 바로가기:
                    </span>
                    {extractedModel.dnaDiffMarkers.slice(0, 24).map((m) => (
                      <button
                        key={m.colIdx0}
                        type="button"
                        onClick={() => jumpToDnaCol(m.colIdx0)}
                        className={`rounded border px-1.5 py-0.5 font-mono text-[10px] font-semibold transition ${
                          focusedDnaCol === m.colIdx0
                            ? 'border-amber-500 bg-amber-200 text-amber-950'
                            : m.isGap
                              ? 'border-amber-200 bg-amber-50 text-amber-800 hover:bg-amber-100'
                              : 'border-sky-200 bg-white text-sky-800 hover:bg-sky-100'
                        }`}
                        title={`${rawAnalysis.chrom}:${m.refPos1.toLocaleString()} (서열 Pos ${m.colIdx0 + 1})`}
                      >
                        Pos {m.colIdx0 + 1} ({m.summary})
                      </button>
                    ))}
                    {extractedModel.dnaDiffMarkers.length > 24 && (
                      <span className="text-[10px] text-slate-500">
                        외 {extractedModel.dnaDiffMarkers.length - 24}건
                      </span>
                    )}
                  </div>
                )}

                <SequenceBlockAlignment
                  title="염기 정렬"
                  unitLabel="bp"
                  refLabel="Reference"
                  refSeq={extractedModel.alignedRefOriented}
                  pageSizeUnits={DNA_PAGE_SIZE}
                  diffOnly={dnaDiffOnly}
                  onToggleDiffOnly={() => {
                    setDnaDiffOnly((v) => !v);
                    setDnaPage(0);
                  }}
                  page={dnaPage}
                  onPageChange={setDnaPage}
                  focusedCol0={focusedDnaCol}
                  formatMatchSymbol={(_colIdx0, rc, sChars) =>
                    sChars.some((sc) => sc !== rc) ? '*' : '|'
                  }
                  formatBlockCoord={(s0, e0) => {
                    const c0 = extractedModel.orderedCols[s0];
                    const c1 = extractedModel.orderedCols[Math.min(e0 - 1, extractedModel.orderedCols.length - 1)];
                    if (!c0 || !c1) return '';
                    return `${rawAnalysis.chrom}:${(c0.refPos + 1).toLocaleString()}–${(c1.refPos + 1).toLocaleString()}`;
                  }}
                  emptyMessage="선택된 구간의 염기 서열이 없습니다."
                  strainRows={extractedModel.strains.map((s) => ({
                    key: s.key,
                    label: s.title,
                    seq: s.alignedOriented,
                    labelClass: s.style.labelClass,
                    classifyChar: (ch, refCh, colIdx0) => {
                      const isLow = Boolean(s.lowDepthFlagsByAlignedCol?.[colIdx0]);
                      if (ch === 'N' && isLow) {
                        return {
                          className: 'bg-amber-200 font-bold text-amber-950',
                          title: `Pos ${colIdx0 + 1}: 저심도/미매핑 (N)`,
                        };
                      }
                      if (ch === '-' && refCh !== '-') {
                        return {
                          className: 'bg-rose-200 font-bold text-rose-950',
                          title: `Pos ${colIdx0 + 1}: 결실 (${refCh} → -)`,
                        };
                      }
                      if (ch !== refCh) {
                        return {
                          className: s.style.diffClass,
                          title: `Pos ${colIdx0 + 1}: ${refCh} → ${ch}`,
                        };
                      }
                      return { className: 'text-slate-600' };
                    },
                  }))}
                />
              </div>

              {/* Full Amino Acid Translation & Protein Comparison Panel */}
              <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
                <div className="mb-3 flex flex-wrap items-center justify-between gap-2">
                  <div className="flex flex-wrap items-center gap-2">
                    <FlaskConical className="h-4 w-4 text-violet-600" />
                    <span className="text-xs font-bold text-slate-800">
                      아미노산 전체 번역 및 단백질 서열 비교 ({extractedModel.refAA.length.toLocaleString()} aa)
                    </span>
                    {selectedFeature && !isExactFeatureBounds && (
                      <span className="rounded bg-amber-100 px-2 py-0.5 text-[10px] font-semibold text-amber-900">
                        주의: 좌표 확장은 시작 프레임(+1)을 이동시킬 수 있습니다
                      </span>
                    )}
                  </div>
                  <div className="flex flex-wrap items-center gap-2 text-xs">
                    <select
                      value={activeTableId}
                      onChange={(e) =>
                        setManualTableId(Number(e.target.value) as TranslationTableId)
                      }
                      className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                    >
                      <option value={1}>Table 1: 표준 핵 유전부호 (chrI–chrXVI)</option>
                      <option value={3}>Table 3: 효모 미토콘드리아 (chrM · TGA=W)</option>
                    </select>
                    <select
                      value={frame}
                      onChange={(e) => setFrame(Number(e.target.value) as 0 | 1 | 2)}
                      className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                    >
                      <option value={0}>Frame +1</option>
                      <option value={1}>Frame +2</option>
                      <option value={2}>Frame +3</option>
                    </select>
                    <button
                      type="button"
                      onClick={() => copyText('protein_faa', extractedModel.proteinMultiFasta)}
                      className="inline-flex items-center gap-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                    >
                      {copiedKey === 'protein_faa' ? (
                        <Check className="h-3 w-3 text-emerald-600" />
                      ) : (
                        <Copy className="h-3 w-3" />
                      )}
                      단백질 FASTA 복사
                    </button>
                    <button
                      type="button"
                      onClick={() =>
                        downloadTextFile(
                          `${rawAnalysis.chrom}_${rawAnalysis.start + 1}-${rawAnalysis.end}_protein.faa`,
                          extractedModel.proteinMultiFasta,
                        )
                      }
                      className="inline-flex items-center gap-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                    >
                      <Download className="h-3 w-3" />
                      .faa 저장
                    </button>
                  </div>
                </div>

                {/* Per-Strain Protein Summary & Clickable AA Mutations */}
                <div className="mb-3 space-y-1.5">
                  <div className="flex flex-wrap items-center gap-2 rounded-lg border border-slate-200 bg-white px-3 py-1.5 text-xs">
                    <span className="font-semibold text-slate-700">Reference 단백질:</span>
                    <span className="font-mono text-slate-800">
                      {extractedModel.refAA.length.toLocaleString()} aa
                      {extractedModel.refAA.endsWith('*') ? ' (종결 코돈 * 포함)' : ''}
                    </span>
                  </div>

                  {extractedModel.strains.map((s) => {
                    const pe = s.proteinEval;
                    const lenDiff = pe.aa.length - extractedModel.refAA.length;
                    return (
                      <div
                        key={s.key}
                        className="space-y-1.5 rounded-lg border border-slate-200 bg-white px-3 py-2 text-xs"
                      >
                        <div className="flex flex-wrap items-center justify-between gap-2">
                          <div className="flex flex-wrap items-center gap-2">
                            <span className={`font-bold ${s.style.badge}`}>{s.title}</span>
                            <span
                              className={`rounded px-1.5 py-0.5 text-[10px] font-bold ${
                                pe.status === 'Identical'
                                  ? 'bg-emerald-100 text-emerald-800'
                                  : pe.status === 'Synonymous'
                                    ? 'bg-teal-100 text-teal-800'
                                    : pe.status === 'Nonsense' || pe.status === 'Frameshift'
                                      ? 'bg-rose-100 text-rose-800'
                                      : 'bg-sky-100 text-sky-900'
                              }`}
                            >
                              {pe.status}
                            </span>
                            <span className="font-mono text-[11px] text-slate-600">
                              {pe.aa.length.toLocaleString()} aa
                              {lenDiff !== 0 ? ` (${lenDiff > 0 ? `+${lenDiff}` : lenDiff} aa)` : ''}
                            </span>
                            <span className="text-[11px] text-slate-600">
                              · 아미노산 치환 {pe.realChangeCount}건
                              {pe.synonymousPositions.length > 0
                                ? ` · 동의 코돈 ${pe.synonymousPositions.length}건`
                                : ''}
                              {pe.ambiguousXCount > 0
                                ? ` · 미매핑/모호 X ${pe.ambiguousXCount}건`
                                : ''}
                            </span>
                          </div>
                        </div>

                        {pe.details.length > 0 && (
                          <div className="flex flex-wrap items-center gap-1">
                            <span className="mr-1 text-[10px] font-semibold text-slate-500">
                              아미노산 변이 (클릭 시 위치 이동):
                            </span>
                            {pe.details.slice(0, 60).map((d: AaChangeDetail) => (
                              <button
                                key={`${s.key}_${d.aaPos1}`}
                                type="button"
                                onClick={() => jumpToAaPos(d.aaPos1, true)}
                                className={`rounded border px-1.5 py-0.5 font-mono text-[10px] font-semibold transition ${
                                  focusedAaPos === d.aaPos1 - 1
                                    ? 'border-amber-500 bg-amber-200 text-amber-950'
                                    : d.isStopGain
                                      ? 'border-rose-300 bg-rose-50 text-rose-800 hover:bg-rose-100'
                                      : d.isAmbiguousGap
                                        ? 'border-amber-200 bg-amber-50 text-amber-800 hover:bg-amber-100'
                                        : 'border-sky-200 bg-sky-50 text-sky-900 hover:bg-sky-100'
                                }`}
                                title={`AA ${d.aaPos1}: ${d.refCodon}(${d.refAA}) → ${d.strainCodon}(${d.strainAA}) [DNA 코돈 Pos ${frame + (d.aaPos1 - 1) * 3 + 1}로 이동]`}
                              >
                                {d.label} ({d.refCodon}→{d.strainCodon})
                              </button>
                            ))}
                            {pe.details.length > 60 && (
                              <span className="text-[10px] text-slate-500">
                                외 {pe.details.length - 60}건
                              </span>
                            )}
                          </div>
                        )}
                      </div>
                    );
                  })}
                </div>

                {/* Full AA Block Alignment */}
                <SequenceBlockAlignment
                  title="아미노산 정렬"
                  unitLabel="aa"
                  refLabel="Ref AA"
                  refSeq={extractedModel.refAA}
                  pageSizeUnits={AA_PAGE_SIZE}
                  diffOnly={aaDiffOnly}
                  onToggleDiffOnly={() => {
                    setAaDiffOnly((v) => !v);
                    setAaPage(0);
                  }}
                  page={aaPage}
                  onPageChange={setAaPage}
                  focusedCol0={focusedAaPos}
                  formatMatchSymbol={(colIdx0, rc, sChars) => {
                    if (sChars.some((sc) => sc !== rc)) return '*';
                    const pos1 = colIdx0 + 1;
                    const hasSyn = extractedModel.strains.some((s) =>
                      s.proteinEval.synonymousPositions.includes(pos1),
                    );
                    return hasSyn ? '.' : '|';
                  }}
                  formatBlockCoord={(s0, e0) => {
                    const dnaS = frame + s0 * 3 + 1;
                    const dnaE = Math.min(
                      extractedModel.ungappedRefOriented.length,
                      frame + e0 * 3,
                    );
                    return `DNA ${dnaS.toLocaleString()}–${dnaE.toLocaleString()} nt`;
                  }}
                  emptyMessage="번역 가능한 길이(3 bp 이상)가 아닙니다."
                  strainRows={extractedModel.strains.map((s) => {
                    const synSet = new Set(s.proteinEval.synonymousPositions);
                    const detailMap = new Map(
                      s.proteinEval.details.map((d) => [d.aaPos1 - 1, d]),
                    );
                    return {
                      key: s.key,
                      label: s.title,
                      seq: s.proteinEval.aa,
                      labelClass: s.style.labelClass,
                      classifyChar: (ch, refCh, colIdx0) => {
                        const detail = detailMap.get(colIdx0);
                        if (detail) {
                          const codonNote = `${detail.refCodon}(${detail.refAA}) → ${detail.strainCodon}(${detail.strainAA})`;
                          if (detail.isAmbiguousGap) {
                            return {
                              className: 'bg-amber-200 font-bold text-amber-950',
                              title: `AA ${colIdx0 + 1}: 미매핑/IUPAC X (${codonNote})`,
                            };
                          }
                          if (detail.isStopGain) {
                            return {
                              className: 'bg-rose-200 font-bold text-rose-950',
                              title: `AA ${colIdx0 + 1}: 조기 종결 (${codonNote})`,
                            };
                          }
                          return {
                            className: s.style.diffClass,
                            title: `AA ${colIdx0 + 1}: ${detail.label} (${codonNote})`,
                          };
                        }
                        if (ch !== refCh) {
                          return {
                            className: s.style.diffClass,
                            title: `AA ${colIdx0 + 1}: ${refCh || '-'} → ${ch}`,
                          };
                        }
                        if (synSet.has(colIdx0 + 1)) {
                          const rc = extractedModel.refCodons[colIdx0] ?? '';
                          const sc = s.proteinEval.codons[colIdx0] ?? '';
                          return {
                            className: 'bg-teal-100 text-teal-900',
                            title: `AA ${colIdx0 + 1}: 동의 변이 (${rc}→${sc}, ${ch})`,
                          };
                        }
                        return { className: 'text-slate-600' };
                      },
                    };
                  })}
                />
              </div>
            </>
          )}
        </div>

        {/* Footer Export Actions */}
        {rawAnalysis && extractedModel && (
          <div className="flex flex-wrap items-center justify-between gap-3 border-t border-slate-200 bg-slate-50 px-6 py-3.5">
            <div className="flex flex-wrap items-center gap-2">
              <button
                type="button"
                onClick={() => copyText('multi_fasta', extractedModel.multiFasta)}
                className="inline-flex items-center gap-1.5 rounded-xl bg-sky-600 px-4 py-2 text-xs font-semibold text-white hover:bg-sky-700"
              >
                {copiedKey === 'multi_fasta' ? (
                  <Check className="h-4 w-4" />
                ) : (
                  <Copy className="h-4 w-4" />
                )}
                {extractedModel.strains.length === 0
                  ? 'Reference FASTA 복사'
                  : extractedModel.strains.length === 1
                    ? 'FASTA 전체 복사 (Ref + 균주1)'
                    : `Multi-FASTA 전체 복사 (Ref + 균주 ${extractedModel.strains.length}개)`}
              </button>

              <button
                type="button"
                onClick={() =>
                  downloadTextFile(
                    `${rawAnalysis.chrom}_${rawAnalysis.start + 1}-${rawAnalysis.end}.fasta`,
                    extractedModel.multiFasta,
                  )
                }
                className="inline-flex items-center gap-1.5 rounded-xl border border-slate-300 bg-white px-3.5 py-2 text-xs font-semibold text-slate-700 hover:bg-slate-100"
              >
                <Download className="h-3.5 w-3.5" />
                .fasta 저장
              </button>

              <button
                type="button"
                onClick={() => {
                  const header = hasStrain2
                    ? 'Chrom\tStart_1based\tEnd_1based\tType\tCategory\tRef\tStrain1\tStrain2\tS1_Depth\tS2_Depth'
                    : 'Chrom\tStart_1based\tEnd_1based\tType\tCategory\tRef\tStrain1\tS1_Depth';
                  const rows = activeAnalysis!.variants.map((v) =>
                    hasStrain2
                      ? `${v.chrom}\t${v.start + 1}\t${v.end}\t${v.type}\t${v.category}\t${v.refAllele}\t${v.s1Allele}\t${v.s2Allele}\t${v.s1Depth}\t${v.s2Depth}`
                      : `${v.chrom}\t${v.start + 1}\t${v.end}\t${v.type}\t${v.category}\t${v.refAllele}\t${v.s1Allele}\t${v.s1Depth}`,
                  );
                  downloadTextFile(
                    `${rawAnalysis.chrom}_${rawAnalysis.start + 1}-${rawAnalysis.end}_variants.tsv`,
                    [header, ...rows].join('\n'),
                    'text/tab-separated-values;charset=utf-8',
                  );
                }}
                className="inline-flex items-center gap-1.5 rounded-xl border border-slate-300 bg-white px-3.5 py-2 text-xs font-semibold text-slate-700 hover:bg-slate-100"
              >
                <Download className="h-3.5 w-3.5" />
                구간 변이 .tsv 저장
              </button>
            </div>

            <div className="flex flex-wrap items-center gap-2">
              {onSendToLocalBlast && (
                <button
                  type="button"
                  onClick={() => {
                    const targetSeq =
                      extractedModel.strains.length > 0
                        ? extractedModel.strains[0].cardSeq.replace(/-/g, '')
                        : extractedModel.refCardSeq.replace(/-/g, '');
                    const label =
                      extractedModel.strains.length > 0
                        ? `${extractedModel.strains[0].title}_${extractedModel.coordTag}`
                        : `Ref_${extractedModel.coordTag}`;
                    onSendToLocalBlast(`>${label}\n${targetSeq}`);
                  }}
                  className="inline-flex items-center gap-1.5 rounded-xl border border-indigo-300 bg-indigo-50 px-3.5 py-2 text-xs font-semibold text-indigo-800 hover:bg-indigo-100"
                >
                  로컬 BLAST로 매핑 확인
                </button>
              )}

              {extractedModel.strains.length > 0 && (
                <a
                  href={`https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY=${encodeURIComponent(
                    extractedModel.strains[0].cardSeq.replace(/-/g, '').slice(0, 7000),
                  )}`}
                  target="_blank"
                  rel="noreferrer"
                  className="inline-flex items-center gap-1.5 rounded-xl border border-emerald-300 bg-emerald-50 px-3.5 py-2 text-xs font-semibold text-emerald-800 hover:bg-emerald-100"
                  title="GET URL 길이 제한으로 최대 7,000 bp까지만 자동 입력됩니다"
                >
                  <ExternalLink className="h-3.5 w-3.5" />
                  균주1 서열 NCBI BLAST
                </a>
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
