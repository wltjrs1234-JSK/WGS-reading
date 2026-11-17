import React, { useEffect, useMemo, useState } from 'react';
import {
  AlertCircle,
  Compass,
  Dna,
  Play,
  Radar,
  Square,
  Trash2,
  X,
} from 'lucide-react';
import type { ConsensusOptions } from '../bio/types';
import type {
  BlastHitWithStrains,
  BlastProgressEvent,
  BlastRunPayload,
  BlastSearchResult,
} from '../workers/wgsProtocol';

interface BlastModalProps {
  isOpen: boolean;
  isRunning: boolean;
  progress: BlastProgressEvent | null;
  errorMsg: string | null;
  result: BlastSearchResult | null;
  initialQuery: string;
  referenceLabel: string;
  strainLabels: string[];
  currentChrom: string;
  options: ConsensusOptions;
  onClose: () => void;
  onRunBlast: (payload: BlastRunPayload) => void;
  onCancelBlast: () => void;
  onJumpToHit: (chrom: string, start: number, end: number) => void;
}

const STRAIN_PILL_STYLES = [
  'border-sky-200 bg-sky-50 text-sky-900',
  'border-emerald-200 bg-emerald-50 text-emerald-900',
  'border-violet-200 bg-violet-50 text-violet-900',
  'border-indigo-200 bg-indigo-50 text-indigo-900',
];

function renderAlignedDiffLine(
  seq: string,
  refComparisonSeq: string,
  mismatchClass: string,
): React.ReactNode[] {
  const nodes: React.ReactNode[] = [];
  let runStart = 0;
  let prevState = 0; // 0=match, 1=mismatch, 2=lowDepthN, 3=gap

  const getState = (i: number): number => {
    const ch = seq[i];
    const rc = refComparisonSeq[i] ?? '';
    if (ch === 'N') return 2;
    if (ch === '-') return ch !== rc ? 3 : 0;
    if (ch !== rc) return 1;
    return 0;
  };

  const stateClass = (st: number): string => {
    if (st === 1) return mismatchClass;
    if (st === 2) return 'bg-amber-200 font-bold text-amber-950';
    if (st === 3) return 'bg-rose-200 font-bold text-rose-950';
    return 'text-slate-700';
  };

  if (seq.length === 0) return nodes;
  prevState = getState(0);

  for (let i = 1; i < seq.length; i++) {
    const st = getState(i);
    if (st !== prevState) {
      nodes.push(
        <span key={runStart} className={stateClass(prevState)}>
          {seq.slice(runStart, i)}
        </span>,
      );
      runStart = i;
      prevState = st;
    }
  }
  nodes.push(
    <span key={runStart} className={stateClass(prevState)}>
      {seq.slice(runStart)}
    </span>,
  );
  return nodes;
}

const BlastHitAlignmentDetail: React.FC<{ hit: BlastHitWithStrains }> = ({ hit }) => {
  const blockSize = 60;
  const blocks = useMemo(() => {
    const len = hit.alignedQueryDisplay.length;
    const out: Array<{ start: number; end: number }> = [];
    for (let s = 0; s < len; s += blockSize) {
      out.push({ start: s, end: Math.min(len, s + blockSize) });
    }
    return out;
  }, [hit]);

  return (
    <div className="space-y-2 rounded-xl border border-slate-200 bg-slate-50 p-3.5">
      <div className="flex flex-wrap items-center justify-between gap-2 text-[13px]">
        <div className="flex flex-wrap items-center gap-2">
          <span className="font-bold text-slate-900">
            선택 Hit 정렬 상세 ({hit.chrom}:{(hit.refStart + 1).toLocaleString()}–
            {hit.refEnd.toLocaleString()} · {hit.strand === '+' ? '정방향 (+)' : '역상보 (-)'})
          </span>
          <span className="font-mono text-xs text-slate-700">
            Query {(hit.qStart + 1).toLocaleString()}–{hit.qEnd.toLocaleString()} /{' '}
            {hit.queryLen.toLocaleString()} bp · Ref 일치 {hit.identities}/{hit.alignLen} (
            {hit.identityPct}%) · Mismatch {hit.mismatches} · Gap {hit.gaps}
          </span>
        </div>
        <div className="flex flex-wrap items-center gap-2.5 text-xs text-slate-700">
          <span className="font-semibold text-indigo-700">■ Ref vs Query 차이</span>
          <span className="font-semibold text-sky-700">■ 균주 vs Query 차이</span>
          <span className="font-semibold text-amber-700">■ 저심도(N)</span>
          <span className="font-semibold text-rose-700">■ 결실(-)</span>
        </div>
      </div>

      <div className="max-h-80 space-y-3 overflow-x-auto overflow-y-auto rounded-lg border border-slate-200 bg-white p-3.5 font-mono text-[14px] leading-[22px] tracking-wide">
        {blocks.slice(0, 120).map((b) => {
          const qSlice = hit.alignedQueryDisplay.slice(b.start, b.end);
          const rSlice = hit.alignedRefDisplay.slice(b.start, b.end);
          const midSlice = qSlice
            .split('')
            .map((qc, i) => (qc === rSlice[i] && qc !== '-' ? '|' : ' '))
            .join('');

          return (
            <div key={b.start} className="space-y-1 border-b border-slate-100 pb-2.5 last:border-b-0 last:pb-0">
              <div className="text-xs font-medium text-slate-600">
                Col {(b.start + 1).toLocaleString()}–{b.end.toLocaleString()}
              </div>
              <div className="flex items-baseline whitespace-nowrap">
                <span className="mr-2.5 inline-block w-36 shrink-0 truncate border-r border-slate-200 pr-2.5 text-[13px] font-bold text-indigo-700">
                  Query
                </span>
                <span className="font-medium text-slate-900">{qSlice}</span>
              </div>
              <div className="flex items-baseline whitespace-nowrap text-slate-400 select-none">
                <span className="mr-2.5 inline-block w-36 shrink-0 truncate border-r border-slate-200 pr-2.5 font-sans text-xs font-medium text-slate-500">
                  Ref 일치
                </span>
                <span>{midSlice}</span>
              </div>
              <div className="flex items-baseline whitespace-nowrap">
                <span className="mr-2.5 inline-block w-36 shrink-0 truncate border-r border-slate-200 pr-2.5 text-[13px] font-bold text-emerald-700">
                  Reference
                </span>
                {renderAlignedDiffLine(
                  rSlice,
                  qSlice,
                  'bg-indigo-200 font-bold text-indigo-950',
                )}
              </div>
              {hit.strains.map((st, idx) => {
                const sSlice = st.alignedStrainDisplay.slice(b.start, b.end);
                const labelClr = idx % 2 === 0 ? 'text-sky-700' : 'text-teal-700';
                const diffClr =
                  idx % 2 === 0
                    ? 'bg-sky-200 font-bold text-sky-950'
                    : 'bg-teal-200 font-bold text-teal-950';
                return (
                  <div key={st.id} className="flex items-baseline whitespace-nowrap">
                    <span
                      className={`mr-2.5 inline-block w-36 shrink-0 truncate border-r border-slate-200 pr-2.5 text-[13px] font-bold ${labelClr}`}
                      title={st.label}
                    >
                      {st.label}
                    </span>
                    {renderAlignedDiffLine(sSlice, qSlice, diffClr)}
                  </div>
                );
              })}
            </div>
          );
        })}
        {blocks.length > 120 && (
          <div className="py-1.5 text-center font-sans text-xs text-slate-600">
            매우 긴 서열 정렬: 처음 7,200 컬럼만 미리 표시했습니다. (전체 상세는 뷰어 이동 후 서열 추출에서 확인 가능)
          </div>
        )}
      </div>
    </div>
  );
};

export const BlastModal: React.FC<BlastModalProps> = ({
  isOpen,
  isRunning,
  progress,
  errorMsg,
  result,
  initialQuery,
  referenceLabel,
  strainLabels,
  currentChrom,
  options,
  onClose,
  onRunBlast,
  onCancelBlast,
  onJumpToHit,
}) => {
  const [rawInput, setRawInput] = useState(initialQuery);
  const [minIdentityPct, setMinIdentityPct] = useState(75);
  const [primerMismatches, setPrimerMismatches] = useState(2);
  const [maxHits, setMaxHits] = useState(25);
  const [foreignScope, setForeignScope] = useState<'off' | 'current_chrom' | 'all_chroms'>('off');

  const [activeQueryIdx, setActiveQueryIdx] = useState(0);
  const [selectedHitIdx, setSelectedHitIdx] = useState(0);

  useEffect(() => {
    if (initialQuery) {
      setRawInput(initialQuery);
    }
  }, [initialQuery]);

  useEffect(() => {
    setActiveQueryIdx(0);
    setSelectedHitIdx(0);
  }, [result]);

  if (!isOpen) return null;

  const currentQueryRes = result?.queries[activeQueryIdx] ?? null;
  const selectedHit = currentQueryRes?.hits[selectedHitIdx] ?? null;

  const triggerSearch = () => {
    onRunBlast({
      rawInput,
      minIdentityPct,
      shortPrimerMaxMismatches: primerMismatches,
      maxHitsPerQuery: maxHits,
      foreignBamScope: strainLabels.length > 0 ? foreignScope : 'off',
      currentChrom,
      options,
    });
  };

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
            <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-xl bg-indigo-100 text-indigo-700">
              <Dna className="h-5 w-5" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-900">
                로컬 DNA BLAST (Reference + BAM 균주 동시 매핑 비교)
              </h2>
              <p className="text-xs text-slate-600">
                Reference: <b className="text-slate-800">{referenceLabel}</b> · 비교 대상 균주:{' '}
                <b className="text-sky-800">
                  {strainLabels.length > 0 ? strainLabels.join(', ') : '없음 (Reference 단독)'}
                </b>{' '}
                · 외부 업로드 없이 브라우저 내 실행
              </p>
            </div>
          </div>
          <button
            type="button"
            onClick={onClose}
            className="rounded-lg p-1.5 text-slate-500 hover:bg-slate-200 hover:text-slate-900"
            title="닫기 (Esc)"
          >
            <X className="h-5 w-5" />
          </button>
        </div>

        {/* Body */}
        <div className="flex-1 space-y-4 overflow-y-auto p-5">
          {/* Query Box */}
          <div className="space-y-2.5 rounded-xl border border-slate-200 bg-slate-50 p-4">
            <div className="flex flex-wrap items-center justify-between gap-2">
              <label className="text-xs font-bold text-slate-800">
                탐색할 DNA / 프라이머 / FASTA 서열 입력 (10 bp ~ 60,000 bp)
              </label>
              <div className="flex items-center gap-2 text-[11px]">
                <button
                  type="button"
                  onClick={() =>
                    setRawInput(
                      '>MET10_T511I_probe (chrVI+)\nATCGATGGTGAGACCATTAACGACGAGTCCAAATTG\n>ACT1_exon1_primer (chrVI-)\nAGGTTGCTGCTTTGGTTATTGATAACGGTTCT',
                    )
                  }
                  className="font-semibold text-indigo-700 hover:underline"
                >
                  예제 서열 채우기
                </button>
                {rawInput && (
                  <button
                    type="button"
                    onClick={() => setRawInput('')}
                    className="inline-flex items-center gap-1 text-slate-500 hover:text-rose-700"
                  >
                    <Trash2 className="h-3 w-3" /> 지우기
                  </button>
                )}
              </div>
            </div>

            <textarea
              rows={4}
              value={rawInput}
              onChange={(e) => setRawInput(e.target.value)}
              placeholder=">query_name&#10;ATGCCAGTTGAGTTTGCTACCAATCCTTTTGG... (또는 순수 염기서열 붙여넣기)"
              className="w-full rounded-lg border border-slate-300 bg-white p-2.5 font-mono text-xs text-slate-900 focus:border-indigo-500 focus:outline-none"
            />

            {/* Options & Action Row */}
            <div className="flex flex-wrap items-center justify-between gap-3">
              <div className="flex flex-wrap items-center gap-2.5 text-xs">
                <label className="flex items-center gap-1.5 text-slate-700">
                  <span>최소 일치율:</span>
                  <select
                    value={minIdentityPct}
                    onChange={(e) => setMinIdentityPct(Number(e.target.value))}
                    className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                  >
                    <option value={70}>70%</option>
                    <option value={75}>75%</option>
                    <option value={80}>80%</option>
                    <option value={90}>90%</option>
                    <option value={95}>95%</option>
                    <option value={100}>100%</option>
                  </select>
                </label>

                <label className="flex items-center gap-1.5 text-slate-700">
                  <span>프라이머 허용 불일치:</span>
                  <select
                    value={primerMismatches}
                    onChange={(e) => setPrimerMismatches(Number(e.target.value))}
                    className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                  >
                    <option value={0}>0 bp</option>
                    <option value={1}>1 bp</option>
                    <option value={2}>2 bp</option>
                    <option value={3}>3 bp</option>
                  </select>
                </label>

                <label className="flex items-center gap-1.5 text-slate-700">
                  <span>최대 Hit:</span>
                  <select
                    value={maxHits}
                    onChange={(e) => setMaxHits(Number(e.target.value))}
                    className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                  >
                    <option value={10}>10개</option>
                    <option value={25}>25개</option>
                    <option value={50}>50개</option>
                  </select>
                </label>

                {strainLabels.length > 0 && (
                  <label className="flex items-center gap-1.5 text-slate-700">
                    <Radar className="h-3.5 w-3.5 text-violet-600" />
                    <span>외래 서열(Soft-clip):</span>
                    <select
                      value={foreignScope}
                      onChange={(e) =>
                        setForeignScope(e.target.value as 'off' | 'current_chrom' | 'all_chroms')
                      }
                      className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                      title="Reference에 없는 카세트/벡터 서열을 BAM Soft-clip·Insertion 리드에서 직접 찾습니다"
                    >
                      <option value="off">끄기 (Ref+BAM 대조만)</option>
                      <option value="current_chrom">
                        현재 염색체 ({currentChrom === 'all' ? '전체' : currentChrom})
                      </option>
                      <option value="all_chroms">전체 염색체 스캔</option>
                    </select>
                  </label>
                )}
              </div>

              <div className="flex items-center gap-2">
                {isRunning ? (
                  <button
                    type="button"
                    onClick={onCancelBlast}
                    className="inline-flex items-center gap-1.5 rounded-xl bg-rose-600 px-4 py-2 text-xs font-semibold text-white hover:bg-rose-700"
                  >
                    <Square className="h-3.5 w-3.5" /> 검색 중지
                  </button>
                ) : (
                  <button
                    type="button"
                    disabled={!rawInput.trim()}
                    onClick={triggerSearch}
                    className="inline-flex items-center gap-1.5 rounded-xl bg-indigo-600 px-5 py-2 text-xs font-semibold text-white shadow-sm hover:bg-indigo-700 disabled:cursor-not-allowed disabled:bg-slate-300 disabled:text-slate-500"
                  >
                    <Play className="h-3.5 w-3.5" /> 로컬 BLAST 실행
                  </button>
                )}
              </div>
            </div>

            {isRunning && progress && (
              <div className="space-y-1 pt-1">
                <div className="h-1.5 overflow-hidden rounded-full bg-slate-200">
                  <div
                    className="h-full bg-indigo-600 transition-all"
                    style={{ width: `${Math.max(4, progress.percent)}%` }}
                  />
                </div>
                <div className="flex justify-between text-[11px] text-slate-600">
                  <span>{progress.message}</span>
                  <span className="font-mono">{progress.percent}%</span>
                </div>
              </div>
            )}
          </div>

          {errorMsg && (
            <div className="flex items-center gap-2 rounded-xl border border-rose-300 bg-rose-50 p-3 text-xs text-rose-900">
              <AlertCircle className="h-4 w-4 shrink-0 text-rose-600" />
              <span>{errorMsg}</span>
            </div>
          )}

          {/* Results */}
          {result && (
            <div className="space-y-4">
              {/* Query Tabs & Runtime Summary */}
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div className="flex flex-wrap items-center gap-1.5">
                  {result.queries.map((q, idx) => (
                    <button
                      key={`${q.queryId}_${idx}`}
                      type="button"
                      onClick={() => {
                        setActiveQueryIdx(idx);
                        setSelectedHitIdx(0);
                      }}
                      className={`rounded-lg border px-3 py-1.5 text-xs font-semibold transition ${
                        activeQueryIdx === idx
                          ? 'border-indigo-500 bg-indigo-600 text-white'
                          : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-100'
                      }`}
                    >
                      {q.queryId} ({q.queryLength.toLocaleString()} bp) · Ref Hit {q.hits.length}건
                      {q.foreignHits.length > 0 ? ` · 외래 ${q.foreignHits.length}건` : ''}
                    </button>
                  ))}
                </div>
                <span className="font-mono text-[11px] text-slate-500">
                  검색 지놈 {result.totalRefLength.toLocaleString()} bp · 소요 시간{' '}
                  {result.elapsedMs} ms
                </span>
              </div>

              {currentQueryRes && (
                <>
                  {/* Hit Comparison Table */}
                  {currentQueryRes.hits.length === 0 ? (
                    <div className="rounded-xl border border-slate-200 bg-slate-50 p-6 text-center text-xs text-slate-600">
                      Reference에서 일치율 {minIdentityPct}% 이상의 유의미한 매핑 구간을 찾지 못했습니다.
                      {strainLabels.length > 0 && foreignScope === 'off'
                        ? " (Reference에 없는 외래 삽입 서열이라면 '외래 서열(Soft-clip)' 옵션을 켜고 다시 검색해 보세요.)"
                        : ''}
                    </div>
                  ) : (
                    <div className="overflow-x-auto rounded-xl border border-slate-200">
                      <table className="w-full border-collapse text-left text-xs">
                        <thead>
                          <tr className="border-b border-slate-200 bg-slate-100 text-[11px] font-bold text-slate-700">
                            <th className="px-3 py-2">#</th>
                            <th className="px-3 py-2">매핑 위치 (Ref)</th>
                            <th className="px-3 py-2">겹치는 유전자</th>
                            <th className="px-3 py-2">Query 범위 / E-value</th>
                            <th className="px-3 py-2">Reference 일치</th>
                            {result.strainLabels.map((lbl, sIdx) => (
                              <th key={sIdx} className="px-3 py-2">
                                균주 대조: {lbl}
                              </th>
                            ))}
                            <th className="px-3 py-2 text-right">이동</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-200 bg-white">
                          {currentQueryRes.hits.map((hit, hIdx) => {
                            const isSelected = selectedHitIdx === hIdx;
                            return (
                              <tr
                                key={`${hit.chrom}_${hit.refStart}_${hit.strand}_${hIdx}`}
                                onClick={() => setSelectedHitIdx(hIdx)}
                                className={`cursor-pointer transition ${
                                  isSelected ? 'bg-indigo-50/70' : 'hover:bg-slate-50'
                                }`}
                              >
                                <td className="px-3 py-2 font-mono font-semibold text-slate-500">
                                  {hIdx + 1}
                                </td>
                                <td className="px-3 py-2 font-mono font-bold text-slate-900">
                                  {hit.chrom}:{(hit.refStart + 1).toLocaleString()}–
                                  {hit.refEnd.toLocaleString()}{' '}
                                  <span
                                    className={`ml-1 rounded px-1 py-0.5 text-[10px] ${
                                      hit.strand === '+'
                                        ? 'bg-slate-100 text-slate-700'
                                        : 'bg-violet-100 text-violet-800'
                                    }`}
                                  >
                                    {hit.strand}
                                  </span>
                                </td>
                                <td className="px-3 py-2">
                                  {hit.overlappingGenes.length > 0 ? (
                                    <div className="flex flex-wrap gap-1">
                                      {hit.overlappingGenes.slice(0, 3).map((g) => (
                                        <span
                                          key={g.id}
                                          className="rounded bg-blue-50 px-1.5 py-0.5 text-[10px] font-bold text-blue-800"
                                        >
                                          {g.name} ({g.strand})
                                        </span>
                                      ))}
                                    </div>
                                  ) : (
                                    <span className="text-[11px] text-slate-400">Intergenic</span>
                                  )}
                                </td>
                                <td className="px-3 py-2 font-mono text-[11px] text-slate-600">
                                  <div>
                                    {hit.qStart + 1}–{hit.qEnd} ({hit.queryCoveragePct}%)
                                  </div>
                                  <div className="text-[10px] text-slate-400">
                                    Bit {hit.bitScore} · E={hit.eValue}
                                  </div>
                                </td>
                                <td className="px-3 py-2">
                                  <span
                                    className={`inline-block rounded px-1.5 py-0.5 font-mono text-[11px] font-bold ${
                                      hit.identityPct === 100
                                        ? 'bg-emerald-100 text-emerald-800'
                                        : 'bg-amber-100 text-amber-900'
                                    }`}
                                  >
                                    {hit.identityPct}% ({hit.identities}/{hit.alignLen})
                                  </span>
                                </td>
                                {hit.strains.map((st, sIdx) => (
                                  <td key={st.id} className="px-3 py-2">
                                    {!st.matchedChrom || st.coveragePct === 0 ? (
                                      <span className="inline-block rounded border border-rose-200 bg-rose-50 px-1.5 py-0.5 text-[10px] font-bold text-rose-800">
                                        {st.matchedChrom ? '커버리지 0% (결실 가능)' : '염색체 미매핑'}
                                      </span>
                                    ) : (
                                      <div
                                        className={`space-y-0.5 rounded border px-2 py-1 text-[10px] ${
                                          STRAIN_PILL_STYLES[sIdx % STRAIN_PILL_STYLES.length]
                                        }`}
                                      >
                                        <div className="font-mono font-bold">
                                          Query 일치 {st.identityVsQueryPct}% ({st.identitiesVsQuery}/
                                          {st.comparedCols})
                                        </div>
                                        <div className="text-slate-600">
                                          커버 {st.coveragePct}% · {st.meanDepth}x
                                          {st.snpCountVsRef > 0 ? ` · Ref대비 변이 ${st.snpCountVsRef}bp` : ''}
                                          {st.softClipReads > 0 ? ` · Clip ${st.softClipReads}개` : ''}
                                        </div>
                                      </div>
                                    )}
                                  </td>
                                ))}
                                <td className="px-3 py-2 text-right">
                                  <button
                                    type="button"
                                    onClick={(e) => {
                                      e.stopPropagation();
                                      onJumpToHit(hit.chrom, hit.refStart, hit.refEnd);
                                    }}
                                    className="inline-flex items-center gap-1 rounded-lg bg-sky-600 px-2.5 py-1 text-[11px] font-semibold text-white hover:bg-sky-700"
                                  >
                                    <Compass className="h-3 w-3" /> 이동
                                  </button>
                                </td>
                              </tr>
                            );
                          })}
                        </tbody>
                      </table>
                    </div>
                  )}

                  {/* Selected Hit Sequence Alignment Detail */}
                  {selectedHit && <BlastHitAlignmentDetail hit={selectedHit} />}

                  {/* Foreign / Soft-clip Read Matches */}
                  {currentQueryRes.foreignHits.length > 0 && (
                    <div className="space-y-2 rounded-xl border border-violet-200 bg-violet-50/50 p-4">
                      <div className="text-xs font-bold text-violet-900">
                        BAM Soft-clip / 삽입 리드 내 직접 탐색 결과 ({currentQueryRes.foreignHits.length}곳)
                      </div>
                      <div className="space-y-2">
                        {currentQueryRes.foreignHits.map((fh, idx) => (
                          <div
                            key={`${fh.strainId}_${fh.chrom}_${fh.anchorPos}_${idx}`}
                            className="flex flex-wrap items-center justify-between gap-2 rounded-lg border border-violet-200 bg-white p-2.5 text-xs"
                          >
                            <div className="space-y-1">
                              <div className="flex flex-wrap items-center gap-2">
                                <span className="rounded bg-violet-100 px-1.5 py-0.5 text-[10px] font-bold text-violet-800">
                                  {fh.strainLabel}
                                </span>
                                <span className="font-mono font-bold text-slate-900">
                                  {fh.chrom}:{(fh.anchorPos + 1).toLocaleString()} (
                                  {fh.clipSide === 'left_clip'
                                    ? '좌측 Soft-clip 접합'
                                    : fh.clipSide === 'right_clip'
                                      ? '우측 Soft-clip 접합'
                                      : 'Insertion 리드'}
                                  )
                                </span>
                                <span className="font-semibold text-emerald-700">
                                  지지 리드 {fh.supportingReads}개 · 일치 {fh.identityPct}% ({fh.matchLen} bp)
                                </span>
                              </div>
                              <div className="font-mono text-[10px] text-slate-500">
                                대표 리드: {fh.sampleReadName} ({fh.sampleCigar}) ·{' '}
                                {fh.alignedQueryDisplay.slice(0, 60)}
                                {fh.alignedQueryDisplay.length > 60 ? '...' : ''}
                              </div>
                            </div>
                            <button
                              type="button"
                              onClick={() =>
                                onJumpToHit(
                                  fh.chrom,
                                  Math.max(0, fh.anchorPos - 60),
                                  fh.anchorPos + 60,
                                )
                              }
                              className="inline-flex items-center gap-1 rounded-lg bg-violet-600 px-2.5 py-1 text-[11px] font-semibold text-white hover:bg-violet-700"
                            >
                              <Compass className="h-3 w-3" /> 접합부 이동
                            </button>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
