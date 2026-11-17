import React, { useMemo, useState } from 'react';
import {
  Activity,
  AlertCircle,
  Compass,
  Download,
  GitBranch,
  Layers,
  Play,
  Square,
  X,
} from 'lucide-react';
import {
  exportCnvSvToBed,
  exportCnvSvToTsv,
  type ChromCopySummary,
  type CnvPloidyPreference,
  type CnvScanMode,
  type CnvSvAnalysisResult,
} from '../bio/cnvSvAnalyzer';
import { normalizeYeastChromName } from '../bio/yeastAliases';
import type { CnvProgressEvent, CnvSvPayload } from '../workers/wgsProtocol';

interface CnvSvModalProps {
  isOpen: boolean;
  isRunning: boolean;
  progress: CnvProgressEvent | null;
  errorMsg: string | null;
  result: CnvSvAnalysisResult | null;
  strainLabels: string[];
  onClose: () => void;
  onStartAnalyze: (payload: CnvSvPayload) => void;
  onCancelAnalyze: () => void;
  onJumpToRegion: (chrom: string, start: number, end: number) => void;
}

const STATUS_STYLE: Record<ChromCopySummary['status'], { label: string; badge: string; bar: string }> = {
  normal: {
    label: '정상',
    badge: 'bg-emerald-50 text-emerald-800 border-emerald-200',
    bar: 'bg-sky-500',
  },
  gain: {
    label: '복제수 증가 (+)',
    badge: 'bg-rose-100 text-rose-800 border-rose-300 font-bold',
    bar: 'bg-rose-500',
  },
  loss: {
    label: '복제수 감소 (-)',
    badge: 'bg-indigo-100 text-indigo-800 border-indigo-300 font-bold',
    bar: 'bg-indigo-500',
  },
  mitochondrial_low: {
    label: 'mtDNA 극저/Petite',
    badge: 'bg-amber-100 text-amber-900 border-amber-300 font-bold',
    bar: 'bg-amber-500',
  },
  mitochondrial_high: {
    label: 'mtDNA 고복제',
    badge: 'bg-emerald-100 text-emerald-900 border-emerald-300',
    bar: 'bg-emerald-500',
  },
  unmapped: {
    label: '미매핑',
    badge: 'bg-slate-100 text-slate-500 border-slate-200',
    bar: 'bg-slate-300',
  },
};

function downloadText(filename: string, content: string, mime: string) {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export const CnvSvModal: React.FC<CnvSvModalProps> = ({
  isOpen,
  isRunning,
  progress,
  errorMsg,
  result,
  strainLabels,
  onClose,
  onStartAnalyze,
  onCancelAnalyze,
  onJumpToRegion,
}) => {
  const [scanMode, setScanMode] = useState<CnvScanMode>('fast_sample');
  const [ploidyPref, setPloidyPref] = useState<CnvPloidyPreference>('auto');
  const [targetStrainIdx, setTargetStrainIdx] = useState<number | null>(null);
  const [excludeRdnaSv, setExcludeRdnaSv] = useState(true);
  const [activeStrainTab, setActiveStrainTab] = useState(0);
  const [activeSubTab, setActiveSubTab] = useState<'cnv' | 'loh' | 'sv'>('cnv');

  const activeProfile = useMemo(() => {
    if (!result || result.strains.length === 0) return null;
    return result.strains[activeStrainTab] ?? result.strains[0];
  }, [result, activeStrainTab]);

  if (!isOpen) return null;

  const triggerRun = () => {
    setActiveStrainTab(0);
    onStartAnalyze({
      scanMode,
      ploidyPreference: ploidyPref,
      targetStrainIndex: targetStrainIdx,
      excludeRdnaSv,
      minSvSupport: 4,
    });
  };

  const bafMax = activeProfile ? Math.max(1, ...activeProfile.bafHistogram) : 1;

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/50 p-4 backdrop-blur-xs"
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div className="flex max-h-[93vh] w-full max-w-6xl flex-col overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-2xl">
        {/* Header */}
        <div className="flex items-center justify-between gap-3 border-b border-slate-200 bg-slate-50 px-6 py-3.5">
          <div className="flex items-center gap-3">
            <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-xl bg-teal-100 text-teal-700">
              <Layers className="h-5 w-5" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-900">
                염색체 복제수(CNV · 이수성) · LOH · 거대 전좌(SV) 탐지기
              </h2>
              <p className="text-xs text-slate-600">
                <strong className="font-semibold text-slate-800">반수체(1n)·이배체(2n)·사배체(4n)</strong>{' '}
                배수성을 자동 판별해 염색체 이수성, BAF/LOH, 염색체 간 전좌 브레이크포인트를 분석합니다.
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
        <div className="flex flex-1 flex-col space-y-4 overflow-y-auto p-5">
          {/* Controls */}
          <div className="space-y-3 rounded-xl border border-slate-200 bg-slate-50 p-4">
            <div className="grid gap-3 md:grid-cols-3">
              <label className="text-xs">
                <span className="font-semibold text-slate-700">배수성 기준 (Ploidy)</span>
                <select
                  value={String(ploidyPref)}
                  onChange={(e) => {
                    const v = e.target.value;
                    setPloidyPref(v === 'auto' ? 'auto' : (Number(v) as 1 | 2 | 4));
                  }}
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value="auto">BAF 자동 추정 ★ (1n · 2n · 4n 자동 판별)</option>
                  <option value="1">반수체 고정 (1n · 실험실 표준 균주)</option>
                  <option value="2">이배체 고정 (2n · 야생/이배체)</option>
                  <option value="4">사배체 고정 (4n · 배수체/산업 균주)</option>
                </select>
              </label>

              <label className="text-xs">
                <span className="font-semibold text-slate-700">분석 정밀도 모드</span>
                <select
                  value={scanMode}
                  onChange={(e) => setScanMode(e.target.value as CnvScanMode)}
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value="fast_sample">고속 CNV + BAF + 샘플 SV 스캔 (2~4초)</option>
                  <option value="deep_sv">정밀 전체 염색체 전좌(SV) + CNV 스캔</option>
                </select>
              </label>

              <label className="text-xs">
                <span className="font-semibold text-slate-700">대상 균주</span>
                <select
                  value={targetStrainIdx === null ? '' : String(targetStrainIdx)}
                  onChange={(e) =>
                    setTargetStrainIdx(e.target.value === '' ? null : Number(e.target.value))
                  }
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value="">
                    {strainLabels.length > 1
                      ? `전체 균주 (${strainLabels.length}개) 동시 비교`
                      : strainLabels[0] ?? '로드된 균주 없음'}
                  </option>
                  {strainLabels.map((lbl, idx) => (
                    <option key={idx} value={String(idx)}>
                      {lbl} 단독 분석
                    </option>
                  ))}
                </select>
              </label>
            </div>

            <div className="flex flex-wrap items-center justify-between gap-3 pt-1">
              <label className="inline-flex cursor-pointer items-center gap-1.5 text-[11px] font-medium text-slate-700">
                <input
                  type="checkbox"
                  checked={excludeRdnaSv}
                  onChange={(e) => setExcludeRdnaSv(e.target.checked)}
                  className="rounded border-slate-300 text-teal-600 focus:ring-teal-500"
                />
                chrXII rDNA 배열 내부 자가 전좌 노이즈 제외
              </label>

              <div className="flex items-center gap-2">
                {isRunning ? (
                  <button
                    type="button"
                    onClick={onCancelAnalyze}
                    className="inline-flex items-center gap-1.5 rounded-xl bg-rose-600 px-4 py-2 text-xs font-semibold text-white hover:bg-rose-700"
                  >
                    <Square className="h-3.5 w-3.5" /> 분석 중지
                  </button>
                ) : (
                  <button
                    type="button"
                    disabled={strainLabels.length === 0}
                    onClick={triggerRun}
                    className="inline-flex items-center gap-1.5 rounded-xl bg-teal-600 px-5 py-2 text-xs font-semibold text-white shadow-xs hover:bg-teal-700 disabled:cursor-not-allowed disabled:bg-slate-300 disabled:text-slate-500"
                  >
                    <Play className="h-3.5 w-3.5" /> CNV · LOH · SV 분석 시작
                  </button>
                )}
              </div>
            </div>

            {isRunning && progress && (
              <div className="space-y-1 pt-1">
                <div className="h-1.5 overflow-hidden rounded-full bg-slate-200">
                  <div
                    className="h-full bg-teal-500 transition-all"
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

          {result && activeProfile && (
            <div className="space-y-3">
              {/* Strain switcher + Export */}
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div className="flex flex-wrap items-center gap-1.5">
                  {result.strains.map((st, idx) => (
                    <button
                      key={st.strainId}
                      type="button"
                      onClick={() => setActiveStrainTab(idx)}
                      className={`rounded-lg border px-3 py-1.5 text-xs font-semibold transition ${
                        idx === activeStrainTab
                          ? 'border-teal-600 bg-teal-600 text-white'
                          : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-100'
                      }`}
                    >
                      {st.strainLabel} ({st.effectivePloidy}n)
                    </button>
                  ))}
                  <span className="ml-2 font-mono text-[11px] text-slate-500">
                    ({result.elapsedMs} ms)
                  </span>
                </div>

                <div className="flex items-center gap-2">
                  <button
                    type="button"
                    onClick={() =>
                      downloadText(
                        `yeast_cnv_sv_${result.scanMode}.tsv`,
                        exportCnvSvToTsv(result),
                        'text/tab-separated-values;charset=utf-8',
                      )
                    }
                    className="inline-flex items-center gap-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                  >
                    <Download className="h-3.5 w-3.5" /> .tsv 저장
                  </button>
                  <button
                    type="button"
                    onClick={() =>
                      downloadText(
                        `yeast_cnv_sv_${result.scanMode}.bed`,
                        exportCnvSvToBed(result),
                        'text/plain;charset=utf-8',
                      )
                    }
                    className="inline-flex items-center gap-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                  >
                    <Download className="h-3.5 w-3.5" /> .bed 저장
                  </button>
                </div>
              </div>

              {/* Ploidy + BAF Telemetry Card */}
              <div className="grid gap-3 rounded-xl border border-slate-200 bg-slate-50/70 p-3.5 md:grid-cols-4">
                <div className="md:col-span-2 space-y-1">
                  <div className="flex items-center gap-2">
                    <span className="rounded-md bg-teal-100 px-2 py-0.5 text-xs font-bold text-teal-900">
                      적용 배수성: {activeProfile.effectivePloidy}n
                    </span>
                    <span className="rounded-md bg-slate-200 px-2 py-0.5 text-[11px] font-medium text-slate-700">
                      자동 추정: {activeProfile.detectedPloidy}n ({activeProfile.ploidyConfidence})
                    </span>
                  </div>
                  <p className="text-[11px] text-slate-600">{activeProfile.ploidyReason}</p>
                  <div className="flex flex-wrap gap-3 pt-1 font-mono text-[11px] text-slate-700">
                    <span>핵 중앙값 깊이: {activeProfile.genomeMedianDepth}x</span>
                    <span>이형접합: {activeProfile.genomeMedianHetPerKb}/kb</span>
                    <span>mtDNA/핵 비율: {activeProfile.mitoToNuclearRatio}x</span>
                  </div>
                </div>

                {/* BAF 20-bin mini histogram */}
                <div className="md:col-span-2 flex flex-col justify-between">
                  <div className="flex items-center justify-between text-[11px] text-slate-600">
                    <span className="font-semibold">
                      이형접합 BAF 분포 (표본 {activeProfile.totalHetSitesSampled}개)
                    </span>
                    <span className="font-mono text-[10px] text-slate-400">0.25 · 0.50 · 0.75</span>
                  </div>
                  <div className="mt-1.5 flex h-11 items-end gap-0.5 rounded-lg border border-slate-200 bg-white px-1.5 py-1">
                    {activeProfile.bafHistogram.map((cnt, bIdx) => {
                      const hPct = Math.round((cnt / bafMax) * 100);
                      const isHalf = bIdx >= 8 && bIdx <= 11;
                      const isQuarter = (bIdx >= 4 && bIdx <= 6) || (bIdx >= 13 && bIdx <= 15);
                      return (
                        <div
                          key={bIdx}
                          title={`BAF ${(bIdx * 0.05).toFixed(2)}–${((bIdx + 1) * 0.05).toFixed(2)}: ${cnt}건`}
                          className={`flex-1 rounded-t-xs ${
                            isHalf ? 'bg-teal-500' : isQuarter ? 'bg-amber-400' : 'bg-slate-300'
                          }`}
                          style={{ height: `${Math.max(cnt > 0 ? 12 : 4, hPct)}%` }}
                        />
                      );
                    })}
                  </div>
                </div>
              </div>

              {/* Sub-navigation Tabs */}
              <div className="flex items-center gap-1.5 border-b border-slate-200 pb-2">
                <button
                  type="button"
                  onClick={() => setActiveSubTab('cnv')}
                  className={`inline-flex items-center gap-1.5 rounded-lg px-3 py-1.5 text-xs font-semibold ${
                    activeSubTab === 'cnv' ? 'bg-teal-600 text-white' : 'text-slate-600 hover:bg-slate-100'
                  }`}
                >
                  <Activity className="h-3.5 w-3.5" />
                  염색체 복제수 · 이수성 (구간 {activeProfile.cnvSegments.length})
                </button>
                <button
                  type="button"
                  onClick={() => setActiveSubTab('loh')}
                  className={`inline-flex items-center gap-1.5 rounded-lg px-3 py-1.5 text-xs font-semibold ${
                    activeSubTab === 'loh' ? 'bg-teal-600 text-white' : 'text-slate-600 hover:bg-slate-100'
                  }`}
                >
                  LOH 이형접합 소실 ({activeProfile.lohSegments.length})
                </button>
                <button
                  type="button"
                  onClick={() => setActiveSubTab('sv')}
                  className={`inline-flex items-center gap-1.5 rounded-lg px-3 py-1.5 text-xs font-semibold ${
                    activeSubTab === 'sv' ? 'bg-teal-600 text-white' : 'text-slate-600 hover:bg-slate-100'
                  }`}
                >
                  <GitBranch className="h-3.5 w-3.5" />
                  전좌 · 거대 SV ({activeProfile.svEvents.length})
                </button>
              </div>

              {/* Tab Content */}
              {activeSubTab === 'cnv' && (
                <div className="space-y-3">
                  {/* Chromosome Ploidy Bar Chart */}
                  <div className="rounded-xl border border-slate-200 bg-white p-3.5">
                    <div className="mb-2.5 flex items-center justify-between text-[13px]">
                      <span className="font-bold text-slate-900">
                        염색체별 추정 복제수 (기준 {activeProfile.effectivePloidy}.0n)
                      </span>
                      <span className="text-xs text-slate-600">
                        클릭 시 해당 염색체로 이동 · chrM은 핵 대비 비율({activeProfile.mitoToNuclearRatio}x)
                      </span>
                    </div>
                    <div className="grid grid-cols-6 gap-1.5 sm:grid-cols-9 md:grid-cols-[repeat(17,minmax(0,1fr))]">
                      {activeProfile.chromSummaries.map((cs) => {
                        const st = STATUS_STYLE[cs.status];
                        const shortName = normalizeYeastChromName(cs.chrom).replace(/^chr/, '');
                        const isMito = cs.chrom.toLowerCase().includes('chrm');
                        const relHeight = isMito
                          ? Math.min(100, Math.max(8, Math.round(Math.log10(Math.max(1, cs.copyRatio)) * 45)))
                          : Math.min(
                              100,
                              Math.max(8, Math.round((cs.estCopies / (activeProfile.effectivePloidy * 2)) * 100)),
                            );
                        return (
                          <button
                            key={cs.chrom}
                            type="button"
                            onClick={() => {
                              onClose();
                              onJumpToRegion(cs.chrom, 0, Math.min(cs.length, 30_000));
                            }}
                            className="flex flex-col items-center rounded-lg border border-slate-200 bg-slate-50/70 p-1.5 hover:border-teal-400 hover:bg-teal-50/50"
                            title={`${cs.chrom}: ${cs.estCopies} copies (Depth ${cs.medianDepth}x, ratio ${cs.copyRatio})`}
                          >
                            <div className="relative flex h-16 w-full items-end justify-center rounded bg-slate-100 px-1">
                              {!isMito && (
                                <div
                                  className="pointer-events-none absolute inset-x-0 border-t border-dashed border-slate-400"
                                  style={{ bottom: '50%' }}
                                />
                              )}
                              <div
                                className={`w-full rounded-t-xs ${st.bar}`}
                                style={{ height: `${relHeight}%` }}
                              />
                            </div>
                            <span className="mt-1 font-mono text-[11px] font-bold text-slate-900">
                              {shortName}
                            </span>
                            <span className="font-mono text-[11px] text-slate-700">
                              {isMito ? `${cs.copyRatio}x` : `${cs.estCopies}n`}
                            </span>
                          </button>
                        );
                      })}
                    </div>
                  </div>

                  {/* Aneuploidy & Segment Table */}
                  <div className="overflow-hidden rounded-xl border border-slate-200 bg-white">
                    <table className="w-full border-collapse text-left text-[13px]">
                      <thead className="border-b border-slate-200 bg-slate-100 text-xs font-bold text-slate-800">
                        <tr>
                          <th className="px-3 py-2.5">염색체</th>
                          <th className="px-3 py-2.5">판정</th>
                          <th className="px-3 py-2.5">추정 복제수</th>
                          <th className="px-3 py-2.5">핵 기준 배율</th>
                          <th className="px-3 py-2.5">중앙 Depth</th>
                          <th className="px-3 py-2.5">Het/kb</th>
                          <th className="px-3 py-2.5">비고</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-slate-100">
                        {activeProfile.chromSummaries
                          .filter((c) => c.status !== 'normal' || activeProfile.chromSummaries.length <= 20)
                          .map((cs) => {
                            const st = STATUS_STYLE[cs.status];
                            return (
                              <tr key={cs.chrom} className="hover:bg-slate-50">
                                <td className="px-3 py-2 font-mono font-bold text-slate-900">{cs.chrom}</td>
                                <td className="px-3 py-2">
                                  <span className={`rounded-md border px-2 py-0.5 text-xs font-semibold ${st.badge}`}>
                                    {st.label}
                                  </span>
                                </td>
                                <td className="px-3 py-2 font-mono font-semibold text-slate-900">
                                  {cs.estCopies}x
                                </td>
                                <td className="px-3 py-2 font-mono text-slate-700">{cs.copyRatio}x</td>
                                <td className="px-3 py-2 font-mono text-slate-700">{cs.medianDepth}x</td>
                                <td className="px-3 py-2 font-mono text-slate-700">{cs.hetPerKb}</td>
                                <td className="px-3 py-2 text-xs text-slate-700">{cs.note ?? '-'}</td>
                              </tr>
                            );
                          })}
                      </tbody>
                    </table>
                  </div>

                  {activeProfile.cnvSegments.length > 0 && (
                    <div className="space-y-2">
                      <h4 className="text-[13px] font-bold text-slate-900">
                        구간 증폭·결실 세그먼트 ({activeProfile.cnvSegments.length}건)
                      </h4>
                      {activeProfile.cnvSegments.map((seg) => (
                        <div
                          key={seg.id}
                          className="flex flex-wrap items-center justify-between gap-2 rounded-lg border border-slate-200 bg-white px-3.5 py-2.5 text-[13px]"
                        >
                          <div className="flex flex-wrap items-center gap-2">
                            <span
                              className={`rounded-md px-2 py-0.5 text-xs font-bold ${
                                seg.direction === 'amplification'
                                  ? 'bg-rose-100 text-rose-800'
                                  : 'bg-indigo-100 text-indigo-800'
                              }`}
                            >
                              {seg.direction === 'amplification'
                                ? '구간 증폭'
                                : seg.direction === 'homozygous_deletion'
                                  ? '완전 결실'
                                  : '구간 결실'}
                            </span>
                            <span className="font-mono font-bold text-slate-900">
                              {seg.chrom}:{(seg.start + 1).toLocaleString()}–{seg.end.toLocaleString()}
                            </span>
                            <span className="font-mono text-slate-600">
                              ({(seg.lengthBp / 1000).toFixed(1)} kb · 추정 {seg.estCopies}x)
                            </span>
                            {seg.overlappingGenes.length > 0 && (
                              <span className="text-[11px] text-slate-600">
                                유전자: {seg.overlappingGenes.join(', ')}
                              </span>
                            )}
                          </div>
                          <button
                            type="button"
                            onClick={() => {
                              onClose();
                              onJumpToRegion(seg.chrom, seg.start, seg.end);
                            }}
                            className="inline-flex items-center gap-1 rounded border border-teal-200 bg-teal-50 px-2 py-0.5 text-[11px] font-semibold text-teal-800 hover:bg-teal-100"
                          >
                            <Compass className="h-3 w-3" /> 이동
                          </button>
                        </div>
                      ))}
                    </div>
                  )}
                </div>
              )}

              {activeSubTab === 'loh' && (
                <div className="space-y-2">
                  {activeProfile.effectivePloidy === 1 ? (
                    <div className="rounded-xl border border-slate-200 bg-slate-50 p-6 text-center text-xs text-slate-600">
                      현재 선택된 배수성이 <strong>반수체(1n)</strong> 입니다. LOH(이형접합 소실)는 야생{' '}
                      <strong>이배체(2n)</strong> 또는 <strong>사배체(4n)</strong> 균주에서 의미가 있으며,
                      상단 배수성 기준을 2n/4n 또는 자동 추정에 두면 자동 탐지됩니다.
                    </div>
                  ) : activeProfile.lohSegments.length === 0 ? (
                    <div className="rounded-xl border border-emerald-200 bg-emerald-50/70 p-6 text-center text-xs text-emerald-900">
                      연속 30 kb 이상 이어지는 LOH(이형접합 소실) 블록이 발견되지 않았습니다.
                    </div>
                  ) : (
                    activeProfile.lohSegments.map((loh) => (
                      <div
                        key={loh.id}
                        className="flex flex-wrap items-center justify-between gap-2 rounded-xl border border-amber-200 bg-amber-50/40 px-4 py-2.5 text-xs"
                      >
                        <div className="space-y-0.5">
                          <div className="flex flex-wrap items-center gap-2">
                            <span className="rounded bg-amber-100 px-2 py-0.5 text-[10px] font-bold text-amber-900">
                              {loh.lohType === 'cn_neutral_loh' ? 'CN-Neutral LOH' : 'Deletion LOH'}
                            </span>
                            <span className="font-mono font-bold text-slate-900">
                              {loh.chrom}:{(loh.start + 1).toLocaleString()}–{loh.end.toLocaleString()}
                            </span>
                            <span className="font-mono text-slate-600">
                              ({(loh.lengthBp / 1000).toFixed(1)} kb)
                            </span>
                          </div>
                          <div className="text-[11px] text-slate-600">
                            구간 Het: <span className="font-mono">{loh.meanHetPerKb}/kb</span> (게놈 중앙값{' '}
                            <span className="font-mono">{loh.genomeMedianHetPerKb}/kb</span>) · Copy Ratio{' '}
                            <span className="font-mono">{loh.meanCopyRatio}x</span>
                            {loh.overlappingGenes.length > 0 && ` · 유전자: ${loh.overlappingGenes.join(', ')}`}
                          </div>
                        </div>
                        <button
                          type="button"
                          onClick={() => {
                            onClose();
                            onJumpToRegion(loh.chrom, loh.start, loh.end);
                          }}
                          className="inline-flex items-center gap-1 rounded border border-teal-200 bg-white px-2.5 py-1 text-[11px] font-semibold text-teal-800 hover:bg-teal-50"
                        >
                          <Compass className="h-3 w-3" /> 구간 이동
                        </button>
                      </div>
                    ))
                  )}
                </div>
              )}

              {activeSubTab === 'sv' && (
                <div className="space-y-2">
                  {activeProfile.svEvents.length === 0 ? (
                    <div className="rounded-xl border border-slate-200 bg-slate-50 p-6 text-center text-xs text-slate-600">
                      지지 페어 4쌍 이상의 염색체 간 전좌·대결실·역위 클러스터가 발견되지 않았습니다.
                      (더 민감하게 전체 BAM을 탐색하려면 위 모드에서 &apos;정밀 전체 염색체 전좌(SV) 스캔&apos;을 선택하세요.)
                    </div>
                  ) : (
                    activeProfile.svEvents.map((sv) => (
                      <div
                        key={sv.id}
                        className="flex flex-wrap items-center justify-between gap-2 rounded-xl border border-slate-200 bg-white px-4 py-2.5 text-xs"
                      >
                        <div className="space-y-1">
                          <div className="flex flex-wrap items-center gap-2">
                            <span
                              className={`rounded px-2 py-0.5 text-[10px] font-bold ${
                                sv.svType === 'TRANSLOCATION'
                                  ? 'bg-purple-100 text-purple-800'
                                  : sv.svType === 'LARGE_DELETION'
                                    ? 'bg-rose-100 text-rose-800'
                                    : 'bg-sky-100 text-sky-800'
                              }`}
                            >
                              {sv.svType === 'TRANSLOCATION'
                                ? '염색체 간 전좌'
                                : sv.svType === 'LARGE_DELETION'
                                  ? '거대 결실/점프'
                                  : '역위/중복'}
                            </span>
                            {sv.specificity === 'strain_specific' && (
                              <span className="rounded bg-amber-100 px-1.5 py-0.5 text-[10px] font-bold text-amber-900">
                                균주 특이 SV
                              </span>
                            )}
                            <span className="font-mono font-bold text-slate-900">
                              {sv.chrom1}:{(sv.pos1 + 1).toLocaleString()} ↔ {sv.chrom2}:
                              {(sv.pos2 + 1).toLocaleString()}
                            </span>
                          </div>
                          <div className="text-[11px] text-slate-600">
                            지지 증거: Paired <span className="font-mono font-semibold">{sv.supportPairs}</span> ·
                            Split-read(SA) <span className="font-mono font-semibold">{sv.splitReads}</span>
                            {sv.spanBp !== null && ` · 거리 ${(sv.spanBp / 1000).toFixed(1)} kb`}
                            {(sv.genes1.length > 0 || sv.genes2.length > 0) &&
                              ` · 유전자 [${sv.genes1.join(',') || 'Intergenic'}] ↔ [${sv.genes2.join(',') || 'Intergenic'}]`}
                          </div>
                        </div>

                        <div className="flex items-center gap-1.5">
                          <button
                            type="button"
                            onClick={() => {
                              onClose();
                              onJumpToRegion(sv.chrom1, Math.max(0, sv.pos1 - 400), sv.pos1 + 400);
                            }}
                            className="inline-flex items-center gap-1 rounded border border-teal-200 bg-teal-50 px-2 py-1 text-[11px] font-semibold text-teal-800 hover:bg-teal-100"
                          >
                            <Compass className="h-3 w-3" /> BP1 이동
                          </button>
                          <button
                            type="button"
                            onClick={() => {
                              onClose();
                              onJumpToRegion(sv.chrom2, Math.max(0, sv.pos2 - 400), sv.pos2 + 400);
                            }}
                            className="inline-flex items-center gap-1 rounded border border-slate-300 bg-white px-2 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                          >
                            <Compass className="h-3 w-3" /> BP2 이동
                          </button>
                        </div>
                      </div>
                    ))
                  )}
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
