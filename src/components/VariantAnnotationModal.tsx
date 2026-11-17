import React, { useMemo, useState } from 'react';
import {
  AlertTriangle,
  Compass,
  Download,
  FileSpreadsheet,
  Filter,
  Play,
  Radar,
  Search,
  Square,
  X,
} from 'lucide-react';
import type { ContigInfo, VariantImpact, VariantRecord } from '../bio/types';
import {
  exportVariantsToTsv,
  exportVariantsToVcf,
  IMPACT_ORDER,
  summarizeVariantsByGene,
} from '../bio/variantAnnotator';
import type { ScanProgressEvent, VariantScanScope } from '../workers/wgsProtocol';

interface VariantAnnotationModalProps {
  isOpen: boolean;
  isScanning: boolean;
  scanProgress: ScanProgressEvent | null;
  variants: VariantRecord[];
  scanScope: VariantScanScope;
  currentChrom: string;
  refContigs: ContigInfo[];
  s1Label: string;
  s2Label: string | null;
  hasStrains: boolean;
  onClose: () => void;
  onStartScan: (scope: VariantScanScope) => void;
  onCancelScan: () => void;
  onJumpToVariant: (chrom: string, start: number, end: number) => void;
}

const IMPACT_BADGE: Record<VariantImpact, { label: string; cls: string }> = {
  HIGH: {
    label: 'HIGH',
    cls: 'border-rose-300 bg-rose-100 text-rose-800',
  },
  MODERATE: {
    label: 'MODERATE',
    cls: 'border-amber-300 bg-amber-100 text-amber-800',
  },
  LOW: {
    label: 'LOW',
    cls: 'border-sky-200 bg-sky-50 text-sky-800',
  },
  MODIFIER: {
    label: 'MODIFIER',
    cls: 'border-slate-200 bg-slate-100 text-slate-600',
  },
};

const CONSEQUENCE_LABEL: Record<string, string> = {
  nonsense: 'Nonsense (조기 종결)',
  frameshift: 'Frameshift (틀이동)',
  start_lost: 'Start Lost (개시코돈 소실)',
  stop_lost: 'Stop Lost (종결코돈 소실)',
  splice_site: 'Splice Site (±2bp)',
  gene_loss: '유전자 결실/미커버',
  missense: 'Missense (아미노산 치환)',
  inframe_indel: 'In-frame Indel',
  synonymous: 'Synonymous (침묵 변이)',
  ncrna_or_trna: '비코딩 RNA/tRNA',
  promoter_500bp: '프로모터 (-500bp)',
  intron: 'Intron (인론)',
  low_depth_gap: '국소 저심도',
  intergenic: 'Intergenic (유전자간)',
};

function triggerDownload(filename: string, content: string, mime: string) {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export const VariantAnnotationModal: React.FC<VariantAnnotationModalProps> = ({
  isOpen,
  isScanning,
  scanProgress,
  variants,
  scanScope,
  currentChrom,
  refContigs,
  s1Label,
  s2Label,
  hasStrains,
  onClose,
  onStartScan,
  onCancelScan,
  onJumpToVariant,
}) => {
  const [activeTab, setActiveTab] = useState<'genes' | 'variants'>('genes');
  const [impactFilter, setImpactFilter] = useState<'all' | VariantImpact>('all');
  const [categoryFilter, setCategoryFilter] = useState<string>('all');
  const [queryText, setQueryText] = useState('');

  const countsByImpact = useMemo(() => {
    const counts: Record<VariantImpact, number> = { HIGH: 0, MODERATE: 0, LOW: 0, MODIFIER: 0 };
    for (const v of variants) {
      counts[v.annotation?.impact ?? 'MODIFIER']++;
    }
    return counts;
  }, [variants]);

  const geneSummaries = useMemo(() => summarizeVariantsByGene(variants), [variants]);

  const filteredGenes = useMemo(() => {
    const q = queryText.trim().toLowerCase();
    return geneSummaries.filter((g) => {
      if (impactFilter !== 'all' && g.highestImpact !== impactFilter) return false;
      if (categoryFilter !== 'all' && !g.categories.includes(categoryFilter)) return false;
      if (
        q &&
        !g.geneName.toLowerCase().includes(q) &&
        !(g.systematicName?.toLowerCase().includes(q) ?? false) &&
        !g.chrom.toLowerCase().includes(q) &&
        !g.topExamples.some((ex) => ex.toLowerCase().includes(q))
      ) {
        return false;
      }
      return true;
    });
  }, [geneSummaries, impactFilter, categoryFilter, queryText]);

  const filteredVariants = useMemo(() => {
    const q = queryText.trim().toLowerCase();
    const list = variants.filter((v) => {
      const imp = v.annotation?.impact ?? 'MODIFIER';
      if (impactFilter !== 'all' && imp !== impactFilter) return false;
      if (categoryFilter !== 'all' && v.category !== categoryFilter) return false;
      if (q) {
        const hay = `${v.chrom} ${v.annotation?.geneName ?? ''} ${v.annotation?.systematicName ?? ''} ${v.annotation?.hgvsP ?? ''} ${v.annotation?.hgvsC ?? ''}`.toLowerCase();
        if (!hay.includes(q)) return false;
      }
      return true;
    });
    return [...list].sort(
      (a, b) =>
        IMPACT_ORDER[b.annotation?.impact ?? 'MODIFIER'] -
          IMPACT_ORDER[a.annotation?.impact ?? 'MODIFIER'] ||
        a.chrom.localeCompare(b.chrom) ||
        a.start - b.start,
    );
  }, [variants, impactFilter, categoryFilter, queryText]);

  if (!isOpen) return null;

  const progressPct = scanProgress
    ? Math.min(100, Math.round((scanProgress.scannedBp / Math.max(1, scanProgress.totalBp)) * 100))
    : 0;

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/50 p-4 backdrop-blur-xs"
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div className="flex max-h-[93vh] w-full max-w-6xl flex-col overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-2xl">
        {/* Header */}
        <div className="flex flex-wrap items-center justify-between gap-3 border-b border-slate-200 bg-slate-50 px-6 py-3.5">
          <div className="flex items-center gap-3">
            <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-xl bg-rose-100 text-rose-700">
              <AlertTriangle className="h-5 w-5" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-900">
                전유전체 변이 기능 예측 (SnpEff 스타일 코돈·단백질 주석 &amp; VCF 내보내기)
              </h2>
              <p className="text-xs text-slate-600">
                S288C 엑손·CDS 좌표와 유전부호(핵/ChrM)를 결합해{' '}
                <strong className="font-semibold text-slate-800">Nonsense · Frameshift · Missense</strong>{' '}
                아미노산 변화를 자동 판정합니다.
              </p>
            </div>
          </div>

          <div className="flex items-center gap-2">
            {isScanning ? (
              <button
                type="button"
                onClick={onCancelScan}
                className="inline-flex items-center gap-1.5 rounded-xl bg-rose-600 px-3.5 py-1.5 text-xs font-semibold text-white hover:bg-rose-700"
              >
                <Square className="h-3.5 w-3.5" /> 스캔 중지
              </button>
            ) : (
              <>
                {currentChrom !== 'all' && (
                  <button
                    type="button"
                    disabled={!hasStrains}
                    onClick={() => onStartScan('current_chrom')}
                    className="inline-flex items-center gap-1.5 rounded-xl border border-sky-300 bg-sky-50 px-3 py-1.5 text-xs font-semibold text-sky-800 hover:bg-sky-100 disabled:cursor-not-allowed disabled:opacity-50"
                  >
                    <Radar className="h-3.5 w-3.5" /> 현재 염색체 ({currentChrom}) 스캔
                  </button>
                )}
                <button
                  type="button"
                  disabled={!hasStrains}
                  onClick={() => onStartScan('all_chroms')}
                  className="inline-flex items-center gap-1.5 rounded-xl bg-sky-600 px-4 py-1.5 text-xs font-semibold text-white shadow-xs hover:bg-sky-700 disabled:cursor-not-allowed disabled:bg-slate-300 disabled:text-slate-500"
                >
                  <Play className="h-3.5 w-3.5" /> 전체 17개 염색체 주석 스캔
                </button>
              </>
            )}
            <button
              type="button"
              onClick={onClose}
              className="rounded-lg p-1.5 text-slate-500 hover:bg-slate-200 hover:text-slate-900"
              title="닫기 (Esc)"
            >
              <X className="h-5 w-5" />
            </button>
          </div>
        </div>

        {/* Progress Bar */}
        {isScanning && scanProgress && (
          <div className="border-b border-sky-200 bg-sky-50 px-6 py-2">
            <div className="mb-1 flex items-center justify-between text-xs text-sky-900">
              <span className="font-semibold">
                {scanProgress.totalChroms && scanProgress.totalChroms > 1
                  ? `${scanProgress.chrom} 스캔 중 (${scanProgress.chromIndex}/${scanProgress.totalChroms} 염색체)`
                  : `${scanProgress.chrom} 스캔 중`}
              </span>
              <span className="font-mono">
                {progressPct}% · {scanProgress.scannedBp.toLocaleString()} /{' '}
                {scanProgress.totalBp.toLocaleString()} bp · 변이 {scanProgress.variantsSoFar.toLocaleString()}건
              </span>
            </div>
            <div className="h-1.5 overflow-hidden rounded-full bg-sky-200">
              <div
                className="h-full bg-sky-600 transition-all"
                style={{ width: `${Math.max(3, progressPct)}%` }}
              />
            </div>
          </div>
        )}

        {/* Body */}
        <div className="flex flex-1 flex-col space-y-3 overflow-hidden p-5">
          {/* Summary KPI Row */}
          <div className="grid grid-cols-2 gap-2.5 sm:grid-cols-5">
            <button
              type="button"
              onClick={() => setImpactFilter('all')}
              className={`rounded-xl border px-3 py-2 text-left transition ${
                impactFilter === 'all'
                  ? 'border-sky-400 bg-sky-50/70'
                  : 'border-slate-200 bg-slate-50 hover:bg-slate-100'
              }`}
            >
              <div className="text-[11px] font-medium text-slate-500">
                전체 변이 ({scanScope === 'all_chroms' ? '17chr' : currentChrom})
              </div>
              <div className="font-mono text-lg font-bold text-slate-900">
                {variants.length.toLocaleString()}
              </div>
            </button>
            <button
              type="button"
              onClick={() => setImpactFilter((p) => (p === 'HIGH' ? 'all' : 'HIGH'))}
              className={`rounded-xl border px-3 py-2 text-left transition ${
                impactFilter === 'HIGH'
                  ? 'border-rose-400 bg-rose-50'
                  : 'border-rose-200 bg-rose-50/40 hover:bg-rose-50'
              }`}
            >
              <div className="text-[11px] font-semibold text-rose-700">HIGH (종결·틀이동·결실)</div>
              <div className="font-mono text-lg font-bold text-rose-800">
                {countsByImpact.HIGH.toLocaleString()}
              </div>
            </button>
            <button
              type="button"
              onClick={() => setImpactFilter((p) => (p === 'MODERATE' ? 'all' : 'MODERATE'))}
              className={`rounded-xl border px-3 py-2 text-left transition ${
                impactFilter === 'MODERATE'
                  ? 'border-amber-400 bg-amber-50'
                  : 'border-amber-200 bg-amber-50/40 hover:bg-amber-50'
              }`}
            >
              <div className="text-[11px] font-semibold text-amber-700">MODERATE (미스센스 등)</div>
              <div className="font-mono text-lg font-bold text-amber-800">
                {countsByImpact.MODERATE.toLocaleString()}
              </div>
            </button>
            <button
              type="button"
              onClick={() => setImpactFilter((p) => (p === 'LOW' ? 'all' : 'LOW'))}
              className={`rounded-xl border px-3 py-2 text-left transition ${
                impactFilter === 'LOW'
                  ? 'border-sky-400 bg-sky-50'
                  : 'border-slate-200 bg-slate-50 hover:bg-slate-100'
              }`}
            >
              <div className="text-[11px] font-medium text-slate-600">LOW (침묵 변이)</div>
              <div className="font-mono text-lg font-bold text-slate-800">
                {countsByImpact.LOW.toLocaleString()}
              </div>
            </button>
            <div className="rounded-xl border border-slate-200 bg-slate-50 px-3 py-2">
              <div className="text-[11px] font-medium text-slate-500">변이 포함 유전자</div>
              <div className="font-mono text-lg font-bold text-slate-900">
                {geneSummaries.length.toLocaleString()}개
              </div>
            </div>
          </div>

          {/* Filter & Export Bar */}
          <div className="flex flex-wrap items-center justify-between gap-2 rounded-xl border border-slate-200 bg-slate-50 px-3.5 py-2 text-xs">
            <div className="flex flex-wrap items-center gap-2">
              <div className="inline-flex rounded-lg border border-slate-300 bg-white p-0.5">
                <button
                  type="button"
                  onClick={() => setActiveTab('genes')}
                  className={`rounded-md px-2.5 py-1 text-xs font-semibold ${
                    activeTab === 'genes' ? 'bg-sky-600 text-white' : 'text-slate-600 hover:text-slate-900'
                  }`}
                >
                  유전자별 요약 ({filteredGenes.length})
                </button>
                <button
                  type="button"
                  onClick={() => setActiveTab('variants')}
                  className={`rounded-md px-2.5 py-1 text-xs font-semibold ${
                    activeTab === 'variants' ? 'bg-sky-600 text-white' : 'text-slate-600 hover:text-slate-900'
                  }`}
                >
                  전체 변이 목록 ({filteredVariants.length})
                </button>
              </div>

              <div className="relative">
                <Search className="absolute left-2.5 top-2 h-3.5 w-3.5 text-slate-400" />
                <input
                  type="text"
                  value={queryText}
                  onChange={(e) => setQueryText(e.target.value)}
                  placeholder="유전자·HGVS(p.T511I)·염색체 검색"
                  className="w-56 rounded-lg border border-slate-300 bg-white py-1.5 pl-8 pr-2.5 text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                />
              </div>

              {s2Label && (
                <div className="flex items-center gap-1">
                  <Filter className="h-3.5 w-3.5 text-slate-400" />
                  <select
                    value={categoryFilter}
                    onChange={(e) => setCategoryFilter(e.target.value)}
                    className="rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-800"
                  >
                    <option value="all">모든 균주 구분</option>
                    <option value="strain1_only">{s1Label} 특이</option>
                    <option value="strain2_only">{s2Label} 특이</option>
                    <option value="shared">공통 변이</option>
                    <option value="divergent">분기 변이</option>
                  </select>
                </div>
              )}
            </div>

            <div className="flex items-center gap-2">
              <button
                type="button"
                disabled={variants.length === 0}
                onClick={() =>
                  triggerDownload(
                    `yeast_variants_${scanScope}.vcf`,
                    exportVariantsToVcf(variants, refContigs, s1Label, s2Label),
                    'text/plain;charset=utf-8',
                  )
                }
                className="inline-flex items-center gap-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-[11px] font-semibold text-slate-700 hover:bg-slate-100 disabled:cursor-not-allowed disabled:opacity-50"
              >
                <Download className="h-3.5 w-3.5 text-indigo-600" /> .vcf 저장 (VCF 4.2)
              </button>
              <button
                type="button"
                disabled={variants.length === 0}
                onClick={() =>
                  triggerDownload(
                    `yeast_variants_${scanScope}.tsv`,
                    exportVariantsToTsv(variants, s1Label, s2Label),
                    'text/tab-separated-values;charset=utf-8',
                  )
                }
                className="inline-flex items-center gap-1.5 rounded-lg border border-slate-300 bg-white px-3 py-1.5 text-xs font-semibold text-slate-800 hover:bg-slate-100 disabled:cursor-not-allowed disabled:opacity-50"
              >
                <FileSpreadsheet className="h-3.5 w-3.5 text-emerald-600" /> .tsv 저장 (Excel)
              </button>
            </div>
          </div>

          {/* Table Area */}
          <div className="flex-1 overflow-y-auto rounded-xl border border-slate-200 bg-white">
            {variants.length === 0 ? (
              <div className="flex h-56 flex-col items-center justify-center gap-2 p-6 text-center text-sm text-slate-600">
                <p>
                  {hasStrains
                    ? '상단의 [전체 17개 염색체 주석 스캔] 또는 [현재 염색체 스캔]을 실행하면 유전자별 단백질 영향도 분석 결과가 여기에 표시됩니다.'
                    : '먼저 BAM + BAI 파일을 불러온 뒤 변이 주석 스캔을 실행해주세요.'}
                </p>
              </div>
            ) : activeTab === 'genes' ? (
              <table className="w-full border-collapse text-left text-[13px]">
                <thead className="sticky top-0 z-10 border-b border-slate-200 bg-slate-100 text-xs font-bold text-slate-800">
                  <tr>
                    <th className="px-3 py-2.5">최고 영향</th>
                    <th className="px-3 py-2.5">유전자 (Systematic)</th>
                    <th className="px-3 py-2.5">위치</th>
                    <th className="px-3 py-2.5 text-center">HIGH</th>
                    <th className="px-3 py-2.5 text-center">MOD</th>
                    <th className="px-3 py-2.5 text-center">LOW</th>
                    <th className="px-3 py-2.5">주요 아미노산·코딩 변화</th>
                    <th className="px-3 py-2.5 text-right">이동</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {filteredGenes.slice(0, 500).map((g) => {
                    const badge = IMPACT_BADGE[g.highestImpact];
                    return (
                      <tr key={g.geneKey} className="hover:bg-sky-50/50">
                        <td className="px-3 py-2.5">
                          <span
                            className={`inline-block rounded-md border px-2 py-0.5 text-xs font-bold ${badge.cls}`}
                          >
                            {badge.label}
                          </span>
                        </td>
                        <td className="px-3 py-2.5 font-semibold text-slate-900">
                          {g.geneName}
                          {g.systematicName && g.systematicName !== g.geneName && (
                            <span className="ml-1.5 font-mono text-xs font-normal text-slate-600">
                              ({g.systematicName})
                            </span>
                          )}
                          <span className="ml-1 text-xs text-slate-500">({g.strand})</span>
                        </td>
                        <td className="px-3 py-2.5 font-mono text-[13px] text-slate-700">
                          {g.chrom}:{(g.start + 1).toLocaleString()}
                        </td>
                        <td className="px-3 py-2.5 text-center font-mono font-bold text-rose-700">
                          {g.highCount || '-'}
                        </td>
                        <td className="px-3 py-2.5 text-center font-mono font-semibold text-amber-700">
                          {g.moderateCount || '-'}
                        </td>
                        <td className="px-3 py-2.5 text-center font-mono text-slate-700">
                          {g.lowCount || '-'}
                        </td>
                        <td className="max-w-xs truncate px-3 py-2.5 font-mono text-[13px] text-slate-900">
                          {g.topExamples.join(', ') || '-'}
                        </td>
                        <td className="px-3 py-2.5 text-right">
                          <button
                            type="button"
                            onClick={() => onJumpToVariant(g.chrom, g.start, g.end)}
                            className="inline-flex items-center gap-1 rounded-md border border-sky-200 bg-sky-50 px-2.5 py-1 text-xs font-semibold text-sky-800 hover:bg-sky-100"
                          >
                            <Compass className="h-3.5 w-3.5" /> 보기
                          </button>
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            ) : (
              <table className="w-full border-collapse text-left text-[13px]">
                <thead className="sticky top-0 z-10 border-b border-slate-200 bg-slate-100 text-xs font-bold text-slate-800">
                  <tr>
                    <th className="px-3 py-2.5">영향도</th>
                    <th className="px-3 py-2.5">좌표</th>
                    <th className="px-3 py-2.5">유전자</th>
                    <th className="px-3 py-2.5">분류 (Consequence)</th>
                    <th className="px-3 py-2.5">단백질 (HGVS.p)</th>
                    <th className="px-3 py-2.5">염기 (Ref→Alt)</th>
                    <th className="px-3 py-2.5 text-right">이동</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {filteredVariants.slice(0, 600).map((v, idx) => {
                    const a = v.annotation;
                    const badge = IMPACT_BADGE[a?.impact ?? 'MODIFIER'];
                    return (
                      <tr key={`${v.chrom}_${v.start}_${idx}`} className="hover:bg-sky-50/50">
                        <td className="px-3 py-2">
                          <span
                            className={`inline-block rounded-md border px-2 py-0.5 text-xs font-bold ${badge.cls}`}
                          >
                            {badge.label}
                          </span>
                        </td>
                        <td className="px-3 py-2 font-mono text-[13px] text-slate-800">
                          {v.chrom}:{(v.start + 1).toLocaleString()}
                        </td>
                        <td className="px-3 py-2 font-semibold text-slate-900">
                          {a?.geneName ?? 'Intergenic'}
                        </td>
                        <td className="px-3 py-2 text-xs font-medium text-slate-800">
                          {CONSEQUENCE_LABEL[a?.consequence ?? 'intergenic'] ?? a?.consequence}
                        </td>
                        <td className="px-3 py-2 font-mono text-[13px] font-bold text-indigo-800">
                          {a?.hgvsP ?? a?.hgvsC ?? '-'}
                        </td>
                        <td className="max-w-xs truncate px-3 py-2 font-mono text-[13px] text-slate-700">
                          {v.refAllele.slice(0, 10)} → {(v.category === 'strain2_only' ? v.s2Allele : v.s1Allele).slice(0, 10)}
                        </td>
                        <td className="px-3 py-2 text-right">
                          <button
                            type="button"
                            onClick={() => onJumpToVariant(v.chrom, v.start, v.end)}
                            className="inline-flex items-center gap-1 rounded-md border border-sky-200 bg-sky-50 px-2.5 py-1 text-xs font-semibold text-sky-800 hover:bg-sky-100"
                          >
                            <Compass className="h-3.5 w-3.5" /> 보기
                          </button>
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

