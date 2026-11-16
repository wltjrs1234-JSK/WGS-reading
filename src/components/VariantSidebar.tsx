import React, { useMemo, useState } from 'react';
import {
  AlertTriangle,
  CheckCircle2,
  Filter,
  GitCompare,
  Radar,
  Scissors,
  Search,
  Settings2,
  Square,
} from 'lucide-react';
import type { ConsensusOptions, ContigInfo, DiffCategory, GffFeature, VariantRecord } from '../bio/types';
import type { AnnotationOrigin, ScanProgressEvent, TrackSummaryMeta } from '../workers/wgsProtocol';

interface VariantSidebarProps {
  chrom: string;
  windowVariants: VariantRecord[];
  scannedVariants: VariantRecord[];
  isScanning: boolean;
  scanProgress: ScanProgressEvent | null;
  hasGff: boolean;
  annotationOrigin: AnnotationOrigin;
  strains: TrackSummaryMeta[];
  primaryStrainIndex: number;
  secondaryStrainIndex: number | null;
  options: ConsensusOptions;
  contigOverrides: Record<string, Record<string, string>>;
  refContigs: ContigInfo[];
  onStartScan: () => void;
  onCancelScan: () => void;
  onJumpToRange: (chrom: string, start: number, end: number) => void;
  onJumpToChrom: (chrom: string) => void;
  onExtractFeature: (feat: GffFeature) => void;
  onSearchGenes: (query: string) => Promise<GffFeature[]>;
  onChangePrimaryIndex: (idx: number) => void;
  onChangeSecondaryIndex: (idx: number | null) => void;
  onChangeOptions: (next: ConsensusOptions) => void;
  onChangeContigOverride: (strainId: string, refChrom: string, bamChrom: string) => void;
}

const MULTI_CATEGORY_LABELS: Record<DiffCategory, { text: string; badge: string }> = {
  match: { text: '일치', badge: 'bg-slate-200 text-slate-700' },
  strain1_only: { text: '균주1 특이', badge: 'bg-sky-100 text-sky-800 border border-sky-200' },
  strain2_only: { text: '균주2 특이', badge: 'bg-emerald-100 text-emerald-800 border border-emerald-200' },
  shared: { text: '공통 변이', badge: 'bg-orange-100 text-orange-800 border border-orange-200' },
  divergent: { text: '분기 변이', badge: 'bg-rose-100 text-rose-800 border border-rose-200' },
  gap_s1: { text: '균주1 결실', badge: 'bg-slate-200 text-slate-700' },
  gap_s2: { text: '균주2 결실', badge: 'bg-slate-200 text-slate-700' },
  gap_both: { text: '공통 미커버', badge: 'bg-slate-200 text-slate-600' },
};

const SINGLE_CATEGORY_LABELS: Record<DiffCategory, { text: string; badge: string }> = {
  match: { text: '일치', badge: 'bg-slate-200 text-slate-700' },
  strain1_only: { text: '변이', badge: 'bg-sky-100 text-sky-800 border border-sky-200' },
  strain2_only: { text: '변이', badge: 'bg-sky-100 text-sky-800 border border-sky-200' },
  shared: { text: '변이', badge: 'bg-sky-100 text-sky-800 border border-sky-200' },
  divergent: { text: '변이', badge: 'bg-sky-100 text-sky-800 border border-sky-200' },
  gap_s1: { text: '결실/미커버', badge: 'bg-slate-200 text-slate-700' },
  gap_s2: { text: '결실/미커버', badge: 'bg-slate-200 text-slate-700' },
  gap_both: { text: '결실/미커버', badge: 'bg-slate-200 text-slate-600' },
};

export const VariantSidebar: React.FC<VariantSidebarProps> = ({
  chrom,
  windowVariants,
  scannedVariants,
  isScanning,
  scanProgress,
  hasGff,
  annotationOrigin,
  strains,
  primaryStrainIndex,
  secondaryStrainIndex,
  options,
  contigOverrides,
  refContigs,
  onStartScan,
  onCancelScan,
  onJumpToRange,
  onJumpToChrom,
  onExtractFeature,
  onSearchGenes,
  onChangePrimaryIndex,
  onChangeSecondaryIndex,
  onChangeOptions,
  onChangeContigOverride,
}) => {
  const [tab, setTab] = useState<'variants' | 'genes' | 'settings'>('variants');
  const [filterCat, setFilterCat] = useState<string>('all');
  const [geneQuery, setGeneQuery] = useState('');
  const [geneResults, setGeneResults] = useState<GffFeature[]>([]);

  const hasStrain2 =
    secondaryStrainIndex !== null &&
    secondaryStrainIndex >= 0 &&
    secondaryStrainIndex < strains.length;

  const activeVariants = scannedVariants.length > 0 ? scannedVariants : windowVariants;
  const filteredVariants = useMemo(() => {
    if (filterCat === 'all') return activeVariants;
    if (filterCat === 'gaps') {
      return activeVariants.filter(
        (v) => v.category === 'gap_s1' || v.category === 'gap_s2' || v.category === 'gap_both',
      );
    }
    if (!hasStrain2 && (filterCat === 'SNP' || filterCat === 'INS' || filterCat === 'DEL')) {
      return activeVariants.filter((v) => v.type === filterCat);
    }
    return activeVariants.filter((v) => v.category === filterCat);
  }, [activeVariants, filterCat, hasStrain2]);

  const handleGeneSearch = async (val: string) => {
    setGeneQuery(val);
    if (!val.trim()) {
      setGeneResults([]);
      return;
    }
    const hits = await onSearchGenes(val);
    setGeneResults(hits);
  };

  const labelMap = hasStrain2 ? MULTI_CATEGORY_LABELS : SINGLE_CATEGORY_LABELS;

  return (
    <aside className="flex w-84 shrink-0 flex-col border-l border-slate-200 bg-slate-50">
      {/* Tabs */}
      <div className="grid grid-cols-3 border-b border-slate-200 bg-white text-xs font-semibold">
        <button
          type="button"
          onClick={() => setTab('variants')}
          className={`flex items-center justify-center gap-1.5 py-2.5 ${
            tab === 'variants'
              ? 'border-b-2 border-sky-600 text-sky-700'
              : 'text-slate-500 hover:text-slate-900'
          }`}
        >
          <GitCompare className="h-3.5 w-3.5" />
          변이 비교
        </button>
        <button
          type="button"
          onClick={() => setTab('genes')}
          className={`flex items-center justify-center gap-1.5 py-2.5 ${
            tab === 'genes'
              ? 'border-b-2 border-sky-600 text-sky-700'
              : 'text-slate-500 hover:text-slate-900'
          }`}
        >
          <Search className="h-3.5 w-3.5" />
          유전자
        </button>
        <button
          type="button"
          onClick={() => setTab('settings')}
          className={`flex items-center justify-center gap-1.5 py-2.5 ${
            tab === 'settings'
              ? 'border-b-2 border-sky-600 text-sky-700'
              : 'text-slate-500 hover:text-slate-900'
          }`}
        >
          <Settings2 className="h-3.5 w-3.5" />
          매핑·설정
        </button>
      </div>

      {/* Tab 1: Variants */}
      {tab === 'variants' && (
        <div className="flex flex-1 flex-col overflow-hidden">
          <div className="space-y-2.5 border-b border-slate-200 bg-white p-3">
            <div className="flex items-center justify-between gap-2">
              <span className="text-xs font-bold text-slate-800">
                {scannedVariants.length > 0
                  ? `${chrom} 전체 스캔 (${scannedVariants.length}건)`
                  : `현재 화면 변이 (${windowVariants.length}건)`}
              </span>
              {isScanning ? (
                <button
                  type="button"
                  onClick={onCancelScan}
                  className="inline-flex items-center gap-1 rounded-lg bg-rose-100 px-2.5 py-1 text-[11px] font-semibold text-rose-700 hover:bg-rose-200"
                >
                  <Square className="h-3 w-3" /> 중지
                </button>
              ) : (
                <button
                  type="button"
                  onClick={onStartScan}
                  className="inline-flex items-center gap-1 rounded-lg bg-sky-100 px-2.5 py-1 text-[11px] font-semibold text-sky-800 hover:bg-sky-200"
                >
                  <Radar className="h-3.5 w-3.5" /> 염색체 스캔
                </button>
              )}
            </div>

            {isScanning && scanProgress && (
              <div className="space-y-1">
                <div className="h-1.5 overflow-hidden rounded-full bg-slate-200">
                  <div
                    className="h-full bg-sky-600 transition-all"
                    style={{
                      width: `${Math.min(
                        100,
                        Math.round((scanProgress.scannedBp / Math.max(1, scanProgress.totalBp)) * 100),
                      )}%`,
                    }}
                  />
                </div>
                <div className="flex justify-between text-[10px] text-slate-600">
                  <span>
                    {scanProgress.scannedBp.toLocaleString()} / {scanProgress.totalBp.toLocaleString()} bp
                  </span>
                  <span>{scanProgress.variantsSoFar}건 발견</span>
                </div>
              </div>
            )}

            <div className="flex items-center gap-1.5">
              <Filter className="h-3.5 w-3.5 text-slate-500" />
              <select
                value={filterCat}
                onChange={(e) => setFilterCat(e.target.value)}
                className="flex-1 rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-800"
              >
                <option value="all">모든 차이 보기 ({activeVariants.length})</option>
                {hasStrain2 ? (
                  <>
                    <option value="strain1_only">균주 1 특이 변이</option>
                    <option value="strain2_only">균주 2 특이 변이</option>
                    <option value="shared">공통 변이</option>
                    <option value="divergent">서로 다른 다중 변이</option>
                    <option value="gaps">커버리지 결실 (Gap)</option>
                  </>
                ) : (
                  <>
                    <option value="SNP">SNP (단일 염기 치환)</option>
                    <option value="INS">INS (삽입)</option>
                    <option value="DEL">DEL (결실)</option>
                    <option value="gaps">커버리지 결실 (Gap)</option>
                  </>
                )}
              </select>
            </div>
          </div>

          <div className="flex-1 space-y-1.5 overflow-y-auto p-2">
            {filteredVariants.length === 0 ? (
              <div className="py-10 text-center text-xs text-slate-500">
                현재 조건에 해당하는 변이가 없습니다.
              </div>
            ) : (
              filteredVariants.slice(0, 400).map((v, idx) => {
                const badge = labelMap[v.category];
                return (
                  <button
                    key={`${v.start}_${idx}`}
                    type="button"
                    onClick={() => {
                      const pad = Math.max(80, Math.round((v.end - v.start) * 1.5));
                      onJumpToRange(v.chrom, Math.max(0, v.start - pad), v.end + pad);
                    }}
                    className="w-full rounded-lg border border-slate-200 bg-white p-2.5 text-left shadow-2xs transition hover:border-sky-400 hover:bg-sky-50/40"
                  >
                    <div className="flex items-center justify-between gap-2">
                      <span className="font-mono text-xs font-bold text-slate-900">
                        {(v.start + 1).toLocaleString()}
                        {v.end - v.start > 1 ? `–${v.end.toLocaleString()}` : ''}
                      </span>
                      <span className={`rounded px-1.5 py-0.5 text-[10px] font-semibold ${badge.badge}`}>
                        {badge.text} · {v.type}
                      </span>
                    </div>
                    <div className="mt-1 flex items-center gap-2 truncate font-mono text-[11px] text-slate-600">
                      <span>Ref:{v.refAllele.slice(0, 8)}</span>
                      <span>→</span>
                      <span className="font-semibold text-sky-700">
                        {hasStrain2 ? 'A' : '균주'}:{v.s1Allele.slice(0, 8)}
                      </span>
                      {hasStrain2 && (
                        <span className="font-semibold text-emerald-700">
                          B:{v.s2Allele.slice(0, 8)}
                        </span>
                      )}
                    </div>
                  </button>
                );
              })
            )}
          </div>
        </div>
      )}

      {/* Tab 2: Gene Search */}
      {tab === 'genes' && (
        <div className="flex flex-1 flex-col space-y-3 overflow-hidden p-3">
          {!hasGff ? (
            <div className="py-10 text-center text-xs text-slate-500">
              유전자 주석 정보가 없습니다. SGD .gff3 파일을 첨부하거나 S. cerevisiae 레퍼런스 FASTA를 불러오면 유전자 검색이 활성화됩니다.
            </div>
          ) : (
            <>
              <div className="flex items-center justify-between text-[11px] text-slate-600">
                <span>유전자 주석 소스:</span>
                <span className="rounded bg-violet-100 px-2 py-0.5 font-semibold text-violet-800">
                  {annotationOrigin === 'builtin_s288c'
                    ? '내장 S288C R64 (자동)'
                    : annotationOrigin === 'user_gff'
                      ? '사용자 GFF3'
                      : '데모 주석'}
                </span>
              </div>
              <div className="relative">
                <Search className="absolute left-2.5 top-2 h-4 w-4 text-slate-400" />
                <input
                  type="text"
                  value={geneQuery}
                  onChange={(e) => handleGeneSearch(e.target.value)}
                  placeholder="유전자명 검색 (예: URA3, LEU2, ACT1, COX1)"
                  className="w-full rounded-lg border border-slate-300 bg-white py-1.5 pl-8 pr-3 text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                />
              </div>
              <div className="flex-1 space-y-1.5 overflow-y-auto">
                {geneResults.length === 0 && geneQuery.trim() !== '' && (
                  <p className="py-6 text-center text-xs text-slate-500">
                    일치하는 유전자를 찾지 못했습니다.
                  </p>
                )}
                {geneResults.map((f) => (
                  <div
                    key={f.id}
                    className="rounded-lg border border-slate-200 bg-white p-2.5 text-xs shadow-2xs"
                  >
                    <div className="flex items-center justify-between gap-2">
                      <button
                        type="button"
                        onClick={() => onJumpToRange(f.chrom, Math.max(0, f.start - 150), f.end + 150)}
                        className="text-left font-bold text-sky-700 hover:underline"
                      >
                        {f.name}
                        {f.systematicName ? ` (${f.systematicName})` : ''}
                      </button>
                      <button
                        type="button"
                        onClick={() => onExtractFeature(f)}
                        className="inline-flex items-center gap-1 rounded bg-sky-100 px-2 py-0.5 text-[10px] font-semibold text-sky-800 hover:bg-sky-200"
                        title="유전자 전체 서열 바로 추출"
                      >
                        <Scissors className="h-3 w-3" /> 추출
                      </button>
                    </div>
                    <div className="mt-0.5 font-mono text-[10px] text-slate-600">
                      {f.chrom}:{(f.start + 1).toLocaleString()}–{f.end.toLocaleString()} ({f.strand})
                    </div>
                    {f.note && (
                      <p className="mt-1 line-clamp-2 text-[10px] text-slate-500">{f.note}</p>
                    )}
                  </div>
                ))}
              </div>
            </>
          )}
        </div>
      )}

      {/* Tab 3: Settings & All-Chromosome Mapping Matrix */}
      {tab === 'settings' && (
        <div className="flex-1 space-y-4 overflow-y-auto p-3.5 text-xs">
          {strains.length >= 1 && (
            <div className="space-y-2 rounded-xl border border-slate-200 bg-white p-3">
              <div className="font-bold text-slate-800">비교 대상 슬롯 설정</div>
              <div>
                <label className="text-[11px] font-semibold text-sky-700">비교 대상 A (균주 1)</label>
                <select
                  value={primaryStrainIndex}
                  onChange={(e) => onChangePrimaryIndex(Number(e.target.value))}
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                >
                  {strains.map((s, i) => (
                    <option key={s.id} value={i}>
                      {s.label}
                    </option>
                  ))}
                </select>
              </div>
              {strains.length >= 2 && (
                <div>
                  <label className="text-[11px] font-semibold text-emerald-700">
                    비교 대상 B (균주 2)
                  </label>
                  <select
                    value={secondaryStrainIndex === null ? '' : String(secondaryStrainIndex)}
                    onChange={(e) =>
                      onChangeSecondaryIndex(e.target.value === '' ? null : Number(e.target.value))
                    }
                    className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900"
                  >
                    <option value="">(선택 안 함 — 단일 균주만 비교)</option>
                    {strains.map((s, i) => (
                      <option key={s.id} value={i}>
                        {s.label}
                      </option>
                    ))}
                  </select>
                </div>
              )}
            </div>
          )}

          <div className="space-y-2.5 rounded-xl border border-slate-200 bg-white p-3">
            <div className="font-bold text-slate-800">Consensus 판정 파라미터</div>
            <label className="block">
              <span className="text-[11px] text-slate-600">최소 Depth (미만은 N/Gap 처리)</span>
              <input
                type="number"
                min={1}
                max={200}
                value={options.minDepth}
                onChange={(e) =>
                  onChangeOptions({ ...options, minDepth: Math.max(1, Number(e.target.value)) })
                }
                className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2 py-1 text-slate-900"
              />
            </label>
            <label className="block">
              <span className="text-[11px] text-slate-600">최소 MAPQ</span>
              <input
                type="number"
                min={0}
                max={60}
                value={options.minMapQ}
                onChange={(e) =>
                  onChangeOptions({ ...options, minMapQ: Math.max(0, Number(e.target.value)) })
                }
                className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2 py-1 text-slate-900"
              />
            </label>
            <label className="block">
              <span className="text-[11px] text-slate-600">배수성 모드</span>
              <select
                value={options.ploidyMode}
                onChange={(e) =>
                  onChangeOptions({
                    ...options,
                    ploidyMode: e.target.value as 'haploid' | 'diploid',
                  })
                }
                className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2 py-1 text-slate-900"
              >
                <option value="haploid">Haploid (반수체 — 단일 최빈 염기)</option>
                <option value="diploid">Diploid (이배체 — 이형접합 IUPAC 코드)</option>
              </select>
            </label>
          </div>

          {/* All Chromosomes Mapping Matrix */}
          <div className="space-y-2 rounded-xl border border-slate-200 bg-white p-3">
            <div className="font-bold text-slate-800">
              전체 염색체 매핑 현황 ({refContigs.length}개)
            </div>
            {strains.map((s) => (
              <div key={s.id} className="rounded bg-slate-50 p-2 text-[10px] text-slate-600">
                <div className="font-semibold text-slate-800">
                  {s.label} — BAM 헤더 ({s.contigs.length}개):
                </div>
                <div className="mt-0.5 truncate font-mono text-slate-500" title={s.contigs.map((c) => c.name).join(', ')}>
                  {s.contigs.map((c) => c.name).join(', ') || '(없음)'}
                </div>
              </div>
            ))}
            <p className="text-[10px] text-slate-500">
              염색체 이름을 클릭하면 해당 염색체로 바로 이동하며, 미매칭된 경우 BAM 컨티그를 직접 지정할 수 있습니다.
            </p>
            <div className="max-h-80 space-y-2 overflow-y-auto pr-1">
              {refContigs.map((rc) => {
                const isActive = rc.name === chrom;
                return (
                  <div
                    key={rc.name}
                    className={`rounded-lg border p-2 ${
                      isActive ? 'border-sky-400 bg-sky-50/50' : 'border-slate-200 bg-slate-50'
                    }`}
                  >
                    <div className="flex items-center justify-between gap-2">
                      <button
                        type="button"
                        onClick={() => onJumpToChrom(rc.name)}
                        className="font-mono text-xs font-bold text-sky-700 hover:underline"
                      >
                        {rc.name}
                      </button>
                      <span className="font-mono text-[10px] text-slate-500">
                        {rc.length.toLocaleString()} bp
                      </span>
                    </div>
                    <div className="mt-1.5 space-y-1">
                      {strains.map((s) => {
                        const overrideVal = contigOverrides[s.id]?.[rc.name];
                        const resolved = overrideVal || s.mappings[rc.name] || '';
                        const isMatched = Boolean(resolved);
                        return (
                          <div key={s.id} className="flex items-center gap-1.5 text-[10px]">
                            {isMatched ? (
                              <CheckCircle2 className="h-3 w-3 shrink-0 text-emerald-600" />
                            ) : (
                              <AlertTriangle className="h-3 w-3 shrink-0 text-amber-500" />
                            )}
                            <span className="w-16 shrink-0 truncate font-medium text-slate-700">
                              {s.label}
                            </span>
                            <select
                              value={overrideVal ?? s.mappings[rc.name] ?? ''}
                              onChange={(e) => onChangeContigOverride(s.id, rc.name, e.target.value)}
                              className="min-w-0 flex-1 rounded border border-slate-300 bg-white px-1.5 py-0.5 font-mono text-[10px] text-slate-800"
                            >
                              <option value="">(미매칭 — 선택)</option>
                              {s.contigs.map((c) => (
                                <option key={c.name} value={c.name}>
                                  {c.name} ({c.length.toLocaleString()})
                                </option>
                              ))}
                            </select>
                          </div>
                        );
                      })}
                    </div>
                  </div>
                );
              })}
            </div>
          </div>
        </div>
      )}
    </aside>
  );
};
