import React, { useMemo, useState } from 'react';
import {
  AlertCircle,
  Check,
  CheckCircle2,
  Compass,
  Copy,
  Dna,
  Download,
  Play,
  QrCode,
  ShieldAlert,
  Square,
  X,
} from 'lucide-react';
import type { BarcodeCandidate } from '../bio/barcodeScanner';
import type {
  BarcodeProgressEvent,
  BarcodeScanPayload,
  BarcodeScanResult,
  BarcodeScanScope,
} from '../workers/wgsProtocol';

interface BarcodeScannerModalProps {
  isOpen: boolean;
  isScanning: boolean;
  progress: BarcodeProgressEvent | null;
  errorMsg: string | null;
  result: BarcodeScanResult | null;
  strainLabels: string[];
  currentChrom: string;
  onClose: () => void;
  onStartScan: (payload: BarcodeScanPayload) => void;
  onCancelScan: () => void;
  onJumpToBreakpoint: (chrom: string, anchorPos: number) => void;
  onSendToBlast: (fastaText: string) => void;
}

function downloadTsv(filename: string, candidates: BarcodeCandidate[]) {
  const header =
    'Strain\tSpecificity\tCategory\tSignalType\tChrom\tAnchor_1based\tLength_bp\tSupport_Reads\tLocal_Depth\tSupport_Pct\tAgreement_Pct\tMotifs\tExtracted_20mer\tOverlapping_Genes\tAssembled_Seq';
  const rows = candidates.map((c) =>
    [
      c.strainLabel,
      c.specificity,
      c.classification,
      c.signalType,
      c.chrom,
      c.anchorPos + 1,
      c.seqLength,
      c.supportReads,
      c.localDepth,
      c.supportRatioPct,
      c.consensusAgreementPct,
      c.motifHits.map((m) => `${m.name}(${m.strand})`).join(';') || '-',
      c.extractedBarseq20mer ?? '-',
      c.overlappingGenes.join(';') || 'Intergenic',
      c.assembledSeq,
    ].join('\t'),
  );
  const blob = new Blob([[header, ...rows].join('\n')], {
    type: 'text/tab-separated-values;charset=utf-8',
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export const BarcodeScannerModal: React.FC<BarcodeScannerModalProps> = ({
  isOpen,
  isScanning,
  progress,
  errorMsg,
  result,
  strainLabels,
  currentChrom,
  onClose,
  onStartScan,
  onCancelScan,
  onJumpToBreakpoint,
  onSendToBlast,
}) => {
  const [scope, setScope] = useState<BarcodeScanScope>('hotspots');
  const [targetStrainIdx, setTargetStrainIdx] = useState<number | null>(null);
  const [minSupport, setMinSupport] = useState(5);
  const [minLen, setMinLen] = useState(10);
  const [customMotif, setCustomMotif] = useState('');
  const [excludeRepeats, setExcludeRepeats] = useState(true);
  const [uniqueMapOnly, setUniqueMapOnly] = useState(true);
  const [excludeSelfGenome, setExcludeSelfGenome] = useState(true);
  const [specificOnly, setSpecificOnly] = useState(false);
  const [copiedId, setCopiedId] = useState<string | null>(null);

  const visibleCandidates = useMemo(() => {
    if (!result) return [];
    if (!specificOnly || strainLabels.length <= 1) return result.candidates;
    return result.candidates.filter((c) => c.specificity !== 'shared');
  }, [result, specificOnly, strainLabels.length]);

  if (!isOpen) return null;

  const copySeq = async (id: string, text: string) => {
    try {
      await navigator.clipboard.writeText(text);
    } catch {
      const ta = document.createElement('textarea');
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      document.body.removeChild(ta);
    }
    setCopiedId(id);
    setTimeout(() => setCopiedId(null), 1600);
  };

  const triggerScan = () => {
    onStartScan({
      scope,
      currentChrom,
      targetStrainIndex: targetStrainIdx,
      options: {
        minSupportReads: minSupport,
        minClipLen: minLen,
        minInsertLen: Math.max(8, minLen - 2),
        minMapQ: uniqueMapOnly ? 10 : 0,
        excludeRepeatLoci: excludeRepeats,
        excludeEndogenousRefMatches: excludeSelfGenome,
        customMotifSeq: customMotif.trim() || undefined,
      },
    });
  };

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/45 p-4 backdrop-blur-xs"
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div className="flex max-h-[93vh] w-full max-w-5xl flex-col overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-2xl">
        {/* Header */}
        <div className="flex items-center justify-between gap-3 border-b border-slate-200 bg-slate-50 px-6 py-3.5">
          <div className="flex items-center gap-3">
            <div className="flex h-9 w-9 shrink-0 items-center justify-center rounded-xl bg-amber-100 text-amber-700">
              <QrCode className="h-5 w-5" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-900">
                외부 삽입 바코드 · 유전공학 스카(Scar) 탐지 스캐너
              </h2>
              <p className="text-xs text-slate-600">
                표준 바코드(<span className="font-mono">UPTAG/DNTAG · loxP · FRT · MX4/TEF</span>)와{' '}
                <strong className="font-semibold text-slate-800">미지 독자 워터마크 서열</strong>을
                Soft-clip/Insertion 합의 조립으로 자동 탐지합니다 (어댑터·VNTR 반복·rDNA/내인성 노이즈 자동 제외).
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
          {/* Controls */}
          <div className="space-y-3 rounded-xl border border-slate-200 bg-slate-50 p-4">
            <div className="grid gap-3 md:grid-cols-4">
              <label className="text-xs">
                <span className="font-semibold text-slate-700">스캔 범위</span>
                <select
                  value={scope}
                  onChange={(e) => setScope(e.target.value as BarcodeScanScope)}
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value="hotspots">주요 통합 핫스팟 빠른 스캔 (HO·URA3·HIS3 등)</option>
                  <option value="current_chrom">
                    현재 염색체 전체 ({currentChrom === 'all' ? '전체' : currentChrom})
                  </option>
                  <option value="all_chroms">전체 17개 염색체 정밀 스캔</option>
                </select>
              </label>

              <label className="text-xs">
                <span className="font-semibold text-slate-700">분석 균주</span>
                <select
                  value={targetStrainIdx === null ? '' : String(targetStrainIdx)}
                  onChange={(e) =>
                    setTargetStrainIdx(e.target.value === '' ? null : Number(e.target.value))
                  }
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value="">
                    {strainLabels.length > 1
                      ? `전체 균주 (${strainLabels.length}개) 교차 비교 ★`
                      : strainLabels[0] ?? '로드된 균주 없음'}
                  </option>
                  {strainLabels.map((lbl, idx) => (
                    <option key={idx} value={String(idx)}>
                      {lbl} 단독 스캔
                    </option>
                  ))}
                </select>
              </label>

              <label className="text-xs">
                <span className="font-semibold text-slate-700">최소 지지 리드</span>
                <select
                  value={minSupport}
                  onChange={(e) => setMinSupport(Number(e.target.value))}
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value={3}>3 reads 이상 (민감)</option>
                  <option value={5}>5 reads 이상 (기본)</option>
                  <option value={10}>10 reads 이상 (고신뢰)</option>
                  <option value={15}>15 reads 이상 (매우 엄격)</option>
                </select>
              </label>

              <label className="text-xs">
                <span className="font-semibold text-slate-700">최소 클립/삽입 길이</span>
                <select
                  value={minLen}
                  onChange={(e) => setMinLen(Number(e.target.value))}
                  className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-2.5 py-1.5 text-xs text-slate-900"
                >
                  <option value={8}>8 bp 이상</option>
                  <option value={10}>10 bp 이상 (기본)</option>
                  <option value={15}>15 bp 이상</option>
                  <option value={20}>20 bp 이상</option>
                </select>
              </label>
            </div>

            <div className="flex flex-wrap items-center gap-x-5 gap-y-1.5 border-t border-slate-200/80 pt-2.5 text-[11px] text-slate-700">
              <label className="inline-flex cursor-pointer items-center gap-1.5 font-medium">
                <input
                  type="checkbox"
                  checked={excludeRepeats}
                  onChange={(e) => setExcludeRepeats(e.target.checked)}
                  className="rounded border-slate-300 text-amber-600 focus:ring-amber-500"
                />
                반복 영역(rDNA/Ty/텔로미어) 제외
              </label>
              <label className="inline-flex cursor-pointer items-center gap-1.5 font-medium">
                <input
                  type="checkbox"
                  checked={uniqueMapOnly}
                  onChange={(e) => setUniqueMapOnly(e.target.checked)}
                  className="rounded border-slate-300 text-amber-600 focus:ring-amber-500"
                />
                유일 매핑 리드만 (MAPQ≥10)
              </label>
              <label className="inline-flex cursor-pointer items-center gap-1.5 font-medium">
                <input
                  type="checkbox"
                  checked={excludeSelfGenome}
                  onChange={(e) => setExcludeSelfGenome(e.target.checked)}
                  className="rounded border-slate-300 text-amber-600 focus:ring-amber-500"
                />
                sacCer3 내인성 서열 제외
              </label>
            </div>

            <div className="flex flex-wrap items-center justify-between gap-3">
              <div className="flex min-w-0 flex-1 items-center gap-2 text-xs">
                <span className="shrink-0 text-slate-600">추가 의심 서열(선택):</span>
                <input
                  type="text"
                  value={customMotif}
                  onChange={(e) => setCustomMotif(e.target.value)}
                  placeholder="특정 링커·프라이머 서열이 있으면 입력 (미입력 시 자동 탐지)"
                  className="min-w-0 flex-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1 font-mono text-xs text-slate-900 focus:border-amber-500 focus:outline-none"
                />
              </div>

              <div className="flex items-center gap-2">
                {isScanning ? (
                  <button
                    type="button"
                    onClick={onCancelScan}
                    className="inline-flex items-center gap-1.5 rounded-xl bg-rose-600 px-4 py-2 text-xs font-semibold text-white hover:bg-rose-700"
                  >
                    <Square className="h-3.5 w-3.5" /> 스캔 중지
                  </button>
                ) : (
                  <button
                    type="button"
                    disabled={strainLabels.length === 0}
                    onClick={triggerScan}
                    className="inline-flex items-center gap-1.5 rounded-xl bg-amber-600 px-5 py-2 text-xs font-semibold text-white shadow-sm hover:bg-amber-700 disabled:cursor-not-allowed disabled:bg-slate-300 disabled:text-slate-500"
                  >
                    <Play className="h-3.5 w-3.5" /> 바코드 스캔 시작
                  </button>
                )}
              </div>
            </div>

            {strainLabels.length === 0 && (
              <p className="text-[11px] text-rose-700">
                ※ 스캔을 실행하려면 먼저 화면에 BAM + BAI 파일을 1개 이상 불러와 주세요.
              </p>
            )}

            {isScanning && progress && (
              <div className="space-y-1 pt-1">
                <div className="h-1.5 overflow-hidden rounded-full bg-slate-200">
                  <div
                    className="h-full bg-amber-500 transition-all"
                    style={{ width: `${Math.max(4, progress.percent)}%` }}
                  />
                </div>
                <div className="flex justify-between text-[11px] text-slate-600">
                  <span>{progress.message}</span>
                  <span className="font-mono">
                    {progress.percent}% · 후보 {progress.foundSoFar}건
                  </span>
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
            <div className="space-y-3">
              <div className="flex flex-wrap items-center justify-between gap-2 rounded-xl border border-slate-200 bg-white px-4 py-2.5 text-xs">
                <div className="flex flex-wrap items-center gap-2">
                  <span className="font-bold text-slate-900">
                    탐지 결과: 총 {result.candidates.length}건
                  </span>
                  <span className="font-mono text-[11px] text-slate-500">
                    (구간 {result.scannedLociCount}개 ·{' '}
                    {result.scannedBpTotal.toLocaleString()} bp · {result.elapsedMs} ms)
                  </span>
                  {result.filteredStats &&
                    result.filteredStats.repeatLociFiltered +
                      result.filteredStats.tandemOrArtifactFiltered +
                      result.filteredStats.endogenousRefFiltered >
                      0 && (
                      <span
                        className="rounded-md bg-slate-100 px-2 py-0.5 text-[11px] font-medium text-slate-600"
                        title="위양성 억제를 위해 자동 제외된 반복·내인성·어댑터 신호 건수"
                      >
                        노이즈 자동 제외:{' '}
                        {[
                          result.filteredStats.repeatLociFiltered > 0
                            ? `rDNA/반복 ${result.filteredStats.repeatLociFiltered}건`
                            : null,
                          result.filteredStats.tandemOrArtifactFiltered > 0
                            ? `VNTR/어댑터 ${result.filteredStats.tandemOrArtifactFiltered}건`
                            : null,
                          result.filteredStats.endogenousRefFiltered > 0
                            ? `sacCer3 내인성 ${result.filteredStats.endogenousRefFiltered}건`
                            : null,
                        ]
                          .filter(Boolean)
                          .join(' · ')}
                      </span>
                    )}
                </div>

                <div className="flex flex-wrap items-center gap-2">
                  {strainLabels.length > 1 && (
                    <label className="inline-flex cursor-pointer items-center gap-1.5 text-xs font-medium text-slate-700">
                      <input
                        type="checkbox"
                        checked={specificOnly}
                        onChange={(e) => setSpecificOnly(e.target.checked)}
                        className="rounded border-slate-300 text-amber-600"
                      />
                      특정 균주에만 존재하는 후보만 보기
                    </label>
                  )}
                  {visibleCandidates.length > 0 && (
                    <button
                      type="button"
                      onClick={() =>
                        downloadTsv(`yeast_barcode_scan_${result.scope}.tsv`, visibleCandidates)
                      }
                      className="inline-flex items-center gap-1 rounded-lg border border-slate-300 bg-white px-2.5 py-1 text-[11px] font-semibold text-slate-700 hover:bg-slate-100"
                    >
                      <Download className="h-3.5 w-3.5" /> .tsv 저장
                    </button>
                  )}
                </div>
              </div>

              {visibleCandidates.length === 0 ? (
                <div className="flex items-center justify-center gap-2 rounded-xl border border-emerald-200 bg-emerald-50/70 p-6 text-xs font-medium text-emerald-900">
                  <CheckCircle2 className="h-4 w-4 text-emerald-600" />
                  <span>
                    선택된 조건·범위에서 외부 바코드 또는 인위적 비-레퍼런스 삽입 흔적이 발견되지
                    않았습니다.
                    {scope === 'hotspots'
                      ? ' (전체 무작위 위치를 확인하려면 위 범위를 전체 염색체로 바꿔 스캔하세요.)'
                      : ''}
                  </span>
                </div>
              ) : (
                <div className="space-y-2.5">
                  {visibleCandidates.map((cand) => (
                    <div
                      key={cand.id}
                      className={`space-y-2 rounded-xl border p-3.5 ${
                        cand.classification === 'known_barcode_or_scar'
                          ? 'border-amber-300 bg-amber-50/50'
                          : cand.specificity === 'strain_specific'
                            ? 'border-sky-300 bg-sky-50/40'
                            : 'border-slate-200 bg-white'
                      }`}
                    >
                      <div className="flex flex-wrap items-center justify-between gap-2">
                        <div className="flex flex-wrap items-center gap-1.5 text-xs">
                          <span className="rounded bg-slate-800 px-2 py-0.5 font-semibold text-white">
                            {cand.strainLabel}
                          </span>

                          {cand.specificity === 'strain_specific' && (
                            <span className="inline-flex items-center gap-1 rounded-md bg-rose-100 px-2.5 py-0.5 text-xs font-bold text-rose-800">
                              <ShieldAlert className="h-3.5 w-3.5" /> 이 균주 특이 삽입
                            </span>
                          )}
                          {cand.specificity === 'shared' && (
                            <span
                              className="rounded-md bg-slate-200 px-2.5 py-0.5 text-xs font-semibold text-slate-700"
                              title={`공통 검출: ${cand.sharedWithStrains.join(', ')}`}
                            >
                              대조 균주와 공통
                            </span>
                          )}

                          <span
                            className={`rounded-md px-2.5 py-0.5 text-xs font-bold ${
                              cand.classification === 'known_barcode_or_scar'
                                ? 'bg-amber-200 text-amber-950'
                                : cand.classification === 'novel_insertion_barcode'
                                  ? 'bg-sky-200 text-sky-950'
                                  : 'bg-violet-100 text-violet-900'
                            }`}
                          >
                            {cand.classification === 'known_barcode_or_scar'
                              ? '표준 바코드/스카 매칭'
                              : cand.classification === 'novel_insertion_barcode'
                                ? '신규 삽입 바코드(조립)'
                                : '외래 카세트 접합부(Breakpoint)'}
                          </span>

                          <span className="font-mono text-[13px] font-bold text-slate-900">
                            {cand.chrom}:{(cand.anchorPos + 1).toLocaleString()}
                          </span>
                          <span className="text-xs text-slate-600">
                            ({cand.signalType === 'paired_insertion'
                              ? '양방향 클립 병합'
                              : cand.signalType === 'insertion'
                                ? '리드 내부 Insertion'
                                : cand.signalType === 'left_clip'
                                  ? '좌측 Soft-clip'
                                  : '우측 Soft-clip'}
                            )
                          </span>
                        </div>

                        <div className="flex items-center gap-1.5">
                          <button
                            type="button"
                            onClick={() =>
                              onSendToBlast(
                                `>barcode_${cand.strainLabel}_${cand.chrom}_${cand.anchorPos + 1}\n${cand.assembledSeq}`,
                              )
                            }
                            className="inline-flex items-center gap-1.5 rounded-lg border border-indigo-300 bg-indigo-50 px-3 py-1.5 text-xs font-semibold text-indigo-800 hover:bg-indigo-100"
                          >
                            <Dna className="h-3.5 w-3.5" /> 로컬 BLAST 대조
                          </button>
                          <button
                            type="button"
                            onClick={() => onJumpToBreakpoint(cand.chrom, cand.anchorPos)}
                            className="inline-flex items-center gap-1.5 rounded-lg bg-sky-600 px-3 py-1.5 text-xs font-semibold text-white hover:bg-sky-700"
                          >
                            <Compass className="h-3.5 w-3.5" /> 뷰어 위치 보기
                          </button>
                        </div>
                      </div>

                      <div className="flex flex-wrap items-center gap-3 text-xs text-slate-700">
                        <span>
                          조립 길이: <b className="font-mono text-slate-900">{cand.seqLength} bp</b>
                        </span>
                        <span>
                          지지 리드:{' '}
                          <b className="font-mono text-slate-900">
                            {cand.supportReads}/{cand.localDepth} reads ({cand.supportRatioPct}%)
                          </b>
                        </span>
                        <span>
                          서열 합의도: <b className="font-mono text-emerald-700">{cand.consensusAgreementPct}%</b>
                        </span>
                        <span>
                          인근 유전자:{' '}
                          <b className="text-blue-800">
                            {cand.overlappingGenes.join(', ') || 'Intergenic'}
                          </b>
                        </span>
                      </div>

                      {(cand.motifHits.length > 0 || cand.extractedBarseq20mer) && (
                        <div className="flex flex-wrap items-center gap-1.5">
                          {cand.motifHits.map((mh, idx) => (
                            <span
                              key={`${mh.motifId}_${idx}`}
                              className="rounded-md border border-amber-300 bg-white px-2.5 py-0.5 text-xs font-bold text-amber-900"
                            >
                              시그니처: {mh.name} ({mh.strand} · {mh.matchedBases} bp{' '}
                              {mh.isExact ? '정확 일치' : '부분 일치'})
                            </span>
                          ))}
                          {cand.extractedBarseq20mer && (
                            <span className="rounded-md border border-emerald-400 bg-emerald-100 px-2.5 py-1 font-mono text-[13px] font-bold text-emerald-950">
                              추출된 20-mer 바코드: {cand.extractedBarseq20mer}
                            </span>
                          )}
                        </div>
                      )}

                      <div className="flex items-center justify-between gap-2 rounded-lg border border-slate-200 bg-white px-3 py-2">
                        <code className="min-w-0 flex-1 overflow-x-auto whitespace-nowrap font-mono text-[14px] font-medium leading-6 text-slate-900">
                          {cand.assembledSeq}
                        </code>
                        <button
                          type="button"
                          onClick={() => copySeq(cand.id, cand.assembledSeq)}
                          className="inline-flex shrink-0 items-center gap-1 rounded-md border border-slate-300 bg-slate-50 px-2.5 py-1 text-xs font-semibold text-slate-800 hover:bg-slate-100"
                        >
                          {copiedId === cand.id ? (
                            <Check className="h-3 w-3 text-emerald-600" />
                          ) : (
                            <Copy className="h-3 w-3" />
                          )}
                          복사
                        </button>
                      </div>
                    </div>
                  ))}
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
