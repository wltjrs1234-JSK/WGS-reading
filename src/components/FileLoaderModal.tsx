import React, { useRef, useState } from 'react';
import {
  AlertCircle,
  CheckCircle2,
  Dna,
  FileCode2,
  FolderOpen,
  Layers,
  Play,
  Plus,
  Trash2,
  Upload,
  X,
} from 'lucide-react';
import { formatBytes, matchDroppedFiles, type PairedStrainFiles } from '../utils/fileMatcher';
import type { LoadLocalFilesPayload, StrainFileInput } from '../workers/wgsProtocol';

interface FileLoaderModalProps {
  isOpen: boolean;
  mode: 'tracks' | 'genome';
  activeGenomeLabel: string;
  isLoading: boolean;
  loadingMsg?: string | null;
  elapsedSec?: number;
  errorMsg: string | null;
  onClose: () => void;
  onCancelLoad?: () => void;
  onAddTracksOnly: (strains: StrainFileInput[]) => Promise<boolean>;
  onLoadCustomGenome: (payload: LoadLocalFilesPayload) => Promise<boolean>;
  onSelectSacCer3: () => void;
  onLoadDemo: () => void;
}

export const FileLoaderModal: React.FC<FileLoaderModalProps> = ({
  isOpen,
  mode,
  activeGenomeLabel,
  isLoading,
  loadingMsg,
  elapsedSec = 0,
  errorMsg,
  onClose,
  onCancelLoad,
  onAddTracksOnly,
  onLoadCustomGenome,
  onSelectSacCer3,
  onLoadDemo,
}) => {
  const [fastaFile, setFastaFile] = useState<File | undefined>();
  const [faiFile, setFaiFile] = useState<File | undefined>();
  const [gffFile, setGffFile] = useState<File | undefined>();
  const [strains, setStrains] = useState<PairedStrainFiles[]>([]);
  const [unpairedBams, setUnpairedBams] = useState<File[]>([]);
  const [unpairedIndexes, setUnpairedIndexes] = useState<File[]>([]);
  const [isDragging, setIsDragging] = useState(false);

  const bulkInputRef = useRef<HTMLInputElement>(null);
  const bamPairInputRef = useRef<HTMLInputElement>(null);

  if (!isOpen) return null;

  const mergeFiles = (incoming: File[]) => {
    const allExisting: File[] = [];
    if (fastaFile) allExisting.push(fastaFile);
    if (faiFile) allExisting.push(faiFile);
    if (gffFile) allExisting.push(gffFile);
    for (const s of strains) {
      allExisting.push(s.bamFile, s.indexFile);
    }
    allExisting.push(...unpairedBams, ...unpairedIndexes);

    const byName = new Map<string, File>();
    for (const f of allExisting) byName.set(f.name, f);
    for (const f of incoming) byName.set(f.name, f);

    const matched = matchDroppedFiles([...byName.values()]);
    if (matched.fastaFile) setFastaFile(matched.fastaFile);
    if (matched.faiFile) setFaiFile(matched.faiFile);
    if (matched.gffFile) setGffFile(matched.gffFile);

    const prevLabelMap = new Map(strains.map((s) => [s.bamFile.name, s.label]));
    setStrains(
      matched.pairedStrains.map((s) => ({
        ...s,
        label: prevLabelMap.get(s.bamFile.name) ?? s.label,
      })),
    );
    setUnpairedBams(matched.unpairedBams);
    setUnpairedIndexes(matched.unpairedIndexes);
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
    const files = Array.from(e.dataTransfer.files);
    if (files.length > 0) mergeFiles(files);
  };

  const updateStrainLabel = (id: string, label: string) => {
    setStrains((prev) => prev.map((s) => (s.id === id ? { ...s, label } : s)));
  };

  const removeStrain = (id: string) => {
    setStrains((prev) => prev.filter((s) => s.id !== id));
  };

  const handleReset = () => {
    setFastaFile(undefined);
    setFaiFile(undefined);
    setGffFile(undefined);
    setStrains([]);
    setUnpairedBams([]);
    setUnpairedIndexes([]);
  };

  const isCustomGenomeReady = Boolean(fastaFile);
  const isTracksReady = strains.length >= 1;
  const canSubmit = !isLoading && (isTracksReady || isCustomGenomeReady);

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/45 p-4 backdrop-blur-xs"
      onDragOver={(e) => {
        if (e.dataTransfer.types.includes('Files')) {
          e.preventDefault();
          e.stopPropagation();
        }
      }}
      onDrop={(e) => {
        e.preventDefault();
        e.stopPropagation();
        setIsDragging(false);
        const files = Array.from(e.dataTransfer.files);
        if (files.length > 0) mergeFiles(files);
      }}
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div className="flex max-h-[92vh] w-full max-w-3xl flex-col overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-2xl">
        {/* Header */}
        <div className="flex items-center justify-between border-b border-slate-200 bg-slate-50 px-6 py-4">
          <div className="flex items-center gap-3">
            <div className="flex h-10 w-10 items-center justify-center rounded-xl bg-sky-100 text-sky-700">
              <Dna className="h-5 w-5" />
            </div>
            <div>
              <h2 className="text-base font-bold text-slate-900">
                {mode === 'tracks'
                  ? `로컬 BAM / BAI 트랙 추가 (현재 지놈: ${activeGenomeLabel})`
                  : '사용자 Reference FASTA · GFF3 · BAM/BAI 불러오기'}
              </h2>
              <p className="text-xs text-slate-600">
                {mode === 'tracks'
                  ? 'sacCer3가 기본 탑재되어 있어 BAM + BAI 파일만 선택하시면 17개 전체 염색체에 바로 매핑됩니다.'
                  : '사내 PC의 Reference FASTA와 BAM/BAI를 서버 전송 없이 로컬 브라우저에서 직접 읽습니다.'}
              </p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="rounded-lg p-1.5 text-slate-500 hover:bg-slate-200 hover:text-slate-900"
          >
            <X className="h-5 w-5" />
          </button>
        </div>

        {/* Body */}
        <div className="flex-1 space-y-4 overflow-y-auto p-6">
          {/* Bulk Dropzone */}
          <div
            onDragOver={(e) => {
              e.preventDefault();
              setIsDragging(true);
            }}
            onDragLeave={() => setIsDragging(false)}
            onDrop={handleDrop}
            onClick={() => bulkInputRef.current?.click()}
            className={`cursor-pointer rounded-xl border-2 border-dashed p-5 text-center transition ${
              isDragging
                ? 'border-sky-500 bg-sky-50'
                : 'border-slate-300 bg-slate-50 hover:border-sky-400 hover:bg-sky-50/40'
            }`}
          >
            <input
              ref={bulkInputRef}
              type="file"
              multiple
              className="hidden"
              onChange={(e) => {
                if (e.target.files) mergeFiles(Array.from(e.target.files));
                e.target.value = '';
              }}
            />
            <Upload className="mx-auto mb-2 h-7 w-7 text-sky-600" />
            <p className="text-sm font-semibold text-slate-800">
              {mode === 'tracks'
                ? '여기에 균주 BAM + BAI 파일을 드래그 앤 드롭하세요 (1개 이상)'
                : '여기에 FASTA, GFF3, 균주별 BAM + BAI 파일을 한꺼번에 드래그 앤 드롭하세요'}
            </p>
            <p className="mt-1 text-xs text-slate-600">
              또는 클릭해 여러 파일 선택 (.bam + .bai/.csi 자동 페어링)
            </p>
          </div>

          {(unpairedBams.length > 0 || unpairedIndexes.length > 0) && (
            <div className="rounded-xl border border-amber-300 bg-amber-50 p-3 text-xs text-amber-900">
              <div className="flex items-start gap-2">
                <AlertCircle className="mt-0.5 h-4 w-4 shrink-0 text-amber-600" />
                <div className="space-y-1">
                  {unpairedBams.length > 0 && (
                    <p>
                      인덱스(.bai/.csi)가 없는 BAM:{' '}
                      <span className="font-mono font-semibold">
                        {unpairedBams.map((f) => f.name).join(', ')}
                      </span>{' '}
                      — 동일한 이름의 .bai 파일을 함께 선택해주세요.
                    </p>
                  )}
                  {unpairedIndexes.length > 0 && (
                    <p>
                      대응 BAM이 없는 인덱스:{' '}
                      <span className="font-mono">
                        {unpairedIndexes.map((f) => f.name).join(', ')}
                      </span>
                    </p>
                  )}
                </div>
              </div>
            </div>
          )}

          {isLoading && (
            <div className="flex items-center justify-between gap-2 rounded-xl border border-amber-300 bg-amber-50 p-3 text-xs font-medium text-amber-900">
              <div className="flex items-center gap-2 min-w-0">
                <span className="inline-block h-2.5 w-2.5 shrink-0 animate-pulse rounded-full bg-amber-500" />
                <span className="truncate">
                  {loadingMsg || '파일 읽는 중...'}
                  {elapsedSec >= 2 ? ` (${elapsedSec}초)` : ''}
                </span>
              </div>
              {onCancelLoad && (
                <button
                  type="button"
                  onClick={onCancelLoad}
                  className="shrink-0 rounded-lg border border-amber-400 bg-white px-2.5 py-1 text-[11px] font-semibold text-amber-900 hover:bg-amber-100"
                >
                  중단
                </button>
              )}
            </div>
          )}

          {errorMsg && (
            <div className="rounded-xl border border-rose-300 bg-rose-50 p-3 text-xs text-rose-900">
              <div className="flex items-center gap-2">
                <AlertCircle className="h-4 w-4 shrink-0 text-rose-600" />
                <span>{errorMsg}</span>
              </div>
            </div>
          )}

          {/* Strain BAM+BAI List */}
          <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
            <div className="mb-3 flex items-center justify-between">
              <span className="flex items-center gap-1.5 text-xs font-bold uppercase tracking-wide text-slate-800">
                <Layers className="h-3.5 w-3.5 text-sky-600" />
                추가할 균주 BAM + BAI ({strains.length}개 매칭됨)
              </span>
              <button
                type="button"
                onClick={() => bamPairInputRef.current?.click()}
                className="inline-flex items-center gap-1 rounded-lg border border-sky-200 bg-sky-50 px-2.5 py-1 text-xs font-semibold text-sky-700 hover:bg-sky-100"
              >
                <Plus className="h-3.5 w-3.5" /> BAM/BAI 선택
              </button>
              <input
                ref={bamPairInputRef}
                type="file"
                accept=".bam,.bai,.csi"
                multiple
                className="hidden"
                onChange={(e) => {
                  if (e.target.files) mergeFiles(Array.from(e.target.files));
                  e.target.value = '';
                }}
              />
            </div>

            {strains.length === 0 ? (
              <p className="py-3 text-center text-xs text-slate-500">
                1개 BAM만 추가하면 sacCer3 레퍼런스 단독 비교, 2개 이상 추가하면 균주 간 비교가 됩니다.
              </p>
            ) : (
              <div className="space-y-2">
                {strains.map((s, idx) => (
                  <div
                    key={s.id}
                    className="flex flex-wrap items-center justify-between gap-2 rounded-lg border border-slate-200 bg-white px-3 py-2"
                  >
                    <div className="flex min-w-0 flex-1 items-center gap-2.5">
                      <span className="flex h-6 w-6 shrink-0 items-center justify-center rounded-md bg-sky-100 text-xs font-bold text-sky-700">
                        {idx + 1}
                      </span>
                      <input
                        type="text"
                        value={s.label}
                        onChange={(e) => updateStrainLabel(s.id, e.target.value)}
                        className="w-40 rounded border border-slate-300 bg-white px-2 py-1 text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                      />
                      <div className="min-w-0 flex-1 truncate text-xs text-slate-600">
                        <span className="font-mono font-medium text-slate-800">{s.bamFile.name}</span>{' '}
                        <span>({formatBytes(s.bamFile.size)})</span>
                        <span className="mx-1.5 text-slate-400">+</span>
                        <span className="font-mono font-medium text-emerald-700">
                          {s.indexFile.name}
                        </span>
                      </div>
                    </div>
                    <button
                      type="button"
                      onClick={() => removeStrain(s.id)}
                      className="rounded p-1 text-slate-400 hover:bg-rose-50 hover:text-rose-600"
                    >
                      <Trash2 className="h-4 w-4" />
                    </button>
                  </div>
                ))}
              </div>
            )}
          </div>

          {/* Optional Custom Reference FASTA & GFF3 */}
          {mode === 'genome' && (
            <div className="grid gap-3 md:grid-cols-2">
              <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
                <div className="mb-2 flex items-center justify-between">
                  <span className="flex items-center gap-1.5 text-xs font-bold text-slate-800">
                    <FileCode2 className="h-3.5 w-3.5 text-emerald-600" />
                    커스텀 Reference FASTA (필수)
                  </span>
                  <label className="cursor-pointer text-xs font-semibold text-sky-700 hover:underline">
                    선택
                    <input
                      type="file"
                      accept=".fa,.fasta,.fna,.fsa,.gz,.fai"
                      multiple
                      className="hidden"
                      onChange={(e) => {
                        if (e.target.files) mergeFiles(Array.from(e.target.files));
                        e.target.value = '';
                      }}
                    />
                  </label>
                </div>
                {fastaFile ? (
                  <div className="space-y-1 text-xs">
                    <div className="flex items-center gap-1.5 font-medium text-emerald-800">
                      <CheckCircle2 className="h-3.5 w-3.5 shrink-0 text-emerald-600" />
                      <span className="truncate font-mono">{fastaFile.name}</span>
                    </div>
                    <div className="text-slate-600">
                      {faiFile ? `인덱스: ${faiFile.name}` : '자동 인덱싱 지원'}
                    </div>
                  </div>
                ) : (
                  <p className="text-xs text-slate-500">
                    S288C가 아닌 별도 균주(CEN.PK 등) FASTA를 쓸 때만 선택하세요.
                  </p>
                )}
              </div>

              <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
                <div className="mb-2 flex items-center justify-between">
                  <span className="flex items-center gap-1.5 text-xs font-bold text-slate-800">
                    <FolderOpen className="h-3.5 w-3.5 text-violet-600" />
                    GFF3 주석 파일 (선택)
                  </span>
                  <label className="cursor-pointer text-xs font-semibold text-sky-700 hover:underline">
                    선택
                    <input
                      type="file"
                      accept=".gff,.gff3,.gtf,.gz"
                      className="hidden"
                      onChange={(e) => {
                        if (e.target.files) mergeFiles(Array.from(e.target.files));
                        e.target.value = '';
                      }}
                    />
                  </label>
                </div>
                {gffFile ? (
                  <div className="flex items-center gap-1.5 text-xs font-medium text-violet-800">
                    <CheckCircle2 className="h-3.5 w-3.5 shrink-0 text-violet-600" />
                    <span className="truncate font-mono">{gffFile.name}</span>
                  </div>
                ) : (
                  <p className="text-xs text-slate-500">미첨부 시 S288C 내장 주석 자동 사용</p>
                )}
              </div>
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="flex flex-wrap items-center justify-between gap-3 border-t border-slate-200 bg-slate-50 px-6 py-4">
          <div className="flex items-center gap-2">
            {mode === 'genome' && (
              <>
                <button
                  type="button"
                  disabled={isLoading}
                  onClick={() => {
                    handleReset();
                    onSelectSacCer3();
                  }}
                  className="rounded-xl border border-sky-300 bg-sky-50 px-3 py-2 text-xs font-semibold text-sky-800 hover:bg-sky-100"
                >
                  sacCer3 기본 지놈으로 전환
                </button>
                <button
                  type="button"
                  disabled={isLoading}
                  onClick={() => {
                    handleReset();
                    onLoadDemo();
                  }}
                  className="inline-flex items-center gap-1 rounded-xl border border-emerald-300 bg-emerald-50 px-3 py-2 text-xs font-semibold text-emerald-800 hover:bg-emerald-100"
                >
                  <Play className="h-3.5 w-3.5" />
                  데모 데이터 (17chr)
                </button>
              </>
            )}
            {(fastaFile || strains.length > 0) && (
              <button
                type="button"
                onClick={handleReset}
                className="rounded-xl px-2.5 py-2 text-xs text-slate-600 hover:bg-slate-200"
              >
                목록 초기화
              </button>
            )}
          </div>

          <div className="flex items-center gap-2">
            <button
              type="button"
              onClick={onClose}
              className="rounded-xl border border-slate-300 bg-white px-4 py-2 text-xs font-medium text-slate-700 hover:bg-slate-100"
            >
              취소
            </button>
            <button
              type="button"
              disabled={!canSubmit}
              onClick={async () => {
                if (fastaFile) {
                  const ok = await onLoadCustomGenome({
                    fastaFile,
                    faiFile,
                    gffFile,
                    strains: strains.map((s) => ({
                      id: s.id,
                      label: s.label,
                      bamFile: s.bamFile,
                      indexFile: s.indexFile,
                    })),
                  });
                  if (ok) handleReset();
                } else if (strains.length > 0) {
                  const ok = await onAddTracksOnly(
                    strains.map((s) => ({
                      id: s.id,
                      label: s.label,
                      bamFile: s.bamFile,
                      indexFile: s.indexFile,
                    })),
                  );
                  if (ok) handleReset();
                }
              }}
              className="inline-flex items-center gap-2 rounded-xl bg-sky-600 px-5 py-2 text-xs font-semibold text-white shadow-sm hover:bg-sky-700 disabled:cursor-not-allowed disabled:bg-slate-300 disabled:text-slate-500"
            >
              {isLoading
                ? '불러오는 중...'
                : fastaFile
                  ? '커스텀 지놈 적용'
                  : `트랙 ${strains.length}개 추가하기`}
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};
