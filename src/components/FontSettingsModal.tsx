import React, { useEffect, useMemo, useState } from 'react';
import { Check, Laptop, RotateCcw, Search, Sparkles, Type, X } from 'lucide-react';
import {
  DEFAULT_FONT_SETTINGS,
  FONT_COMBO_PRESETS,
  MONO_FONT_PRESETS,
  UI_FONT_PRESETS,
  preloadAllPresetStylesheets,
  resolveFontStacks,
  sanitizeFontName,
  scanLocalSystemFonts,
  type FontSettings,
  type UiScalePct,
} from '../utils/fontSettings';

interface FontSettingsModalProps {
  isOpen: boolean;
  settings: FontSettings;
  onChangeSettings: (next: FontSettings) => void;
  onClose: () => void;
}

const SCALE_OPTIONS: { pct: UiScalePct; label: string; desc: string }[] = [
  { pct: 95, label: '95% (조밀)', desc: '한 화면에 더 많은 정보' },
  { pct: 100, label: '100% (표준)', desc: '기본 권장 크기' },
  { pct: 108, label: '108% (크게)', desc: '라벨·표가 더 큼직함' },
  { pct: 115, label: '115% (특대)', desc: '고해상도·노안 편안함' },
];

export const FontSettingsModal: React.FC<FontSettingsModalProps> = ({
  isOpen,
  settings,
  onChangeSettings,
  onClose,
}) => {
  const [uiCustomDraft, setUiCustomDraft] = useState(settings.uiCustomFont);
  const [monoCustomDraft, setMonoCustomDraft] = useState(settings.monoCustomFont);

  const [localFonts, setLocalFonts] = useState<string[]>([]);
  const [localScanMsg, setLocalScanMsg] = useState<string | null>(null);
  const [isScanningLocal, setIsScanningLocal] = useState(false);
  const [localFilter, setLocalFilter] = useState('');

  useEffect(() => {
    setUiCustomDraft(settings.uiCustomFont);
    setMonoCustomDraft(settings.monoCustomFont);
  }, [settings.uiCustomFont, settings.monoCustomFont]);

  useEffect(() => {
    if (isOpen) {
      preloadAllPresetStylesheets();
    }
  }, [isOpen]);

  const activeStacks = useMemo(() => resolveFontStacks(settings), [settings]);

  const filteredLocalFonts = useMemo(() => {
    const q = localFilter.trim().toLowerCase();
    if (!q) return localFonts;
    return localFonts.filter((f) => f.toLowerCase().includes(q));
  }, [localFonts, localFilter]);

  if (!isOpen) return null;

  const handleScanLocal = async () => {
    setIsScanningLocal(true);
    setLocalScanMsg(null);
    const res = await scanLocalSystemFonts();
    setLocalFonts(res.families);
    if (res.message) {
      setLocalScanMsg(res.message);
    } else if (res.families.length > 0) {
      setLocalScanMsg(`내 PC에서 ${res.families.length.toLocaleString()}개의 설치 글꼴을 불러왔습니다.`);
    }
    setIsScanningLocal(false);
  };

  const applyCustomUi = (name: string) => {
    const cleaned = sanitizeFontName(name);
    if (!cleaned) return;
    setUiCustomDraft(cleaned);
    onChangeSettings({
      ...settings,
      uiPresetId: 'custom',
      uiCustomFont: cleaned,
    });
  };

  const applyCustomMono = (name: string) => {
    const cleaned = sanitizeFontName(name);
    if (!cleaned) return;
    setMonoCustomDraft(cleaned);
    onChangeSettings({
      ...settings,
      monoPresetId: 'custom',
      monoCustomFont: cleaned,
    });
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/45 p-4 backdrop-blur-xs">
      <div className="flex max-h-[94vh] w-full max-w-5xl flex-col overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-2xl">
        {/* Header */}
        <div className="flex items-center justify-between border-b border-slate-200 bg-slate-50 px-6 py-4">
          <div className="flex items-center gap-3">
            <div className="flex h-10 w-10 items-center justify-center rounded-xl bg-sky-100 text-sky-700">
              <Type className="h-5 w-5" />
            </div>
            <div>
              <h2 className="text-lg font-bold text-slate-900">
                전체 화면 글꼴 및 가독성 설정
              </h2>
              <p className="text-xs text-slate-600">
                메인 유전체 뷰어(Canvas), 메뉴, 표, 서열 정렬 모달 전체에 즉시 반영되며 브라우저에 자동 저장됩니다.
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
        <div className="flex-1 space-y-5 overflow-y-auto p-6">
          {/* Quick Recommended Combos & Scale */}
          <div className="grid grid-cols-1 gap-4 lg:grid-cols-12">
            <div className="rounded-xl border border-slate-200 bg-slate-50/70 p-4 lg:col-span-8">
              <div className="mb-2.5 flex items-center gap-1.5 text-xs font-bold text-slate-800">
                <Sparkles className="h-4 w-4 text-sky-600" />
                원클릭 추천 글꼴 조합 (UI + 서열)
              </div>
              <div className="grid grid-cols-2 gap-2 sm:grid-cols-3">
                {FONT_COMBO_PRESETS.map((combo) => {
                  const isActive =
                    settings.uiPresetId === combo.uiPresetId &&
                    settings.monoPresetId === combo.monoPresetId;
                  return (
                    <button
                      key={combo.id}
                      type="button"
                      onClick={() =>
                        onChangeSettings({
                          ...settings,
                          uiPresetId: combo.uiPresetId,
                          monoPresetId: combo.monoPresetId,
                        })
                      }
                      className={`rounded-lg border p-2.5 text-left transition ${
                        isActive
                          ? 'border-sky-500 bg-sky-50 ring-2 ring-sky-200'
                          : 'border-slate-200 bg-white hover:border-sky-300 hover:bg-slate-50'
                      }`}
                    >
                      <div className="flex items-center justify-between gap-1 text-xs font-bold text-slate-900">
                        <span>{combo.name}</span>
                        {isActive && <Check className="h-3.5 w-3.5 text-sky-600" />}
                      </div>
                      <div className="mt-0.5 truncate text-[11px] text-slate-600">
                        {combo.desc}
                      </div>
                    </button>
                  );
                })}
              </div>
            </div>

            <div className="rounded-xl border border-slate-200 bg-slate-50/70 p-4 lg:col-span-4">
              <div className="mb-2.5 text-xs font-bold text-slate-800">
                전체 글자 크기 배율 (화면 + 트랙)
              </div>
              <div className="grid grid-cols-2 gap-2">
                {SCALE_OPTIONS.map((opt) => {
                  const active = settings.uiScalePct === opt.pct;
                  return (
                    <button
                      key={opt.pct}
                      type="button"
                      onClick={() =>
                        onChangeSettings({
                          ...settings,
                          uiScalePct: opt.pct,
                        })
                      }
                      className={`rounded-lg border px-2.5 py-2 text-left transition ${
                        active
                          ? 'border-sky-500 bg-sky-600 text-white'
                          : 'border-slate-200 bg-white text-slate-800 hover:bg-slate-100'
                      }`}
                    >
                      <div className="text-xs font-bold">{opt.label}</div>
                      <div
                        className={`text-[11px] ${active ? 'text-sky-100' : 'text-slate-500'}`}
                      >
                        {opt.desc}
                      </div>
                    </button>
                  );
                })}
              </div>
            </div>
          </div>

          {/* Live Combined Preview */}
          <div className="rounded-xl border border-indigo-200 bg-indigo-50/40 p-4">
            <div className="mb-2 flex flex-wrap items-center justify-between gap-2 text-xs font-bold text-indigo-900">
              <span>현재 선택 글꼴 실시간 미리보기</span>
              <span className="font-normal text-indigo-700">
                UI: {settings.uiPresetId === 'custom' ? settings.uiCustomFont : UI_FONT_PRESETS.find((p) => p.id === settings.uiPresetId)?.label}
                {' · '}
                서열: {settings.monoPresetId === 'custom' ? settings.monoCustomFont : MONO_FONT_PRESETS.find((p) => p.id === settings.monoPresetId)?.label}
              </span>
            </div>
            <div className="flex flex-col gap-2.5 rounded-lg border border-slate-200 bg-white p-3.5 sm:flex-row sm:items-center sm:justify-between">
              <div style={{ fontFamily: activeStacks.uiFontFamily }} className="space-y-1">
                <div className="flex flex-wrap items-center gap-2">
                  <span className="rounded bg-sky-600 px-2 py-0.5 text-xs font-semibold text-white">
                    NPT1 (YOR209C)
                  </span>
                  <span className="text-sm font-bold text-slate-900">
                    117_WT_rmdup · 효모 유전체 변이 비교 및 단백질 번역
                  </span>
                </div>
                <div className="text-xs text-slate-600">
                  전체 글자 배율 {settings.uiScalePct}% · 메뉴·표·캔버스 트랙 라벨에 실시간 반영됩니다.
                </div>
              </div>
              <div
                style={{ fontFamily: activeStacks.monoFontFamily }}
                className="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2 text-sm text-slate-900"
              >
                <div className="text-xs font-semibold text-slate-500">chrXV:736,437–737,726</div>
                <div className="mt-0.5 tracking-wider">
                  <span className="mr-3">ATGTCAGAAC</span>
                  <span className="mr-3">TTGGATCCTT</span>
                  <span className="rounded bg-sky-200 px-1 font-bold text-sky-950">G</span>
                  <span>TCGT</span>
                </div>
              </div>
            </div>
          </div>

          {/* Detailed Two-Column Selector: UI Sans & Sequence Mono */}
          <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
            {/* Left Column: UI / Korean Sans */}
            <div className="space-y-3 rounded-xl border border-slate-200 bg-slate-50/50 p-4">
              <div className="flex items-center justify-between">
                <h3 className="text-sm font-bold text-slate-900">
                  1. 전체 UI · 한글 글꼴 (메뉴 · 유전자명 · 표)
                </h3>
              </div>

              <div className="grid max-h-72 grid-cols-1 gap-2 overflow-y-auto pr-1">
                {UI_FONT_PRESETS.map((preset) => {
                  const selected = settings.uiPresetId === preset.id;
                  return (
                    <button
                      key={preset.id}
                      type="button"
                      onClick={() =>
                        onChangeSettings({
                          ...settings,
                          uiPresetId: preset.id,
                        })
                      }
                      className={`rounded-lg border p-3 text-left transition ${
                        selected
                          ? 'border-sky-500 bg-white ring-2 ring-sky-200'
                          : 'border-slate-200 bg-white hover:border-sky-300'
                      }`}
                    >
                      <div className="flex items-center justify-between gap-2">
                        <span className="text-xs font-bold text-slate-900">{preset.label}</span>
                        <div className="flex items-center gap-1.5">
                          <span className="rounded bg-slate-100 px-1.5 py-0.5 text-[11px] font-medium text-slate-600">
                            {preset.sourceBadge}
                          </span>
                          {selected && <Check className="h-4 w-4 text-sky-600" />}
                        </div>
                      </div>
                      <div className="mt-0.5 text-[11px] text-slate-500">{preset.subtitle}</div>
                      <div
                        style={{ fontFamily: preset.cssFamily }}
                        className="mt-2 rounded bg-slate-50 px-2.5 py-1.5 text-sm font-medium text-slate-900"
                      >
                        {preset.sampleText}
                      </div>
                    </button>
                  );
                })}
              </div>

              {/* Direct Custom UI Font Input */}
              <div className="rounded-lg border border-slate-200 bg-white p-3">
                <label className="mb-1.5 block text-xs font-semibold text-slate-700">
                  UI 글꼴 직접 입력 (내 PC에 설치된 폰트명)
                </label>
                <div className="flex gap-2">
                  <input
                    type="text"
                    value={uiCustomDraft}
                    onChange={(e) => setUiCustomDraft(e.target.value)}
                    onKeyDown={(e) => {
                      if (e.key === 'Enter') {
                        e.preventDefault();
                        applyCustomUi(uiCustomDraft);
                      }
                    }}
                    placeholder="예: Malgun Gothic, 나눔바른고딕, Arial"
                    className="min-w-0 flex-1 rounded-lg border border-slate-300 px-2.5 py-1.5 text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                  />
                  <button
                    type="button"
                    onClick={() => applyCustomUi(uiCustomDraft)}
                    className="rounded-lg bg-slate-800 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-900"
                  >
                    적용
                  </button>
                </div>
                {settings.uiPresetId === 'custom' && settings.uiCustomFont && (
                  <p className="mt-1 text-[11px] font-semibold text-sky-700">
                    현재 사용자 지정 글꼴 적용 중: {settings.uiCustomFont}
                  </p>
                )}
              </div>
            </div>

            {/* Right Column: Sequence / Coordinate Monospace */}
            <div className="space-y-3 rounded-xl border border-slate-200 bg-slate-50/50 p-4">
              <div className="flex items-center justify-between">
                <h3 className="text-sm font-bold text-slate-900">
                  2. 서열 · 좌표 고정폭 글꼴 (DNA · 아미노산 · 좌표)
                </h3>
              </div>

              <div className="grid max-h-72 grid-cols-1 gap-2 overflow-y-auto pr-1">
                {MONO_FONT_PRESETS.map((preset) => {
                  const selected = settings.monoPresetId === preset.id;
                  return (
                    <button
                      key={preset.id}
                      type="button"
                      onClick={() =>
                        onChangeSettings({
                          ...settings,
                          monoPresetId: preset.id,
                        })
                      }
                      className={`rounded-lg border p-3 text-left transition ${
                        selected
                          ? 'border-sky-500 bg-white ring-2 ring-sky-200'
                          : 'border-slate-200 bg-white hover:border-sky-300'
                      }`}
                    >
                      <div className="flex items-center justify-between gap-2">
                        <span className="text-xs font-bold text-slate-900">{preset.label}</span>
                        <div className="flex items-center gap-1.5">
                          <span className="rounded bg-slate-100 px-1.5 py-0.5 text-[11px] font-medium text-slate-600">
                            {preset.sourceBadge}
                          </span>
                          {selected && <Check className="h-4 w-4 text-sky-600" />}
                        </div>
                      </div>
                      <div className="mt-0.5 text-[11px] text-slate-500">{preset.subtitle}</div>
                      <div
                        style={{ fontFamily: preset.cssFamily }}
                        className="mt-2 rounded bg-slate-50 px-2.5 py-1.5 text-sm font-medium tracking-wide text-slate-900"
                      >
                        {preset.sampleText}
                      </div>
                    </button>
                  );
                })}
              </div>

              {/* Direct Custom Mono Font Input */}
              <div className="rounded-lg border border-slate-200 bg-white p-3">
                <label className="mb-1.5 block text-xs font-semibold text-slate-700">
                  고정폭 서열 글꼴 직접 입력
                </label>
                <div className="flex gap-2">
                  <input
                    type="text"
                    value={monoCustomDraft}
                    onChange={(e) => setMonoCustomDraft(e.target.value)}
                    onKeyDown={(e) => {
                      if (e.key === 'Enter') {
                        e.preventDefault();
                        applyCustomMono(monoCustomDraft);
                      }
                    }}
                    placeholder="예: Consolas, Cascadia Mono, Courier New"
                    className="min-w-0 flex-1 rounded-lg border border-slate-300 px-2.5 py-1.5 text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                  />
                  <button
                    type="button"
                    onClick={() => applyCustomMono(monoCustomDraft)}
                    className="rounded-lg bg-slate-800 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-900"
                  >
                    적용
                  </button>
                </div>
                {settings.monoPresetId === 'custom' && settings.monoCustomFont && (
                  <p className="mt-1 text-[11px] font-semibold text-sky-700">
                    현재 사용자 지정 서열 글꼴 적용 중: {settings.monoCustomFont}
                  </p>
                )}
              </div>
            </div>
          </div>

          {/* Scan Installed Local System Fonts (Chrome / Edge Local Font Access API) */}
          <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
            <div className="flex flex-wrap items-center justify-between gap-2">
              <div className="flex items-center gap-2">
                <Laptop className="h-4 w-4 text-slate-700" />
                <span className="text-xs font-bold text-slate-900">
                  내 PC에 설치된 폰트에서 직접 고르기 (Windows/Mac 로컬 폰트)
                </span>
              </div>
              <button
                type="button"
                disabled={isScanningLocal}
                onClick={() => void handleScanLocal()}
                className="inline-flex items-center gap-1.5 rounded-lg border border-sky-300 bg-white px-3 py-1.5 text-xs font-semibold text-sky-800 hover:bg-sky-50 disabled:opacity-50"
              >
                <Laptop className="h-3.5 w-3.5" />
                {isScanningLocal ? '폰트 스캔 중...' : '내 PC 설치 폰트 스캔'}
              </button>
            </div>

            {localScanMsg && (
              <p className="mt-2 text-xs text-slate-600">{localScanMsg}</p>
            )}

            {localFonts.length > 0 && (
              <div className="mt-3 space-y-2.5">
                <div className="relative max-w-sm">
                  <Search className="absolute left-2.5 top-2 h-3.5 w-3.5 text-slate-400" />
                  <input
                    type="text"
                    value={localFilter}
                    onChange={(e) => setLocalFilter(e.target.value)}
                    placeholder="설치된 폰트명 검색 (예: Noto, 맑은, Consolas)"
                    className="w-full rounded-lg border border-slate-300 bg-white py-1.5 pl-8 pr-3 text-xs text-slate-900 focus:border-sky-500 focus:outline-none"
                  />
                </div>
                <div className="grid max-h-48 grid-cols-1 gap-2 overflow-y-auto sm:grid-cols-2 md:grid-cols-3">
                  {filteredLocalFonts.slice(0, 120).map((fam) => (
                    <div
                      key={fam}
                      className="flex items-center justify-between gap-2 rounded-lg border border-slate-200 bg-white px-2.5 py-2"
                    >
                      <span
                        style={{ fontFamily: `"${fam}", sans-serif` }}
                        className="min-w-0 flex-1 truncate text-xs font-medium text-slate-900"
                        title={fam}
                      >
                        {fam}
                      </span>
                      <div className="flex shrink-0 items-center gap-1">
                        <button
                          type="button"
                          onClick={() => applyCustomUi(fam)}
                          className="rounded border border-sky-200 bg-sky-50 px-1.5 py-0.5 text-[11px] font-semibold text-sky-800 hover:bg-sky-100"
                        >
                          UI
                        </button>
                        <button
                          type="button"
                          onClick={() => applyCustomMono(fam)}
                          className="rounded border border-indigo-200 bg-indigo-50 px-1.5 py-0.5 text-[11px] font-semibold text-indigo-800 hover:bg-indigo-100"
                        >
                          서열
                        </button>
                      </div>
                    </div>
                  ))}
                </div>
                {filteredLocalFonts.length > 120 && (
                  <p className="text-[11px] text-slate-500">
                    상위 120개만 표시 중입니다. 검색어로 원하는 폰트 이름을 좁혀보세요.
                  </p>
                )}
              </div>
            )}
          </div>
        </div>

        {/* Footer */}
        <div className="flex flex-wrap items-center justify-between gap-2 border-t border-slate-200 bg-slate-50 px-6 py-3.5">
          <button
            type="button"
            onClick={() => onChangeSettings({ ...DEFAULT_FONT_SETTINGS })}
            className="inline-flex items-center gap-1.5 rounded-xl border border-slate-300 bg-white px-3.5 py-2 text-xs font-semibold text-slate-700 hover:bg-slate-100"
          >
            <RotateCcw className="h-3.5 w-3.5" />
            기본값으로 되돌리기
          </button>

          <button
            type="button"
            onClick={onClose}
            className="rounded-xl bg-sky-600 px-5 py-2 text-xs font-semibold text-white hover:bg-sky-700"
          >
            적용 완료 (닫기)
          </button>
        </div>
      </div>
    </div>
  );
};
