export type UiScalePct = 95 | 100 | 108 | 115;

export interface FontPreset {
  id: string;
  label: string;
  subtitle: string;
  cssFamily: string;
  sampleText: string;
  stylesheetUrl: string | null;
  sourceBadge: '웹폰트' | '로컬/웹폰트' | 'Windows 기본';
}

export interface FontComboPreset {
  id: string;
  name: string;
  desc: string;
  uiPresetId: string;
  monoPresetId: string;
}

export interface FontSettings {
  uiPresetId: string;
  uiCustomFont: string;
  monoPresetId: string;
  monoCustomFont: string;
  uiScalePct: UiScalePct;
}

export interface ResolvedFontStacks {
  uiFontFamily: string;
  monoFontFamily: string;
  canvasDeltaPx: number;
}

export interface LocalFontScanResult {
  status: 'ok' | 'unsupported' | 'denied' | 'empty';
  families: string[];
  message?: string;
}

export const FONT_SETTINGS_STORAGE_KEY = 'yeast_wgs_font_settings_v1';
export const FONT_CHANGE_EVENT = 'wgs-font-change';

const GENERIC_FONT_KEYWORDS = new Set([
  'sans-serif',
  'serif',
  'monospace',
  'cursive',
  'fantasy',
  'system-ui',
  'ui-sans-serif',
  'ui-serif',
  'ui-monospace',
  'ui-rounded',
]);

export const UI_FONT_PRESETS: FontPreset[] = [
  {
    id: 'pretendard',
    label: 'Pretendard (프리텐다드)',
    subtitle: '현대적이고 또렷한 대시보드·UI 추천 표준 글꼴',
    cssFamily: '"Pretendard Variable", "Pretendard", "Noto Sans KR", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://cdn.jsdelivr.net/gh/orioncactus/pretendard@v1.3.9/dist/web/variable/pretendardvariable-dynamic-subset.min.css',
    sourceBadge: '웹폰트',
  },
  {
    id: 'ibm_plex_kr',
    label: 'IBM Plex Sans KR',
    subtitle: '정갈하고 직선미가 살아있는 연구·분석용 고딕',
    cssFamily: '"IBM Plex Sans KR", "Noto Sans KR", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=IBM+Plex+Sans+KR:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'suit',
    label: 'SUIT (수트)',
    subtitle: '작은 화면과 표에서도 숫자와 한글 균형이 우수한 고딕',
    cssFamily: '"SUIT Variable", "SUIT", "Noto Sans KR", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://cdn.jsdelivr.net/gh/sun-typeface/SUIT@2/fonts/variable/woff2/SUIT-Variable.css',
    sourceBadge: '웹폰트',
  },
  {
    id: 'noto_sans_kr',
    label: 'Noto Sans KR (본고딕)',
    subtitle: '안정적이고 차분한 구글·어도비 표준 고딕',
    cssFamily: '"Noto Sans KR", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Noto+Sans+KR:wght@400;500;600;700;800&display=swap',
    sourceBadge: '로컬/웹폰트',
  },
  {
    id: 'nanum_square_neo',
    label: '나눔스퀘어 네오 (NanumSquare Neo)',
    subtitle: '각진 모던 스타일로 제목·트랙 라벨 시인성이 높음',
    cssFamily: '"NanumSquareNeo", "Noto Sans KR", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://cdn.jsdelivr.net/gh/moonspam/NanumSquareNeo@1.0/nanumsquareneo.css',
    sourceBadge: '웹폰트',
  },
  {
    id: 'nanum_gothic',
    label: '나눔고딕 (Nanum Gothic)',
    subtitle: '부드럽고 친숙한 가독성의 클래식 한글 글꼴',
    cssFamily: '"Nanum Gothic", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Nanum+Gothic:wght@400;700;800&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'gowun_dodum',
    label: '고운돋움 (Gowun Dodum)',
    subtitle: '눈이 편안하고 둥근 인상의 장시간 판독용 글꼴',
    cssFamily: '"Gowun Dodum", "Noto Sans KR", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Gowun+Dodum&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'inter_noto',
    label: 'Inter + Noto Sans KR',
    subtitle: '유전자명·영문 숫자 시인성을 극대화한 영문+한글 조합',
    cssFamily: '"Inter", "Noto Sans KR", "Malgun Gothic", sans-serif',
    sampleText: 'YeastWGS · YER160C · NPT1 (254.3x)',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=Noto+Sans+KR:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'malgun',
    label: '맑은 고딕 (Windows 기본)',
    subtitle: '네트워크 없이 즉시 적용되는 Windows 표준 글꼴',
    cssFamily: '"Malgun Gothic", "Segoe UI", sans-serif',
    sampleText: 'YeastWGS · chrV 유전자 S288C 비교 뷰어 (146.2x)',
    stylesheetUrl: null,
    sourceBadge: 'Windows 기본',
  },
];

export const MONO_FONT_PRESETS: FontPreset[] = [
  {
    id: 'jetbrains_mono',
    label: 'JetBrains Mono',
    subtitle: '0/O · 1/l/I · 염기 A/T/G/C 구분이 뛰어난 서열 추천 글꼴',
    cssFamily: '"JetBrains Mono", "Consolas", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'fira_code',
    label: 'Fira Code',
    subtitle: '글자 폭이 시원하고 10-bp 청크 정렬 판독이 편한 고정폭',
    cssFamily: '"Fira Code", "Consolas", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Fira+Code:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'ibm_plex_mono',
    label: 'IBM Plex Mono',
    subtitle: '학술·논문 그림에 어울리는 샤프하고 균일한 고정폭',
    cssFamily: '"IBM Plex Mono", "Consolas", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'source_code_pro',
    label: 'Source Code Pro',
    subtitle: '둥글고 부드러운 윤곽으로 장시간 서열 비교에 적합',
    cssFamily: '"Source Code Pro", "Consolas", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Source+Code+Pro:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'roboto_mono',
    label: 'Roboto Mono',
    subtitle: '깔끔하고 단순한 Google 표준 고정폭 글꼴',
    cssFamily: '"Roboto Mono", "Consolas", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl:
      'https://fonts.googleapis.com/css2?family=Roboto+Mono:wght@400;500;600;700&display=swap',
    sourceBadge: '웹폰트',
  },
  {
    id: 'd2coding',
    label: 'D2Coding (네이버)',
    subtitle: '한글 2폭 + 영문/염기 1폭 비율이 완벽히 맞는 한국어 코딩 폰트',
    cssFamily: '"D2Coding", "Consolas", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl:
      'https://cdn.jsdelivr.net/gh/joungkyun/font-d2coding/d2coding.css',
    sourceBadge: '웹폰트',
  },
  {
    id: 'consolas',
    label: 'Consolas (Windows 기본)',
    subtitle: '네트워크 없이 즉시 적용되는 Windows 내장 고정폭',
    cssFamily: '"Consolas", "Courier New", monospace',
    sampleText: 'ATGTCAGAAC  TTGGATCCTT  491,001 bp',
    stylesheetUrl: null,
    sourceBadge: 'Windows 기본',
  },
];

export const FONT_COMBO_PRESETS: FontComboPreset[] = [
  {
    id: 'combo_standard',
    name: '바이오 표준 (추천)',
    desc: 'Pretendard + JetBrains Mono',
    uiPresetId: 'pretendard',
    monoPresetId: 'jetbrains_mono',
  },
  {
    id: 'combo_research',
    name: '샤프 연구용',
    desc: 'IBM Plex Sans KR + IBM Plex Mono',
    uiPresetId: 'ibm_plex_kr',
    monoPresetId: 'ibm_plex_mono',
  },
  {
    id: 'combo_clean_ui',
    name: '고밀도 대시보드',
    desc: 'SUIT + Fira Code',
    uiPresetId: 'suit',
    monoPresetId: 'fira_code',
  },
  {
    id: 'combo_hangul_coding',
    name: '한글·염기 정렬형',
    desc: '나눔스퀘어 네오 + D2Coding',
    uiPresetId: 'nanum_square_neo',
    monoPresetId: 'd2coding',
  },
  {
    id: 'combo_gene_en',
    name: '영문 유전자 특화',
    desc: 'Inter·Noto + Roboto Mono',
    uiPresetId: 'inter_noto',
    monoPresetId: 'roboto_mono',
  },
  {
    id: 'combo_offline',
    name: '오프라인 Windows 기본',
    desc: '맑은 고딕 + Consolas',
    uiPresetId: 'malgun',
    monoPresetId: 'consolas',
  },
];

export const DEFAULT_FONT_SETTINGS: FontSettings = {
  uiPresetId: 'pretendard',
  uiCustomFont: '',
  monoPresetId: 'jetbrains_mono',
  monoCustomFont: '',
  uiScalePct: 100,
};

export function sanitizeFontName(raw: string): string {
  return raw
    .replace(/[\r\n\t;"'{}<>\\]/g, ' ')
    .replace(/\s+/g, ' ')
    .trim()
    .slice(0, 80);
}

export function quoteCssFontName(name: string): string {
  const cleaned = sanitizeFontName(name);
  if (!cleaned) return '';
  if (GENERIC_FONT_KEYWORDS.has(cleaned.toLowerCase())) {
    return cleaned.toLowerCase();
  }
  return `"${cleaned}"`;
}

export function normalizeScalePct(val: unknown): UiScalePct {
  if (val === 95 || val === 100 || val === 108 || val === 115) return val;
  return 100;
}

export function getCanvasDeltaPx(scalePct: UiScalePct): number {
  switch (scalePct) {
    case 95:
      return -0.5;
    case 108:
      return 1;
    case 115:
      return 2;
    default:
      return 0;
  }
}

export function resolveFontStacks(settings: FontSettings): ResolvedFontStacks {
  const cleanUiCustom = sanitizeFontName(settings.uiCustomFont);
  let uiFontFamily: string;
  if (settings.uiPresetId === 'custom' && cleanUiCustom) {
    uiFontFamily = `${quoteCssFontName(cleanUiCustom)}, "Noto Sans KR", "Malgun Gothic", sans-serif`;
  } else {
    const preset =
      UI_FONT_PRESETS.find((p) => p.id === settings.uiPresetId) ?? UI_FONT_PRESETS[0];
    uiFontFamily = preset.cssFamily;
  }

  const cleanMonoCustom = sanitizeFontName(settings.monoCustomFont);
  let monoFontFamily: string;
  if (settings.monoPresetId === 'custom' && cleanMonoCustom) {
    monoFontFamily = `${quoteCssFontName(cleanMonoCustom)}, "Consolas", monospace`;
  } else {
    const preset =
      MONO_FONT_PRESETS.find((p) => p.id === settings.monoPresetId) ?? MONO_FONT_PRESETS[0];
    monoFontFamily = preset.cssFamily;
  }

  return {
    uiFontFamily,
    monoFontFamily,
    canvasDeltaPx: getCanvasDeltaPx(normalizeScalePct(settings.uiScalePct)),
  };
}

export function loadFontSettings(): FontSettings {
  try {
    if (typeof localStorage === 'undefined') return { ...DEFAULT_FONT_SETTINGS };
    const raw = localStorage.getItem(FONT_SETTINGS_STORAGE_KEY);
    if (!raw) return { ...DEFAULT_FONT_SETTINGS };
    const parsed = JSON.parse(raw) as Partial<FontSettings>;
    const uiPresetId =
      parsed.uiPresetId === 'custom' ||
      UI_FONT_PRESETS.some((p) => p.id === parsed.uiPresetId)
        ? String(parsed.uiPresetId)
        : DEFAULT_FONT_SETTINGS.uiPresetId;
    const monoPresetId =
      parsed.monoPresetId === 'custom' ||
      MONO_FONT_PRESETS.some((p) => p.id === parsed.monoPresetId)
        ? String(parsed.monoPresetId)
        : DEFAULT_FONT_SETTINGS.monoPresetId;

    return {
      uiPresetId,
      uiCustomFont: sanitizeFontName(String(parsed.uiCustomFont ?? '')),
      monoPresetId,
      monoCustomFont: sanitizeFontName(String(parsed.monoCustomFont ?? '')),
      uiScalePct: normalizeScalePct(parsed.uiScalePct),
    };
  } catch {
    return { ...DEFAULT_FONT_SETTINGS };
  }
}

export function saveFontSettings(settings: FontSettings): void {
  try {
    if (typeof localStorage !== 'undefined') {
      localStorage.setItem(FONT_SETTINGS_STORAGE_KEY, JSON.stringify(settings));
    }
  } catch {
    // Ignore storage quota or private browsing errors
  }
}

export function injectFontStylesheet(presetId: string, url: string | null): void {
  if (!url || typeof document === 'undefined' || !document.head) return;
  const selector = `link[data-wgs-font="${presetId}"]`;
  if (document.head.querySelector(selector)) return;

  const link = document.createElement('link');
  link.rel = 'stylesheet';
  link.href = url;
  link.setAttribute('data-wgs-font', presetId);
  link.onload = () => {
    window.dispatchEvent(new CustomEvent(FONT_CHANGE_EVENT));
  };
  link.onerror = () => {
    // Offline fallback stays silent; OS fallback stack remains active
  };
  document.head.appendChild(link);
}

export function preloadAllPresetStylesheets(): void {
  for (const p of UI_FONT_PRESETS) {
    if (p.stylesheetUrl) injectFontStylesheet(`ui_${p.id}`, p.stylesheetUrl);
  }
  for (const p of MONO_FONT_PRESETS) {
    if (p.stylesheetUrl) injectFontStylesheet(`mono_${p.id}`, p.stylesheetUrl);
  }
}

export function applyFontSettings(settings: FontSettings): ResolvedFontStacks {
  const stacks = resolveFontStacks(settings);

  if (typeof document !== 'undefined') {
    const uiPreset = UI_FONT_PRESETS.find((p) => p.id === settings.uiPresetId);
    const monoPreset = MONO_FONT_PRESETS.find((p) => p.id === settings.monoPresetId);

    if (uiPreset?.stylesheetUrl) {
      injectFontStylesheet(`ui_${uiPreset.id}`, uiPreset.stylesheetUrl);
    }
    if (monoPreset?.stylesheetUrl) {
      injectFontStylesheet(`mono_${monoPreset.id}`, monoPreset.stylesheetUrl);
    }

    const root = document.documentElement;
    root.style.setProperty('--wgs-font-sans', stacks.uiFontFamily);
    root.style.setProperty('--wgs-font-mono', stacks.monoFontFamily);
    root.style.setProperty('--font-sans', stacks.uiFontFamily);
    root.style.setProperty('--font-mono', stacks.monoFontFamily);
    root.style.setProperty('--wgs-canvas-delta-px', String(stacks.canvasDeltaPx));
    root.style.fontSize = `${((16 * normalizeScalePct(settings.uiScalePct)) / 100).toFixed(2)}px`;

    if (typeof window !== 'undefined') {
      window.dispatchEvent(new CustomEvent(FONT_CHANGE_EVENT, { detail: stacks }));
    }
    if (document.fonts?.ready) {
      void document.fonts.ready.then(() => {
        if (typeof window !== 'undefined') {
          window.dispatchEvent(new CustomEvent(FONT_CHANGE_EVENT, { detail: stacks }));
        }
      });
    }
  }

  return stacks;
}

export function getActiveCanvasFontConfig(): ResolvedFontStacks {
  if (typeof document === 'undefined') {
    return resolveFontStacks(DEFAULT_FONT_SETTINGS);
  }
  const styles = getComputedStyle(document.documentElement);
  const uiFontFamily =
    styles.getPropertyValue('--wgs-font-sans').trim() ||
    resolveFontStacks(DEFAULT_FONT_SETTINGS).uiFontFamily;
  const monoFontFamily =
    styles.getPropertyValue('--wgs-font-mono').trim() ||
    resolveFontStacks(DEFAULT_FONT_SETTINGS).monoFontFamily;
  const rawDelta = Number(styles.getPropertyValue('--wgs-canvas-delta-px').trim());
  const canvasDeltaPx = Number.isFinite(rawDelta) ? rawDelta : 0;

  return { uiFontFamily, monoFontFamily, canvasDeltaPx };
}

export function formatCanvasFont(
  familyStack: string,
  sizePx: number,
  weight: 'normal' | 'bold' | number = 'normal',
): string {
  const safeSize = Math.max(7.5, Math.round(sizePx * 10) / 10);
  return weight === 'normal'
    ? `${safeSize}px ${familyStack}`
    : `${weight} ${safeSize}px ${familyStack}`;
}

export async function scanLocalSystemFonts(): Promise<LocalFontScanResult> {
  if (typeof window === 'undefined') {
    return { status: 'unsupported', families: [] };
  }
  const queryFn = (window as unknown as { queryLocalFonts?: () => Promise<Array<{ family: string }>> })
    .queryLocalFonts;
  if (typeof queryFn !== 'function') {
    return {
      status: 'unsupported',
      families: [],
      message:
        '이 브라우저는 자동 폰트 스캔(Local Font Access API)을 지원하지 않습니다. 아래 입력창에 Windows 글꼴 이름을 직접 입력해 적용할 수 있습니다.',
    };
  }

  try {
    const fonts = await queryFn();
    const seen = new Set<string>();
    for (const f of fonts) {
      const cleaned = sanitizeFontName(f.family ?? '');
      if (cleaned) seen.add(cleaned);
    }
    const families = Array.from(seen).sort((a, b) => a.localeCompare(b, 'ko'));
    if (families.length === 0) {
      return {
        status: 'empty',
        families: [],
        message: '사용 가능한 로컬 폰트 권한이 허용되지 않았거나 목록이 비어 있습니다. 아래 입력창에서 직접 폰트명을 입력할 수 있습니다.',
      };
    }
    return { status: 'ok', families };
  } catch {
    return {
      status: 'denied',
      families: [],
      message: '브라우저에서 로컬 폰트 접근 권한이 거부되었습니다. 주소창 왼쪽 권한 아이콘에서 허용하거나 아래에 직접 폰트 이름을 입력해주세요.',
    };
  }
}

