import { describe, expect, it } from 'vitest';
import {
  DEFAULT_FONT_SETTINGS,
  formatCanvasFont,
  normalizeScalePct,
  quoteCssFontName,
  resolveFontStacks,
  sanitizeFontName,
} from '../fontSettings';

describe('fontSettings utilities', () => {
  it('sanitizes dangerous characters from custom font names', () => {
    expect(sanitizeFontName('  "Nanum"; {color:red} <script>  ')).toBe(
      'Nanum color:red script',
    );
    expect(sanitizeFontName('맑은\t고딕\nBold')).toBe('맑은 고딕 Bold');
  });

  it('quotes multi-word custom font names while leaving generic families unquoted', () => {
    expect(quoteCssFontName('Malgun Gothic')).toBe('"Malgun Gothic"');
    expect(quoteCssFontName('sans-serif')).toBe('sans-serif');
    expect(quoteCssFontName('MONOSPACE')).toBe('monospace');
  });

  it('resolves curated presets and custom font overrides cleanly', () => {
    const presetStacks = resolveFontStacks({
      uiPresetId: 'ibm_plex_kr',
      uiCustomFont: '',
      monoPresetId: 'fira_code',
      monoCustomFont: '',
      uiScalePct: 108,
    });
    expect(presetStacks.uiFontFamily).toContain('"IBM Plex Sans KR"');
    expect(presetStacks.monoFontFamily).toContain('"Fira Code"');
    expect(presetStacks.canvasDeltaPx).toBe(1);

    const customStacks = resolveFontStacks({
      uiPresetId: 'custom',
      uiCustomFont: 'Gmarket Sans',
      monoPresetId: 'custom',
      monoCustomFont: 'Cascadia Code',
      uiScalePct: 95,
    });
    expect(customStacks.uiFontFamily).toMatch(/^"Gmarket Sans",/);
    expect(customStacks.monoFontFamily).toMatch(/^"Cascadia Code",/);
    expect(customStacks.canvasDeltaPx).toBe(-0.5);
  });

  it('falls back to defaults on unknown preset ids or invalid scale values', () => {
    const fallback = resolveFontStacks({
      ...DEFAULT_FONT_SETTINGS,
      uiPresetId: 'nonexistent',
      monoPresetId: 'nonexistent',
      uiScalePct: normalizeScalePct(999),
    });
    expect(fallback.uiFontFamily).toContain('Pretendard');
    expect(fallback.monoFontFamily).toContain('JetBrains Mono');
    expect(fallback.canvasDeltaPx).toBe(0);
  });

  it('formats Canvas 2D font strings validly', () => {
    expect(formatCanvasFont('"Pretendard", sans-serif', 11, 600)).toBe(
      '600 11px "Pretendard", sans-serif',
    );
    expect(formatCanvasFont('"Consolas", monospace', 5, 'bold')).toBe(
      'bold 7.5px "Consolas", monospace',
    );
  });
});

