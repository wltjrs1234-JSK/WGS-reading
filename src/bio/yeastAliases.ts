import type { ContigInfo } from './types';

/**
 * Saccharomyces cerevisiae (S288C R64) 염색체 이름 정규화 및 별칭 매핑 유틸리티
 *
 * SGD, UCSC, Ensembl, NCBI RefSeq/GenBank 등 다양한 출처의 표기 및
 * 파이프 구분 헤더(ref|NC_001133|, gb|BK006935.2|), 0-패딩(chr01) 등을 모두 처리합니다.
 */

export interface YeastChromDef {
  canonical: string;
  roman: string;
  arabic: string;
  refseq: string;
  genbank: string;
  standardLength: number;
  isMito: boolean;
  aliases: string[];
}

export const YEAST_CHROMOSOMES: YeastChromDef[] = [
  { canonical: 'chrI', roman: 'I', arabic: '1', refseq: 'NC_001133.9', genbank: 'BK006935.2', standardLength: 230218, isMito: false, aliases: ['chri', 'chr1', 'chr01', 'i', '1', '01', 'nc_001133', 'nc_001133.9', 'bk006935', 'bk006935.2', 'chri_s288c'] },
  { canonical: 'chrII', roman: 'II', arabic: '2', refseq: 'NC_001134.8', genbank: 'BK006936.2', standardLength: 813184, isMito: false, aliases: ['chrii', 'chr2', 'chr02', 'ii', '2', '02', 'nc_001134', 'nc_001134.8', 'bk006936', 'bk006936.2'] },
  { canonical: 'chrIII', roman: 'III', arabic: '3', refseq: 'NC_001135.5', genbank: 'BK006937.2', standardLength: 316620, isMito: false, aliases: ['chriii', 'chr3', 'chr03', 'iii', '3', '03', 'nc_001135', 'nc_001135.5', 'bk006937', 'bk006937.2'] },
  { canonical: 'chrIV', roman: 'IV', arabic: '4', refseq: 'NC_001136.10', genbank: 'BK006938.2', standardLength: 1531933, isMito: false, aliases: ['chriv', 'chr4', 'chr04', 'iv', '4', '04', 'nc_001136', 'nc_001136.10', 'bk006938', 'bk006938.2'] },
  { canonical: 'chrV', roman: 'V', arabic: '5', refseq: 'NC_001137.3', genbank: 'BK006939.2', standardLength: 576874, isMito: false, aliases: ['chrv', 'chr5', 'chr05', 'v', '5', '05', 'nc_001137', 'nc_001137.3', 'bk006939', 'bk006939.2'] },
  { canonical: 'chrVI', roman: 'VI', arabic: '6', refseq: 'NC_001138.5', genbank: 'BK006940.2', standardLength: 270161, isMito: false, aliases: ['chrvi', 'chr6', 'chr06', 'vi', '6', '06', 'nc_001138', 'nc_001138.5', 'bk006940', 'bk006940.2'] },
  { canonical: 'chrVII', roman: 'VII', arabic: '7', refseq: 'NC_001139.9', genbank: 'BK006941.2', standardLength: 1090940, isMito: false, aliases: ['chrvii', 'chr7', 'chr07', 'vii', '7', '07', 'nc_001139', 'nc_001139.9', 'bk006941', 'bk006941.2'] },
  { canonical: 'chrVIII', roman: 'VIII', arabic: '8', refseq: 'NC_001140.6', genbank: 'BK006934.2', standardLength: 562643, isMito: false, aliases: ['chrviii', 'chr8', 'chr08', 'viii', '8', '08', 'nc_001140', 'nc_001140.6', 'bk006934', 'bk006934.2'] },
  { canonical: 'chrIX', roman: 'IX', arabic: '9', refseq: 'NC_001141.2', genbank: 'BK006942.2', standardLength: 439888, isMito: false, aliases: ['chrix', 'chr9', 'chr09', 'ix', '9', '09', 'nc_001141', 'nc_001141.2', 'bk006942', 'bk006942.2'] },
  { canonical: 'chrX', roman: 'X', arabic: '10', refseq: 'NC_001142.9', genbank: 'BK006943.2', standardLength: 745751, isMito: false, aliases: ['chrx', 'chr10', 'x', '10', 'nc_001142', 'nc_001142.9', 'bk006943', 'bk006943.2'] },
  { canonical: 'chrXI', roman: 'XI', arabic: '11', refseq: 'NC_001143.9', genbank: 'BK006944.2', standardLength: 666816, isMito: false, aliases: ['chrxi', 'chr11', 'xi', '11', 'nc_001143', 'nc_001143.9', 'bk006944', 'bk006944.2'] },
  { canonical: 'chrXII', roman: 'XII', arabic: '12', refseq: 'NC_001144.5', genbank: 'BK006945.2', standardLength: 1078177, isMito: false, aliases: ['chrxii', 'chr12', 'xii', '12', 'nc_001144', 'nc_001144.5', 'bk006945', 'bk006945.2'] },
  { canonical: 'chrXIII', roman: 'XIII', arabic: '13', refseq: 'NC_001145.3', genbank: 'BK006946.2', standardLength: 924431, isMito: false, aliases: ['chrxiii', 'chr13', 'xiii', '13', 'nc_001145', 'nc_001145.3', 'bk006946', 'bk006946.2'] },
  { canonical: 'chrXIV', roman: 'XIV', arabic: '14', refseq: 'NC_001146.8', genbank: 'BK006947.3', standardLength: 784333, isMito: false, aliases: ['chrxiv', 'chr14', 'xiv', '14', 'nc_001146', 'nc_001146.8', 'bk006947', 'bk006947.3'] },
  { canonical: 'chrXV', roman: 'XV', arabic: '15', refseq: 'NC_001147.6', genbank: 'BK006948.2', standardLength: 1091291, isMito: false, aliases: ['chrxv', 'chr15', 'xv', '15', 'nc_001147', 'nc_001147.6', 'bk006948', 'bk006948.2'] },
  { canonical: 'chrXVI', roman: 'XVI', arabic: '16', refseq: 'NC_001148.4', genbank: 'BK006949.2', standardLength: 948066, isMito: false, aliases: ['chrxvi', 'chr16', 'xvi', '16', 'nc_001148', 'nc_001148.4', 'bk006949', 'bk006949.2'] },
  { canonical: 'chrM', roman: 'Mito', arabic: 'MT', refseq: 'NC_001224.1', genbank: 'AJ011856.1', standardLength: 85779, isMito: true, aliases: ['chrm', 'chrmt', 'mt', 'mito', 'mitochondrion', 'chr17', '17', 'chrmito', 'nc_001224', 'nc_001224.1', 'kp263414', 'kp263414.1', 'aj011856', 'aj011856.1'] },
  { canonical: '2-micron', roman: '2-micron', arabic: '2u', refseq: 'NC_001398.1', genbank: 'J01347.1', standardLength: 6318, isMito: false, aliases: ['2-micron', '2micron', '2_micron', '2u', '2um', 'nc_001398', 'nc_001398.1', 'j01347', 'j01347.1'] },
];

const aliasLookup = new Map<string, YeastChromDef>();
for (const def of YEAST_CHROMOSOMES) {
  aliasLookup.set(def.canonical.toLowerCase(), def);
  for (const a of def.aliases) {
    aliasLookup.set(a.toLowerCase(), def);
  }
}

function extractCandidateTokens(rawName: string): string[] {
  const trimmed = rawName.trim().toLowerCase();
  const tokens = new Set<string>();
  const addToken = (t: string) => {
    const clean = t.trim();
    if (!clean) return;
    tokens.add(clean);
    tokens.add(clean.replace(/\.\d+$/, ''));
    // chromosome_I / chromosome1 / chr_I / chr_1 / sc_chr1 등 접두어 정리
    const noPrefix = clean
      .replace(/^chromosome[_-]?/, 'chr')
      .replace(/^chr[_-]/, 'chr')
      .replace(/^(s288c|saccer\d*|sc)[_-]/, '');
    tokens.add(noPrefix);
    tokens.add(noPrefix.replace(/\.\d+$/, ''));
    // chr01 -> chr1, 01 -> 1
    tokens.add(noPrefix.replace(/^(chr)?0+(\d+)$/, '$1$2'));
    if (noPrefix === 'mitochondria' || noPrefix === 'chr_mt') {
      tokens.add('chrmt');
    }
  };

  addToken(trimmed);
  if (trimmed.includes('|')) {
    for (const part of trimmed.split('|')) {
      addToken(part);
    }
  }
  return [...tokens];
}

export function findYeastDef(name: string): YeastChromDef | undefined {
  for (const tok of extractCandidateTokens(name)) {
    const found = aliasLookup.get(tok);
    if (found) return found;
  }
  return undefined;
}

/**
 * 표준 효모 염색체 이름(예: 'chrI', 'chrM') 반환. 목록에 없으면 원본 유지.
 */
export function normalizeYeastChromName(name: string): string {
  return findYeastDef(name)?.canonical ?? name.trim();
}

/**
 * 미토콘드리아 염색체 여부 판별 (NCBI Table 3 코돈표 적용 대상)
 */
export function isMitochondrialChrom(name: string): boolean {
  const def = findYeastDef(name);
  if (def) return def.isMito;
  return /^(chr)?m(t|ito|itochondrion)?$/i.test(name.trim());
}

/**
 * S. cerevisiae 컨티그 목록과 하나라도 일치하는지 확인 (내장 S288C 주석 자동 활성화 판단용)
 */
export function looksLikeYeastGenome(contigs: ContigInfo[]): boolean {
  let matchedCount = 0;
  for (const c of contigs) {
    const def = findYeastDef(c.name);
    if (def) {
      matchedCount++;
      continue;
    }
    if (YEAST_CHROMOSOMES.some((d) => Math.abs(d.standardLength - c.length) <= 5)) {
      matchedCount++;
    }
  }
  return matchedCount >= 1;
}

/**
 * 기준 이름(targetChrom)과 일치하는 트랙 내부 실제 컨티그 이름을 자동 매칭
 * 1) 수동 지정(manualOverrides)
 * 2) 정확 일치 / 대소문자 무시 일치
 * 3) S. cerevisiae 별칭 정규화 일치 (NCBI/GenBank/Ensembl/SGD/0-패딩)
 * 4) 선택적 길이 기반 유일 일치 (이름이 완전히 다른 커스텀 컨티그 대비)
 */
export function resolveMatchingRefName(
  targetChrom: string,
  availableRefNames: string[],
  manualOverrides?: Record<string, string>,
  targetLength?: number,
  availableContigs?: ContigInfo[],
): string | null {
  if (manualOverrides && manualOverrides[targetChrom]) {
    const mapped = manualOverrides[targetChrom];
    if (availableRefNames.includes(mapped)) return mapped;
  }

  if (availableRefNames.includes(targetChrom)) return targetChrom;

  const targetLower = targetChrom.trim().toLowerCase();
  const caseInsensitive = availableRefNames.find((n) => n.trim().toLowerCase() === targetLower);
  if (caseInsensitive) return caseInsensitive;

  const targetDef = findYeastDef(targetChrom);
  if (targetDef) {
    for (const candidate of availableRefNames) {
      const candDef = findYeastDef(candidate);
      if (candDef && candDef.canonical === targetDef.canonical) {
        return candidate;
      }
    }
  }

  // 이름이 완전히 다른 경우(예: contig_04 vs chrIV), 컨티그 길이가 일치(±50bp 이내)하는 유일 후보 연결
  if (targetLength && targetLength > 0 && availableContigs && availableContigs.length > 0) {
    const exactLenMatches = availableContigs.filter((c) => c.length === targetLength);
    if (exactLenMatches.length === 1) {
      return exactLenMatches[0].name;
    }
    const nearLenMatches = availableContigs.filter(
      (c) => Math.abs(c.length - targetLength) <= 50,
    );
    if (nearLenMatches.length === 1) {
      return nearLenMatches[0].name;
    }
  }

  return null;
}

const CANONICAL_ORDER = new Map<string, number>(
  YEAST_CHROMOSOMES.map((def, idx) => [def.canonical, idx]),
);

/**
 * 효모 염색체를 생물학적 표준 순서(chrI ~ chrXVI, chrM, 2-micron)로 정렬
 */
export function sortYeastContigs(contigs: ContigInfo[]): ContigInfo[] {
  return [...contigs].sort((a, b) => {
    const ca = normalizeYeastChromName(a.name);
    const cb = normalizeYeastChromName(b.name);
    const ia = CANONICAL_ORDER.get(ca) ?? 999;
    const ib = CANONICAL_ORDER.get(cb) ?? 999;
    if (ia !== ib) return ia - ib;
    return a.name.localeCompare(b.name, undefined, { numeric: true });
  });
}
