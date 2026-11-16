import type { AlignedRead, CigarOp, CigarOperator, SoftClipInfo } from './types';

/**
 * 정규표현식 할당 없이 CIGAR 문자열을 순회하며 (op, len)을 전달합니다.
 */
export function forEachCigarOp(
  cigar: string,
  cb: (op: CigarOperator, len: number) => void,
): void {
  if (!cigar || cigar === '*') return;
  let num = 0;
  for (let i = 0; i < cigar.length; i++) {
    const ch = cigar.charCodeAt(i);
    if (ch >= 48 && ch <= 57) {
      num = num * 10 + (ch - 48);
    } else if (num > 0) {
      cb(cigar[i] as CigarOperator, num);
      num = 0;
    }
  }
}

export function parseCigar(cigar: string): CigarOp[] {
  if (!cigar || cigar === '*') return [];
  const ops: CigarOp[] = [];
  forEachCigarOp(cigar, (op, len) => {
    ops.push({ len, op });
  });
  return ops;
}

export function getSoftClipInfo(ops: CigarOp[]): SoftClipInfo {
  let left = 0;
  let right = 0;
  for (let i = 0; i < ops.length; i++) {
    const op = ops[i];
    if (op.op === 'S') {
      if (i === 0 || (i === 1 && ops[0].op === 'H')) {
        left = op.len;
      } else {
        right = op.len;
      }
    }
  }
  return { left, right };
}

export interface ReadSegment {
  kind: 'match' | 'del' | 'skip';
  start: number; // 0-based ref start
  end: number;   // 0-based ref end
}

export interface ReadInsertionMarker {
  anchorPos: number; // 0-based ref pos immediately before insertion (or read.start - 1)
  seq: string;
}

export interface ExpandedRead {
  segments: ReadSegment[];
  insertions: ReadInsertionMarker[];
  bases: Array<{ refPos: number; base: string; qual?: number }>;
  deletions: number[]; // refPos positions deleted by 'D'
  softClips: SoftClipInfo;
}

/**
 * CIGAR 연산을 레퍼런스 좌표계로 전개합니다.
 * rangeStart/rangeEnd가 주어지면 해당 윈도우 내부 정보만 수집해 성능을 최적화합니다.
 */
export function expandReadAlignment(
  read: AlignedRead,
  rangeStart = 0,
  rangeEnd = Number.MAX_SAFE_INTEGER,
): ExpandedRead {
  const ops = parseCigar(read.cigar);
  const seq = read.seq && read.seq !== '*' ? read.seq.toUpperCase() : '';
  const qual = read.qual;

  let refPos = read.start;
  let queryPos = 0;

  const segments: ReadSegment[] = [];
  const insertions: ReadInsertionMarker[] = [];
  const bases: Array<{ refPos: number; base: string; qual?: number }> = [];
  const deletions: number[] = [];

  for (const { op, len } of ops) {
    switch (op) {
      case 'M':
      case '=':
      case 'X': {
        const segStart = refPos;
        const segEnd = refPos + len;
        segments.push({ kind: 'match', start: segStart, end: segEnd });

        const overlapStart = Math.max(segStart, rangeStart);
        const overlapEnd = Math.min(segEnd, rangeEnd);
        for (let rp = overlapStart; rp < overlapEnd; rp++) {
          const qp = queryPos + (rp - segStart);
          const base = qp < seq.length ? seq[qp] : 'N';
          const q = qual && qp < qual.length ? qual[qp] : undefined;
          bases.push({ refPos: rp, base, qual: q });
        }

        refPos += len;
        queryPos += len;
        break;
      }
      case 'I': {
        const anchorPos = refPos - 1;
        if (anchorPos >= rangeStart - 1 && anchorPos < rangeEnd) {
          const insSeq =
            queryPos + len <= seq.length
              ? seq.slice(queryPos, queryPos + len)
              : 'N'.repeat(len);
          insertions.push({ anchorPos, seq: insSeq });
        }
        queryPos += len;
        break;
      }
      case 'D': {
        const segStart = refPos;
        const segEnd = refPos + len;
        segments.push({ kind: 'del', start: segStart, end: segEnd });

        const overlapStart = Math.max(segStart, rangeStart);
        const overlapEnd = Math.min(segEnd, rangeEnd);
        for (let rp = overlapStart; rp < overlapEnd; rp++) {
          deletions.push(rp);
        }
        refPos += len;
        break;
      }
      case 'N': {
        segments.push({ kind: 'skip', start: refPos, end: refPos + len });
        refPos += len;
        break;
      }
      case 'S': {
        queryPos += len;
        break;
      }
      case 'H':
      case 'P': {
        break;
      }
    }
  }

  return {
    segments,
    insertions,
    bases,
    deletions,
    softClips: getSoftClipInfo(ops),
  };
}

