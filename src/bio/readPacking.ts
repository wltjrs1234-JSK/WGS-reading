import type { AlignedRead, BaseCounts, ConsensusPositionCall } from './types';

export const MAX_VISIBLE_PILEUP_ROWS = 28;

export interface PackedRead extends AlignedRead {
  row: number;
}

export interface AlleleTick {
  pos: number;
  depth: number;
  counts: BaseCounts;
}

function heapPush(heap: number[], val: number): void {
  heap.push(val);
  let i = heap.length - 1;
  while (i > 0) {
    const p = (i - 1) >> 1;
    if (heap[p] <= val) break;
    heap[i] = heap[p];
    i = p;
  }
  heap[i] = val;
}

function heapReplaceTop(heap: number[], val: number): void {
  let i = 0;
  const n = heap.length;
  while (true) {
    const l = (i << 1) + 1;
    if (l >= n) break;
    const r = l + 1;
    const smallest = r < n && heap[r] < heap[l] ? r : l;
    if (val <= heap[smallest]) break;
    heap[i] = heap[smallest];
    i = smallest;
  }
  heap[i] = val;
}

/**
 * 상단 가시 행(0 .. maxVisibleRows - 1)은 IGV와 동일하게 위에서부터 빈 레인을 우선 배정하고,
 * 가시 영역을 초과하는 고심도 리드는 Min-Heap으로 O(N log D) 시간에 총 심도만 집계합니다.
 */
export function packReadsWithHeap(
  reads: AlignedRead[],
  maxVisibleRows = MAX_VISIBLE_PILEUP_ROWS,
): { visiblePacked: PackedRead[]; totalRows: number } {
  const sorted = reads
    .filter((r) => (r.flags & 0x704) === 0)
    .sort((a, b) => a.start - b.start || a.end - b.end);

  const visibleEnds: number[] = [];
  const overflowHeap: number[] = [];
  const visiblePacked: PackedRead[] = [];

  for (const r of sorted) {
    let placedVisible = -1;
    for (let i = 0; i < visibleEnds.length; i++) {
      if (visibleEnds[i] + 2 <= r.start) {
        placedVisible = i;
        visibleEnds[i] = r.end;
        break;
      }
    }

    if (placedVisible !== -1) {
      visiblePacked.push({ ...r, row: placedVisible });
      continue;
    }

    if (visibleEnds.length < maxVisibleRows) {
      const row = visibleEnds.length;
      visibleEnds.push(r.end);
      visiblePacked.push({ ...r, row });
      continue;
    }

    // 가시 행(28줄)을 넘는 리드는 화면 렌더링 전송에서 생략하고 총 행 수(totalRows)만 O(log D) 추적
    if (overflowHeap.length > 0 && overflowHeap[0] + 2 <= r.start) {
      heapReplaceTop(overflowHeap, r.end);
    } else {
      heapPush(overflowHeap, r.end);
    }
  }

  return {
    visiblePacked,
    totalRows: visibleEnds.length + overflowHeap.length,
  };
}

/**
 * Worker에서 이미 계산된 포지션별 Consensus BaseCounts로부터 커버리지 변이 틱을 추출합니다.
 */
export function deriveAlleleTicksFromCalls(calls: ConsensusPositionCall[]): AlleleTick[] {
  const ticks: AlleleTick[] = [];
  for (const c of calls) {
    const cnt = c.counts;
    const depth = cnt.A + cnt.C + cnt.G + cnt.T + cnt.N;
    if (depth < 2) continue;
    const refBase = c.refBase;
    const nonRef =
      depth -
      (refBase === 'A'
        ? cnt.A
        : refBase === 'C'
          ? cnt.C
          : refBase === 'G'
            ? cnt.G
            : refBase === 'T'
              ? cnt.T
              : 0);
    if (nonRef >= 2 && nonRef / depth >= 0.2) {
      ticks.push({ pos: c.pos, depth, counts: cnt });
    }
  }
  return ticks;
}

