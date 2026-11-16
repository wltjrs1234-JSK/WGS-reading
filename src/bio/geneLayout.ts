import type { GffFeature } from './types';

export interface PackedGene {
  feature: GffFeature;
  row: number;
  drawLabel: boolean;
  labelWidthPx: number;
}

/**
 * IGV RefSeq/SGD 유전자 트랙에서 인접하거나 겹치는 유전자 막대·라벨이 충돌하지 않도록
 * 다중 레인(행)에 배치한다.
 * 좌표 0을 기준으로 한 픽셀 간격으로 계산되므로 좌우 팬/스크롤 중 행 번호가 흔들리지 않는다.
 */
export function packGeneFeatures(
  features: GffFeature[],
  pxPerBp: number,
  measureLabelWidth: (name: string) => number = (name) => name.length * 6.5,
): { packed: PackedGene[]; totalRows: number } {
  if (features.length === 0) {
    return { packed: [], totalRows: 0 };
  }

  const seen = new Set<string>();
  const unique: GffFeature[] = [];
  for (const feat of features) {
    const key = `${feat.chrom}:${feat.start}-${feat.end}:${feat.strand}:${feat.name}`;
    if (seen.has(key)) continue;
    seen.add(key);
    unique.push(feat);
  }

  const sorted = unique.sort(
    (a, b) =>
      a.start - b.start ||
      (b.end - b.start) - (a.end - a.start) ||
      a.name.localeCompare(b.name),
  );

  const includeLabelBounds = pxPerBp >= 0.03;
  const maxRows = includeLabelBounds ? 16 : 6;
  const rowEnds: number[] = [];
  const packed: PackedGene[] = [];

  for (const feat of sorted) {
    const barStart = feat.start * pxPerBp;
    const barEnd = Math.max(barStart + 4, feat.end * pxPerBp);
    const barMid = (barStart + barEnd) / 2;
    const labelWidthPx = measureLabelWidth(feat.name);

    let occStart: number;
    let occEnd: number;
    let drawLabel: boolean;

    if (includeLabelBounds) {
      occStart = Math.min(barStart, barMid - labelWidthPx / 2) - 7;
      occEnd = Math.max(barEnd, barMid + labelWidthPx / 2) + 7;
      drawLabel = true;
    } else {
      occStart = barStart - 2;
      occEnd = barEnd + 2;
      drawLabel = barEnd - barStart >= 22;
    }

    let placedRow = -1;
    for (let r = 0; r < rowEnds.length; r++) {
      if (rowEnds[r] <= occStart) {
        placedRow = r;
        rowEnds[r] = occEnd;
        break;
      }
    }

    if (placedRow === -1) {
      if (rowEnds.length < maxRows) {
        placedRow = rowEnds.length;
        rowEnds.push(occEnd);
      } else {
        // 희소 조감 뷰 등에서 상한 초과 시 가장 먼저 비는 레인 선택
        let earliestRow = 0;
        for (let r = 1; r < rowEnds.length; r++) {
          if (rowEnds[r] < rowEnds[earliestRow]) earliestRow = r;
        }
        placedRow = earliestRow;
        rowEnds[earliestRow] = Math.max(rowEnds[earliestRow], occEnd);
      }
    }

    packed.push({
      feature: feat,
      row: placedRow,
      drawLabel,
      labelWidthPx,
    });
  }

  return { packed, totalRows: rowEnds.length };
}
