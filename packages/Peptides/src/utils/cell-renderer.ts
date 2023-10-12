import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import * as types from './types';
import {PositionStats, MonomerPositionStats, CLUSTER_TYPE} from '../model';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';

export function renderCellSelection(canvasContext: CanvasRenderingContext2D, bound: DG.Rect): void {
  canvasContext.strokeStyle = '#000';
  canvasContext.lineWidth = 1;
  canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
}

/** A function that sets amino acid residue cell renderer to the specified column */
export function setMonomerRenderer(col: DG.Column, alphabet: string): void {
  col.semType = C.SEM_TYPES.MONOMER;
  col.setTag('cell.renderer', C.SEM_TYPES.MONOMER);
  col.tags[C.TAGS.ALPHABET] = alphabet;
}

export function renderMutationCliffCell(canvasContext: CanvasRenderingContext2D, currentMonomer: string,
  currentPosition: string, monomerPositionStats: MonomerPositionStats, bound: DG.Rect,
  mutationCliffsSelection: types.Selection, substitutionsInfo: types.MutationCliffs | null = null,
  _twoColorMode: boolean = false, renderNums: boolean = true): void {
  const positionStats = monomerPositionStats[currentPosition];
  const pVal = positionStats![currentMonomer]!.pValue;
  const currentMeanDifference = positionStats![currentMonomer]!.meanDifference;

  // Transform p-value to increase intensity for smaller values and decrease for larger values
  const maxPValComplement = 1 - positionStats!.general.maxPValue;
  const minPValComplement = 1 - positionStats!.general.minPValue;
  const pValCentering = Math.min(maxPValComplement, minPValComplement);
  const centeredMaxPValComplement = maxPValComplement - pValCentering;
  const centeredMinPValComplement = minPValComplement - pValCentering;
  const centeredPValLimit = Math.max(centeredMaxPValComplement, centeredMinPValComplement);
  const pValComplement = pVal === null ? 0 : 1 - pVal - pValCentering;

  const coef = DG.Color.toHtml(pVal === null ? DG.Color.lightLightGray :
    DG.Color.scaleColor(currentMeanDifference >= 0 ? pValComplement : -pValComplement, -centeredPValLimit, centeredPValLimit));

  const maxMeanDifference = Math.max(Math.abs(monomerPositionStats.general.minMeanDifference), monomerPositionStats.general.maxMeanDifference);
  const rCoef = Math.abs(currentMeanDifference) / maxMeanDifference;
  const maxRadius = 0.9 * (bound.width > bound.height ? bound.height : bound.width) / 2;
  const radius = Math.floor(maxRadius * rCoef);

  const midX = bound.x + bound.width / 2;
  const midY = bound.y + bound.height / 2;
  canvasContext.beginPath();
  canvasContext.fillStyle = coef;
  canvasContext.arc(midX, midY, radius < 3 || pVal === null ? 3 : radius, 0, Math.PI * 2, true);
  canvasContext.closePath();
  canvasContext.fill();

  if (renderNums) {
    const substitutions = substitutionsInfo?.get(currentMonomer)?.get(currentPosition)?.entries() ?? null;
    if (substitutions !== null) {
      canvasContext.textBaseline = 'middle';
      canvasContext.textAlign = 'center';
      canvasContext.fillStyle = DG.Color.toHtml(DG.Color.black);
      canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
      canvasContext.shadowBlur = 5;
      canvasContext.shadowColor = DG.Color.toHtml(DG.Color.white);
      const uniqueValues = new Set<number>();

      for (const [key, value] of substitutions) {
        uniqueValues.add(key);
        for (const val of value)
          uniqueValues.add(val);
      }
      if (uniqueValues.size !== 0)
        canvasContext.fillText(uniqueValues.size.toString(), midX, midY);
    }
  }

  const monomerSelection = mutationCliffsSelection[currentPosition];
  if (monomerSelection && monomerSelection.includes(currentMonomer))
    renderCellSelection(canvasContext, bound);
}

export function renderInvaraintMapCell(canvasContext: CanvasRenderingContext2D, currentMonomer: string,
  currentPosition: string, invariantMapSelection: types.Selection, cellValue: number, bound: DG.Rect,
  color: number): void {
  //FIXME: This is a hack, because `color` value sometimes comes incomplete. E.g. we found that here `color` value is
  // 255 and its contrast color would be black, which is not visible on blue (color code) background. The full number
  // is actually 4278190335.
  color = DG.Color.fromHtml(DG.Color.toHtml(color));
  canvasContext.fillStyle = DG.Color.toHtml(color);
  canvasContext.fillRect(bound.x, bound.y, bound.width, bound.height);
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(color));
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  const monomerSelection = invariantMapSelection[currentPosition];
  if (monomerSelection && monomerSelection.includes(currentMonomer))
    renderCellSelection(canvasContext, bound);
}

export function renderLogoSummaryCell(canvasContext: CanvasRenderingContext2D, cellValue: string,
  clusterSelection: types.Selection, bound: DG.Rect): void {
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = '#000';
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  if (clusterSelection[CLUSTER_TYPE.CUSTOM].includes(cellValue) || clusterSelection[CLUSTER_TYPE.ORIGINAL].includes(cellValue))
    renderCellSelection(canvasContext, bound);
}


export function drawLogoInBounds(ctx: CanvasRenderingContext2D, bounds: DG.Rect, stats: PositionStats, position: string,
  sortedOrder: string[], rowCount: number, cp: SeqPalette, monomerSelectionStats: { [monomer: string]: number } = {},
  drawOptions: types.DrawOptions = {}): { [monomer: string]: DG.Rect } {
  const pr = window.devicePixelRatio;
  drawOptions.symbolStyle ??= '16px Roboto, Roboto Local, sans-serif';
  drawOptions.upperLetterHeight ??= 12.2;
  drawOptions.upperLetterAscent ??= 0.25;
  drawOptions.marginVertical ??= 2;
  drawOptions.marginHorizontal ??= 2;
  drawOptions.selectionWidth ??= 1;
  drawOptions.textHeight ??= 13;
  drawOptions.headerStyle ??= `bold ${drawOptions.textHeight * pr}px Roboto, Roboto Local, sans-serif`;

  const totalSpace = (sortedOrder.length - 1) * drawOptions.upperLetterAscent; // Total space between letters
  const barHeight = (bounds.height - 2 * drawOptions.marginVertical - totalSpace - 1.25 * drawOptions.textHeight) * pr;
  const leftShift = drawOptions.marginHorizontal * 2;
  const barWidth = (bounds.width - (leftShift + drawOptions.marginHorizontal)) * pr;
  const xStart = (bounds.x + leftShift) * pr;
  const selectionWidth = Math.min(drawOptions.selectionWidth * pr, 1);
  const xSelection = (bounds.x + 1) * pr;
  let currentY = (bounds.y + drawOptions.marginVertical) * pr;

  const monomerBounds: { [monomer: string]: DG.Rect } = {};
  for (const monomer of sortedOrder) {
    const monomerHeight = barHeight * (stats[monomer]!.count / rowCount);
    const selectionHeight = barHeight * ((monomerSelectionStats[monomer] ?? 0) / rowCount);
    const currentBound = new DG.Rect(xStart / pr, currentY / pr, barWidth / pr, monomerHeight / pr);
    monomerBounds[monomer] = currentBound;

    ctx.resetTransform();
    if (monomer !== '-' && monomer !== '') {
      const monomerTxt = monomerToShort(monomer, 5);
      const mTm: TextMetrics = ctx.measureText(monomerTxt);

      // Filling selection
      ctx.lineWidth = selectionWidth;
      ctx.line(xSelection, currentY, xSelection, currentY + selectionHeight, DG.Color.rowSelection);

      ctx.fillStyle = cp.get(monomer) ?? cp.get('other');
      ctx.textAlign = 'left';
      ctx.textBaseline = 'top';
      ctx.font = drawOptions.symbolStyle;
      // Hacks to scale uppercase characters to target rectangle
      const widthTransform = barWidth / mTm.width;
      const heightTransfrom = monomerHeight / drawOptions.upperLetterHeight;
      ctx.setTransform(widthTransform, 0, 0, heightTransfrom, xStart, currentY);
      ctx.fillText(monomerTxt, 0, 0);
    }
    currentY += monomerHeight + drawOptions.upperLetterAscent * pr;
  }

  // Drawing column header
  ctx.resetTransform();
  ctx.fillStyle = DG.Color.toHtml(DG.Color.black);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.font = drawOptions.headerStyle;
  ctx.fillText(position, (bounds.x + bounds.width / 2) * pr, (bounds.y + bounds.height - drawOptions.textHeight) * pr);

  return monomerBounds;
}
