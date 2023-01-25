import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import * as types from './types';
import {PositionStats, SummaryStats, MonomerPositionStats} from '../model';
import {monomerToShort, SeqPalette} from '@datagrok-libraries/bio';

function renderCellSelection(canvasContext: CanvasRenderingContext2D, bound: DG.Rect): void {
  canvasContext.strokeStyle = '#000';
  canvasContext.lineWidth = 1;
  canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
}

/** A function that sets amino acid residue cell renderer to the specified column */
export function setAARRenderer(col: DG.Column, alphabet: string): void {
  col.semType = C.SEM_TYPES.MONOMER;
  col.setTag('cell.renderer', C.SEM_TYPES.MONOMER);
  col.tags[C.TAGS.ALPHABET] = alphabet;
}

export function renderMutationCliffCell(canvasContext: CanvasRenderingContext2D, currentAAR: string,
  currentPosition: string, monomerPositionStats: MonomerPositionStats, bound: DG.Rect,
  mutationCliffsSelection: types.PositionToAARList, substitutionsInfo: types.SubstitutionsInfo,
  twoColorMode: boolean = false): void {
  const positionStats = monomerPositionStats[currentPosition];
  const pVal: number = positionStats[currentAAR].pValue;
  const currentMeanDiff = positionStats[currentAAR].meanDifference;

  let coef: string;
  const isMeanDeltaNegative = currentMeanDiff < 0;
  if (pVal < 0.01)
    coef = isMeanDeltaNegative && twoColorMode ? '#FF7900' : '#299617';
  else if (pVal < 0.05)
    coef = isMeanDeltaNegative && twoColorMode ? '#FFA500' : '#32CD32';
  else if (pVal < 0.1)
    coef = isMeanDeltaNegative && twoColorMode ? '#FBCEB1' : '#98FF98';
  else
    coef = DG.Color.toHtml(DG.Color.lightLightGray);


  const minMeanDifference = twoColorMode ? 0 : monomerPositionStats.general.minMeanDifference;
  const maxMeanDifference = twoColorMode ?
    Math.max(Math.abs(monomerPositionStats.general.minMeanDifference), monomerPositionStats.general.maxMeanDifference) :
    monomerPositionStats.general.maxMeanDifference;
  const currentMeanDifference = twoColorMode ? Math.abs(currentMeanDiff) : currentMeanDiff;

  const rCoef = (currentMeanDifference - minMeanDifference) / (maxMeanDifference - minMeanDifference);

  const maxRadius = 0.9 * (bound.width > bound.height ? bound.height : bound.width) / 2;
  const radius = Math.floor(maxRadius * rCoef);

  const midX = bound.x + bound.width / 2;
  const midY = bound.y + bound.height / 2;
  canvasContext.beginPath();
  canvasContext.fillStyle = coef;
  canvasContext.arc(midX, midY, radius < 3 ? 3 : radius, 0, Math.PI * 2, true);
  canvasContext.closePath();

  canvasContext.fill();
  if (substitutionsInfo.size > 0) {
    canvasContext.textBaseline = 'middle';
    canvasContext.textAlign = 'center';
    canvasContext.fillStyle = DG.Color.toHtml(DG.Color.black);
    canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
    canvasContext.shadowBlur = 5;
    canvasContext.shadowColor = DG.Color.toHtml(DG.Color.white);
    let substValue = 0;
    substitutionsInfo.get(currentAAR)?.get(currentPosition)?.forEach((idxs) => substValue += idxs.length);
    if (substValue && substValue != 0)
      canvasContext.fillText(substValue.toString(), midX, midY);
  }

  const aarSelection = mutationCliffsSelection[currentPosition];
  if (aarSelection && aarSelection.includes(currentAAR))
    renderCellSelection(canvasContext, bound);
}

export function renderInvaraintMapCell(canvasContext: CanvasRenderingContext2D, currentAAR: string,
  currentPosition: string, invariantMapSelection: types.PositionToAARList, cellValue: number, bound: DG.Rect,
  color: number): void {
  canvasContext.fillStyle = DG.Color.toHtml(color);
  canvasContext.fillRect(bound.x, bound.y, bound.width, bound.height);
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = '#000';
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  const aarSelection = invariantMapSelection[currentPosition];
  if (aarSelection && aarSelection.includes(currentAAR))
    renderCellSelection(canvasContext, bound);
}

export function renderLogoSummaryCell(canvasContext: CanvasRenderingContext2D, cellValue: string, cellRawData: number,
  clusterSelection: number[], bound: DG.Rect): void {
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = '#000';
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  if (clusterSelection.includes(cellRawData))
    renderCellSelection(canvasContext, bound);
}


export function drawLogoInBounds(ctx: CanvasRenderingContext2D, bounds: DG.Rect, stats: PositionStats,
  sortedOrder: string[], rowCount: number, cp: SeqPalette, monomerSelectionStats: { [monomer: string]: number } = {},
  drawOptions: types.DrawOptions = {}): { [monomer: string]: DG.Rect } {
  drawOptions.fontStyle ??= '16px Roboto, Roboto Local, sans-serif';
  drawOptions.upperLetterHeight ??= 12.2;
  drawOptions.upperLetterAscent ??= 0.25;
  drawOptions.marginVertical ??= 5;
  drawOptions.marginHorizontal ??= 5;

  const pr = window.devicePixelRatio;
  const totalSpaceBetweenLetters = (sortedOrder.length - 1) * drawOptions.upperLetterAscent;
  const barHeight = (bounds.height - 2 * drawOptions.marginVertical - totalSpaceBetweenLetters) * pr;
  const leftShift = drawOptions.marginHorizontal * 2;
  const barWidth = (bounds.width - leftShift * 2) * pr;
  const xStart = (bounds.x + leftShift) * pr;
  const selectionWidth = 4 * pr;
  const xSelection = (bounds.x + 3) * pr;
  let currentY = (bounds.y + drawOptions.marginVertical) * pr;

  const monomerBounds: { [monomer: string]: DG.Rect } = {};
  for (const monomer of sortedOrder) {
    const monomerHeight = barHeight * (stats[monomer].count / rowCount);
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
      ctx.font = drawOptions.fontStyle;
      // Hacks to scale uppercase characters to target rectangle
      const widthTransform = barWidth / mTm.width;
      const heightTransfrom = monomerHeight / drawOptions.upperLetterHeight;
      ctx.setTransform(widthTransform, 0, 0, heightTransfrom, xStart, currentY);
      ctx.fillText(monomerTxt, 0, 0);
    }
    currentY += monomerHeight + drawOptions.upperLetterAscent * pr;
  }

  return monomerBounds;
}
