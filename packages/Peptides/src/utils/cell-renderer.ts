import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import * as types from './types';
import * as bio from '@datagrok-libraries/bio';

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
  currentPosition: string, statsDf: DG.DataFrame, mdCol: DG.Column<number>, bound: DG.Rect, cellValue: number,
  mutationCliffsSelection: types.PositionToAARList, substitutionsInfo: types.SubstitutionsInfo,
  twoColorMode: boolean = false): void {
  const queryAAR = `${C.COLUMNS_NAMES.MONOMER} = ${currentAAR}`;
  const query = `${queryAAR} and ${C.COLUMNS_NAMES.POSITION} = ${currentPosition}`;
  const pVal: number = statsDf
    .groupBy([C.COLUMNS_NAMES.P_VALUE])
    .where(query)
    .aggregate()
    .get(C.COLUMNS_NAMES.P_VALUE, 0);

  let coef: string;
  const variant = cellValue < 0;
  if (pVal < 0.01)
    coef = variant && twoColorMode ? '#FF7900' : '#299617';
  else if (pVal < 0.05)
    coef = variant && twoColorMode ? '#FFA500' : '#32CD32';
  else if (pVal < 0.1)
    coef = variant && twoColorMode ? '#FBCEB1' : '#98FF98';
  else
    coef = DG.Color.toHtml(DG.Color.lightLightGray);


  const chooseMin = (): number => twoColorMode ? 0 : mdCol.min;
  const chooseMax = (): number => twoColorMode ? Math.max(Math.abs(mdCol.min), mdCol.max) : mdCol.max;
  const chooseCurrent = (): any => twoColorMode ? Math.abs(cellValue) : cellValue;

  const rCoef = (chooseCurrent() - chooseMin()) / (chooseMax() - chooseMin());

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
    // canvasContext.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(DG.Color.fromHtml(coef)));
    canvasContext.fillStyle = DG.Color.toHtml(DG.Color.black);
    canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
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
  currentPosition: string, invariantMapSelection: types.PositionToAARList, cellValue: number, bound: DG.Rect): void {
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


export function drawLogoInBounds(ctx: CanvasRenderingContext2D, bounds: DG.Rect, statsInfo: types.StatsInfo,
  rowCount: number, cp: bio.SeqPalette, monomerSelectionStats: {[monomer: string]: number} = {},
  drawOptions: types.DrawOptions = {}): {[monomer: string]: DG.Rect} {
  drawOptions.fontStyle ??= '16px Roboto, Roboto Local, sans-serif';
  drawOptions.upperLetterHeight ??= 12.2;
  drawOptions.upperLetterAscent ??= 0.25;
  drawOptions.marginVertical ??= 5;
  drawOptions.marginHorizontal ??= 5;

  const totalSpaceBetweenLetters = (statsInfo.orderedIndexes.length - 1) * drawOptions.upperLetterAscent;
  const barHeight = bounds.height - 2 * drawOptions.marginVertical - totalSpaceBetweenLetters;
  const leftShift = drawOptions.marginHorizontal * 2;
  const barWidth = bounds.width - leftShift * 2;
  const xStart = bounds.x + leftShift;
  const selectionWidth = 4;
  const xSelection = bounds.x + 3;
  let currentY = bounds.y + drawOptions.marginVertical;

  const monomerBounds: {[monomer: string]: DG.Rect} = {};
  for (const index of statsInfo.orderedIndexes) {
    const monomer = statsInfo.monomerCol.get(index)!;
    const monomerHeight = barHeight * (statsInfo.countCol.get(index)! / rowCount);
    const selectionHeight = barHeight * ((monomerSelectionStats[monomer] ?? 0) / rowCount);
    const currentBound = new DG.Rect(xStart, currentY, barWidth, monomerHeight);
    monomerBounds[monomer] = currentBound;

    ctx.resetTransform();
    if (monomer !== '-') {
      const monomerTxt = bio.monomerToShort(monomer, 5);
      const mTm: TextMetrics = ctx.measureText(monomerTxt);

      // Filling selection
      ctx.lineWidth = selectionWidth;
      ctx.line(xSelection, currentY, xSelection, currentY + selectionHeight, DG.Color.rowSelection);

      ctx.fillStyle = cp.get(monomer) ?? cp.get('other');
      ctx.textAlign = 'left';
      ctx.textBaseline = 'top';
      ctx.font = drawOptions.fontStyle;
      // Hacks to scale uppercase characters to target rectangle
      ctx.setTransform(barWidth / mTm.width, 0, 0, monomerHeight / drawOptions.upperLetterHeight, xStart, currentY);
      ctx.fillText(monomerTxt, 0, 0);
    }
    currentY += monomerHeight + drawOptions.upperLetterAscent;
  }

  return monomerBounds;
}
