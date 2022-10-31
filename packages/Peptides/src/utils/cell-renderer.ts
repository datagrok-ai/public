import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import * as C from './constants';
import * as types from './types';

function renderCellSelection(canvasContext: CanvasRenderingContext2D, bound: DG.Rect): void {
  canvasContext.strokeStyle = '#000';
  canvasContext.lineWidth = 1;
  canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
}

/** A function that sets amino acid residue cell renderer to the specified column */
export function setAARRenderer(col: DG.Column, alphabet: string, grid: DG.Grid, timeout: number = 500): void {
  col.semType = C.SEM_TYPES.MONOMER;
  col.setTag('cell.renderer', C.SEM_TYPES.MONOMER);
  col.tags[C.TAGS.ALPHABET] = alphabet;
  // setTimeout(() => grid.columns.byName(col.name)!.width = 60, timeout);
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
    canvasContext.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(DG.Color.fromHtml(coef)));
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

export function renderLogoSummaryCell(canvasContext: CanvasRenderingContext2D, cellValue: number,
  clusterSelection: number[], bound: DG.Rect): void {
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = '#000';
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  if (clusterSelection.includes(cellValue))
    renderCellSelection(canvasContext, bound);
}

export function renderBarchart(ctx: CanvasRenderingContext2D, col: DG.Column, monomerColStats: types.MonomerColStats,
  bounds: DG.Rect, max: number): types.BarCoordinates {
  let sum = col.length - (monomerColStats['-']?.count ?? 0);
  const colorPalette = bio.getPaletteByType(col.tags[C.TAGS.ALPHABET]);
  const name = col.name;
  const colNameSize = ctx.measureText(name);
  const margin = 0.2;
  const innerMargin = 0.02;
  const selectLineRatio = 0.1;
  const fontSize = 11;

  const xMargin = bounds.x + bounds.width * margin;
  const yMargin = bounds.y + bounds.height * margin / 4;
  const wMargin = bounds.width - bounds.width * margin * 2;
  const hMargin = bounds.height - bounds.height * margin;
  const barWidth = 10;
  ctx.fillStyle = 'black';
  ctx.textBaseline = 'top';
  ctx.font = `${hMargin * margin / 2}px`;
  ctx.fillText(name, xMargin + (wMargin - colNameSize.width) / 2, yMargin + hMargin + hMargin * margin / 4);


  const barCoordinates: types.BarCoordinates = {};

  const xStart = xMargin + (wMargin - barWidth) / 2;
  for (const [monomer, monomerStats] of Object.entries(monomerColStats)) {
    if (monomer == '-')
      continue;

    const count = monomerStats.count;
    const sBarHeight = hMargin * count / max;
    const gapSize = sBarHeight * innerMargin;
    const verticalShift = (max - sum) / max;
    const textSize = ctx.measureText(monomer);
    const subBarHeight = sBarHeight - gapSize;
    const yStart = yMargin + hMargin * verticalShift + gapSize / 2;
    barCoordinates[monomer] = new DG.Rect(xStart, yStart, barWidth, subBarHeight);

    const color = colorPalette.get(monomer);
    ctx.strokeStyle = color;
    ctx.fillStyle = color;

    if (textSize.width <= subBarHeight) {
      if (color != bio.SeqPaletteBase.undefinedColor)
        ctx.fillRect(xStart, yStart, barWidth, subBarHeight);
      else {
        ctx.strokeRect(xStart + 0.5, yStart, barWidth - 1, subBarHeight);
        barCoordinates[monomer].x -= 0.5;
        barCoordinates[monomer].width -= 1;
      }

      const leftMargin = (wMargin - (monomer.length > 1 ? fontSize : textSize.width - 8)) / 2;
      const absX = xMargin + leftMargin;
      const absY = yStart + subBarHeight / 2 + (monomer.length == 1 ? 4 : 0);
      const origTransform = ctx.getTransform();

      if (monomer.length > 1) {
        ctx.translate(absX, absY);
        ctx.rotate(Math.PI / 2);
        ctx.translate(-absX, -absY);
      }

      ctx.fillStyle = 'black';
      ctx.font = `${fontSize}px monospace`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText(monomer, absX, absY);
      ctx.setTransform(origTransform);
    } else
      ctx.fillRect(xStart, yStart, barWidth, subBarHeight);

    const selectedCount = monomerStats.selected;
    if (selectedCount) {
      ctx.fillStyle = 'rgb(255,165,0)';
      ctx.fillRect(xStart - wMargin * selectLineRatio * 2, yStart,
        barWidth * selectLineRatio, hMargin * selectedCount / max - gapSize);
    }

    sum -= count;
  }

  return barCoordinates;
}
