import {SeqPaletteBase} from '@datagrok-libraries/bio/src/seq-palettes';
import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import {getPalleteByType} from './misc';
import * as types from './types';

/**
 * A function to expand column size based on its contents.
 *
 * @export
 * @param {DG.Column} col Column to expand.
 * @param {DG.Grid} grid Grid containing colum for expansion.
 * @param {(cellVal: string) => number} cellRenderSize An anonymous function that calculates cell value length.
 * @param {number} [textSizeMult=10] Text size muliplier.
 * @param {number} [minSize=30] Minimal column width.
 * @param {number} [maxSize=650] Maximum column width.
 * @param {number} [timeout=500] Timeout value.
 */
export function expandColumn(col: DG.Column, grid: DG.Grid, cellRenderSize: (cellVal: string) => number,
  textSizeMult = 10, minSize = 30, maxSize = 650, timeout = 500): void {
  let maxLen = 0;
  col.categories.forEach((ent: string) => {
    const len = cellRenderSize(ent);
    if (len > maxLen)
      maxLen = len;
  });
  setTimeout(() => {
      grid.columns.byName(col.name)!.width = Math.min(Math.max(maxLen * textSizeMult, minSize), maxSize);
  },
  timeout);
}

/**
 * A function that sets amino acid residue to the specified column.
 *
 * @export
 * @param {DG.Column} col Column to set renderer for.
 * @param {(DG.Grid | null)} [grid=null] Grid that contains the col column.
 * @param {boolean} [grouping=false] Is grouping enabled.
 */
export function setAARRenderer(col: DG.Column, alphabet: string, grid?: DG.Grid): void {
  col.semType = C.SEM_TYPES.MONOMER;
  col.setTag('cell.renderer', C.SEM_TYPES.MONOMER);
  col.tags[C.TAGS.ALPHABET] = alphabet;

  if (grid)
    expandColumn(col, grid, (ent) => measureAAR(ent));
}

/**
 * A function to measure amino acid residue
 *
 * @export
 * @param {string} s Amino acid residue string.
 * @return {number} Amino acid residue size.
 */
export function measureAAR(s: string): number {
  const end = s.lastIndexOf(')');
  const beg = s.indexOf('(');
  return end == beg ? s.length : s.length - (end - beg) + 1;
}

export function renderSARCell(canvasContext: CanvasRenderingContext2D, currentAAR: string, currentPosition: string,
  statsDf: DG.DataFrame, twoColorMode: boolean, mdCol: DG.Column<number>, bound: DG.Rect, cellValue: number,
  currentSelection: types.SelectionObject, substitutionsInfo: types.SubstitutionsInfo | null): void {
  const queryAAR = `${C.COLUMNS_NAMES.AMINO_ACID_RESIDUE} = ${currentAAR}`;
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
  if (substitutionsInfo) {
    canvasContext.textBaseline = 'middle';
    canvasContext.textAlign = 'center';
    canvasContext.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(DG.Color.fromHtml(coef)));
    canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
    let substValue = 0;
    substitutionsInfo.get(currentAAR)?.get(currentPosition)?.forEach((idxs) => substValue += idxs.length);
    if (substValue && substValue != 0)
      canvasContext.fillText(substValue.toString(), midX, midY);
  }

  //TODO: frame based on currentSelection
  const aarSelection = currentSelection[currentPosition];
  if (aarSelection && aarSelection.includes(currentAAR)) {
    canvasContext.strokeStyle = '#000';
    canvasContext.lineWidth = 1;
    canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
  }
}

export function renderBarchart(ctx: CanvasRenderingContext2D, col: DG.Column, monomerColStats: types.MonomerColStats,
  bounds: DG.Rect, max: number): types.BarCoordinates {
  let sum = col.length - (monomerColStats['-']?.count ?? 0);
  const colorPalette = getPalleteByType(col.tags[C.TAGS.ALPHABET]);
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
      if (color != SeqPaletteBase.undefinedColor)
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
