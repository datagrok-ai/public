import * as DG from 'datagrok-api/dg';

import * as C from './constants';
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
  col.semType = C.SEM_TYPES.AMINO_ACIDS;
  col.setTag('cell.renderer', C.SEM_TYPES.AMINO_ACIDS);
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
