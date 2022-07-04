import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import * as types from './types';
import {getSeparator} from './misc';

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
export function expandColumn(
  col: DG.Column, grid: DG.Grid, cellRenderSize: (cellVal: string) => number,
  textSizeMult = 10, minSize = 30, maxSize = 650, timeout = 500,
): void {
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
export function setAARRenderer(col: DG.Column, grid: DG.Grid | null = null, grouping = false): void {
  col.semType = C.SEM_TYPES.AMINO_ACIDS;
  col.setTag('cell.renderer', C.SEM_TYPES.AMINO_ACIDS);
  if (grouping)
    col.setTag('groups', `${grouping}`);

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

/**
 * A function that prints a string aligned to left or centered.
 *
 * @param {number} x x coordinate.
 * @param {number} y y coordinate.
 * @param {number} w Width.
 * @param {number} h Height.
 * @param {CanvasRenderingContext2D} g Canvas rendering context.
 * @param {string} s String to print.
 * @param {string} [color=ChemPalette.undefinedColor] String color.
 * @param {number} [pivot=0] Pirvot.
 * @param {boolean} [left=false] Is left aligned.
 * @param {boolean} [hideMod=false] Hide amino acid redidue modifications.
 * @param {number} [transparencyRate=0.0] Transparency rate where 1.0 is fully transparent
 * @return {number} x coordinate to start printing at.
 */
function printLeftOrCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color = ChemPalette.undefinedColor,
  pivot: number = 0, left = false, hideMod = false, transparencyRate: number = 1.0,
): number {
  g.textAlign = 'start';
  let colorPart = pivot == -1 ? s.substring(0) : s.substring(0, pivot);
  if (colorPart.length == 1)
    colorPart = colorPart.toUpperCase();

  if (colorPart.length >= 3) {
    if (colorPart.substring(0, 3) in ChemPalette.AAFullNames)
      colorPart = ChemPalette.AAFullNames[s.substring(0, 3)] + colorPart.substring(3);
    else if (colorPart.substring(1, 4) in ChemPalette.AAFullNames)
      colorPart = colorPart[0] + ChemPalette.AAFullNames[s.substring(1, 4)] + colorPart.substring(4);
  }
  let grayPart = pivot == -1 ? '' : s.substring(pivot);
  if (hideMod) {
    let end = colorPart.lastIndexOf(')');
    let beg = colorPart.indexOf('(');
    if (beg > -1 && end > -1 && end - beg > 2)
      colorPart = colorPart.substring(0, beg) + '(+)' + colorPart.substring(end + 1);


    end = grayPart.lastIndexOf(')');
    beg = grayPart.indexOf('(');
    if (beg > -1 && end > -1 && end - beg > 2)
      grayPart = grayPart.substring(0, beg) + '(+)' + grayPart.substring(end + 1);
  }
  const textSize = g.measureText(colorPart + grayPart);
  const indent = 5;

  const colorTextSize = g.measureText(colorPart);
  const dy = (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2;

  function draw(dx1: number, dx2: number): void {
    g.fillStyle = color;
    g.globalAlpha = transparencyRate;
    g.fillText(colorPart, x + dx1, y + dy);
    g.fillStyle = ChemPalette.undefinedColor;
    g.fillText(grayPart, x + dx2, y + dy);
  }


  if (left || textSize.width > w) {
    draw(indent, indent + colorTextSize.width);
    return x + colorTextSize.width + g.measureText(grayPart).width;
  } else {
    const dx = (w - textSize.width) / 2;
    draw(dx, dx + colorTextSize.width);
    return x + dx + colorTextSize.width;
  }
}

export class AminoAcidsCellRenderer extends DG.GridCellRenderer {
  chemPalette: ChemPalette | null;

  get name(): string {return 'aminoAcidsCR';}

  get cellType(): string {return C.SEM_TYPES.AMINO_ACIDS;}

  get defaultHeight(): number {return 15;}

  get defaultWidth(): number {return 30;}

  constructor() {
    super();
    this.chemPalette = null;
  }

  /**
     * Cell renderer function.
     *
     * @param {CanvasRenderingContext2D} g Canvas rendering context.
     * @param {number} x x coordinate on the canvas.
     * @param {number} y y coordinate on the canvas.
     * @param {number} w width of the cell.
     * @param {number} h height of the cell.
     * @param {DG.GridCell} gridCell Grid cell.
     * @param {DG.GridCellStyle} cellStyle Cell style.
     */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle): void {
    y -= 2;
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.font = `12px monospace`;
    g.textBaseline = 'top';
    const s: string = gridCell.cell.value ? gridCell.cell.value : '-';
    let [color, outerS, innerS, pivot] = ChemPalette.getColorAAPivot(s);
    if (innerS)
      outerS = s;

    printLeftOrCentered(x, y, w, h, g, outerS, color, pivot, false, true);
    g.restore();
  }
}

export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
  get name(): string {return 'alignedSequenceCR';}

  get cellType(): string {return C.SEM_TYPES.ALIGNED_SEQUENCE;}

  get defaultHeight(): number {return 30;}

  get defaultWidth(): number {return 230;}

  /**
   * Cell renderer function.
   *
   * @param {CanvasRenderingContext2D} g Canvas rendering context.
   * @param {number} x x coordinate on the canvas.
   * @param {number} y y coordinate on the canvas.
   * @param {number} w width of the cell.
   * @param {number} h height of the cell.
   * @param {DG.GridCell} gridCell Grid cell.
   * @param {DG.GridCellStyle} cellStyle Cell style.
   * @memberof AlignedSequenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle,
  ): void {
    const grid = gridCell.grid;
    const cell = gridCell.cell;
    w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.font = '12px monospace';
    g.textBaseline = 'top';
    const s: string = cell.value ?? '';

    //TODO: can this be replaced/merged with splitSequence?
    const subParts = s.split(getSeparator(cell.column));
    const [text, simplified] = processSequence(subParts);
    const textSize = g.measureText(text.join(''));
    let x1 = Math.max(x, x + (w - textSize.width) / 2);

    subParts.forEach((amino, index) => {
      let [color, outerAmino,, pivot] = ChemPalette.getColorAAPivot(amino);
      g.fillStyle = ChemPalette.undefinedColor;
      if (index + 1 < subParts.length) {
        const gap = simplified ? '' : ' ';
        outerAmino += `${outerAmino ? '' : '-'}${gap}`;
      }
      x1 = printLeftOrCentered(x1, y, w, h, g, outerAmino, color, pivot, true);
    });

    g.restore();
  }
}

export function processSequence(subParts: string[]): [string[], boolean] {
  const simplified = !subParts.some((amino, index) =>
    amino.length > 1 &&
    index != 0 &&
    index != subParts.length - 1);

  const text: string[] = [];
  const gap = simplified ? '' : ' ';
  subParts.forEach((amino: string, index) => {
    if (index < subParts.length)
      amino += `${amino ? '' : '-'}${gap}`;

    text.push(amino);
  });
  return [text, simplified];
}

export class AlignedSequenceDifferenceCellRenderer extends DG.GridCellRenderer {
  get name(): string {return 'alignedSequenceDifferenceCR';}

  get cellType(): string {return C.SEM_TYPES.ALIGNED_SEQUENCE_DIFFERENCE;}

  get defaultHeight(): number {return 30;}

  get defaultWidth(): number {return 230;}

  /**
   * Cell renderer function.
   *
   * @param {CanvasRenderingContext2D} g Canvas rendering context.
   * @param {number} x x coordinate on the canvas.
   * @param {number} y y coordinate on the canvas.
   * @param {number} w width of the cell.
   * @param {number} h height of the cell.
   * @param {DG.GridCell} gridCell Grid cell.
   * @param {DG.GridCellStyle} cellStyle Cell style.
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle): void {
    const grid = gridCell.grid;
    const cell = gridCell.cell;

    w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.font = '12px monospace';
    g.textBaseline = 'top';
    const s: string = cell.value ?? '';

    //TODO: can this be replaced/merged with splitSequence?
    const [s1, s2] = s.split('#');
    const separator = getSeparator(gridCell.tableColumn!);
    const subParts1 = s1.split(separator);
    const subParts2 = s2.split(separator);
    const [text] = processSequence(subParts1);
    const textSize = g.measureText(text.join(''));
    let updatedX = Math.max(x, x + (w - (textSize.width + subParts1.length * 4)) / 2);
    // 28 is the height of the two substitutions on top of each other + space
    const updatedY = Math.max(y, y + (h - 28) / 2);

    let amino2;
    let updatedAmino1: string;
    let updatedAmino2: string;
    subParts1.forEach((amino1: string, index) => {
      amino2 = subParts2[index];
      const [color1, amino1Outer, amino1Inner, pivot1] = ChemPalette.getColorAAPivot(amino1);
      const [color2, amino2Outer, amino2Inner, pivot2] = ChemPalette.getColorAAPivot(amino2);

      updatedAmino1 = amino1Outer + (amino1Inner !== '' ? '(' + amino1Inner + ')' : '');
      updatedAmino1 = updatedAmino1 === '' ? '-' : updatedAmino1;

      if (amino1 != amino2) {
        updatedAmino2 = amino2Outer + (amino2Inner !== '' ? '(' + amino2Inner + ')' : '');
        updatedAmino2 = updatedAmino2 === '' ? '-' : updatedAmino2;

        const vShift = 7;
        const subX0 = printLeftOrCentered(updatedX, updatedY - vShift, w, h, g, updatedAmino1, color1, pivot1, true);
        const subX1 = printLeftOrCentered(updatedX, updatedY + vShift, w, h, g, updatedAmino2, color2, pivot2, true);
        updatedX = Math.max(subX1, subX0);
      } else
        updatedX = printLeftOrCentered(updatedX, updatedY, w, h, g, updatedAmino1, color1, pivot1, true, true, 0.5);
      updatedX += 4;
    });
    g.restore();
  }
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
