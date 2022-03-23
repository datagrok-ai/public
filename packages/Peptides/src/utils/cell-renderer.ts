import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

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
) {
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
export function setAARRenderer(col: DG.Column, grid: DG.Grid | null = null, grouping = false) {
  col.semType = 'aminoAcids';
  col.setTag('cell.renderer', 'aminoAcids');
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
) {
  g.textAlign = 'start';
  let colorPart = pivot == -1 ? s.substring(0) : s.substring(0, pivot);
  if (colorPart.length == 1)
    colorPart = colorPart.toUpperCase();

  if (colorPart.length >= 3) {
    if (colorPart.substring(0, 3) in ChemPalette.AAFullNames)
      colorPart = ChemPalette.AAFullNames[s.substring(0, 3)] + colorPart.substr(3);
    else if (colorPart.substring(1, 4) in ChemPalette.AAFullNames)
      colorPart = colorPart[0] + ChemPalette.AAFullNames[s.substring(1, 4)] + colorPart.substr(4);
  }
  let grayPart = pivot == -1 ? '' : s.substr(pivot);
  if (hideMod) {
    let end = colorPart.lastIndexOf(')');
    let beg = colorPart.indexOf('(');
    if (beg > -1 && end > -1 && end - beg > 2)
      colorPart = colorPart.substr(0, beg) + '(+)' + colorPart.substr(end + 1);


    end = grayPart.lastIndexOf(')');
    beg = grayPart.indexOf('(');
    if (beg > -1 && end > -1 && end - beg > 2)
      grayPart = grayPart.substr(0, beg) + '(+)' + grayPart.substr(end + 1);
  }
  const textSize = g.measureText(colorPart + grayPart);
  const indent = 5;

  const colorTextSize = g.measureText(colorPart);

  function draw(dx1: number, dx2: number) {
    g.fillStyle = color;
    g.globalAlpha = transparencyRate;
    g.fillText(
      colorPart,
      x + dx1,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    g.fillStyle = ChemPalette.undefinedColor;
    g.fillText(
      grayPart,
      x + dx2,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
  }


  if (left || textSize.width > w) {
    draw(indent, indent + colorTextSize.width);
    return x + colorTextSize.width + g.measureText(grayPart).width;
  } else {
    draw((w - textSize.width) / 2, (w - textSize.width) / 2 + colorTextSize.width);
    return x + (w - textSize.width) / 2 + colorTextSize.width;
  }
}

/**
 * Amino acid residue cell renderer.
 *
 * @export
 * @class AminoAcidsCellRenderer
 * @extends {DG.GridCellRenderer}
 */
export class AminoAcidsCellRenderer extends DG.GridCellRenderer {
  chemPalette: ChemPalette | null;

  /**
     * Renderer name.
     *
     * @readonly
     * @memberof AminoAcidsCellRenderer
     */
  get name() {
    return 'aminoAcidsCR';
  }

  /**
     * Cell type.
     *
     * @readonly
     * @memberof AminoAcidsCellRenderer
     */
  get cellType() {
    return 'aminoAcids';
  }

  /**
     * Cell height.
     *
     * @readonly
     * @memberof AminoAcidsCellRenderer
     */
  get defaultHeight() {
    return 15;
  }

  /**
     * Cell width.
     *
     * @readonly
     * @memberof AminoAcidsCellRenderer
     */
  get defaultWidth() {
    return 30;
  }

  /**
     * Constructor.
     */
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
    cellStyle: DG.GridCellStyle) {
    // this.chemPalette ??= new ChemPalette('grok', gridCell.tableColumn?.getTag('groups') ? true : false);

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

/**
 * Aligned sequence cell renderer.
 *
 * @export
 * @class AlignedSequenceCellRenderer
 * @extends {DG.GridCellRenderer}
 */
export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
  /**
   * Renderer name.
   *
   * @readonly
   * @memberof AlignedSequenceCellRenderer
   */
  get name() {
    return 'alignedSequenceCR';
  }

  /**
   * Cell type.
   *
   * @readonly
   * @memberof AlignedSequenceCellRenderer
   */
  get cellType() {
    return 'alignedSequence';
  }

  /**
   * Cell height.
   *
   * @readonly
   * @memberof AlignedSequenceCellRenderer
   */
  get defaultHeight() {
    return 30;
  }

  /**
   * Cell width.
   *
   * @readonly
   * @memberof AlignedSequenceCellRenderer
   */
  get defaultWidth() {
    return 230;
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
   * @memberof AlignedSequenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle,
  ) {
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
    const subParts = s.split('-');
    const [text, simplified] = processSequence(subParts);
    const textSize = g.measureText(text.join(''));
    let x_1 = Math.max(x, x + (w - textSize.width) / 2);

    subParts.forEach((amino, index) => {
      let [color, outerAmino,, pivot] = ChemPalette.getColorAAPivot(amino);
      g.fillStyle = ChemPalette.undefinedColor;
      if (index + 1 < subParts.length) {
        const gap = simplified ? '' : ' ';
        outerAmino += `${outerAmino ? '' : '-'}${gap}`;
      }
      x_1 = printLeftOrCentered(x_1, y, w, h, g, outerAmino, color, pivot, true);
    });

    g.restore();
  }
}

/**
 * Function for sequence processing.
 *
 * @export
 * @param {string[]} subParts Sequence subparts.
 * @return {[string[], boolean]}
 */
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

/**
 * Aligned sequence difference cell renderer.
 *
 * @export
 * @class AlignedSequenceDifferenceCellRenderer
 * @extends {DG.GridCellRenderer}
 */
export class AlignedSequenceDifferenceCellRenderer extends DG.GridCellRenderer {
  /**
   * Renderer name.
   *
   * @readonly
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  get name() {
    return 'alignedSequenceDifferenceCR';
  }

  /**
   * Cell type.
   *
   * @readonly
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  get cellType() {
    return 'alignedSequenceDifference';
  }

  /**
   * Cell height.
   *
   * @readonly
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  get defaultHeight() {
    return 30;
  }

  /**
   * Cell width.
   *
   * @readonly
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  get defaultWidth() {
    return 230;
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
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle,
  ) {
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
    const subParts1 = s1.split('-');
    const subParts2 = s2.split('-');
    const [text] = processSequence(subParts1);
    const textSize = g.measureText(text.join(''));
    x = Math.max(x, x + (w - textSize.width) / 2);

    subParts1.forEach((amino1: string, index) => {
      let amino2 = subParts2[index];
      const [color1, amino1Outer, amino1Inner, pivot1] = ChemPalette.getColorAAPivot(amino1);
      const [color2, amino2Outer, amino2Inner, pivot2] = ChemPalette.getColorAAPivot(amino2);

      if (amino1 != amino2) {
        const verticalShift = 7;

        amino1 = amino1Outer + (amino1Inner !== '' ? '(' + amino1Inner + ')' : '');
        amino2 = amino2Outer + (amino2Inner !== '' ? '(' + amino2Inner + ')' : '');
        amino1 = amino1 === '' ? '-' : amino1;
        amino2 = amino2 === '' ? '-' : amino2;

        const x1 = printLeftOrCentered(x, y - verticalShift, w, h, g, amino1, color1, pivot1, true);
        x = printLeftOrCentered(x, y + verticalShift, w, h, g, amino2, color2, pivot2, true);
        x = Math.max(x, x1) + 4;
      } else
        x = printLeftOrCentered(x, y, w, h, g, amino1 ? amino1 : '-', color1, pivot1, true, true, 0.5) + 4;
    });
    g.restore();
  }
}
