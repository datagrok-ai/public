import * as C from "./constants";
import {getSeparator} from "./misc";
import {ChemPalette} from "./chem-palette";
import * as DG from 'datagrok-api/dg';
import {AminoacidsPalettes} from "@datagrok-libraries/bio/src/aminoacids";
import {NucleotidesPalettes} from "@datagrok-libraries/bio/src/nucleotides";
import {UnknownSeqPalettes} from "@datagrok-libraries/bio/src/unknown";
import {SplitterFunc, WebLogo} from "@datagrok-libraries/bio/src/viewers/web-logo";
import {SeqPalette} from "@datagrok-libraries/bio/src/seq-palettes";
import {HelmSequenceCellRenderer} from "../../../Helm/src/cell-renderer";

function getPalleteByType(paletteType: string): SeqPalette  {
  switch (paletteType) {
    case 'PT':
      return  AminoacidsPalettes.GrokGroups;
    case 'NT':
      return  NucleotidesPalettes.Chromatogram
      // other
    default:
      return UnknownSeqPalettes.Color;
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
  let colorPart = s.substring(0);
  const textSize = g.measureText(colorPart);
  const indent = 5;

  const colorTextSize = g.measureText(colorPart);
  const dy = (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2;

  function draw(dx1: number, dx2: number): void {
    g.fillStyle = color;
    g.globalAlpha = transparencyRate;
    g.fillText(colorPart, x + dx1, y + dy);
  }


  if (left || textSize.width > w) {
    draw(indent, indent + colorTextSize.width);
    return x + colorTextSize.width;
  } else {
    const dx = (w - textSize.width) / 2;
    draw(dx, dx + colorTextSize.width);
    return x + dx + colorTextSize.width;
  }
}
export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
  constructor() {
    super();
  }

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
    const [type, subtype, paletteType] =  gridCell.cell.column.getTag(DG.TAGS.UNITS).split(":");
    w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.font = '12px monospace';
    g.textBaseline = 'top';
    const s: string = cell.value ?? '';

    //TODO: can this be replaced/merged with splitSequence?
    const units = gridCell.cell.column.getTag(DG.TAGS.UNITS);

    const palette = getPalleteByType(paletteType);

    const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, gridCell.cell.column.getTag('separator') );// splitter,

    const subParts:string[] =  splitterFunc(cell.value);
    console.log(subParts);

    const textSize = g.measureText(subParts.join(''));
    let x1 = Math.max(x, x + (w - textSize.width) / 2);

    subParts.forEach((amino, index) => {
      let [color, outerAmino,, pivot] = ChemPalette.getColorAAPivot(amino);
      color = palette.get(amino);
      g.fillStyle = ChemPalette.undefinedColor;
      x1 = printLeftOrCentered(x1, y, w, h, g, amino, color, pivot, true);
    });

    g.restore();
  }
}
