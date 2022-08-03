import * as C from './constants';
import * as DG from 'datagrok-api/dg';
import {AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {UnknownSeqPalette, UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {SplitterFunc, WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import * as ui from 'datagrok-api/ui';

export const lru = new DG.LruCache<any, any>();
const undefinedColor = 'rgb(100,100,100)';
const grayColor = '#808080';

function getPalleteByType(paletteType: string): SeqPalette {
  switch (paletteType) {
  case 'PT':
    return AminoacidsPalettes.GrokGroups;
  case 'NT':
    return NucleotidesPalettes.Chromatogram;
  case 'DNA':
    return NucleotidesPalettes.Chromatogram;
  case 'RNA':
    return NucleotidesPalettes.Chromatogram;
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
 * @param {string} [color=undefinedColor] String color.
 * @param {number} [pivot=0] Pirvot.
 * @param {boolean} [left=false] Is left aligned.
 * @param {number} [transparencyRate=0.0] Transparency rate where 1.0 is fully transparent
 * @param {string} [separator=''] Is separator for sequence.
 * @param {boolean} [last=false] Is checker if element last or not.
 * @return {number} x coordinate to start printing at.
 */
function printLeftOrCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color = undefinedColor,
  pivot: number = 0, left = false, transparencyRate: number = 1.0,
  separator: string = '', last: boolean = false, drawStyle: string = 'classic', maxWord: string = ''): number {
  g.textAlign = 'start';
  const colorPart = s.substring(0);
  let grayPart = last ? '' : separator;

  let textSize = g.measureText(colorPart + grayPart);
  const indent = 5;

  let colorTextSize = g.measureText(colorPart);
  if (drawStyle === 'msa') {
    colorTextSize = g.measureText(maxWord);
    textSize = g.measureText(maxWord + grayPart);
  }
  const dy = (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2;

  function draw(dx1: number, dx2: number): void {
    g.fillStyle = color;
    g.globalAlpha = transparencyRate;
    g.fillText(colorPart, x + dx1, y + dy);
    g.fillStyle = grayColor;
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

function findMonomers(helmString: string) {
  //@ts-ignore
  const types = Object.keys(org.helm.webeditor.monomerTypeList());
  const monomers: any = [];
  const monomer_names: any = [];
  for (var i = 0; i < types.length; i++) {
    //@ts-ignore
    monomers.push(new scil.helm.Monomers.getMonomerSet(types[i]));
    Object.keys(monomers[i]).forEach(k => {
      monomer_names.push(monomers[i][k].id);
    });
  }
  const split_string = WebLogo.splitterAsHelm(helmString);
  return new Set(split_string.filter(val => !monomer_names.includes(val)));
}

export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'macromoleculeSequence'; }

  get cellType(): string { return C.SEM_TYPES.MACROMOLECULE; }

  get defaultHeight(): number { return 30; }

  get defaultWidth(): number { return 230; }

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
    cellStyle: DG.GridCellStyle
  ): void {
    const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
    const cell = gridCell.cell;
    const tag = gridCell.cell.column.getTag(DG.TAGS.UNITS);
    if (tag === 'HELM') {
      const monomers = findMonomers(cell.value);
      if (monomers.size == 0) {
        const host = ui.div([], {style: {width: `${w}px`, height: `${h}px`}});
        host.setAttribute('dataformat', 'helm');
        host.setAttribute('data', gridCell.cell.value);
        gridCell.element = host;
        //@ts-ignore
        const canvas = new JSDraw2.Editor(host, {width: w, height: h, skin: 'w8', viewonly: true});
        const formula = canvas.getFormula(true);
        if (!formula) {
          gridCell.element = ui.divText(gridCell.cell.value, {style: {color: 'red'}});
        }
        const molWeight = Math.round(canvas.getMolWeight() * 100) / 100;
        const coef = Math.round(canvas.getExtinctionCoefficient(true) * 100) / 100;
        const molfile = canvas.getMolfile();
        const result = formula + ', ' + molWeight + ', ' + coef + ', ' + molfile;
        lru.set(gridCell.cell.value, result);
        return;
      }
      if (monomers.size > 0) {
        w = grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
        g.save();
        g.beginPath();
        g.rect(x, y, w, h);
        g.clip();
        g.font = '12px monospace';
        g.textBaseline = 'top';
        let x1 = x;
        const s: string = cell.value ?? '';
        let subParts: string[] = WebLogo.splitterAsHelm(s);
        subParts.forEach((amino, index) => {
          let color = monomers.has(amino) ? 'red' : grayColor;
          g.fillStyle = undefinedColor;
          let last = index === subParts.length - 1;
          x1 = printLeftOrCentered(x1, y, w, h, g, amino, color, 0, true, 1.0, '/', last);
        });
        g.restore();
        return;
      }
    } else {
      const [type, subtype, paletteType] = gridCell.cell.column.getTag(DG.TAGS.UNITS).split(':');
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

      const separator = gridCell.cell.column.getTag('separator') ?? '';
      const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, gridCell.cell.column.getTag('separator'));

      const subParts: string[] = splitterFunc(cell.value);
      let x1 = x;
      let color = undefinedColor;
      // get max length word in subParts
      let tagUnits = gridCell.cell.column.getTag(DG.TAGS.UNITS);
      let maxLength = 0;
      let maxWord = '';
      let drawStyle = 'classic';
      if (tagUnits.includes('MSA')) {
        subParts.forEach(part => {
          if (part.length > maxLength) {
            maxLength = part.length;
            maxWord = part;
            drawStyle = 'msa';
          }
        });
      }
      subParts.forEach((amino, index) => {
        color = palette.get(amino);
        g.fillStyle = undefinedColor;
        let last = index === subParts.length - 1;
        x1 = printLeftOrCentered(x1, y, w, h, g, amino, color, 0, true, 1.0, separator, last, drawStyle, maxWord);
      });

      g.restore();
      return;
    }
  }
}

export class MonomerCellRenderer extends DG.GridCellRenderer {
  get name(): string {return 'MonomerCR';}

  get cellType(): string {return C.SEM_TYPES.MONOMER;}

  get defaultHeight(): number {return 15;}

  get defaultWidth(): number {return 30;}

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

    const palette = getPalleteByType(gridCell.tableColumn!.tags[C.TAGS.ALPHABET]);
    const s: string = gridCell.cell.value ? gridCell.cell.value : '-';
    const color = palette.get(s);

    printLeftOrCentered(x, y, w, h, g, s, color, 0, false);
    g.restore();
  }
}

export class MacromoleculeDifferenceCellRenderer extends DG.GridCellRenderer {
  get name(): string {return 'MacromoleculeDifferenceCR';}

  get cellType(): string {return C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;}

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
    const separator = gridCell.tableColumn!.tags[C.TAGS.SEPARATOR];
    const units: string = gridCell.tableColumn!.tags[DG.TAGS.UNITS];
    const splitter = WebLogo.getSplitter(units, separator);
    const subParts1 = splitter(s1);
    const subParts2 = splitter(s2);
    const [text] = processSequence(subParts1);
    const textSize = g.measureText(text.join(''));
    let updatedX = Math.max(x, x + (w - (textSize.width + subParts1.length * 4)) / 2);
    // 28 is the height of the two substitutions on top of each other + space
    const updatedY = Math.max(y, y + (h - 28) / 2);

    let palette: SeqPalette = UnknownSeqPalettes.Color;
    if (units != 'HELM')
      palette = getPalleteByType(units.substring(units.length - 2));

    const vShift = 7;
    for (let i = 0; i < subParts1.length; i++) {
      const amino1 = subParts1[i];
      const amino2 = subParts2[i];
      const color1 = palette.get(amino1);

      if (amino1 != amino2) {
        const color2 = palette.get(amino2);
        const subX0 = printLeftOrCentered(updatedX, updatedY - vShift, w, h, g, amino1, color1, 0, true);
        const subX1 = printLeftOrCentered(updatedX, updatedY + vShift, w, h, g, amino2, color2, 0, true);
        updatedX = Math.max(subX1, subX0);
      } else
        updatedX = printLeftOrCentered(updatedX, updatedY, w, h, g, amino1, color1, 0, true, 0.5);
      updatedX += 4;
    }
    g.restore();
  }
}
