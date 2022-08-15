import * as C from './constants';
import * as DG from 'datagrok-api/dg';
import {AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {UnknownSeqPalette, UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {SplitterFunc, WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import * as ui from 'datagrok-api/ui';

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
export function printLeftOrCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color = undefinedColor,
  pivot: number = 0, left = false, transparencyRate: number = 1.0,
  separator: string = '', last: boolean = false, drawStyle: string = 'classic', maxWord: any = {}, maxWordIdx: number = 0, gridCell: any = {}): number {
  g.textAlign = 'start';
  const colorPart = s.substring(0);
  let grayPart = last ? '' : separator;
  if (drawStyle === 'msa') {
    grayPart = '';
  }

  let textSize: any = g.measureText(colorPart + grayPart);
  const indent = 5;

  let maxColorTextSize = g.measureText(colorPart).width;
  let colorTextSize = g.measureText(colorPart).width;
  const dy = (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2;
  textSize = textSize.width;
  if (drawStyle === 'msa') {
    maxColorTextSize = maxWord[maxWordIdx];
    textSize = maxWord[maxWordIdx];
    if (maxColorTextSize > maxWord) {
      maxWord[maxWordIdx] = maxColorTextSize;
      gridCell.cell.column.temp = maxWord;
    }
    if (maxWordIdx > (maxWord['maxIndex'] ?? 0)) {
      maxWord['maxIndex'] = maxWordIdx;
      gridCell.cell.column.temp = maxWord;
    }
  }

  function draw(dx1: number, dx2: number): void {
    g.fillStyle = color;
    g.globalAlpha = transparencyRate;
    if (drawStyle === 'classic') {
      g.fillText(colorPart, x + dx1, y + dy);
      g.fillStyle = grayColor;
      g.fillText(grayPart, x + dx2, y + dy);
    }
    if (drawStyle === 'msa') {
      g.fillStyle = color;
      g.fillText(colorPart, x + dx1 + ((maxWord[maxWordIdx] - colorTextSize) / 2), y + dy);
    }
  }

  if (left || textSize > w) {
    draw(indent, indent + maxColorTextSize);
    return x + maxColorTextSize + g.measureText(grayPart).width;

  } else {
    const dx = (w - textSize) / 2;
    draw(dx, dx + maxColorTextSize);
    return x + dx + maxColorTextSize;
  }
}


export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'sequence'; }

  get cellType(): string { return 'sequence'; }

  get defaultHeight(): number { return 30; }

  get defaultWidth(): number { return 230; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.cell.column.getTag('aligned').includes('MSA')) {
      let val = gridCell.cell.value;
      let maxLengthWordsSum = gridCell.cell.column.temp['sum'];
      let maxIndex = gridCell.cell.column.temp['maxIndex'];
      //@ts-ignore
      let argsX = e.layerX;
      let left = 0;
      let right = maxIndex;
      let found = false;
      maxLengthWordsSum[maxIndex + 1] = argsX + 1;
      let mid = 0;
      if (argsX > maxLengthWordsSum[0]) {
        while (!found) {
          mid = Math.floor((right + left) / 2);
          if (argsX >= maxLengthWordsSum[mid] && argsX <= maxLengthWordsSum[mid + 1]) {
            left = mid;
            found = true;
          }
          if (argsX < maxLengthWordsSum[mid]) {
            right = mid - 1;
          }
          if (argsX > maxLengthWordsSum[mid + 1]) {
            left = mid + 1;
          }
          if (left  == right) {
            found = true;
          }
        }
      }
      const separator = gridCell.cell.column.getTag('separator') ?? '';
      const splitterFunc: SplitterFunc = WebLogo.getSplitter('separator', separator);
      const subParts: string[] = splitterFunc(gridCell.cell.value);
      ui.tooltip.show(ui.div(subParts[left]), e.x + 16, e.y + 16);
    }
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
    cellStyle: DG.GridCellStyle
  ): void {
    const grid = gridCell.gridRow !== -1 ? gridCell.grid : undefined;
    const cell = gridCell.cell;
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
    const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, separator);

    const columns = gridCell.cell.column.categories;
    let monomerToShortFunction: (amino: string, maxLengthOfMonomer: number) => string = WebLogo.monomerToShort;
    let maxLengthOfMonomer = 8;

    let maxLengthWords = {};
    if (gridCell.cell.column.getTag('.calculatedCellRender') !== 'exist') {
      for (let i = 0; i < columns.length; i++) {
        let subParts: string[] = splitterFunc(columns[i]);
        subParts.forEach((amino, index) => {
          //@ts-ignore
          let textSizeWidth = g.measureText(monomerToShortFunction(amino, maxLengthOfMonomer));
          //@ts-ignore
          if (textSizeWidth.width > (maxLengthWords[index] ?? 0)) {
            //@ts-ignore
            maxLengthWords[index] = textSizeWidth.width;
          }
          //@ts-ignore
          if (index > (maxLengthWords['maxIndex'] ?? 0)) {
            //@ts-ignore
            maxLengthWords['maxIndex'] = index;
          }
        });
      }
      let maxLengthWordSum = {};
      //@ts-ignore
      maxLengthWordSum[-1] = 0;
      //@ts-ignore
      for (let i = 0; i <= maxLengthWords['maxIndex']; i++) {
        //@ts-ignore
        maxLengthWordSum[i] = maxLengthWordSum[i - 1] + maxLengthWords[i];
      }
      //@ts-ignore
      maxLengthWords['sum'] = maxLengthWordSum;
      gridCell.cell.column.temp = maxLengthWords;
      gridCell.cell.column.setTag('.calculatedCellRender', 'exist');
    } else {
      maxLengthWords = gridCell.cell.column.temp;
    }

    const subParts: string[] = splitterFunc(cell.value);
    let x1 = x;
    let color = undefinedColor;
    let drawStyle = 'classic';
    if (gridCell.cell.column.getTag('aligned').includes('MSA')) {
      drawStyle = 'msa';
    }
    subParts.forEach((amino, index) => {
      color = palette.get(amino);
      g.fillStyle = undefinedColor;
      let last = index === subParts.length - 1;
      x1 = printLeftOrCentered(x1, y, w, h, g, monomerToShortFunction(amino, maxLengthOfMonomer), color, 0, true, 1.0, separator, last, drawStyle, maxLengthWords, index, gridCell);
    });

    g.restore();
    return;
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
