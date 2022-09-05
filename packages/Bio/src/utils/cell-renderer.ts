import * as C from './constants';
import * as DG from 'datagrok-api/dg';
import {AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {UnknownSeqPalette, UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {SplitterFunc, WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import * as ui from 'datagrok-api/ui';
import {printLeftOrCentered, DrawStyle} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

const undefinedColor = 'rgb(100,100,100)';
const monomerToShortFunction: (amino: string, maxLengthOfMonomer: number) => string = WebLogo.monomerToShort;
const gapRenderer = 5;


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


export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'sequence'; }

  get cellType(): string { return 'sequence'; }

  get defaultHeight(): number { return 30; }

  get defaultWidth(): number { return 230; }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.cell.column.getTag(UnitsHandler.TAGS.aligned) !== 'SEQ.MSA') {
      return;
    }
    const maxLengthWordsSum = gridCell.cell.column.temp['bio-sum-maxLengthWords'];
    const maxIndex = gridCell.cell.column.temp['bio-maxIndex'];
    //@ts-ignore
    const argsX = e.layerX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCell.bounds.x);
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
        } else if (argsX < maxLengthWordsSum[mid]) {
          right = mid - 1;
        } else if (argsX > maxLengthWordsSum[mid + 1]) {
          left = mid + 1;
        }
        if (left == right) {
          found = true;
        }
      }
    }
    left = (argsX >= maxLengthWordsSum[left]) ? left + 1 : left;
    const separator = gridCell.cell.column.getTag('separator') ?? '';
    const splitterFunc: SplitterFunc = WebLogo.getSplitter('separator', separator);
    const subParts: string[] = splitterFunc(gridCell.cell.value);
    (((subParts[left]?.length ?? 0) > 0)) ? ui.tooltip.show(ui.div(subParts[left]), e.x + 16, e.y + 16) : ui.tooltip.hide();
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
    const minDistanceRenderer = 50;
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
    const splitLimit = gridCell.bounds.width / 5;
    const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, separator, gridCell.bounds.width / 5);


    const maxLengthOfMonomer = 8;

    let maxLengthWords: any = {};
    if (gridCell.cell.column.getTag('.calculatedCellRender') !== splitLimit.toString()) {
      let samples = 0;
      while (samples < Math.min(gridCell.cell.column.length, 100)) {
        let column = gridCell.cell.column.get(samples);
        let subParts: string[] = splitterFunc(column);
        subParts.forEach((amino, index) => {
          let textSize = monomerToShortFunction(amino, maxLengthOfMonomer).length * 7 + gapRenderer;
          if (textSize > (maxLengthWords[index] ?? 0)) {
            maxLengthWords[index] = textSize;
          }
          if (index > (maxLengthWords['bio-maxIndex'] ?? 0)) {
            maxLengthWords['bio-maxIndex'] = index;
          }
        });
        samples += 1;
      }
      let minLength = 3 * 7;
      for (let i = 0; i <= maxLengthWords['bio-maxIndex']; i++) {
        if (maxLengthWords[i] < minLength) {
          maxLengthWords[i] = minLength;
        }
      }
      let maxLengthWordSum: any = {};
      maxLengthWordSum[0] = maxLengthWords[0];
      for (let i = 1; i <= maxLengthWords['bio-maxIndex']; i++) {
        maxLengthWordSum[i] = maxLengthWordSum[i - 1] + maxLengthWords[i];
      }
      gridCell.cell.column.temp = {
        'bio-sum-maxLengthWords': maxLengthWordSum,
        'bio-maxIndex': maxLengthWords['bio-maxIndex'],
        'bio-maxLengthWords': maxLengthWords
      };
      gridCell.cell.column.setTag('.calculatedCellRender', splitLimit.toString());
    } else {
      maxLengthWords = gridCell.cell.column.temp['bio-maxLengthWords'];
    }

    const subParts: string[] = splitterFunc(cell.value);
    let x1 = x;
    let color = undefinedColor;
    let drawStyle = DrawStyle.classic;
    if (gridCell.cell.column.getTag('aligned').includes('MSA') && gridCell.cell.column.getTag('units') === 'separator') {
      drawStyle = DrawStyle.MSA;
    }
    subParts.every((amino, index) => {
      color = palette.get(amino);
      g.fillStyle = undefinedColor;
      let last = index === subParts.length - 1;
      x1 = printLeftOrCentered(x1, y, w, h, g, monomerToShortFunction(amino, maxLengthOfMonomer), color, 0, true, 1.0, separator, last, drawStyle, maxLengthWords, index, gridCell);
      if (x1 - minDistanceRenderer - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCell.bounds.x) > gridCell.bounds.width) {
        return false;
      }
      return true;
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
