import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {printLeftOrCentered, DrawStyle} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import * as C from './constants';
import {
  ALIGNMENT,
  getPaletteByType,
  getSplitter,
  monomerToShort,
  SeqPalette,
  SplitterFunc,
  TAGS as bioTAGS,
  UnknownSeqPalettes
} from '@datagrok-libraries/bio';

const undefinedColor = 'rgb(100,100,100)';
const monomerToShortFunction: (amino: string, maxLengthOfMonomer: number) => string = monomerToShort;

function getUpdatedWidth(grid: DG.Grid | null, g: CanvasRenderingContext2D, x: number, w: number): number {
  return grid ? Math.min(grid.canvas.width - x, w) : g.canvas.width - x;
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

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    gridCell.cell.column.temp['current-word'] = gridCell.cell.value;
    gridCell.grid.invalidate();
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.cell.column.getTag(bioTAGS.aligned) !== ALIGNMENT.SEQ_MSA)
      return;

    const maxLengthWordsSum = gridCell.cell.column.temp['bio-sum-maxLengthWords'];
    const maxIndex = gridCell.cell.column.temp['bio-maxIndex'];
    const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCell.bounds.x);
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
        if (left == right)
          found = true;
      }
    }
    left = (argsX >= maxLengthWordsSum[left]) ? left + 1 : left;
    const separator = gridCell.cell.column.getTag('separator') ?? '';
    const splitterFunc: SplitterFunc = getSplitter('separator', separator);
    const subParts: string[] = splitterFunc(gridCell.cell.value);
    (((subParts[left]?.length ?? 0) > 0)) ?
      ui.tooltip.show(ui.div(subParts[left]), e.x + 16, e.y + 16) : ui.tooltip.hide();
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
  ) {
    const grid = gridCell.gridRow !== -1 ? gridCell.grid : null;
    const cell = gridCell.cell;
    const paletteType = gridCell.cell.column.getTag(C.TAGS.ALPHABET);
    const minDistanceRenderer = 50;
    w = getUpdatedWidth(grid, g, x, w);
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.font = '12px monospace';
    g.textBaseline = 'top';

    //TODO: can this be replaced/merged with splitSequence?
    const units = gridCell.cell.column.getTag(DG.TAGS.UNITS);

    const palette = getPaletteByType(paletteType);

    const separator = gridCell.cell.column.getTag('separator') ?? '';
    const splitLimit = gridCell.bounds.width / 5;
    const splitterFunc: SplitterFunc = getSplitter(units, separator, splitLimit);
    const referenceSequence: string[] = splitterFunc(((gridCell.cell.column?.temp['reference-sequence'] != null) && (gridCell.cell.column?.temp['reference-sequence'] != ''))
      ? gridCell.cell.column.temp['reference-sequence'] : gridCell.cell.column.temp['current-word'] ?? '');
    const monomerWidth = (gridCell.cell.column?.temp['monomer-width'] != null) ? gridCell.cell.column.temp['monomer-width'] : 'short';
    let gapRenderer = 5;

    let maxIndex = 0;

    let maxLengthOfMonomer = 8;

    if (monomerWidth === 'short') {
      gapRenderer = 12;
      maxLengthOfMonomer = 1;
    }

    let maxLengthWords: any = {};
    if (gridCell.cell.column.getTag('.calculatedCellRender') !== splitLimit.toString()) {
      let samples = 0;
      while (samples < Math.min(gridCell.cell.column.length, 100)) {
        const column = gridCell.cell.column.get(samples);
        const subParts: string[] = splitterFunc(column);
        subParts.forEach((amino, index) => {
          const textSize = monomerToShortFunction(amino, maxLengthOfMonomer).length * 7 + gapRenderer;
          if (textSize > (maxLengthWords[index] ?? 0))
            maxLengthWords[index] = textSize;
          if (index > maxIndex) {
            maxIndex = index;
          }
        });
        samples += 1;
      }
      const minLength = 3 * 7;
      for (let i = 0; i <= maxIndex; i++) {
        if (maxLengthWords[i] < minLength) {
          maxLengthWords[i] = minLength;
        }
        const maxLengthWordSum: any = {};
        maxLengthWordSum[0] = maxLengthWords[0];
        for (let i = 1; i <= maxIndex; i++) {
          maxLengthWordSum[i] = maxLengthWordSum[i - 1] + maxLengthWords[i];
        }
        gridCell.cell.column.temp['bio-sum-maxLengthWords'] = maxLengthWordSum;
        gridCell.cell.column.temp['bio-maxIndex'] = maxIndex;
        gridCell.cell.column.temp['bio-maxLengthWords'] = maxLengthWords;
        gridCell.cell.column.setTag('.calculatedCellRender', splitLimit.toString());
      }
    } else {
      maxLengthWords = gridCell.cell.column.temp['bio-maxLengthWords'];
    }

    const subParts: string[] = splitterFunc(cell.value);
    let x1 = x;
    let color = undefinedColor;
    let drawStyle = DrawStyle.classic;
    if (gridCell.cell.column.getTag('aligned').includes('MSA') && gridCell.cell.column.getTag('units') === 'separator')
      drawStyle = DrawStyle.MSA;

    subParts.every((amino, index) => {
      color = palette.get(amino);
      g.fillStyle = undefinedColor;
      const last = index === subParts.length - 1;
      x1 = printLeftOrCentered(x1, y, w, h, g, amino, color, 0, true, 1.0, separator, last, drawStyle, maxLengthWords, index, gridCell, referenceSequence, maxLengthOfMonomer);
      return x1 - minDistanceRenderer - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCell.bounds.x) <= gridCell.bounds.width;
    });

    g.restore();
    return;
  }
}

export class MonomerCellRenderer extends DG.GridCellRenderer {
  get name(): string { return C.SEM_TYPES.MONOMER; }

  get cellType(): string { return C.SEM_TYPES.MONOMER; }

  get defaultHeight(): number { return 15; }

  get defaultWidth(): number { return 30; }

  /**
   * Cell renderer function.
   *
   * @param {CanvasRenderingContext2D} g Canvas rendering context.
   * @param {number} x x coordinate on the canvas.
   * @param {number} y y coordinate on the canvas.
   * @param {number} w width of the cell.
   * @param {number} h height of the cell.
   * @param {DG.GridCell} gridCell Grid cell.
   * @param {DG.GridCellStyle} _cellStyle Cell style.
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    _cellStyle: DG.GridCellStyle): void {
    g.font = `12px monospace`;
    g.textBaseline = 'middle';
    g.textAlign = 'center';

    const palette = getPaletteByType(gridCell.cell.column.getTag(C.TAGS.ALPHABET));
    const s: string = gridCell.cell.value;
    if (!s)
      return;
    const color = palette.get(s);

    g.fillStyle = color;
    g.fillText(monomerToShort(s, 3), x + (w / 2), y + (h / 2), w);
  }
}

export class MacromoleculeDifferenceCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'MacromoleculeDifferenceCR'; }

  get cellType(): string { return C.SEM_TYPES.MACROMOLECULE_DIFFERENCE; }

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
   * @param {DG.GridCellStyle} _cellStyle Cell style.
   * @memberof AlignedSequenceDifferenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    _cellStyle: DG.GridCellStyle): void {
    const grid = gridCell.grid;
    const cell = gridCell.cell;
    const s: string = cell.value ?? '';
    const separator = gridCell.tableColumn!.tags[C.TAGS.SEPARATOR];
    const units: string = gridCell.tableColumn!.tags[DG.TAGS.UNITS];
    w = getUpdatedWidth(grid, g, x, w);
    //TODO: can this be replaced/merged with splitSequence?
    const [s1, s2] = s.split('#');
    const splitter = getSplitter(units, separator);
    const subParts1 = splitter(s1);
    const subParts2 = splitter(s2);
    drawMoleculeDifferenceOnCanvas(g, x, y, w, h, subParts1, subParts2, units);
  }
}

export function drawMoleculeDifferenceOnCanvas(
  g: CanvasRenderingContext2D,
  x: number,
  y: number,
  w: number,
  h: number,
  subParts1: string [],
  subParts2: string [],
  units: string,
  fullStringLength?: boolean,
  molDifferences?: { [key: number]: HTMLCanvasElement }
): void {
  if (subParts1.length !== subParts2.length) {
    const emptyMonomersArray = new Array<string>(Math.abs(subParts1.length - subParts2.length)).fill('');
    subParts1.length > subParts2.length ?
      subParts2 = subParts2.concat(emptyMonomersArray) : subParts1 = subParts1.concat(emptyMonomersArray);
  }
  const textSize1 = g.measureText(processSequence(subParts1).join(''));
  const textSize2 = g.measureText(processSequence(subParts2).join(''));
  const textWidth = Math.max(textSize1.width, textSize2.width);
  if (fullStringLength) {
    w = textWidth + subParts1.length * 4;
    g.canvas.width = textWidth + subParts1.length * 4;
  }
  let updatedX = Math.max(x, x + (w - (textWidth + subParts1.length * 4)) / 2);
  // 28 is the height of the two substitutions on top of each other + space
  const updatedY = Math.max(y, y + (h - 28) / 2);

  g.save();
  g.beginPath();
  g.rect(x, y, fullStringLength ? textWidth + subParts1.length * 4 : w, h);
  g.clip();
  g.font = '12px monospace';
  g.textBaseline = 'top';

  let palette: SeqPalette = UnknownSeqPalettes.Color;
  if (units != 'HELM')
    palette = getPaletteByType(units.substring(units.length - 2));

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
      if (molDifferences)
        molDifferences[i] = createDifferenceCanvas(amino1, amino2, color1, color2, updatedY, vShift, h);
    } else { updatedX = printLeftOrCentered(updatedX, updatedY, w, h, g, amino1, color1, 0, true, 0.5); }
    updatedX += 4;
  }
  g.restore();
}

function createDifferenceCanvas(
  amino1: string,
  amino2: string,
  color1: string,
  color2: string,
  y: number,
  shift: number,
  h: number): HTMLCanvasElement {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d')!;
  context.font = '12px monospace';
  const width1 = context.measureText(processSequence([amino1]).join('')).width;
  const width2 = context.measureText(processSequence([amino2]).join('')).width;
  const width = Math.max(width1, width2);
  canvas.height = h;
  canvas.width = width + 4;
  context.font = '12px monospace';
  context.textBaseline = 'top';
  printLeftOrCentered(0, y - shift, width, h, context, amino1, color1, 0, true);
  printLeftOrCentered(0, y + shift, width, h, context, amino2, color2, 0, true);
  return canvas;
}
