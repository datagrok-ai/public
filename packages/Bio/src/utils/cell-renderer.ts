import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {_package} from '../package';
import {printLeftOrCentered, DrawStyle} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import * as C from './constants';
import {
  ALIGNMENT,
  getPaletteByType,
  getSplitter,
  getSplitterForColumn,
  monomerToShort,
  MonomerToShortFunc,
  NOTATION,
  SplitterFunc,
  TAGS as bioTAGS,
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

const enum tempTAGS {
  referenceSequence = 'reference-sequence',
  currentWord = 'current-word',
  monomerWidth = 'monomer-width',
  bioSeqCol = 'bio-seqCol',
}

const enum rndrTAGS {
  calculatedCellRender = '.calculatedCellRender',
}

type TempType = { [tagName: string]: any };

const undefinedColor = 'rgb(100,100,100)';
const monomerToShortFunction: MonomerToShortFunc = monomerToShort;

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

class SeqColForTemp {
  private readonly _uh: UnitsHandler;
  private readonly _splitter: SplitterFunc;
  private _monomerLengthList: number[][] | null = null;

  private _updated: boolean = false;
  public get updated(): boolean { return this._updated; }

  constructor(
    public col: DG.Column<string>,
    private monomerToShort: MonomerToShortFunc,
    private monomerLabelLengthLimit: number,
    private monomerCharWidth: number = 7,
    private monomerSeparatorWidth: number = 12
  ) {
    this._uh = new UnitsHandler(col);

    this._splitter = getSplitterForColumn(this.col);

    col.dataFrame.onDataChanged.subscribe(() => {
      this._monomerLengthList = null;
    });
  }

  public getCellMonomerLengths(rowIdx: number): number[] {
    if (this._uh.isMsa())
      return this.getCellMonomerLengthsForSeqMsa();
    else
      return this.getCellMonomerLengthsForSeq(rowIdx);
  }

  private getCellMonomerLengthsForSeq(rowIdx: number): number[] {
    if (this._monomerLengthList === null) {
      this._monomerLengthList = new Array(this.col.length).fill(null);
      this._updated = true;
    }

    let res: number[] = this._monomerLengthList[rowIdx];
    if (res === null) {
      const seqMonList: string[] = this.getSeqMonList(rowIdx);
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqMonList.length);

      for (const [seqMonI, seqMonLabel] of seqMonList.entries()) {
        const seqMonWidth: number = this.monomerSeparatorWidth +
          this.monomerToShort(seqMonLabel, this.monomerLabelLengthLimit).length * this.monomerCharWidth;
        res[seqMonI] = seqMonWidth;
      }
      this._updated = true;
    }

    return res;
  }

  private getCellMonomerLengthsForSeqMsa(): number[] {
    if (this._monomerLengthList === null) {
      this._monomerLengthList = new Array(1).fill(null);
      this._updated = true;
    }
    let res: number[] | null = this._monomerLengthList[0];
    if (res === null) {
      res = this._monomerLengthList[0] = new Array(0);
      for (let seqIdx = 0; seqIdx < Math.max(this.col.length, 100); seqIdx++) {
        const seqMonList: string[] = this.getSeqMonList(seqIdx);
        if (seqMonList.length > res.length)
          res.push(...new Array<number>(seqMonList.length - res.length).fill(0));

        for (const [seqMonI, seqMonLabel] of seqMonList.entries()) {
          const seqMonWidth: number = this.monomerSeparatorWidth +
            this.monomerToShort(seqMonLabel, this.monomerLabelLengthLimit).length * this.monomerCharWidth;
          res[seqMonI] = Math.max(res[seqMonI] ?? 0, seqMonWidth);
        }
      }
      this._updated = true;
    }
    return res; // first (and single) row of data
  }

  /** Returns seq position for pointer x */
  public getPosition(rowIdx: number, x: number): number | null {
    const monomerMaxLengths: number[] = this.getCellMonomerLengths(rowIdx);
    const seq: string = this.col.get(rowIdx)!;
    const seqMonList: string[] = this._splitter(seq);

    let left: number | null = null;
    let right = seqMonList.length;
    let found = false;
    let mid = 0;
    if (x > monomerMaxLengths[0]) {
      while (!found) {
        mid = Math.floor((right + (left ?? 0)) / 2);
        if (x >= monomerMaxLengths[mid] && x <= monomerMaxLengths[mid + 1]) {
          left = mid;
          found = true;
        } else if (x < monomerMaxLengths[mid]) {
          right = mid - 1;
        } else if (x > monomerMaxLengths[mid + 1]) {
          left = mid + 1;
        }
        if (left == right)
          found = true;
      }
    }
    return left;
  }

  getSeqMonList(rowIdx: number): string[] {
    const seq: string | null = this.col.get(rowIdx);
    return seq ? this._splitter(seq) : [];
  }
}

export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'sequence'; }

  get cellType(): string { return 'sequence'; }

  get defaultHeight(): number { return 30; }

  get defaultWidth(): number { return 230; }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    const colTemp: TempType = gridCell.cell.column.temp;
    colTemp[tempTAGS.currentWord] = gridCell.cell.value;
    gridCell.grid.invalidate();
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    // if (gridCell.cell.column.getTag(bioTAGS.aligned) !== ALIGNMENT.SEQ_MSA)
    //   return;

    const tableCol: DG.Column = gridCell.cell.column;
    const tableColTemp: TempType = tableCol.temp;
    const seqColTemp: SeqColForTemp = tableCol.temp[tempTAGS.bioSeqCol];
    if (!seqColTemp) return; // Can do nothing without precalculated data

    const gridCellBounds: DG.Rect = gridCell.bounds;
    // const value: any = gridCell.cell.value;
    //
    // const maxLengthWords: number[] = seqColTemp.getCellMonomerLengths(gridCell.tableRowIndex!);
    // const maxLengthWordsSum: number[] = new Array<number>(maxLengthWords.length).fill(0);
    // for (let posI: number = 1; posI < maxLengthWords.length; posI++)
    //   maxLengthWordsSum[posI] = maxLengthWordsSum[posI - 1] + maxLengthWords[posI];
    // const maxIndex = maxLengthWords.length;

    const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCellBounds.x);
    const left: number | null = seqColTemp.getPosition(gridCell.tableRowIndex!, argsX);

    const seqMonList: string[] = seqColTemp.getSeqMonList(gridCell.tableRowIndex!);
    if (left !== null && left < seqMonList.length) {
      const monomer: string = seqMonList[left];
      ui.tooltip.show(ui.divV([monomer, `left: ${left}`, `argsX: ${argsX}`]), e.x + 16, e.y + 16);
    } else {
      ui.tooltip.hide();
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
    let gapRenderer = 5;
    let maxLengthOfMonomer = 8;

    // TODO: Store temp data to GridColumn
    // Now the renderer requires data frame table Column underlying GridColumn
    const tableCol: DG.Column = gridCell.cell.column;
    const tableColTemp: TempType = tableCol.temp;

    // Cell renderer settings
    const tempMonomerWidth: string | null = tableColTemp[tempTAGS.monomerWidth];
    const monomerWidth: string = (tempMonomerWidth != null) ? tempMonomerWidth : 'short';
    if (monomerWidth === 'short') {
      gapRenderer = 12;
      maxLengthOfMonomer = 1;
    }

    let seqColTemp: SeqColForTemp = tableCol.temp[tempTAGS.bioSeqCol];
    if (!seqColTemp)
      seqColTemp = new SeqColForTemp(tableCol, monomerToShortFunction, maxLengthOfMonomer);

    const maxLengthWords: number[] = seqColTemp.getCellMonomerLengths(gridCell.tableRowIndex!);
    const maxLengthWordsSum: number[] = new Array<number>(maxLengthWords.length).fill(0);
    for (let posI: number = 1; posI < maxLengthWords.length; posI++)
      maxLengthWordsSum[posI] = maxLengthWordsSum[posI - 1] + maxLengthWords[posI];
    const maxIndex = maxLengthWords.length;

    // Store updated seqColTemp to the col temp
    if (seqColTemp.updated) tableColTemp[tempTAGS.bioSeqCol] = seqColTemp;

    g.save();
    try {
      const grid = gridCell.gridRow !== -1 ? gridCell.grid : null;
      const value: any = gridCell.cell.value;
      const paletteType = tableCol.getTag(bioTAGS.alphabet);
      const minDistanceRenderer = 50;
      w = getUpdatedWidth(grid, g, x, w);
      g.beginPath();
      g.rect(x, y, w, h);
      g.clip();
      g.font = '12px monospace';
      g.textBaseline = 'top';

      //TODO: can this be replaced/merged with splitSequence?
      const units = tableCol.getTag(DG.TAGS.UNITS);
      const aligned: string = tableCol.getTag(bioTAGS.aligned);

      const palette = getPaletteByType(paletteType);

      const separator = tableCol.getTag(bioTAGS.separator) ?? '';
      const splitLimit = w / 5;
      const splitterFunc: SplitterFunc = getSplitter(units, separator, splitLimit);

      const tempReferenceSequence: string | null = tableColTemp[tempTAGS.referenceSequence];
      const tempCurrentWord: string | null = tableColTemp[tempTAGS.currentWord];
      const referenceSequence: string[] = splitterFunc(
        ((tempReferenceSequence != null) && (tempReferenceSequence != '')) ?
          tempReferenceSequence : tempCurrentWord ?? '');

      // let maxLengthWords: { [pos: number]: number } = {};
      // if (tableCol.getTag(rndrTAGS.calculatedCellRender) !== splitLimit.toString()) {
      //   let sampleCount = 0;
      //   while (sampleCount < Math.min(tableCol.length, 100)) {
      //     const rowIdx: number = sampleCount;
      //     const column = tableCol.get(rowIdx);
      //     const subParts: string[] = splitterFunc(column);
      //     for (const [index, amino] of subParts.entries()) {
      //       const textSize = monomerToShortFunction(amino, maxLengthOfMonomer).length * 7 + gapRenderer;
      //       if (textSize > (maxLengthWords[index] ?? 0))
      //         maxLengthWords[index] = textSize;
      //       if (index > maxIndex) maxIndex = index;
      //     }
      //     sampleCount += 1;
      //   }
      //   const minLength = 3 * 7;
      //   for (let i = 0; i <= maxIndex; i++) {
      //     if (maxLengthWords[i] < minLength) maxLengthWords[i] = minLength;
      //     const maxLengthWordSum: { [pos: number]: number } = {};
      //     maxLengthWordSum[0] = maxLengthWords[0];
      //     for (let i = 1; i <= maxIndex; i++) maxLengthWordSum[i] = maxLengthWordSum[i - 1] + maxLengthWords[i];
      //     colTemp[tempTAGS.bioSumMaxLengthWords] = maxLengthWordSum;
      //     colTemp[tempTAGS.bioMaxIndex] = maxIndex;
      //     colTemp[tempTAGS.bioMaxLengthWords] = maxLengthWords;
      //     tableCol.setTag(rndrTAGS.calculatedCellRender, splitLimit.toString());
      //   }
      // } else {
      //   maxLengthWords = colTemp[tempTAGS.bioMaxLengthWords];
      // }

      const subParts: string[] = splitterFunc(value);
      let x1 = x;
      let color = undefinedColor;
      let drawStyle = DrawStyle.classic;

      if (aligned && aligned.includes('MSA') && units == NOTATION.SEPARATOR)
        drawStyle = DrawStyle.MSA;

      // subParts.every((amino, index) => {
      // });
      for (const [index, amino] of subParts.entries()) {
        color = palette.get(amino);
        g.fillStyle = undefinedColor;
        const last = index === subParts.length - 1;
        x1 = printLeftOrCentered(x1, y, w, h,
          g, amino, color, 0, true, 1.0, separator, last, drawStyle,
          maxLengthWords, index, gridCell, referenceSequence, maxLengthOfMonomer);
        if (minDistanceRenderer > w) break;
      }
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : !!err ? err.toString() : 'Error \'undefined\'';
      _package.logger.error(`Bio: MacromoleculeSequenceCellRenderer.render() error: ${errMsg}`);
      //throw err; // Do not throw to prevent disabling renderer
    } finally {
      g.restore();
    }
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

    const palette = getPaletteByType(gridCell.cell.column.getTag(bioTAGS.alphabet));
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
    const separator = gridCell.tableColumn!.tags[bioTAGS.separator];
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
