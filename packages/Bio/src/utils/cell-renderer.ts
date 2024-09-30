/* eslint-disable max-lines */
/* eslint-disable max-params */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {printLeftOrCentered, DrawStyle, TAGS as mmcrTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {getUpdatedWidth, MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {
  getPaletteByType,
  monomerToShort,
  MonomerToShortFunc,
  NOTATION,
  SplitterFunc,
  TAGS as bioTAGS,
  ALPHABET,
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {getSplitter} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {GapOriginals} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqMonomer} from '@datagrok-libraries/bio/src/helm/types';
import {execMonomerHoverLinks} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';
import {getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {MmcrTemps, rendererSettingsChangedState, tempTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import * as C from './constants';

import {_package, getMonomerLib} from '../package';
import {max} from 'rxjs/operators';

type TempType = { [tagName: string]: any };

export const monomerToShortFunction: MonomerToShortFunc = monomerToShort;

export function processSequence(subParts: string[]): [string[], boolean] {
  const simplified = !wu.enumerate(subParts).some(([amino, index]) =>
    amino.length > 1 &&
    index != 0 &&
    index != subParts.length - 1);

  const text: string[] = [];
  const gap = simplified ? '' : ' ';
  for (const [amino, index] of wu.enumerate(subParts)) {
    let aminoRes = amino;
    if (index < subParts.length)
      aminoRes += `${amino ? '' : '-'}${gap}`;

    text.push(aminoRes);
  }
  return [text, simplified];
}

type RendererGridCellTemp = {
  [MmcrTemps.monomerPlacer]: MonomerPlacer
}

export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'sequence'; }

  get cellType(): string { return 'sequence'; }

  get defaultHeight(): number | null { return 30; }

  get defaultWidth(): number | null { return 230; }

  onClick(gridCell: DG.GridCell, _e: MouseEvent): void {
    const colTemp: TempType = gridCell.cell.column.temp;
    colTemp[tempTAGS.currentWord] = gridCell.cell.value;
    gridCell.grid.invalidate();
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.tableRowIndex == null) return;

    // if (gridCell.cell.column.getTag(bioTAGS.aligned) !== ALIGNMENT.SEQ_MSA)
    //   return;

    const [_gridCol, tableCol, temp] =
      getGridCellColTemp<string, MonomerPlacer>(gridCell);
    const seqColTemp: MonomerPlacer = temp.rendererBack;
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
    const left: number | null = seqColTemp.getPosition(gridCell.tableRowIndex!, argsX, gridCellBounds.width);

    const seqSS = SeqHandler.forColumn(tableCol)
      .getSplitted(gridCell.tableRowIndex!);
    if (left !== null && left < seqSS.length) {
      const sh = SeqHandler.forColumn(tableCol);
      const alphabet = sh.alphabet ?? ALPHABET.UN;
      const seqMonomer = {
        position: left,
        biotype: alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA,
        symbol: seqSS.getCanonical(left),
      } as ISeqMonomer;
      const tooltipElements: HTMLElement[] = [];
      let monomerDiv = seqColTemp._monomerStructureMap[seqMonomer.symbol];
      if (!monomerDiv || true) {
        monomerDiv = seqColTemp._monomerStructureMap[seqMonomer.symbol] = (() => {
          const lib: IMonomerLib | null = getMonomerLib();
          return lib ? lib.getTooltip(seqMonomer.biotype, seqMonomer.symbol) : ui.divText('Monomer library is not available');
        })();
      }
      tooltipElements.push(monomerDiv);
      ui.tooltip.show(ui.divV(tooltipElements), e.x + 16, e.y + 16);

      execMonomerHoverLinks(gridCell, seqMonomer);
    } else {
      //
      ui.tooltip.hide();
      execMonomerHoverLinks(gridCell, null);
    }
  }

  override onMouseLeave(gridCell: DG.GridCell, e: MouseEvent) {
    execMonomerHoverLinks(gridCell, null);
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
   * @param {DG.GridCellStyle} _cellStyle Cell style.
   * @memberof AlignedSequenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    _cellStyle: DG.GridCellStyle
  ): void {
    const logPrefix: string = 'MacromoleculeSequenceCellRenderer.render()';

    const [gridCol, tableCol, _temp] =
      getGridCellColTemp<string, MonomerPlacer>(gridCell);
    if (!tableCol) return;
    const tableColTemp: TempType = tableCol.temp;

    let gapLength = 0;
    const msaGapLength = 8;

    // Cell renderer settings
    let maxLengthOfMonomer: number = (_package.properties ? _package.properties.maxMonomerLength : 4) ?? 50;
    if (mmcrTAGS.maxMonomerLength in tableCol.tags) {
      const v = parseInt(tableCol.getTag(mmcrTAGS.maxMonomerLength));
      maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
    }
    if (MmcrTemps.maxMonomerLength in tableColTemp) {
      const v = tableColTemp[MmcrTemps.maxMonomerLength];
      maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
    }

    const [_gc, _tc, temp] =
      getGridCellColTemp<string, MonomerPlacer>(gridCell);
    let seqColTemp: MonomerPlacer = temp.rendererBack;
    if (!seqColTemp) {
      seqColTemp = temp.rendererBack = new MonomerPlacer(gridCol, tableCol, _package.logger, maxLengthOfMonomer,
        () => {
          const sh = SeqHandler.forColumn(tableCol);
          return {
            seqHandler: sh,
            monomerCharWidth: 7, separatorWidth: !sh.isMsa() ? gapLength : msaGapLength,
            monomerToShort: monomerToShortFunction,
          };
        });
    }

    if (
      tableCol.temp[MmcrTemps.rendererSettingsChanged] === rendererSettingsChangedState.true ||
      seqColTemp.monomerLengthLimit != maxLengthOfMonomer
    ) {
      gapLength = tableColTemp[MmcrTemps.gapLength] as number ?? gapLength;
      // this event means that the mm renderer settings have changed,
      // particularly monomer representation and max width.
      seqColTemp.setMonomerLengthLimit(maxLengthOfMonomer);
      seqColTemp.setSeparatorWidth(seqColTemp.isMsa() ? msaGapLength : gapLength);
      tableCol.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.false;
    }

    seqColTemp.render(g, x, y, w, h, gridCell, _cellStyle);
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
    const dpr = window.devicePixelRatio;
    const grid = gridCell.grid;
    const cell = gridCell.cell;
    const tableCol = gridCell.tableColumn as DG.Column<string>;
    const s: string = cell.value ?? '';
    const separator = tableCol.tags[bioTAGS.separator];
    const units: string = tableCol.meta.units!;
    w = getUpdatedWidth(grid, g, x, w, dpr);
    //TODO: can this be replaced/merged with splitSequence?
    const [s1, s2] = s.split('#');
    const splitter = getSplitter(units, separator);
    const s1SS = splitter(s1);
    const s2SS = splitter(s2);
    const subParts1 = wu.count(0).take(s1SS.length).map((posIdx) => s1SS.getCanonical(posIdx)).toArray();
    const subParts2 = wu.count(0).take(s2SS.length).map((posIdx) => s2SS.getCanonical(posIdx)).toArray();
    drawMoleculeDifferenceOnCanvas(g, x, y, w, h, subParts1, subParts2, units);
  }
}

export function drawMoleculeDifferenceOnCanvas(
  g: CanvasRenderingContext2D,
  x: number,
  y: number,
  w: number,
  h: number,
  subParts1: string[],
  subParts2: string[],
  units: string,
  fullStringLength?: boolean,
  molDifferences?: { [key: number]: HTMLCanvasElement },
): void {
  if (subParts1.length !== subParts2.length) {
    const sequences: IComparedSequences = fillShorterSequence(subParts1, subParts2);
    subParts1 = sequences.subParts1;
    subParts2 = sequences.subParts2;
  }
  const textSize1 = g.measureText(processSequence(subParts1).join(''));
  const textSize2 = g.measureText(processSequence(subParts2).join(''));
  const textWidth = Math.max(textSize1.width, textSize2.width);
  if (fullStringLength) {
    w = textWidth + subParts1.length * 4;
    g.canvas.width = textWidth + subParts1.length * 4;
  }
  let updatedX = Math.max(x, x + (w - (textWidth + subParts1.length * 4)) / 2) + 5;
  // 28 is the height of the two substitutions on top of each other + space
  const updatedY = Math.max(y, y + (h - 28) / 2);

  g.save();
  g.beginPath();
  g.rect(x, y, fullStringLength ? textWidth + subParts1.length * 4 : w, h);
  g.clip();
  g.font = '12px monospace';
  g.textBaseline = 'top';

  let palette: SeqPalette = UnknownSeqPalettes.Color;
  if (units !== 'HELM')
    palette = getPaletteByType(units.substring(units.length - 2));

  const vShift = 7;
  for (let i = 0; i < subParts1.length; i++) {
    const amino1 = subParts1[i];
    const amino2 = subParts2[i];
    const color1 = palette.get(amino1);

    if (amino1 != amino2) {
      const color2 = palette.get(amino2);
      const subX0 = printLeftOrCentered(g, amino1, updatedX, updatedY - vShift, w, h,
        {color: color1, pivot: 0, left: true});
      const subX1 = printLeftOrCentered(g, amino2, updatedX, updatedY + vShift, w, h,
        {color: color2, pivot: 0, left: true});
      updatedX = Math.max(subX1, subX0);
      if (molDifferences)
        molDifferences[i] = createDifferenceCanvas(amino1, amino2, color1, color2, updatedY, vShift, h);
    } else {
      //
      updatedX = printLeftOrCentered(g, amino1, updatedX, updatedY, w, h,
        {color: color1, pivot: 0, left: true, transparencyRate: 0.5});
    }
    updatedX += 4;
  }
  g.restore();
}

interface IComparedSequences {
  subParts1: string[];
  subParts2: string[];
}

function createDifferenceCanvas(amino1: string, amino2: string, color1: string, color2: string,
  y: number, shift: number, h: number
): HTMLCanvasElement {
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
  printLeftOrCentered(context, amino1, 0, y - shift, width, h, {color: color1, pivot: 0, left: true});
  printLeftOrCentered(context, amino2, 0, y + shift, width, h, {color: color2, pivot: 0, left: true});
  return canvas;
}

function fillShorterSequence(subParts1: string[], subParts2: string[]): IComparedSequences {
  let numIdenticalStart = 0;
  let numIdenticalEnd = 0;
  const longerSeq = subParts1.length > subParts2.length ? subParts1 : subParts2;
  const shorterSeq = subParts1.length > subParts2.length ? subParts2 : subParts1;

  for (let i = 0; i < shorterSeq.length; i++) {
    if (longerSeq[i] === shorterSeq[i])
      numIdenticalStart++;
  }

  const lengthDiff = longerSeq.length - shorterSeq.length;
  for (let i = longerSeq.length - 1; i > lengthDiff; i--) {
    if (longerSeq[i] === shorterSeq[i - lengthDiff])
      numIdenticalEnd++;
  }

  const emptyMonomersArray = new Array<string>(Math.abs(subParts1.length - subParts2.length))
    .fill(GapOriginals[NOTATION.FASTA]);

  function concatWithEmptyVals(subparts: string[]): string[] {
    return numIdenticalStart > numIdenticalEnd ?
      subparts.concat(emptyMonomersArray) : emptyMonomersArray.concat(subparts);
  }

  subParts1.length > subParts2.length ?
    subParts2 = concatWithEmptyVals(wu(subParts2).toArray()) : subParts1 = concatWithEmptyVals(wu(subParts1).toArray());
  return {subParts1: subParts1, subParts2: subParts2};
}
