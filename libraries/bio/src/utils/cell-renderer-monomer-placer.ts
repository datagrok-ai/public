import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {getSeqHelper, ISeqHelper} from './seq-helper';
import {ALPHABET, MonomerToShortFunc, NOTATION, SplitterFunc, TAGS as bioTAGS} from './macromolecule';
import {ISeqSplitted} from './macromolecule/types';
import {CellRendererBackBase, getGridCellColTemp} from './cell-renderer-back-base';
import {ILogger} from './logger';
import {getMonomerLibHelper} from '../monomer-works/monomer-utils';
import {errInfo} from './err-info';
import {DrawStyle, printLeftOrCentered, TAGS as mmcrTAGS, PrintOptions} from './cell-renderer';
import {MmcrTemps, rendererSettingsChangedState, tempTAGS} from './cell-renderer-consts';
import {IMonomerLibBase} from '../types/index';
import {HelmTypes} from '../helm/consts';
import {ISeqMonomer} from '../helm/types';
import {execMonomerHoverLinks} from '../monomer-works/monomer-hover';

type MonomerPlacerProps = {
  separatorWidth: number,
  monomerToShort: MonomerToShortFunc,
  font: string, fontCharWidth: number,
};

export const undefinedColor = 'rgb(100,100,100)';

export function hitBounds(bounds: number[], x: number): number | null {
  let iterationCount: number = 100;
  let leftI: number = 0;
  let rightI = bounds.length - 1;
  let midI;
  while (leftI <= rightI) {
    midI = Math.floor((rightI + leftI) / 2);
    if (bounds[midI] <= x && x < bounds[midI + 1])
      return midI;
    else if (x < bounds[midI])
      rightI = midI - 1;
    else /* if (bounds[midI + 1] <= x) */
      leftI = midI + 1;

    if (--iterationCount <= 0) {
      throw new Error(`Get position for pointer x = ${x} searching has not converged ` +
        `on ${JSON.stringify(bounds)}. `);
    }
  }
  return null;
}

export class MonomerPlacer extends CellRendererBackBase<string> {
  private colWidth: number = 0;
  private _monomerLengthList: number[][] | null = null;

  // width of separator symbol
  private get separatorWidth() { return this.props?.fontCharWidth ? this.props?.fontCharWidth : 5; }
  public props: MonomerPlacerProps;
  private _processedRows: DG.BitSet; // rows for which monomer lengths were processed
  private _processedMaxVisibleSeqLength: number = 0;

  public _monomerLengthMap: { [key: string]: TextMetrics } = {}; // caches the lengths to save time on g.measureText
  public _monomerStructureMap: { [key: string]: HTMLElement } = {}; // caches the atomic structures of monomers

  private seqHelper!: ISeqHelper;

  private sysMonomerLib: IMonomerLibBase | null = null;

  /** View is required to subscribe and handle for data frame changes */
  constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>,
    logger: ILogger,
    public monomerLengthLimit: number,
    private readonly propsProvider: () => MonomerPlacerProps,
  ) {
    super(gridCol, tableCol, logger);
    this.props = this.propsProvider();
    this._processedRows = DG.BitSet.create(this.tableCol!.length);

    if (this.gridCol) {
      this.subs.push(this.gridCol.grid.onAfterDrawContent.subscribe(() => {
        this._onRendered.next();
      }));
    }

    if (this.tableCol && this.gridCol) {
      this.subs.push(this.tableCol.dataFrame.onCurrentRowChanged.subscribe(() => {
        const df = this.tableCol.dataFrame;
        const grid = this.gridCol!.grid;
        if (df.currentRowIdx === -1) {
          this.tableCol.temp[tempTAGS.referenceSequence] = null;
          this.tableCol.temp[tempTAGS.currentWord] = null;
          this.invalidateGrid();
        }
      }));
    }
  }

  public async init(): Promise<void> {
    await Promise.all([
      (async () => {
        this.seqHelper = await getSeqHelper();
        this.invalidateGrid();
      })(),
      (async () => {
        const libHelper = await getMonomerLibHelper();
        this.sysMonomerLib = libHelper.getMonomerLib();
      })(),
    ]);

    this.subs.push(this.sysMonomerLib!.onChanged.subscribe(() => {
      this.reset();
    }));

    this.reset();
  }

  public static getFontSettings(tableCol?: DG.Column) {
    const fontFamily = 'monospace';
    let fontSize = 12;
    if (tableCol && tableCol.temp[MmcrTemps.fontSize] &&
        typeof tableCol.temp[MmcrTemps.fontSize] === 'number' && !isNaN(tableCol.temp[MmcrTemps.fontSize]))
      fontSize = Math.max(tableCol.temp[MmcrTemps.fontSize], 1);
    // monospace ratio
    const fontWidth = fontSize * 0.6;
    return {font: `${fontSize}px ${fontFamily}`, fontWidth};
  }

  override toLog(): string {
    return `MonomerPlacer<${this.viewerId}>`;
  }

  protected getMonomerLib(): IMonomerLibBase | null {
    return this.tableCol.temp[MmcrTemps.overriddenLibrary] ?? this.sysMonomerLib;
  }

  protected override reset(): void {
    if (this.propsProvider)
      this.props = this.propsProvider();
    this._processedRows = DG.BitSet.create(this.tableCol!.length);
    this._monomerLengthList = null;
    this._monomerLengthMap = {};
    this._monomerStructureMap = {};
    super.reset();
    this.invalidateGrid();
  }

  protected invalidateGrid(): void {
    if (this.gridCol && this.gridCol.dart) this.gridCol.grid?.invalidate();
  }

  /** Returns monomers lengths of the {@link rowIdx} and cumulative sums for borders, monomer places */
  public getCellMonomerLengths(rowIdx: number, newWidth: number): [number[], number[]] {
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    if (this.colWidth < newWidth) {
      this.colWidth = newWidth;
      this.dirty = true;
    }
    if (this.dirty) {
      try { this.reset(); } catch (err) {
        const [errMsg, errStack] = errInfo(err);
        this.logger.error(errMsg, undefined, errStack);
      }
    }

    const res: number[] = sh.isMsa() ? this.getCellMonomerLengthsForSeqMsa() :
      this.getCellMonomerLengthsForSeq(rowIdx);

    const resSum: number[] = new Array<number>(res.length + 1);
    resSum[0] = 5; // padding
    for (let pos: number = 1; pos < resSum.length; pos++)
      resSum[pos] = resSum[pos - 1] + res[pos - 1];
    return [res, resSum];
  }

  private getCellMonomerLengthsForSeq(rowIdx: number): number[] {
    const logPrefix = `${this.toLog()}.getCellMonomerLengthsForSeqMsa()`;
    // this.logger.debug(`${logPrefix}, start`);

    if (this._monomerLengthList === null)
      this._monomerLengthList = new Array(this.tableCol.length).fill(null);

    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const maxVisibleSeqLength: number = Math.ceil(this.colWidth / minMonWidth);
    const seqSS: ISeqSplitted = sh.getSplitted(rowIdx);
    const visibleSeqLength: number = Math.min(maxVisibleSeqLength, seqSS.length);

    let res: number[] = this._monomerLengthList[rowIdx];
    if (res === null || res.length < visibleSeqLength) {
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqSS.length);

      let seqWidth: number = 0;
      for (let seqMonI = 0; seqMonI < visibleSeqLength; ++seqMonI) {
        const seqMonLabel = seqSS.getOriginal(seqMonI);
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.monomerLengthLimit);
        const separatorWidth = sh.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
        const seqMonWidth: number = separatorWidth + shortMon.length * this.props.fontCharWidth;
        res[seqMonI] = seqMonWidth;
        seqWidth += seqMonWidth;
        if (seqWidth > this.colWidth) break;
      }
    }
    return res;
  }

  private getCellMonomerLengthsForSeqMsa(): number[] {
    const logPrefix = `${this.toLog()}.getCellMonomerLengthsForSeqMsa()`;
    // this.logger.debug(`${logPrefix}, start`);
    if (this._monomerLengthList === null)
      this._monomerLengthList = new Array(1).fill(null);
    this._monomerLengthList[0] ??= new Array(0);
    const res = this._monomerLengthList[0];
    // const startIdx = Math.max(Math.floor((this.grid?.vertScroll.min ?? 0) - 10), 0);
    // const endIdx = Math.min(Math.ceil((this.grid?.vertScroll.max ?? 0) + 10), this.col.length);

    const {startIdx, endIdx} = (() => {
      try {
        const grid: DG.Grid | null = this.gridCol && this.gridCol.dart ? this.gridCol.grid : null;
        if (grid && grid.dart) {
          return {
            startIdx: Math.max(Math.floor((grid?.vertScroll.min ?? 0) - 10), 0),
            endIdx: Math.min(Math.ceil((grid?.vertScroll.max ?? 0) + 10), this.tableCol.length)
          };
        } else return {startIdx: 0, endIdx: Math.min(this.tableCol.length, 10)};
      } catch (_e) {
        return {startIdx: 0, endIdx: Math.min(this.tableCol.length, 10)};
      }
    })();

    const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
    const maxVisibleSeqLength: number = Math.ceil(this.colWidth / minMonWidth);
    for (let seqIdx = startIdx; seqIdx < endIdx; seqIdx++) {
      if (this._processedRows.get(seqIdx) && maxVisibleSeqLength <= this._processedMaxVisibleSeqLength)
        continue;
      const sh = this.seqHelper.getSeqHandler(this.tableCol);
      const seqSS: ISeqSplitted = sh.getSplitted(seqIdx, maxVisibleSeqLength);

      const visibleSeqLength: number = Math.min(maxVisibleSeqLength, seqSS.length);
      if (visibleSeqLength > res.length)
        res.push(...new Array<number>(visibleSeqLength - res.length).fill(minMonWidth));

      let seqWidth: number = 0;
      for (let seqMonI = 0; seqMonI < visibleSeqLength; ++seqMonI) {
        const seqMonLabel: string = seqSS.getOriginal(seqMonI);
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.monomerLengthLimit);
        const seqMonWidth: number = this.props.separatorWidth + shortMon.length * this.props.fontCharWidth;
        res[seqMonI] = Math.max(res[seqMonI] ?? 0, seqMonWidth);
        seqWidth += seqMonWidth;
        if (seqWidth >= this.colWidth) break;
      }
      this._processedMaxVisibleSeqLength = Math.max(this._processedMaxVisibleSeqLength, maxVisibleSeqLength);
    }
    return res; // first (and single) row of data
  }

  /** Returns seq position for pointer x */
  public getPosition(rowIdx: number, x: number, width: number): number | null {
    const [_monomerMaxLengthList, monomerMaxLengthSumList]: [number[], number[]] =
      this.getCellMonomerLengths(rowIdx, width);
    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const seqSS = sh.getSplitted(rowIdx);
    if (seqSS.length === 0) return null;
    return hitBounds(monomerMaxLengthSumList, x);
  }

  public setMonomerLengthLimit(value: number): void {
    if (this.monomerLengthLimit != value) {
      this.monomerLengthLimit = value;
      this.dirty = true;
    }
  }

  public setSeparatorWidth(value: number): void {
    if (this.separatorWidth != value) {
      this.props.separatorWidth = value;
      this.dirty = true;
    }
  }

  private padding: number = 5;

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle
  ) {
    if (!this.seqHelper) return;
    const gridCol = this.gridCol;
    const tableCol = this.tableCol;
    const dpr = window.devicePixelRatio;

    const logPrefix = `${this.toLog()}.render()`;
    this.logger.debug(`${logPrefix}, start`);

    // Cell renderer settings
    let maxLengthOfMonomer: number = this.monomerLengthLimit;
    if (mmcrTAGS.maxMonomerLength in tableCol.tags) {
      const v = parseInt(tableCol.getTag(mmcrTAGS.maxMonomerLength));
      maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
    }
    if (MmcrTemps.maxMonomerLength in tableCol.temp) {
      const v = tableCol.temp[MmcrTemps.maxMonomerLength];
      maxLengthOfMonomer = !isNaN(v) && v ? v : 50;
    }

    g.save();
    try {
      const sh = this.seqHelper.getSeqHandler(tableCol);

      if (
        tableCol.temp[MmcrTemps.rendererSettingsChanged] === rendererSettingsChangedState.true ||
        this.monomerLengthLimit != maxLengthOfMonomer
      ) {
        let gapLength = 0;
        const msaGapLength = 8;
        gapLength = tableCol.temp[MmcrTemps.gapLength] as number ?? gapLength;
        // this event means that the mm renderer settings have changed,
        // particularly monomer representation and max width.
        this.setMonomerLengthLimit(maxLengthOfMonomer);
        this.setSeparatorWidth(sh.isMsa() ? msaGapLength : gapLength);
        tableCol.temp[MmcrTemps.rendererSettingsChanged] = rendererSettingsChangedState.false;
      }

      const [maxLengthWords, maxLengthWordsSum]: [number[], number[]] =
        this.getCellMonomerLengths(gridCell.tableRowIndex!, w);
      const _maxIndex = maxLengthWords.length;

      const value: any = gridCell.cell.value;
      const rowIdx = gridCell.cell.rowIndex;
      const paletteType = tableCol.getTag(bioTAGS.alphabet);
      const minDistanceRenderer = 50;
      w = getUpdatedWidth(gridCol?.grid, g, x, w, dpr);
      g.beginPath();
      g.rect(x + this.padding, y + this.padding, w - this.padding - 1, h - this.padding * 2);
      g.clip();
      g.font = this.props?.font ?? '12px monospace';
      g.textBaseline = 'top';

      //TODO: can this be replaced/merged with splitSequence?
      const units = tableCol.meta.units;
      const aligned: string = tableCol.getTag(bioTAGS.aligned);

      const separator = tableCol.getTag(bioTAGS.separator) ?? '';
      const minMonWidth = this.props.separatorWidth + 1 * this.props.fontCharWidth;
      const splitLimit = Math.ceil(w / minMonWidth);

      const tempReferenceSequence: string | null = tableCol.temp[tempTAGS.referenceSequence];
      const tempCurrentWord: string | null = this.tableCol.temp[tempTAGS.currentWord];
      if (tempCurrentWord && tableCol?.dataFrame?.currentRowIdx === -1)
        this.tableCol.temp[tempTAGS.currentWord] = null;

      const referenceSequence: string[] = (() => {
        // @ts-ignore
        const splitterFunc: SplitterFunc = sh.getSplitter(splitLimit);
        const seqSS = splitterFunc(
          ((tempReferenceSequence != null) && (tempReferenceSequence != '')) ?
            tempReferenceSequence : tempCurrentWord ?? '');
        return wu.count(0).take(seqSS.length).map((posIdx) => seqSS.getOriginal(posIdx)).toArray();
      })();

      const subParts: ISeqSplitted = sh.getSplitted(rowIdx);
      let drawStyle = DrawStyle.classic;

      if (aligned && aligned.includes('MSA') && units == NOTATION.SEPARATOR)
        drawStyle = DrawStyle.MSA;

      const visibleSeqLength = Math.min(subParts.length, splitLimit);
      for (let posIdx: number = 0; posIdx < visibleSeqLength; ++posIdx) {
        const om: string = subParts.getOriginal(posIdx);
        const cm: string = subParts.getCanonical(posIdx);

        let color = undefinedColor;
        const monomerLib = this.getMonomerLib();
        if (monomerLib) {
          const biotype = sh.defaultBiotype;
          //this.logger.debug(`${logPrefix}, biotype: ${biotype}, amino: ${amino}`);
          color = monomerLib.getMonomerTextColor(biotype, cm);
        }
        g.fillStyle = undefinedColor;
        const last = posIdx === subParts.length - 1;
        /*x1 = */
        const opts: Partial<PrintOptions> = {
          color: color, pivot: 0, left: true, transparencyRate: 1.0, separator: separator, last: last,
          drawStyle: drawStyle, maxWord: maxLengthWordsSum, wordIdx: posIdx, gridCell: gridCell,
          referenceSequence: referenceSequence, maxLengthOfMonomer: maxLengthOfMonomer,
          monomerTextSizeMap: this._monomerLengthMap, logger: this.logger
        };
        printLeftOrCentered(g, om, x + this.padding, y, w, h, opts);
        if (minDistanceRenderer > w) break;
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(errMsg, undefined, errStack);
      this.errors.push(err);
      //throw err; // Do not throw to prevent disabling renderer
    } finally {
      g.restore();
    }
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const logPrefix = `${this.toLog()}.onMouseMove()`;
    if (!this.seqHelper || gridCell.tableRowIndex == null) return;

    // if (gridCell.cell.column.getTag(bioTAGS.aligned) !== ALIGNMENT.SEQ_MSA)
    //   return;

    const gridCellBounds: DG.Rect = gridCell.bounds;
    // const value: any = gridCell.cell.value;
    //
    // const maxLengthWords: number[] = seqColTemp.getCellMonomerLengths(gridCell.tableRowIndex!);
    // const maxLengthWordsSum: number[] = new Array<number>(maxLengthWords.length).fill(0);
    // for (let posI: number = 1; posI < maxLengthWords.length; posI++)
    //   maxLengthWordsSum[posI] = maxLengthWordsSum[posI - 1] + maxLengthWords[posI];
    // const maxIndex = maxLengthWords.length;
    const argsX = e.offsetX - gridCell.gridColumn.left + (gridCell.gridColumn.left - gridCellBounds.x);
    const left: number | null = this.getPosition(gridCell.tableRowIndex!, argsX, gridCellBounds.width);
    this.logger.debug(`${logPrefix}, start, argsX: ${argsX}, left: ${left}`);

    const sh = this.seqHelper.getSeqHandler(this.tableCol);
    const seqSS = sh.getSplitted(gridCell.tableRowIndex!);
    if (left !== null && left < seqSS.length) {
      const alphabet = sh.alphabet ?? ALPHABET.UN;
      const seqMonomer = {
        position: left,
        biotype: alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA,
        symbol: seqSS.getCanonical(left),
      } as ISeqMonomer;
      const tooltipElements: HTMLElement[] = [];
      let monomerDiv = this._monomerStructureMap[seqMonomer.symbol];
      if (!monomerDiv || true) {
        const monomerLib = this.getMonomerLib();
        monomerDiv = this._monomerStructureMap[seqMonomer.symbol] = (() => {
          return monomerLib ? monomerLib.getTooltip(seqMonomer.biotype, seqMonomer.symbol) :
            ui.divText('Monomer library is not available');
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
}

export function getUpdatedWidth(
  grid: DG.Grid | null | undefined, g: CanvasRenderingContext2D, x: number, w: number, dpr: number
): number {
  return !!grid ? Math.max(Math.min(grid.canvas.width / dpr - x, w)) : Math.max(g.canvas.width / dpr - x, 0);
}
