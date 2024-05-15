import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import wu from 'wu';

import {SeqHandler} from './seq-handler';
import {MonomerToShortFunc} from './macromolecule';
import {IMonomerLib} from '../types';
import {SeqSplittedBase} from './macromolecule/types';
import {CellRendererBackBase} from './cell-renderer-back-base';
import {ILogger} from './logger';
import {getMonomerLibHelper} from '../monomer-works/monomer-utils';

type MonomerPlacerProps = {
  seqHandler: SeqHandler,
  monomerCharWidth: number, separatorWidth: number,
  monomerToShort: MonomerToShortFunc, monomerLengthLimit: number,
};

export function hitBounds(bounds: number[], x: number): number | null {
  let iterationCount: number = 100;
  let leftI: number = 0;
  let rightI = bounds.length - 1;
  let midI;
  while (leftI <= rightI) {
    midI = Math.floor((rightI + leftI) / 2);
    const logVal: {} = {x, leftI, left: bounds[leftI], midI, mid: bounds[midI], rightI, right: bounds[rightI]};
    console.debug(`bio lib: hitTestBounds(), ${JSON.stringify(logVal)}`);
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
  private _monomerLengthList: number[][] | null = null;

  // width of separator symbol
  private separatorWidth = 5;
  public props: MonomerPlacerProps;
  private _processedRows: DG.BitSet; // rows for which monomer lengths were processed
  private _processedMaxVisibleSeqLength: number = 0;

  private _updated: boolean = false;
  public get updated(): boolean { return this._updated; }

  public _monomerLengthMap: { [key: string]: TextMetrics } = {}; // caches the lengths to save time on g.measureText
  public _monomerStructureMap: { [key: string]: HTMLElement } = {}; // caches the atomic structures of monomers

  /** View is required to subscribe and handle for data frame changes */
  constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>,
    logger: ILogger,
    private readonly propsProvider: () => MonomerPlacerProps
  ) {
    super(gridCol, tableCol, logger);
    this.props = this.propsProvider();
    this._processedRows = DG.BitSet.create(this.tableCol!.length);

    if (this.gridCol) {
      this.subs.push(this.gridCol.grid.onAfterDrawContent.subscribe(() => {
        this._onRendered.next();
      }));
    }

    getMonomerLibHelper().then((libHelper) => {
      if (this.destroyed) return;
      const monomerLib = libHelper.getBioLib();
      this.subs.push(monomerLib.onChanged.subscribe(() => {
        this.reset();
        this.gridCol?.grid?.invalidate();
      }));
    });
  }

  protected override reset(): void {
    if (this.propsProvider) this.props = this.propsProvider();
    this._processedRows = DG.BitSet.create(this.tableCol!.length);
    this._monomerLengthList = null;
    this._monomerLengthMap = {};
    this._monomerStructureMap = {};
  }

  /** Returns monomers lengths of the {@link rowIdx} and cumulative sums for borders, monomer places */
  public getCellMonomerLengths(rowIdx: number, width: number): [number[], number[]] {
    const res: number[] = this.props.seqHandler.isMsa() ? this.getCellMonomerLengthsForSeqMsa(width) :
      this.getCellMonomerLengthsForSeq(rowIdx, width);

    const resSum: number[] = new Array<number>(res.length + 1);
    resSum[0] = 5; // padding
    for (let pos: number = 1; pos < resSum.length; pos++)
      resSum[pos] = resSum[pos - 1] + res[pos - 1];
    return [res, resSum];
  }

  private getCellMonomerLengthsForSeq(rowIdx: number, width: number): number[] {
    if (this._monomerLengthList === null) {
      this._monomerLengthList = new Array(this.tableCol.length).fill(null);
      this._updated = true;
    }

    const minMonWidth = this.props.separatorWidth + 1 * this.props.monomerCharWidth;
    const maxVisibleSeqLength: number = Math.ceil(width / minMonWidth);
    const seqMonList: SeqSplittedBase = SeqHandler.forColumn(this.tableCol).getSplitted(rowIdx).originals;
    const visibleSeqLength: number = Math.min(maxVisibleSeqLength, seqMonList.length);

    let res: number[] = this._monomerLengthList[rowIdx];
    if (res === null || res.length < visibleSeqLength) {
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqMonList.length);

      let seqWidth: number = 0;
      for (let seqMonI = 0; seqMonI < visibleSeqLength; ++seqMonI) {
        const seqMonLabel = seqMonList[seqMonI];
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.props.monomerLengthLimit);
        const separatorWidth = this.props.seqHandler.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
        const seqMonWidth: number = separatorWidth + shortMon.length * this.props.monomerCharWidth;
        res[seqMonI] = seqMonWidth;
        seqWidth += seqMonWidth;
        if (seqWidth > width) break;
      }
      this._updated = true;
    }
    return res;
  }

  private getCellMonomerLengthsForSeqMsa(width: number): number[] {
    if (this._monomerLengthList === null) {
      this._monomerLengthList = new Array(1).fill(null);
      this._updated = true;
    }
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

    const minMonWidth = this.props.separatorWidth + 1 * this.props.monomerCharWidth;
    const maxVisibleSeqLength: number = Math.ceil(width / minMonWidth);
    for (let seqIdx = startIdx; seqIdx < endIdx; seqIdx++) {
      if (this._processedRows.get(seqIdx) && maxVisibleSeqLength <= this._processedMaxVisibleSeqLength)
        continue;
      const seqMonList: SeqSplittedBase = SeqHandler.forColumn(this.tableCol)
        .getSplitted(seqIdx, maxVisibleSeqLength).originals;

      const visibleSeqLength: number = Math.min(maxVisibleSeqLength, seqMonList.length);
      if (visibleSeqLength > res.length)
        res.push(...new Array<number>(visibleSeqLength - res.length).fill(minMonWidth));

      let seqWidth: number = 0;
      for (let seqMonI = 0; seqMonI < visibleSeqLength; ++seqMonI) {
        const seqMonLabel: string = seqMonList[seqMonI];
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.props.monomerLengthLimit);
        const seqMonWidth: number = this.props.separatorWidth + shortMon.length * this.props.monomerCharWidth;
        res[seqMonI] = Math.max(res[seqMonI] ?? 0, seqMonWidth);
        seqWidth += seqMonWidth;
        if (seqWidth >= width) break;
      }
      this._updated = true;
      this._processedMaxVisibleSeqLength = Math.max(this._processedMaxVisibleSeqLength, maxVisibleSeqLength);
    }
    return res; // first (and single) row of data
  }

  /** Returns seq position for pointer x */
  public getPosition(rowIdx: number, x: number, width: number): number | null {
    const [_monomerMaxLengthList, monomerMaxLengthSumList]: [number[], number[]] =
      this.getCellMonomerLengths(rowIdx, width);
    const sh = SeqHandler.forColumn(this.tableCol);
    const seqMonList: string[] = wu(sh.getSplitted(rowIdx).originals).toArray();
    if (seqMonList.length === 0) return null;
    return hitBounds(monomerMaxLengthSumList, x);
  }

  public setMonomerLengthLimit(limit: number): void {
    this.props.monomerLengthLimit = limit;
    this._updated = true;
  }

  public setSeparatorWidth(width: number): void {
    this.props.separatorWidth = width;
    this._updated = true;
  }

  public isMsa(): boolean {
    return this.props.seqHandler.isMsa();
  }
}
