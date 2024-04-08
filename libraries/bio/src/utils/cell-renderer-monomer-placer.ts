import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import wu from 'wu';

import {SeqHandler} from './seq-handler';
import {MonomerToShortFunc, ALPHABET} from './macromolecule';
import {IMonomerLib, Monomer} from '../types';
import {HELM_POLYMER_TYPE} from './const';
import {SeqSplittedBase} from './macromolecule/types';

type MonomerPlacerProps = {
  seqHandler: SeqHandler,
  monomerLib?: IMonomerLib,
  monomerCharWidth: number, separatorWidth: number,
  monomerToShort: MonomerToShortFunc, monomerLengthLimit: number,
};

const polymerTypeMap = {
  [ALPHABET.DNA]: HELM_POLYMER_TYPE.RNA,
  [ALPHABET.RNA]: HELM_POLYMER_TYPE.RNA,
  [ALPHABET.PT]: HELM_POLYMER_TYPE.PEPTIDE,
  [ALPHABET.UN]: HELM_POLYMER_TYPE.PEPTIDE,
};

export class MonomerPlacer {
  private _monomerLengthList: number[][] | null = null;

  // width of separator symbol
  private separatorWidth = 5;
  private props: MonomerPlacerProps;
  private _rowsProcessed: DG.BitSet; // rows for which monomer lengths were processed

  private _updated: boolean = false;
  public get updated(): boolean { return this._updated; }

  public _monomerLengthMap: { [key: string]: TextMetrics } = {}; // caches the lengths to save time on g.measureText
  public _monomerStructureMap: { [key: string]: HTMLDivElement } = {}; // caches the atomic structures of monomers
  private readonly subs: Unsubscribable[] = [];

  /** View is required to subscribe and handle for data frame changes */
  constructor(
    public readonly grid: DG.Grid | null,
    public readonly col: DG.Column<string>,
    private readonly propsProvider: () => MonomerPlacerProps
  ) {
    this.props = this.propsProvider();
    this._rowsProcessed = DG.BitSet.create(this.col.length);
    if (this.grid) {
      // Changes handling is available only in with a view
      this.subs.push(col.dataFrame.onDataChanged.subscribe(() => {
        try {
          this.props = this.propsProvider();
          this._monomerLengthList = null;
          this._rowsProcessed = DG.BitSet.create(this.col.length);
        } catch (e) {
          console.error(e);
        }
      }));
      this.subs.push(grok.events.onViewRemoved.subscribe((view: DG.View) => {
        try {
          if (this.grid?.view?.id === view.id) this.destroy();
        } catch (e) {
          console.error(e);
        }
      }));
    }
  }

  private destroy() {
    for (const sub of this.subs) sub.unsubscribe();
  }

  /** Returns monomers lengths of the {@link rowIdx} and cumulative sums for borders, monomer places */
  public getCellMonomerLengths(rowIdx: number): [number[], number[]] {
    const res: number[] = this.props.seqHandler.isMsa() ? this.getCellMonomerLengthsForSeqMsa() :
      this.getCellMonomerLengthsForSeq(rowIdx);

    const resSum: number[] = new Array<number>(res.length + 1);
    resSum[0] = 5; // padding
    for (let pos: number = 1; pos < resSum.length; pos++)
      resSum[pos] = resSum[pos - 1] + res[pos - 1];
    return [res, resSum];
  }

  private getCellMonomerLengthsForSeq(rowIdx: number): number[] {
    if (this._monomerLengthList === null) {
      this._monomerLengthList = new Array(this.col.length).fill(null);
      this._updated = true;
    }

    let res: number[] = this._monomerLengthList[rowIdx];
    if (res === null) {
      const seqMonList: SeqSplittedBase = SeqHandler.forColumn(this.col).getSplitted(rowIdx).originals;
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqMonList.length);

      for (const [seqMonLabel, seqMonI] of wu.enumerate(seqMonList)) {
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.props.monomerLengthLimit);
        const separatorWidth = this.props.seqHandler.isSeparator() ? this.separatorWidth : this.props.separatorWidth;
        const seqMonWidth: number = separatorWidth + shortMon.length * this.props.monomerCharWidth;
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
    this._monomerLengthList[0] ??= new Array(0);
    const res = this._monomerLengthList[0];
    // const startIdx = Math.max(Math.floor((this.grid?.vertScroll.min ?? 0) - 10), 0);
    // const endIdx = Math.min(Math.ceil((this.grid?.vertScroll.max ?? 0) + 10), this.col.length);

    const {startIdx, endIdx} = (() => {
      try {
        if (this.grid && this.grid.dart) {
          return {
            startIdx: Math.max(Math.floor((this.grid?.vertScroll.min ?? 0) - 10), 0),
            endIdx: Math.min(Math.ceil((this.grid?.vertScroll.max ?? 0) + 10), this.col.length)
          };
        } else return {startIdx: 0, endIdx: Math.min(this.col.length, 10)};
      } catch (_e) {
        return {startIdx: 0, endIdx: Math.min(this.col.length, 10)};
      }
    })();

    for (let seqIdx = startIdx; seqIdx < endIdx; seqIdx++) {
      if (this._rowsProcessed.get(seqIdx))
        continue;
      const seqMonList: SeqSplittedBase = SeqHandler.forColumn(this.col).getSplitted(seqIdx).originals;
      if (seqMonList.length > res.length)
        res.push(...new Array<number>(seqMonList.length - res.length).fill(0));

      for (const [seqMonLabel, seqMonI] of wu.enumerate(seqMonList)) {
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.props.monomerLengthLimit);
        const seqMonWidth: number = this.props.separatorWidth + shortMon.length * this.props.monomerCharWidth;
        res[seqMonI] = Math.max(res[seqMonI] ?? 0, seqMonWidth);
      }
      this._updated = true;
    }
    return res; // first (and single) row of data
  }

  /** Returns seq position for pointer x */
  public getPosition(rowIdx: number, x: number): number | null {
    const [_monomerMaxLengthList, monomerMaxLengthSumList]: [number[], number[]] = this.getCellMonomerLengths(rowIdx);
    const sh = SeqHandler.forColumn(this.col);
    const seqMonList: string[] = wu(sh.getSplitted(rowIdx).originals).toArray();
    if (seqMonList.length === 0) return null;

    let iterationCount: number = 100;
    let left: number | null = null;
    let right = seqMonList.length;
    let found = false;
    let mid = 0;
    if (monomerMaxLengthSumList[0] <= x && x < monomerMaxLengthSumList.slice(-1)[0]) {
      while (!found) {
        mid = Math.floor((right + (left ?? 0)) / 2);
        if (x >= monomerMaxLengthSumList[mid] && x <= monomerMaxLengthSumList[mid + 1]) {
          left = mid;
          found = true;
        } else if (x < monomerMaxLengthSumList[mid])
          right = mid - 1;
        else if (x > monomerMaxLengthSumList[mid + 1])
          left = mid + 1;

        if (left == right)
          found = true;

        if (--iterationCount <= 0) {
          throw new Error(`Get position for pointer x = ${x} searching has not converged ` +
            `on ${JSON.stringify(monomerMaxLengthSumList)}. `);
        }
      }
    }
    return left;
  }

  public getMonomer(symbol: string): Monomer | null {
    const alphabet = this.props.seqHandler.alphabet ?? ALPHABET.UN;
    const polymerType = polymerTypeMap[alphabet as ALPHABET];
    return this.props.monomerLib?.getMonomer(polymerType, symbol) ?? null;
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
