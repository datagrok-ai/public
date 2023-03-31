import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {UnitsHandler} from './units-handler';
import {getSplitterForColumn, MonomerToShortFunc, NOTATION, SplitterFunc} from './macromolecule';

type MonomerPlacerProps = {
  unitsHandler: UnitsHandler,
  monomerCharWidth: number, separatorWidth: number,
  monomerToShort: MonomerToShortFunc, monomerLengthLimit: number,
};

export class MonomerPlacer {
  private readonly _splitter: SplitterFunc;
  private _monomerLengthList: number[][] | null = null;

  private props: MonomerPlacerProps;

  private _updated: boolean = false;
  public get updated(): boolean { return this._updated; }

  private readonly subs: Unsubscribable[] = [];

  /** View is required to subscribe and handle for data frame changes */
  constructor(
    public readonly view: DG.View | null,
    public readonly col: DG.Column<string>,
    private readonly propsProvider: () => MonomerPlacerProps
  ) {
    this._splitter = getSplitterForColumn(this.col);
    this.props = this.propsProvider();

    if (this.view) {
      // Changes handling is available only in with a view
      this.subs.push(col.dataFrame.onDataChanged.subscribe(() => {
        this.props = this.propsProvider();
        this._monomerLengthList = null;
      }));
      this.subs.push(grok.events.onViewRemoved.subscribe((view: DG.View) => {
        if (this.view!.id === view.id) this.destroy();
      }));
    }
  }

  private destroy() {
    for (const sub of this.subs) sub.unsubscribe();
  }

  /** Returns monomers lengths of the {@link rowIdx} and cumulative sums for borders, monomer places */
  public getCellMonomerLengths(rowIdx: number): [number[], number[]] {
    const res: number[] = this.props.unitsHandler.isMsa() ? this.getCellMonomerLengthsForSeqMsa() :
      this.getCellMonomerLengthsForSeq(rowIdx);

    const resSum: number[] = new Array<number>(res.length + 1);
    resSum[0] = 0;
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
      const seqMonList: string[] = this.getSeqMonList(rowIdx);
      res = this._monomerLengthList[rowIdx] = new Array<number>(seqMonList.length);

      for (const [seqMonI, seqMonLabel] of seqMonList.entries()) {
        const shortMon: string = this.props.monomerToShort(seqMonLabel, this.props.monomerLengthLimit);
        const seqMonWidth: number = this.props.separatorWidth + shortMon.length * this.props.monomerCharWidth;
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
          const shortMon: string = this.props.monomerToShort(seqMonLabel, this.props.monomerLengthLimit);
          const seqMonWidth: number = this.props.separatorWidth + shortMon.length * this.props.monomerCharWidth;
          res[seqMonI] = Math.max(res[seqMonI] ?? 0, seqMonWidth);
        }
      }
      this._updated = true;
    }
    return res; // first (and single) row of data
  }

  /** Returns seq position for pointer x */
  public getPosition(rowIdx: number, x: number): number | null {
    const [monomerMaxLengthList, monomerMaxLengthSumList]: [number[], number[]] = this.getCellMonomerLengths(rowIdx);
    const seq: string = this.col.get(rowIdx)!;
    const seqMonList: string[] = this._splitter(seq);

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
        } else if (x < monomerMaxLengthSumList[mid]) {
          right = mid - 1;
        } else if (x > monomerMaxLengthSumList[mid + 1]) {
          left = mid + 1;
        }
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

  getSeqMonList(rowIdx: number): string[] {
    const seq: string | null = this.col.get(rowIdx);
    return seq ? this._splitter(seq) : [];
  }
}
