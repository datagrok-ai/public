import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';

import {getParts, parseHelm} from './utils';
import {getMonomerLib} from './package';

export const enum Temps {
  helmMonomerPlacer = 'bio-helmMonomerPlacer',
}

export class HelmMonomerPlacer {
  private _allPartsList: (string[] | null)[] | null = null;
  private _lengthsList: (number[] | null) [] | null = null;

  private monomerLib: IMonomerLib;

  public monomerCharWidth: number = 7;
  public leftPadding: number = 5;

  constructor(public readonly col: DG.Column<string>) {
    this.col.dataFrame.onDataChanged.subscribe();
    this.monomerLib = getMonomerLib();
    this.monomerLib.onChanged.subscribe(this.monomerLibOnChanged.bind(this));
  }

  public skipCell(rowIdx: number): void {
    if (this._allPartsList === null)
      this._allPartsList = new Array<string[] | null>(this.col.length).fill(null);
    if (this._lengthsList === null)
      this._lengthsList = new Array<number[] | null>(this.col.length).fill(null);

    this._allPartsList[rowIdx] = [];
    this._lengthsList[rowIdx] = [];
  }

  /** @param rowIdx Row index of the table {@link DG.DataFrame}, HelmMonomerPlacer is {@link DG.Column} based */
  public getCellAllPartsLengths(rowIdx: number): [string[], number[], number[]] {
    if (this._allPartsList === null)
      this._allPartsList = new Array<string[] | null>(this.col.length).fill(null);
    if (this._lengthsList === null)
      this._lengthsList = new Array<number[] | null>(this.col.length).fill(null);

    const [allParts, lengths] = this.getCellMonomerLengthsForSeq(rowIdx);

    const sumLengths: number[] = new Array<number>(lengths.length + 1);
    sumLengths[0] = this.leftPadding; // padding
    for (let pos: number = 1; pos < sumLengths.length; pos++)
      sumLengths[pos] = sumLengths[pos - 1] + lengths[pos - 1];
    return [allParts, lengths, sumLengths];
  }

  private getCellMonomerLengthsForSeq(rowIdx: number): [string[], number[]] {
    let allParts: string[] | null = this._allPartsList![rowIdx];
    if (allParts === null)
      allParts = this._allPartsList![rowIdx] = this.getAllParts(rowIdx);

    let lengths: number[] | null = this._lengthsList![rowIdx];
    if (lengths === null) {
      lengths = this._lengthsList![rowIdx] = new Array<number>(allParts.length);
      for (const [part, partI] of wu.enumerate(allParts)) {
        const partWidth: number = part.length * this.monomerCharWidth;
        lengths[partI] = partWidth;
      }
    }

    return [allParts, lengths];
  }

  getAllParts(rowIdx: number): string[] {
    const seq: string | null = this.col.get(rowIdx);
    const monomerList: string[] = seq ? parseHelm(seq) : [];
    return seq ? getParts(monomerList, seq) : [];
  }

  getMonomer(monomerSymbol: any): Monomer | null {
    let res: Monomer | null = null;
    for (const polymerType of this.monomerLib.getPolymerTypes()) {
      res = this.monomerLib.getMonomer(polymerType, monomerSymbol);
      if (res) break;
    }
    return res;
  }

  // -- Handle events --

  private monomerLibOnChanged(_value: any): void {
    this._lengthsList = null;
    this._allPartsList = null;
    // TODO: Invalidate all grids of this.col.dataFrame
  }

  public static getOrCreate(col: DG.Column<string>): HelmMonomerPlacer {
    if (!(Temps.helmMonomerPlacer in col.temp)) col.temp[Temps.helmMonomerPlacer] = new HelmMonomerPlacer(col);
    return col.temp[Temps.helmMonomerPlacer];
  }
}
