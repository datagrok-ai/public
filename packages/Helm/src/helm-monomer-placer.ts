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

export interface IEditor {
  get m(): IEditorMol;

  redraw(): void;
}

export interface IEditorMol {
  get atoms(): IEditorMolAtom[];
}

export interface IEditorMolAtom {
  get p(): IEditorPoint;

  get elem(): string;
}

export interface IEditorPoint {
  get x(): number;

  get y(): number;
}

export interface ISeqMonomer {
  symbol: string,
  polymerType?: string
}

export class HelmMonomerPlacer {
  private readonly _allPartsList: (string[] | null)[];
  private readonly _lengthsList: (number[] | null)[];
  private readonly _editorList: (IEditor | null)[];

  public readonly monomerLib: IMonomerLib;

  public monomerCharWidth: number = 7;
  public leftPadding: number = 5;
  public monomerTextSizeMap: { [p: string]: TextMetrics } = {};

  constructor(public readonly col: DG.Column<string>) {
    this.col.dataFrame.onDataChanged.subscribe();
    this.monomerLib = getMonomerLib();
    this.monomerLib.onChanged.subscribe(this.monomerLibOnChanged.bind(this));

    this._allPartsList = new Array<string[] | null>(this.col.length).fill(null);
    this._lengthsList = new Array<number[] | null>(this.col.length).fill(null);
    this._editorList = new Array<IEditor | null>(this.col.length).fill(null);
  }

  /** Skips cell for the fallback rendering */
  public skipCell(rowIdx: number): void {
    this._allPartsList[rowIdx] = [];
    this._lengthsList[rowIdx] = [];
  }

  /** @param rowIdx Row index of the table {@link DG.DataFrame}, HelmMonomerPlacer is {@link DG.Column} based */
  public getCellAllPartsLengths(rowIdx: number): [string[], number[], number[]] {
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

  getMonomer(monomer: ISeqMonomer): Monomer | null {
    let res: Monomer | null = null;
    if (monomer.polymerType)
      res = this.monomerLib.getMonomer(monomer.polymerType, monomer.symbol);
    else {
      for (const polymerType of this.monomerLib.getPolymerTypes()) {
        res = this.monomerLib.getMonomer(polymerType, monomer.symbol);
        if (res) break;
      }
    }
    return res;
  }

  // -- Handle events --

  private monomerLibOnChanged(_value: any): void {
    this._lengthsList.fill(null);
    this._allPartsList.fill(null);
    this._editorList.fill(null);
    // TODO: Invalidate all grids of this.col.dataFrame
  }

  public static getOrCreate(col: DG.Column<string>): HelmMonomerPlacer {
    if (!(Temps.helmMonomerPlacer in col.temp)) col.temp[Temps.helmMonomerPlacer] = new HelmMonomerPlacer(col);
    return col.temp[Temps.helmMonomerPlacer];
  }

  public setEditor(tableRowIndex: number, editor: IEditor) {
    this._editorList![tableRowIndex] = editor;
  }

  public getEditor(tableRowIndex: number): IEditor | null {
    return this._editorList![tableRowIndex];
  }
}
