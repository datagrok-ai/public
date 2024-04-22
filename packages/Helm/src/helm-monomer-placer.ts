import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {CellRendererBackBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';

import {getParts, parseHelm} from './utils';

import {_package, getMonomerLib} from './package';

export interface IEditor {
  get m(): IEditorMol;

  resize(width: number, height: number): void;
  setData(data: string, format: string): void;
}

export interface IEditorMol {
  get atoms(): IEditorMolAtom[];

  clone(selectedOnly: boolean): IEditorMol;
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
  polymerType: string
  symbol: string,
}

export class HelmMonomerPlacer extends CellRendererBackBase<string> {
  public readonly monomerLib: IMonomerLib;

  private _allPartsList: (string[] | null)[];
  private _lengthsList: (number[] | null)[];
  private _editorMolList: (IEditorMol | null)[];

  public monomerCharWidth: number = 7;
  public leftPadding: number = 5;
  public monomerTextSizeMap: { [p: string]: TextMetrics } = {};

  public constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>
  ) {
    super(gridCol, tableCol, _package.logger);
    this.monomerLib = getMonomerLib();
    this.subs.push(this.monomerLib.onChanged.subscribe(this.monomerLibOnChanged.bind(this)));
  }

  protected override reset(): void {
    this._allPartsList = new Array<string[] | null>(this.tableCol.length).fill(null);
    this._lengthsList = new Array<number[] | null>(this.tableCol.length).fill(null);
    this._editorMolList = new Array<IEditorMol | null>(this.tableCol.length).fill(null);
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
    if (allParts === null) {
      const seq = this.tableCol.get(rowIdx);
      allParts = this._allPartsList![rowIdx] = getAllParts(seq);
    }

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
    this.reset();
    this.invalidateGrid();
  }

  public setEditorMol(tableRowIndex: number, editorMol: IEditorMol) {
    this._editorMolList![tableRowIndex] = editorMol.clone(false);
  }

  public getEditorMol(tableRowIndex: number): IEditorMol | null {
    return this._editorMolList![tableRowIndex];
  }
}

export function getAllParts(seq: string | null): string[] {
  const monomerList: string[] = !!seq ? parseHelm(seq) : [];
  return !!seq ? getParts(monomerList, seq) : [];
}
