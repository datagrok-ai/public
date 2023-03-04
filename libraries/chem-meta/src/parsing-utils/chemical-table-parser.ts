/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type AtomAndBondCounts = {
  atomCount: number,
  bondCount: number
}

type CoordinateArrays = {
  x: number[],
  y: number[],
  z: number[],
}

/** Base singleton for Molfile and Mol2 parser/handler */
export abstract class ChemicalTableParserBase {
  protected constructor(file: string) {
    this.reset(file);
  };

  private static instance: ChemicalTableParserBase;

  public static getInstance<T extends ChemicalTableParserBase>(file: string): T {
    if (!this.instance)
      this.instance = new (ChemicalTableParserBase as any)(file); // a workaround to define an abstract singleton
    return this.instance as T;
  }

  public reset(file: string) {
    this.file = file.replaceAll('\r', '');
  }

  protected file!: string;
  /** Index running along the string/file being parsed  */
  protected _atomCount?: number;
  protected _bondCount?: number;
  /** The array of X, Y, Z arrays for atomic coordinates */
  protected atomCoordinates?: CoordinateArrays;
  protected _atomTypes?: string[];

  protected abstract parseAtomAndBondCounts(): AtomAndBondCounts;
  protected abstract getCountsLineIdx(): number;
  /** Get idx of the first line of the atom block containing atom data  */
  protected abstract getAtomBlockIdx(): number;
  /** Get idx of the first line of the bond block containing bond data */
  protected abstract getBondBlockIdx(): number;
  /** Shift idx from the beginning of the line to X coordinate */
  protected abstract shiftIdxToXColumn(lineStartIdx: number): number;
  /** Shift idx from the beginning of the line to atom type column */
  protected abstract shiftIdxToAtomType(lineStartIdx: number): number;
  /** Parse atom type at idx */
  protected abstract parseAtomType(idx: number): string;

  protected setAtomAndBondCounts(): void {
    const {atomCount, bondCount} = this.parseAtomAndBondCounts();
    this._atomCount = atomCount;
    this._bondCount = bondCount;
  }

  /** Gets the idx of the next column relatively to this._currentIdx */
  protected getNextColumnIdx(idx: number): number {
    // skip non-whitespace, if necessary
    while (!this.isWhitespace(idx))
      ++idx;
    // skip whitespace
    while (this.isWhitespace(idx))
      ++idx;
    return idx;
  }

  /** Shift idx from beginning of the specified line to the specified column  */
  protected shiftToSpecifiedColumn(lineStartIdx: number, columnNumber: number) {
    let idx = lineStartIdx;
    for (let i = 0; i < columnNumber; i++)
      idx = this.getNextColumnIdx(idx);
    return idx;
  }

  protected parseAtomTypes(): string[] {
    const atomTypes = new Array<string>(this.atomCount);
    let idx = this.getAtomBlockIdx();
    const atomCount = this.atomCount;
    for (let i = 0; i < atomCount; i++) {
      idx = this.shiftIdxToAtomType(idx);
      atomTypes[idx] = this.parseAtomType(idx);
      idx = this.getNextLineIdx(idx);
    }
    return atomTypes;
  };

  protected parseAtomCoordinates(): CoordinateArrays {
    const x = new Array<number>(this.atomCount);
    const y = new Array<number>(this.atomCount);
    const z = new Array<number>(this.atomCount);
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < this.atomCount; i++) {
      idx = this.shiftIdxToXColumn(idx);
      for (const item of [x, y, z]) {
        item[i] = this.parseFloatValue(idx);
        idx = this.getNextColumnIdx(idx);
      }
    }
    return {x: x, y: y, z: z};
  }

  /** Check if a character is whitespace including '\t'  */
  protected isWhitespace(idx: number): boolean {
    return /\s/.test(this.file.at(idx)!);
  }

  /** Get index of the next line starting from idx  */
  protected getNextLineIdx(idx: number): number {
    if (this.file.at(idx) !== '\n')
      return this.file.indexOf('\n', idx) + 1;
    else
      return this.file.indexOf('\n', idx + 1) + 1;
  }

  /** Get a float value in the current column (at idx) */
  protected parseFloatValue(idx: number): number {
    return this.parseNumericValue(parseFloat, idx);
  }

  /** Get an int value in the current column (at idx) */
  protected parseIntValue(idx: number): number {
    return this.parseNumericValue(parseInt, idx);
  }

  /** Parse a numeric value depending on the functional argument  */
  protected parseNumericValue(
    parserFunction: (str: string) => number,
    idx: number
  ): number {
    let end = idx + 1;
    while (!this.isWhitespace(end))
      ++end;
    const value = parserFunction(this.file.substring(idx, end));
    return value;
  }

  /** Number of atoms in a molecule  */
  get atomCount(): number {
    if (this._atomCount === undefined)
      this.setAtomAndBondCounts();
    return this._atomCount!;
  }

  /** Number of bonds in a molecule  */
  get bondCount(): number {
    if (this._bondCount === undefined)
      this.setAtomAndBondCounts();
    return this._bondCount!;
  }

  /** X coordinates of all atoms in a molecule  */
  get x(): number[] {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates.x;
  };

  /** Y coordinates of all atoms in a molecule  */
  get y(): number[] {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates!.y;
  };

  /** Z coordinates of all atoms in a molecule  */
  get z(): number[] {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates!.z;
  };

  get atomTypes(): string[] {
    if (this._atomTypes === undefined)
      this._atomTypes = this.parseAtomTypes();
    return this._atomTypes;
  }
}
