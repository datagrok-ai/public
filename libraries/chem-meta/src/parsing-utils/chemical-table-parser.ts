export interface ChemicalTableParser {
  init(file: string): void;
  atomCount: number;
  bondCount: number;
  x: Float32Array;
  y: Float32Array;
  z: Float32Array;
}

export type AtomAndBondCounts = {
  atomCount: number,
  bondCount: number
}

type CoordinateArrays = {
  x: Float32Array,
  y: Float32Array,
  z: Float32Array,
}

/** Base singleton for Molfile or Mol2 parser/handler */
export abstract class ChemicalTableParserBase implements ChemicalTableParser {
  protected constructor(file: string) {
    this.init(file);
  };

  // Public members

  // public abstract static createInstance(file: string): ChemicalTableParserBase;

  public init(file: string): void {
    this.file = file.replace(/\r/g, '');
    this._atomCount = undefined;
    this._atomTypes = undefined;
    this._bondCount = undefined;
    this.atomCoordinates = undefined;
  }

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
  get x(): Float32Array {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates.x;
  };

  /** Y coordinates of all atoms in a molecule  */
  get y(): Float32Array {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates!.y;
  };

  /** Z coordinates of all atoms in a molecule  */
  get z(): Float32Array {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates!.z;
  };

  get atomTypes(): string[] {
    if (this._atomTypes === undefined)
      this._atomTypes = this.parseAtomTypes();
    return this._atomTypes;
  }

  // Protected members

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
    const atomCount = this.atomCount;
    const atomTypes = new Array<string>(atomCount);
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < atomCount; i++) {
      idx = this.shiftIdxToAtomType(idx);
      atomTypes[i] = this.parseAtomType(idx);
      idx = this.getNextLineIdx(idx);
    }
    return atomTypes;
  };

  protected parseAtomCoordinates(): CoordinateArrays {
    const x = new Float32Array(this.atomCount);
    const y = new Float32Array(this.atomCount);
    const z = new Float32Array(this.atomCount);
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < this.atomCount; i++) {
      idx = this.shiftIdxToXColumn(idx);
      for (const item of [x, y, z]) {
        item[i] = this.parseFloatValue(idx);
        idx = this.getNextColumnIdx(idx);
      }
      idx = this.getNextLineIdx(idx);
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
      return idx + 1;
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

  protected static instance: ChemicalTableParserBase;
}
