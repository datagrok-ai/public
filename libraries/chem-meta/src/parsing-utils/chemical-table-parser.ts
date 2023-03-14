export interface ChemicalTableParser {
  init(file: string): void;
  atomCount: number;
  bondCount: number;
  x: Float32Array;
  y: Float32Array;
  z: Float32Array;
  pairsOfBondedAtoms: Uint16Array[];
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

  protected static instance: ChemicalTableParserBase;

  protected file!: string;
  protected _atomCount?: number;
  protected _bondCount?: number;
  protected xyzAtomCoordinates?: CoordinateArrays;
  protected _atomTypes?: string[];
  protected _pairsOfBondedAtoms?: Uint16Array[];

  public init(file: string): void {
    this.file = file.replaceAll('\r', '');
    this._atomCount = undefined;
    this._atomTypes = undefined;
    this._bondCount = undefined;
    this.xyzAtomCoordinates = undefined;
    this._pairsOfBondedAtoms = undefined;
  }

  /** Total number of atoms in a molecule */
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
    if (this.xyzAtomCoordinates === undefined)
      this.xyzAtomCoordinates = this.parseAtomCoordinates();
    return this.xyzAtomCoordinates.x;
  };

  /** Y coordinates of all atoms in a molecule  */
  get y(): Float32Array {
    if (this.xyzAtomCoordinates === undefined)
      this.xyzAtomCoordinates = this.parseAtomCoordinates();
    return this.xyzAtomCoordinates!.y;
  };

  /** Z coordinates of all atoms in a molecule  */
  get z(): Float32Array {
    if (this.xyzAtomCoordinates === undefined)
      this.xyzAtomCoordinates = this.parseAtomCoordinates();
    return this.xyzAtomCoordinates!.z;
  };

  get atomTypes(): string[] {
    if (this._atomTypes === undefined)
      this._atomTypes = this.parseAtomTypes();
    return this._atomTypes;
  }

  get pairsOfBondedAtoms(): Uint16Array[] {
    if (this._pairsOfBondedAtoms === undefined)
      this._pairsOfBondedAtoms = this.parseBondedAtomPairs();
    return this._pairsOfBondedAtoms!;
  }

  protected abstract parseAtomAndBondCounts(): AtomAndBondCounts;
  protected abstract getCountsLineIdx(): number;
  protected abstract getAtomBlockIdx(): number;
  protected abstract getBondBlockIdx(): number;
  protected abstract shiftIdxToXColumn(lineStartIdx: number): number;
  protected abstract shiftIdxToAtomType(lineStartIdx: number): number;
  protected abstract shiftIdxToBondedAtomsPair(lineStartIdx: number): number;
  protected abstract parseAtomType(idx: number): string;

  protected setAtomAndBondCounts(): void {
    const {atomCount, bondCount} = this.parseAtomAndBondCounts();
    this._atomCount = atomCount;
    this._bondCount = bondCount;
  }

  protected getNextColumnIdx(idx: number): number {
    // skip non-whitespace, if necessary
    while (!this.isWhitespace(idx))
      ++idx;
    // skip whitespace
    while (this.isWhitespace(idx))
      ++idx;
    return idx;
  }

  protected shiftIdxToSpecifiedColumn(lineStartIdx: number, columnNumber: number) {
    let idx = lineStartIdx;
    const numberOfJumps = this.isWhitespace(idx) ? columnNumber : columnNumber - 1;
    for (let i = 0; i < numberOfJumps; i++)
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

  protected parseBondedAtomPairs(): Uint16Array[] {
    const bondedAtomPairs = new Array<Uint16Array>(this.bondCount);
    let idx = this.getBondBlockIdx();
    idx = this.getNextLineIdx(idx);
    for (let i = 0; i < this.bondCount; i++) {
      const pair = new Uint16Array(2);
      idx = this.shiftIdxToBondedAtomsPair(idx);
      pair[0] = this.parseIntValue(idx);
      idx = this.getNextColumnIdx(idx);
      pair[1] = this.parseIntValue(idx);
      bondedAtomPairs[i] = pair;
    }
    return bondedAtomPairs;
  }

  protected isWhitespace(idx: number): boolean {
    return /\s/.test(this.file.at(idx)!);
  }

  protected getNextLineIdx(idx: number): number {
    if (this.file.at(idx) !== '\n')
      return this.file.indexOf('\n', idx) + 1;
    else
      return idx + 1;
  }

  protected parseFloatValue(idxOfNumber: number): number {
    return this.parseNumericValue(parseFloat, idxOfNumber);
  }

  protected parseIntValue(idxOfNumber: number): number {
    return this.parseNumericValue(parseInt, idxOfNumber);
  }

  protected parseNumericValue(
    parserFunction: (str: string) => number,
    idxOfNumber: number
  ): number {
    let end = idxOfNumber + 1;
    while (!this.isWhitespace(end))
      ++end;
    const value = parserFunction(this.file.substring(idxOfNumber, end));
    return value;
  }
}
