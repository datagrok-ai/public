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
export abstract class ChemicalTableParserBase {
  constructor(fileContent: string) {
    this.init(fileContent);
  };

  protected fileContent!: string;
  protected _atomCount?: number;
  protected _bondCount?: number;
  protected xyzAtomCoordinates?: CoordinateArrays;
  protected _atomTypes?: string[];
  protected _pairsOfBondedAtoms?: Uint16Array[];
  protected _bondTypes?: Uint16Array;

  protected init(fileContent: string): void {
    this.fileContent = fileContent.replace(/\r/g, '');
    this._atomCount = undefined;
    this._atomTypes = undefined;
    this._bondCount = undefined;
    this._bondTypes = undefined;
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
    this.xyzAtomCoordinates ??= this.parseAtomCoordinates();
    return this.xyzAtomCoordinates.x;
  };

  /** Y coordinates of all atoms in a molecule  */
  get y(): Float32Array {
    this.xyzAtomCoordinates ??= this.parseAtomCoordinates();
    return this.xyzAtomCoordinates!.y;
  };

  /** Z coordinates of all atoms in a molecule  */
  get z(): Float32Array {
    this.xyzAtomCoordinates ??= this.parseAtomCoordinates();
    return this.xyzAtomCoordinates!.z;
  };

  get atomTypes(): string[] {
    this._atomTypes ??= this.parseAtomTypes();
    return this._atomTypes;
  }

  get pairsOfBondedAtoms(): Uint16Array[] {
    this._pairsOfBondedAtoms ??= this.parseBondedAtomPairs();
    return this._pairsOfBondedAtoms!;
  }

  get bondTypes(): Uint16Array {
    this._bondTypes ??= this.parseBondTypes();
    return this._bondTypes!;
  }

  protected abstract parseAtomAndBondCounts(): AtomAndBondCounts;
  protected abstract getCountsLineIdx(): number;
  protected abstract getAtomBlockIdx(): number;
  protected abstract getBondBlockIdx(): number;
  protected abstract shiftIdxToXColumn(lineStartIdx: number): number;
  protected abstract shiftIdxToAtomType(lineStartIdx: number): number;
  protected abstract shiftIdxToBondedAtomsPair(lineStartIdx: number): number;
  protected abstract shiftIdxToBondType(lineStartIdx: number): number;
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
    for (let i = 0; i < this.bondCount; i++) {
      idx = this.shiftIdxToBondedAtomsPair(idx);
      const pair = new Uint16Array(2);
      pair[0] = this.parseIntValue(idx);
      idx = this.getNextColumnIdx(idx);
      pair[1] = this.parseIntValue(idx);
      bondedAtomPairs[i] = pair;
      idx = this.getNextLineIdx(idx);
    }
    return bondedAtomPairs;
  }

  protected parseBondTypes(): Uint16Array {
    const bondCount = this.bondCount;
    const bondTypes = new Uint16Array(bondCount);
    let idx = this.getBondBlockIdx();
    for (let i = 0; i < bondCount; i++) {
      idx = this.shiftIdxToBondType(idx);
      bondTypes[i] = this.parseIntValue(idx);
      idx = this.getNextLineIdx(idx);
    }
    return bondTypes;
  };

  protected isWhitespace(idx: number): boolean {
    return /\s/.test(this.fileContent.at(idx)!);
  }

  protected getNextLineIdx(idx: number): number {
    if (this.fileContent.at(idx) !== '\n')
      return this.fileContent.indexOf('\n', idx) + 1;
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
    const value = parserFunction(this.fileContent.substring(idxOfNumber, end));
    return value;
  }
}
