/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Interface for parsers of Molfile and Mol2 formats */
interface ChemicalTableParser {
  X: Float32Array;
  Y: Float32Array;
  Z: Float32Array;
  AtomTypes: string[];
}

export type AtomAndBondCounts = {
  atomCount: number,
  bondCount: number
}

/** Base class for Molfile and Mol2 parsers/handlers */
export abstract class AbstractChemicalTableParser implements ChemicalTableParser {
  constructor(file: string) {
    this.file = file.replaceAll('\r', '');
  };

  protected readonly file: string;

  /** Index running along the string/file being parsed  */
  protected currentIdx: number = 0;
  protected _atomCount?: number;
  protected _bondCount?: number;
  /** The array of X, Y, Z arrays for atomic coordinates */
  protected atomCoordinates?: Float32Array[];
  protected atomTypes?: string[];

  protected abstract parseAtomAndBondCounts(): AtomAndBondCounts;
  protected abstract parseAtomCoordinates(): Float32Array[];
  protected abstract parseAtomTypes(): string[];
  /** Get idx of the first line of the atom block  */
  protected abstract getAtomBlockIdx(): number;
  /** Get idx of the first line of the bond block  */
  protected abstract getBondBlockIdx(): number;

  protected setAtomAndBondCounts(): void {
    const {atomCount, bondCount} = this.parseAtomAndBondCounts();
    this._atomCount = atomCount;
    this._bondCount = bondCount;
  }

  /** Gets the idx of the next column relatively to this._currentIdx  */
  protected getNextColumnIdx(): number {
    let idx = this.currentIdx;
    // skip non-whitespace, if necessary
    while (!this.isWhitespace(idx))
      ++idx;
    // skip whitespace
    while (this.isWhitespace(idx))
      ++idx;
    return idx;
  }

  /** Check if a character is whitespace including '\t'  */
  protected isWhitespace(idx: number): boolean {
    return /\s/.test(this.file.at(idx)!);
  }

  protected jumpToNextLine(): void {
    this.currentIdx = this.getIdxOfNextLine(this.currentIdx);
  }

  /** Get index of the next line starting from idx  */
  protected getIdxOfNextLine(idx: number): number {
    if (this.file.at(idx) !== '\n')
      return this.file.indexOf('\n', idx) + 1;
    else
      return this.file.indexOf('\n', idx + 1) + 1;
  }

  /** Get a float value in the current column */
  protected parseFloatValue(): number {
    return this.parseNumericValue(parseFloat);
  }

  /** Get an int value in the current column */
  protected parseIntValue(): number {
    return this.parseNumericValue(parseInt);
  }

  /** Parse a numeric value depending on the functional argument  */
  protected parseNumericValue(parserFunction: (str: string) => number): number {
    let end = this.currentIdx + 1;
    while (!this.isWhitespace(end))
      ++end;
    const value = parserFunction(this.file.substring(this.currentIdx, end));
    return value;
  }

  protected get atomCount(): number {
    if (this._atomCount === undefined)
      this.setAtomAndBondCounts();
    return this._atomCount!;
  }

  get X(): Float32Array {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates![0];
  };

  get Y(): Float32Array {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates![1];
  };

  get Z(): Float32Array {
    if (this.atomCoordinates === undefined)
      this.atomCoordinates = this.parseAtomCoordinates();
    return this.atomCoordinates![2];
  };

  get AtomTypes(): string[] {
    if (this.atomTypes === undefined)
      this.atomTypes = this.parseAtomTypes();
    return this.atomTypes;
  }
}
