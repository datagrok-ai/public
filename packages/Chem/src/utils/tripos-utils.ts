import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

const TRIPOS_MOLECULE_LINE = '@<TRIPOS>MOLECULE';
const TRIPOS_ATOM_LINE = '@<TRIPOS>ATOM';
const TRIPOS_BOND_LINE = '@<TRIPOS>BOND';

const HYDROGEN = 'H';

type Atoms = {
  atomTypes: string[],
  x: number[],
  y: number[],
  z: number[],
  charge: number[]
}

type Bonds = {
  bondTypes: number[],
  atomPairs: number[][]
}

/** Class for parsing a TRIPOS string into an array of MDL MOLFILEs. Currently,
 * only the atoms/bonds record types are supported */
class TriposParser {
  /** The content of mol2 file  */
  private _str: string;
  /** Index of the molecule being currently parsed, always points to some actual
   * value*/
  private _molIdx: number;

  /** Running idx within the current molecule */
  private _beginIdx: number;

  /** Running idx within the current molecule */
  private _endIdx: number;

  /** Get index of the next molecule block  */
  private getNextMolIdx(): number {
    const begin = (this._molIdx === 0) ? 0 : this._molIdx + 1;
    return this._str.indexOf(TRIPOS_MOLECULE_LINE, begin);
  }

  /** Returns the idx of the atom block for current molecule, or throws if
   * none */
  private getCurrentAtomBlockIdx(): number {
    const value = this._str.indexOf(TRIPOS_ATOM_LINE, this._molIdx);
    if (value !== -1)
      return value;
    else
      throw new Error('No next valid atom block in mol2');
  }

  /** Returns the idx of the bond block for current molecule, or throws if
   * none */
  private getCurrentBondBlockIdx(): number {
    const value = this._str.indexOf(TRIPOS_BOND_LINE, this._molIdx);
    if (value !== -1)
      return value;
    else
      throw new Error('No next valid bond block in mol2');
  }

  private getIdxOfNextLine(): number {
    if (this._str.at(this._beginIdx) !== '\n')
      return this._str.indexOf('\n', this._beginIdx) + 1;
    else
      return this._str.indexOf('\n', this._beginIdx + 1) + 1;
  }

  /** Returns 'true' if there is a next Tripos molecule to parse */
  public next(): boolean {
    return this.getNextMolIdx() !== -1;
  }

  public getNextMolFile(): string {
    const atoms = this.parseAtoms();
    const bonds = this.parseBonds();
  }

  /** Extract data from TRIPOS atom block  */
  private parseAtoms(): Atoms {
    const atomTypesArray: string[] = [];
    const x: number[] = [];
    const y: number[] = [];
    const z: number[] = [];
    const charge: number[] = [];

    // position at the first atom line
    this._beginIdx = this.getCurrentAtomBlockIdx();
    const endOfAtomBlockIdx = this.getCurrentBondBlockIdx();
    while (this._beginIdx !== endOfAtomBlockIdx) {
      // jump to next line
      this._beginIdx = this.getIdxOfNextLine();
      const atomType = this.parseAtomType();
      if (atomType !== HYDROGEN) {
        atomTypesArray.push(atomType);
        x.push(this.getFloatValue());
        y.push(this.getFloatValue());
        z.push(this.getFloatValue());
        charge.push(this.getFloatValue());
      }
    }
    return {
      atomTypes: atomTypesArray,
      x: x,
      y: y,
      z: z,
      charge: charge,
    };
  }

  private parseAtomType(): string {
    this.jumpToNextColumn(); // used for side effect only: shift to atom type
    this.jumpToNextColumn();
  }

  /** Get a float value and shift running idx values */
  private getFloatValue(): number {
    this.jumpToNextColumn();
  }

  /** Get a float value and shift running idx values */
  private getIntValue(): number {
    this.jumpToNextColumn();
  }

  /** Jumps to the next column relatively to this._beginIdx  */
  private jumpToNextColumn(): void {
    this._beginIdx = this.getNextColumnIdx();
  }

  /** Gets the idx of the next column relatively to this._beginIdx  */
  private getNextColumnIdx(): number {
    let idx = this._beginIdx;
    // skip non-whitespace, if necessary
    while (!this.isWhitespace(idx))
      ++idx;
    // skip whitespace
    while (this.isWhitespace(idx))
      ++idx;
    return idx;
  }

  private isWhitespace(idx: number): bool {
    return this._str.at(idx) === ' ';
  }

  // todo: remove as unnecessary?
  /** Get the total number of molecules in a mol2 file */
  public getNumberOfMols(): number {
    const initIdx = this._molIdx;
    this._molIdx = 0;
    let numOfMols = 0;
    while (this.next()) {
      numOfMols++;
      this._molIdx = this.getNextMolIdx();
    }
    this._molIdx = initIdx;

    return numOfMols;
  }

  constructor(str: string) {
    this._str = str;

    this._molIdx = 0;
    this._molIdx = this.getNextMolIdx();
    if (!this.next()) // check for existence of molecular data
      throw new Error('The mol2 file does not contain a molecule');

    // dummy values
    this._beginIdx = 0;
    this._endIdx = 0;
  }
}

export function _importTripos(bytes: Uint8Array): DG.DataFrame[] {
  const str = new TextDecoder().decode(bytes);

  const parser = new TriposParser(str);

  // const molfileArray = new Array<string>(parser.getNumberOfMols());
  // todo: consider the above alternative
  const molfileArray: string[] = [];
  
  while (parser.next())
    molfileArray.push(parser.getNextMolFile());

  const molCol = DG.Column.fromStrings('molecules', molfileArray);

  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.setTag(DG.TAGS.UNITS, DG.UNITS.Molecule.MOLBLOCK);

  return [DG.DataFrame.fromColumns([molCol])];
}
