import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

const TRIPOS_MOLECULE_LINE = '@<TRIPOS>MOLECULE';
const TRIPOS_ATOM_LINE = '@<TRIPOS>ATOM';
const TRIPOS_BOND_LINE = '@<TRIPOS>BOND';

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
  /** Index of the molecule being currently parsed */
  private _idx: number;

  /** Get index of the next molecule block  */
  private getNextMolIdx(): number {
    return this._str.indexOf(TRIPOS_MOLECULE_LINE, this._idx);
  }

  /** Returns 'true' if there is a next Tripos molecule to parse  */
  public next(): boolean {
    return this.getNextMolIdx() !== -1;
  }

  public getNextMolFile(): string {
    const atoms = this.parseAtoms();


    this._idx = molEndIdx;
  }

  /** Extract data from TRIPOS atom block  */
  private parseAtoms(): Atoms {
    // position at the first atom line
    let begin = this._str.indexOf(TRIPOS_ATOM_LINE, molBeginIdx);
    begin = this._str.indexOf('\n', begin) + 1;

  }

  // todo: remove as unnecessary?
  /** Get the total number of molecules in a mol2 file */
  public getNumberOfMols(): number {
    const initIdx = this._idx;
    this._idx = 0;
    let numOfMols = 0;
    while (this.next()) {
      numOfMols++;
      this._idx = this.getNextMolIdx();
    }
    this._idx = initIdx;

    return numOfMols;
  }

  constructor(str: string) {
    this._str = str;
    this._idx = this.getNextMolIdx();
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
