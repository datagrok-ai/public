// import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
// import * as ui from 'datagrok-api/ui';

// todo: import common constants from to-atomic-level
const V3K_HEADER_FIRST_LINE = '\nDatagrok macromolecule handler\n\n';
const V3K_HEADER_SECOND_LINE = '  0  0  0  0  0  0            999 V3000\n';
const V3K_BEGIN_CTAB_BLOCK = 'M  V30 BEGIN CTAB\n';
const V3K_END_CTAB_BLOCK = 'M  V30 END CTAB\n';
const V3K_BEGIN_COUNTS_LINE = 'M  V30 COUNTS ';
const V3K_COUNTS_LINE_ENDING = ' 0 0 0\n';
const V3K_BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM\n';
const V3K_END_ATOM_BLOCK = 'M  V30 END ATOM\n';
const V3K_BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND\n';
const V3K_END_BOND_BLOCK = 'M  V30 END BOND\n';
// const V3K_BOND_CONFIG = ' CFG=';
const V3K_BEGIN_DATA_LINE = 'M  V30 ';
const V3K_END = 'M  END\n';

const V3K_CHARGE = ' CHG=';

const TRIPOS_MOLECULE_LINE = '@<TRIPOS>MOLECULE';
const TRIPOS_ATOM_LINE = '@<TRIPOS>ATOM';
const TRIPOS_BOND_LINE = '@<TRIPOS>BOND';

const HYDROGEN = 'H';

type Atoms = {
  atomTypes: string[],
  x: number[],
  y: number[],
  z: number[],
  charge: number[],
  initialIdx: number[],
}

type Bonds = {
  bondTypes: number[],
  atomPairs: number[][]
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

/** Class for parsing a TRIPOS string into an array of MDL MOLFILEs. Currently,
 * only the atoms/bonds record types are supported */
class TriposParser {
  /** The content of mol2 file  */
  private _str: string;
  /** Index of the molecule being currently parsed, always points to some actual
   * value*/
  private _currentMolIdx: number;

  /** Running idx within the current molecule */
  private _currentIdx: number; // todo: rename

  /** Set of hydrogen indices  */
  private _hydrogenIndices: Set<number>;

  constructor(str: string) {
    this._str = str;

    this._currentMolIdx = 0;
    //this._currentMolIdx = this.getNextMolIdx();
    if (!this.next()) // check for existence of molecular data
      throw new Error('The mol2 file does not contain a molecule');

    this._currentIdx = 0;
    this._hydrogenIndices = new Set<number>();
  }

  /** Get index of the next molecule block  */
  private getNextMolIdx(): number {
    const begin = (this._currentIdx === 0) ? 0 : this._currentIdx + 1;
    const idx = this._str.indexOf(TRIPOS_MOLECULE_LINE, begin);
    this._currentIdx = idx;
    return idx;
  }

  /** Returns the idx of the atom block for current molecule, or throws if
   * none */
  private getCurrentAtomBlockIdx(): number {
    const value = this._str.indexOf(TRIPOS_ATOM_LINE, this._currentIdx);
    if (value !== -1)
      return value;
    else
      throw new Error('No next valid atom block in mol2');
  }

  /** Returns the idx of the bond block for current molecule, or throws if
   * none */
  private getCurrentBondBlockIdx(): number {
    const value = this._str.indexOf(TRIPOS_BOND_LINE, this._currentIdx);
    if (value !== -1)
      return value;
    else
      throw new Error('No next valid bond block in mol2');
  }

  private getIdxOfNextLine(): number {
    if (this._str.at(this._currentIdx) !== '\n')
      return this._str.indexOf('\n', this._currentIdx) + 1;
    else
      return this._str.indexOf('\n', this._currentIdx + 1) + 1;
  }

  /** Returns 'true' if there is a next Tripos molecule to parse */
  // todo: delete as unnecessary
  public next(): boolean {
    return this.getNextMolIdx() !== -1;
  }

  public getNextMolFile(): string {
    const {atomCount, bondCount} = this.parseAtomAndBondCounts();
    const atoms = this.parseAtoms(atomCount);
    const bonds = this.parseBonds(bondCount, atoms);

    const molfile = this.getMolFile(atoms, bonds);

    //this._currentMolIdx = this.getNextMolIdx();

    return molfile;
  }

  private parseAtomAndBondCounts(): {atomCount: number, bondCount: number} {
    // position at the molecule line
   // this._currentIdx = this._currentMolIdx;
    for (let i = 0; i < 2; ++i)
      this.jumpToNextLine();

    const atomCount = this.parseIntValue();
    this.jumpToNextColumn();
    const bondCount = this.parseIntValue();

    return {atomCount: atomCount, bondCount: bondCount};
  }

  private jumpToNextLine(): void {
    this._currentIdx = this.getIdxOfNextLine();
  }

  /** Construct molfile from the parsed atom and bond data  */
  private getMolFile(atoms: Atoms, bonds: Bonds): string {
    const atomCount = atoms.atomTypes.length;
    const bondCount = bonds.bondTypes.length;

    const atomType = atoms.atomTypes;
    const x = atoms.x;
    const y = atoms.y;
    const z = atoms.z;
    const charge = atoms.charge;

    const bondType = bonds.bondTypes;
    const atomPair = bonds.atomPairs;
    // const bondKwargs = molGraph.bonds.kwargs;
    // const bondConfig = molGraph.bonds.bondConfiguration;

    const molfileCountsLine = V3K_BEGIN_COUNTS_LINE + atomCount + ' ' + bondCount + V3K_COUNTS_LINE_ENDING;

    // atom block
    let molfileAtomBlock = '';
    for (let i = 0; i < atomCount; ++i) {
      const atomIdx = i + 1;
      // const coordinate = [x[i].toString(), y[i].toString(), z[i].toString()];
      const chargeKwarg = (charge[i] === 0) ? '' :
        V3K_CHARGE + charge[i].toString();

      const atomLine = V3K_BEGIN_DATA_LINE + atomIdx + ' ' + atomType[i] + ' ' +
        x[i].toString() + ' ' + y[i].toString() + ' ' + z[i].toString() +
        ' 0' + chargeKwarg + '\n';
      molfileAtomBlock += atomLine;
    }

    // bond block
    let molfileBondBlock = '';
    for (let i = 0; i < bondCount; ++i) {
      const bondIdx = i + 1;
      const firstAtom = atomPair[i][0];
      const secondAtom = atomPair[i][1];
      // const kwargs = bondKwargs.has(i) ? ' ' + bondKwargs.get(i) : '';
      // const bondCfg = bondConfig.has(i) ? ' CFG=' + bondConfig.get(i) : '';
      // const bondLine = V3K_BEGIN_DATA_LINE + bondIdx + ' ' + bondType[i] + ' ' +
      //   firstAtom + ' ' + secondAtom + bondCfg + kwargs + '\n';
      const bondLine = V3K_BEGIN_DATA_LINE + bondIdx + ' ' + bondType[i] + ' ' +
        firstAtom + ' ' + secondAtom + '\n';
      molfileBondBlock += bondLine;
    }

    const molfileParts = [
      V3K_HEADER_FIRST_LINE,
      V3K_HEADER_SECOND_LINE,
      V3K_BEGIN_CTAB_BLOCK,
      molfileCountsLine,
      V3K_BEGIN_ATOM_BLOCK,
      molfileAtomBlock,
      V3K_END_ATOM_BLOCK,
      V3K_BEGIN_BOND_BLOCK,
      molfileBondBlock,
      V3K_END_BOND_BLOCK,
      V3K_END_CTAB_BLOCK,
      V3K_END,
    ];
    const resultingMolfile = molfileParts.join('');
    // console.log(resultingMolfile);

    return resultingMolfile;
  }

  /** Extract data from TRIPOS atom block  */
  private parseAtoms(atomCount: number): Atoms {
    const atomTypesArray: string[] = [];
    const x: number[] = [];
    const y: number[] = [];
    const z: number[] = [];
    const charge: number[] = [];
    const initialIdx: number[] = [];

    // position at the first atom line
    this._currentIdx = this.getCurrentAtomBlockIdx();
    if (this._currentIdx === -1)
      throw new Error('No valid atom block');

    for (let i = 0; i < atomCount; ++i) {
      this.jumpToNextLine();

      const atomType = this.parseAtomType();
      if (atomType !== HYDROGEN) {
        initialIdx.push(i + 1);

        atomTypesArray.push(atomType);

        this.jumpToNextColumn(); // jump to the X coordinates column
        x.push(this.parseFloatValue());

        this.jumpToNextColumn(); // jump to the Y coordinates column
        y.push(this.parseFloatValue());

        this.jumpToNextColumn(); // jump to the Z coordinates column
        z.push(this.parseFloatValue());

        for (let j = 0; j < 4; ++j) // jump to the last column with charge value
          this.jumpToNextColumn();
        charge.push(Math.round(this.parseFloatValue())); // charge values are integer in MDL molfiles
        // todo: molfile charges are limited to the span from -15 to 15,
        // condider the corresponding exception
      } else
        this._hydrogenIndices.add(i + 1);
    }
    return {
      atomTypes: atomTypesArray,
      x: x,
      y: y,
      z: z,
      charge: charge,
      initialIdx: initialIdx,
    };
  }

  /** Extract data from TRIPOS bond block  */
  private parseBonds(bondCount: number, atoms: Atoms): Bonds {
    const bondTypes: number[] = [];
    const atomPairs: number[][] = [];

    // position at the first bond line
    this._currentIdx = this.getCurrentBondBlockIdx();
    if (this._currentIdx === -1)
      throw new Error('No valid bond block');

    for (let i = 0; i < bondCount; ++i) {
      this.jumpToNextLine();

      // jump to the 2nd column
      for (let j = 0; j < 2; ++j)
        this.jumpToNextColumn();

      let pairHasHydrogen = false;

      const pairOfAtoms: number[] = [];
      for (let j = 0; j < 2; ++j) {
        const atomIdx = this.parseIntValue();
        if (this._hydrogenIndices.has(atomIdx)) {
          pairHasHydrogen = true;
          break; // we don't need atom pairs with hydrogens
        }
        pairOfAtoms.push(atomIdx);
        this.jumpToNextColumn();
      }
      
      if (pairHasHydrogen)
        continue; // we don't need atom pairs with hydrogens
      else {
        atomPairs.push(pairOfAtoms);

        // no need to jump to this column, already here
        bondTypes.push(this.parseIntValue());
      }
    }

    const bonds = {
      bondTypes: bondTypes,
      atomPairs: atomPairs,
    };

    this.renumberBonds(atoms, bonds);

    return bonds;
  }

  /** Renumber bonds taking into account the ommitted hydrogens  */
  private renumberBonds(atoms: Atoms, bonds: Bonds) {
    /** Maps the initial atom indices to the ones after hydrogen elimination */
    const map = new Map<number, number>();
    for (let i = 0; i < atoms.initialIdx.length; ++i)
      map.set(atoms.initialIdx[i], i + 1);

    for (let i = 0; i < bonds.atomPairs.length; ++i) {
      for (let j = 0; j < 2; ++j)
        bonds.atomPairs[i][j] = map.get(bonds.atomPairs[i][j])!;
    }
  }

  // private atBondLine(): boolean {
  //   const end = this._str.indexOf('\n', this._currentIdx);
  //   return this._str.substring(this._currentIdx, end) === TRIPOS_BOND_LINE;
  // }

  /** Jumps to atom type column and parses it */
  private parseAtomType(): string {
    // jump to the 2nd column
    for (let i = 0; i < 2; ++i)
      this.jumpToNextColumn();

    // let end = this._endIdx + 1;
    let end = this._currentIdx + 1;
    while (!this.isWhitespace(end))
      ++end;

    const atomType = this._str.substring(this._currentIdx, end);

    return atomType;
  }

  /** Get a float value in the current column */
  private parseFloatValue(): number {
    return this.parseNumericValue(parseFloat);
  }

  /** Get an int value in the current column */
  private parseIntValue(): number {
    return this.parseNumericValue(parseInt);
  }

  /** Parse a numeric value depending on the functional argument  */
  private parseNumericValue(parserFunction: (str: string) => number): number {
    let end = this._currentIdx + 1;
    while (!this.isWhitespace(end))
      ++end;
    const value = parserFunction(this._str.substring(this._currentIdx, end));
    return value;
  }

  /** Jumps to the next column relatively to this._currentIdx  */
  private jumpToNextColumn(): void {
    this._currentIdx = this.getNextColumnIdx();
  }

  /** Gets the idx of the next column relatively to this._currentIdx  */
  private getNextColumnIdx(): number {
    let idx = this._currentIdx;
    // skip non-whitespace, if necessary
    while (!this.isWhitespace(idx))
      ++idx;
    // skip whitespace
    while (this.isWhitespace(idx))
      ++idx;
    return idx;
  }

  private isWhitespace(idx: number): boolean {
    return this._str[idx] === ' ';
  }

  // todo: remove as unnecessary?
  /** Get the total number of molecules in a mol2 file */
  public getNumberOfMols(): number {
    const initIdx = this._currentMolIdx;
    this._currentMolIdx = 0;
    let numOfMols = 0;
    while (this.next()) {
      numOfMols++;
      this._currentMolIdx = this.getNextMolIdx();
    }
    this._currentMolIdx = initIdx;

    return numOfMols;
  }
}
