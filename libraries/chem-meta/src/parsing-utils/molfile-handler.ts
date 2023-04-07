import {ChemicalTableParser, ChemicalTableParserBase, AtomAndBondCounts} from './chemical-table-parser';
import {V2K_CONST} from '../formats/molfile-v2k-const';
import {V3K_CONST} from '../formats/molfile-v3k-const';

const enum MOLFILE_VERSION {
  V2000 = 'V2000',
  V3000 = 'V3000',
}

export class MolfileHandler extends ChemicalTableParserBase implements ChemicalTableParser {
  constructor(molfile: string) {
    super(molfile);
    this.reset(molfile);
  }

  /** Init/reset the state of the handler for a new molfile */
  public reset(molfile: string) {
    super.reset(molfile);

    const molfileVersion = MolfileHandler.determineMolfileVersion(this.fileContent);
    const isV2K = (molfileVersion === MOLFILE_VERSION.V2000);

    this.parseAtomAndBondCounts = isV2K ? this.parseAtomAndBondCountsV2K : this.parseAtomAndBondCountsV3K;

    this.getAtomBlockIdx = isV2K ? this.getAtomBlockIdxV2K : this.getAtomBlockIdxV3K;

    this.getCountsLineIdx = isV2K ? this.getCountsLineV2KIdx : this.getCountsLineV3KIdx;

    this.shiftIdxToXColumn = isV2K ? this.shiftIdxToXColumnV2K : this.shiftIdxToXColumnV3K;

    this.shiftIdxToAtomType = isV2K ? this.shiftIdxToAtomTypeV2K :
      this.shiftIdxToAtomTypeV3K;

    this.getBondBlockIdx = isV2K ? this.getBondBlockIdxV2K : this.getBondBlockIdxV3K;

    this.shiftIdxToBondedAtomsPair = isV2K ? this.shiftIdxToBondedAtomsPairV2K : this.shiftIdxToBondedAtomsPairV3K;

    this.shiftIdxToBondType = isV2K ? this.shiftIdxToBondTypeV2K : this.shiftIdxToBondTypeV3K;
  }

  public static createInstance(file: string): MolfileHandler {
    if (!this.instance)
      this.instance = new MolfileHandler(file);
    return this.instance as MolfileHandler;
  }

  protected parseAtomAndBondCounts!: () => AtomAndBondCounts;
  protected getCountsLineIdx!: () => number;
  protected getAtomBlockIdx!: () => number;
  protected shiftIdxToXColumn!: (lineStartIdx: number) => number;
  protected shiftIdxToAtomType!: (lineStartIdx: number) => number;
  protected getBondBlockIdx!: () => number;
  protected shiftIdxToBondedAtomsPair!: (lineStartIdx: number) => number;
  protected shiftIdxToBondType!: (lineStartIdx: number) => number;

  protected parseAtomType(idx: number): string {
    let begin = idx;
    let end = begin;
    if (this.isQuote(begin)) {
      end = this.getNextIdenticalChar(begin);
      begin++;
    } else {
      end = this.fileContent.indexOf(' ', end);
    }
    return this.fileContent.substring(begin, end);
  }

  /** Determine whether the file is V2000/V3000, or throw */
  private static determineMolfileVersion(molfile: string): MOLFILE_VERSION {
    if (MolfileHandler.validateV3K(molfile))
      return MOLFILE_VERSION.V3000;
    else if (MolfileHandler.validateV2K(molfile))
      return MOLFILE_VERSION.V2000;
    else
      throw new Error('Malformed molfile');
  }

  private shiftIdxToAtomTypeV2K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.ATOM_TYPE_COL);
  }

  private shiftIdxToAtomTypeV3K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.ATOM_TYPE_COL);
  }

  private getCountsLineV2KIdx(): number {
    let idx = 0;
    for (let i = 0; i < V2K_CONST.NUM_OF_HEADER_LINES; ++i)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getCountsLineV3KIdx(): number {
    return this.fileContent.indexOf(V3K_CONST.BEGIN_COUNTS_LINE);
  }

  private getAtomBlockIdxV2K(): number {
    let idx = this.getCountsLineIdx();
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getAtomBlockIdxV3K(): number {
    let idx = this.fileContent.indexOf(V3K_CONST.BEGIN_ATOM_BLOCK);
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  private shiftIdxToXColumnV2K(lineStartIdx: number): number {
    return this.getNextColumnIdx(lineStartIdx);
  }

  private isQuote(idx: number): boolean {
    return this.fileContent.at(idx) === '\"' || this.fileContent.at(idx) === '\'';
  }

  private getNextIdenticalChar(idx: number): number {
    const sym = this.fileContent.at(idx);
    if (sym)
      return this.fileContent.indexOf(sym, idx + 1);
    else
      return -1;
  }

  private shiftIdxToXColumnV3K(lineStartIdx: number): number {
    let idx = this.shiftIdxToAtomType(lineStartIdx);
    if (this.isQuote(idx)) {
      idx = this.getNextIdenticalChar(idx);
      idx = this.getNextColumnIdx(idx);
      return idx;
    }
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.X_COL);
  }

  private shiftIdxToBondedAtomsPairV2K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.FIRST_BONDED_ATOM_COL);
  }

  private shiftIdxToBondedAtomsPairV3K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.FIRST_BONDED_ATOM_COL);
  }

  private shiftIdxToBondTypeV2K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.BOND_TYPE_COL);
  }

  private shiftIdxToBondTypeV3K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.BOND_TYPE_COL);
  }

  private getBondBlockIdxV2K(): number {
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < this.atomCount; i++)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getBondBlockIdxV3K(): number {
    return this.getNextLineIdx(this.fileContent.indexOf(V3K_CONST.BEGIN_BOND_BLOCK));
  }

  private static validateV3K(molfile: string): boolean {
    return (molfile.indexOf(V3K_CONST.HEADER) !== -1 &&
    molfile.indexOf(V3K_CONST.END) !== -1);
  }

  private static validateV2K(molfile: string): boolean {
    return (molfile.indexOf(V2K_CONST.HEADER) !== -1 &&
    molfile.indexOf(V2K_CONST.END) !== -1);
  }

  private parseAtomAndBondCountsV2K(): AtomAndBondCounts {
    let begin = this.getCountsLineIdx();
    let end = begin + V2K_CONST.NUM_OF_COUNTS_DIGITS;
    const atomCount = parseInt(this.fileContent.substring(begin, end));
    begin = end;
    end += V2K_CONST.NUM_OF_COUNTS_DIGITS;
    const bondCount = parseInt(this.fileContent.substring(begin, end));
    return {atomCount: atomCount, bondCount: bondCount};
  };

  private parseAtomAndBondCountsV3K(): AtomAndBondCounts {
    // parse atom count
    let begin = this.fileContent.indexOf(V3K_CONST.BEGIN_COUNTS_LINE) + V3K_CONST.COUNTS_SHIFT;
    let end = this.fileContent.indexOf(' ', begin + 1);
    const numOfAtoms = parseInt(this.fileContent.substring(begin, end));

    // parse bond count
    begin = end + 1;
    end = this.fileContent.indexOf(' ', begin + 1);
    const numOfBonds = parseInt(this.fileContent.substring(begin, end));

    return {atomCount: numOfAtoms, bondCount: numOfBonds};
  }
}
