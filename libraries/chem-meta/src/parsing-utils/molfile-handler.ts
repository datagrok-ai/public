import {ChemicalTableParser, ChemicalTableParserBase, AtomAndBondCounts} from './chemical-table-parser';
import {MOLFILE_VERSION, V2K_CONST, V3K_CONST} from './molfile-parsing-const';

export class MolfileHandler extends ChemicalTableParserBase implements ChemicalTableParser {
  constructor(molfile: string) {
    super(molfile);
    this.init(molfile);
  }

  /** Init/reset the state of the handler for a new molfile */
  public init(molfile: string) {
    super.init(molfile);

    const molfileVersion = MolfileHandler.determineMolfileVersion(this.file);
    const isV2K = (molfileVersion === MOLFILE_VERSION.V2000);

    this.parseAtomAndBondCounts = isV2K ? this.parseAtomAndBondCountsV2K : this.parseAtomAndBondCountsV3K;

    this.getAtomBlockIdx = isV2K ? this.getAtomBlockIdxV2K : this.getAtomBlockIdxV3K;

    this.getCountsLineIdx = isV2K ? this.getCountsLineV2KIdx : this.getCountsLineV3KIdx;

    this.shiftIdxToXColumn = isV2K ? this.shiftIdxToXColumnV2K : this.shiftIdxToXColumnV3K;

    this.shiftIdxToAtomType = isV2K ? this.shiftIdxToAtomTypeV2K :
      this.shiftIdxToAtomTypeV3K;

    this.getBondBlockIdx = isV2K ? this.getBondBlockIdxV2K : this.getBondBlockIdxV3K;

    this.shiftIdxToBondedAtomsPair = isV2K ? this.shiftIdxToBondedAtomsPairV2K : this.shiftIdxToBondedAtomsPairV3K;
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
  protected shiftIdxToAtomType!: (idx: number) => number;
  protected getBondBlockIdx!: () => number;
  protected shiftIdxToBondedAtomsPair!: (lineStartIdx: number) => number;

  protected parseAtomType(idx: number): string {
    const begin = idx;
    const end = this.file.indexOf(' ', begin);
    return this.file.substring(begin, end);
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
    return this.file.indexOf(V3K_CONST.BEGIN_COUNTS_LINE);
  }

  private getAtomBlockIdxV2K(): number {
    let idx = this.getCountsLineIdx();
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getAtomBlockIdxV3K(): number {
    let idx = this.file.indexOf(V3K_CONST.BEGIN_ATOM_BLOCK);
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  private shiftIdxToXColumnV2K(lineStartIdx: number): number {
    return this.getNextColumnIdx(lineStartIdx);
  }

  private shiftIdxToXColumnV3K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.X_COL);
  }

  private shiftIdxToBondedAtomsPairV2K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.FIRST_BONDED_ATOM_COL);
  }

  private shiftIdxToBondedAtomsPairV3K(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.FIRST_BONDED_ATOM_COL);
  }
  private getBondBlockIdxV2K(): number {
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < this.atomCount; i++)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getBondBlockIdxV3K(): number {
    return this.getNextLineIdx(this.file.indexOf(V3K_CONST.BEGIN_BOND_BLOCK));
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
    const atomCount = parseInt(this.file.substring(begin, end));
    begin = end;
    end += V2K_CONST.NUM_OF_COUNTS_DIGITS;
    const bondCount = parseInt(this.file.substring(begin, end));
    return {atomCount: atomCount, bondCount: bondCount};
  };

  private parseAtomAndBondCountsV3K(): AtomAndBondCounts {
    // parse atom count
    let begin = this.file.indexOf(V3K_CONST.BEGIN_COUNTS_LINE) + V3K_CONST.COUNTS_SHIFT;
    let end = this.file.indexOf(' ', begin + 1);
    const numOfAtoms = parseInt(this.file.substring(begin, end));

    // parse bond count
    begin = end + 1;
    end = this.file.indexOf(' ', begin + 1);
    const numOfBonds = parseInt(this.file.substring(begin, end));

    return {atomCount: numOfAtoms, bondCount: numOfBonds};
  }
}
