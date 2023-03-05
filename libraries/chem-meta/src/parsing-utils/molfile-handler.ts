/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ChemicalTableParser, ChemicalTableParserBase, AtomAndBondCounts} from './chemical-table-parser';
import {MOLFILE_VERSION, V2K, V3K} from './molfile-parsing-const';

/** Class for handling MDL Molfiles V2000 and V3000 */
export class MolfileHandler extends ChemicalTableParserBase implements ChemicalTableParser {
  constructor(molfile: string) {
    super(molfile);
    this.init(molfile);
  }

  // Public members

  /** Init/reset the state of the handler for a new molfile  */
  public init(molfile: string) {
    super.init(molfile);

    const molfileVersion = MolfileHandler.determineMolfileVersion(molfile);
    const isV2K = (molfileVersion === MOLFILE_VERSION.V2000);

    this.parseAtomAndBondCounts = isV2K ? this.parseAtomAndBondCountsV2K : this.parseAtomAndBondCountsV3K;

    this.getAtomBlockIdx = isV2K ? this.getAtomBlockIdxV2K : this.getAtomBlockIdxV3K;

    this.getCountsLineIdx = isV2K ? this.getCountsLineV2KIdx : this.getCountsLineV3KIdx;

    this.shiftIdxToXColumn = isV2K ? this.shiftIdxToXColumnV2K : this.shiftIdxToXColumnV3K;

    this.shiftIdxToAtomType = isV2K ? this.shiftIdxToAtomTypeV2K :
      this.shiftIdxToAtomTypeV3K;

    this.getBondBlockIdx = isV2K ? this.getBondBlockIdxV2K : this.getBondBlockIdxV3K;
  }

  // Protected members

  protected parseAtomAndBondCounts!: () => AtomAndBondCounts;
  protected getCountsLineIdx!: () => number;
  protected getAtomBlockIdx!: () => number;
  protected shiftIdxToXColumn!: (lineStartIdx: number) => number;
  protected shiftIdxToAtomType!: (idx: number) => number;
  protected getBondBlockIdx!: () => number;

  // Private members

  /** Determine whether the file is V2000/V3000, or throw */
  private static determineMolfileVersion(molfile: string): MOLFILE_VERSION {
    if (MolfileHandler.validateV3K(molfile))
      return MOLFILE_VERSION.V3000;
    else if (MolfileHandler.validateV2K(molfile))
      return MOLFILE_VERSION.V3000;
    else
      throw new Error('Malformed molfile');
  }

  private shiftIdxToAtomTypeV2K(lineStartIdx: number): number {
    return this.shiftToSpecifiedColumn(lineStartIdx, 4);
  }

  private shiftIdxToAtomTypeV3K(lineStartIdx: number): number {
    return this.shiftToSpecifiedColumn(lineStartIdx, 3);
  }

  private getCountsLineV2KIdx(): number {
    let idx = 0;
    // 3 is for the number of header lines in molfiles
    for (let i = 0; i < 3; ++i)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getCountsLineV3KIdx(): number {
    return this.file.indexOf(V3K.BEGIN_COUNTS_LINE);
  }

  private getAtomBlockIdxV2K(): number {
    let idx = this.getCountsLineIdx();
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getAtomBlockIdxV3K(): number {
    let idx = this.file.indexOf(V3K.BEGIN_ATOM_BLOCK);
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  private shiftIdxToXColumnV2K(lineStartIdx: number): number {
    return this.getNextLineIdx(lineStartIdx);
  }

  private shiftIdxToXColumnV3K(lineStartIdx: number): number {
    return this.shiftToSpecifiedColumn(lineStartIdx, 5);
  }

  private getBondBlockIdxV2K(): number {
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < this.atomCount; i++)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  private getBondBlockIdxV3K(): number {
    return this.getNextLineIdx(this.file.indexOf(V3K.BEGIN_BOND_BLOCK));
  }

  // todo: devise a more reliable validation
  private static validateV3K(molfile: string): boolean {
    if (molfile.indexOf(V3K.HEADER) !== -1 &&
    molfile.indexOf(V3K.END) !== -1)
      return true;
    return false;
  }

  // todo: devise a more reliable validation
  private static validateV2K(molfile: string): boolean {
    if (molfile.indexOf(V2K.HEADER) !== -1 &&
    molfile.indexOf(V2K.END) !== -1)
      return true;
    return false;
  }

  private parseAtomAndBondCountsV2K(): AtomAndBondCounts {
    let begin = this.getCountsLineIdx();
    // 3 is the # of digits allocated for atom/bond counts
    let end = begin + 3;
    const atomCount = parseInt(this.file.substring(begin, end));
    begin = end;
    end += 3;
    const bondCount = parseInt(this.file.substring(begin, end));
    return {atomCount: atomCount, bondCount: bondCount};
  };

  private parseAtomAndBondCountsV3K(): AtomAndBondCounts {
    // parse atom count
    let begin = this.file.indexOf(V3K.BEGIN_COUNTS_LINE) + V3K.COUNTS_SHIFT;
    let end = this.file.indexOf(' ', begin + 1);
    const numOfAtoms = parseInt(this.file.substring(begin, end));

    // parse bond count
    begin = end + 1;
    end = this.file.indexOf(' ', begin + 1);
    const numOfBonds = parseInt(this.file.substring(begin, end));

    return {atomCount: numOfAtoms, bondCount: numOfBonds};
  }
}
