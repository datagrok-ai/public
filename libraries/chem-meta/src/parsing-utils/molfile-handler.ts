/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ChemicalTableParserBase, AtomAndBondCounts} from './chemical-table-parser';
import {MOLFILE_VERSION, V2K, V3K} from './molfile-parsing-const';

/** Class for handling MDL Molfiles V2000 and V3000 */
class MolfileHandler extends ChemicalTableParserBase {
  constructor(molfile: string) {
    super(molfile);
    this.molfileVersion = MolfileHandler.determineMolfileVersion(molfile);
    this.parseAtomAndBondCounts = (this.molfileVersion === MOLFILE_VERSION.V2000) ? this.parseAtomAndBondCountsV2K :
      this.parseAtomAndBondCountsV3K;
    this.getAtomBlockIdx = (this.molfileVersion === MOLFILE_VERSION.V2000) ? this.getAtomBlockIdxV2K :
      this.getAtomBlockIdxV3K;
    this.getCountsLineIdx = (this.molfileVersion === MOLFILE_VERSION.V2000) ? this.getCountsLineV2KIdx :
      this.getCountsLineV3KIdx;
    this.shiftIdxToXColumn = (this.molfileVersion === MOLFILE_VERSION.V2000) ? this.shiftIdxToXColumnV2K :
      this.shiftIdxToXColumnV3K;
    this.shiftIdxToAtomType = (this.molfileVersion === MOLFILE_VERSION.V2000) ? this.shiftIdxToAtomTypeV2K :
      this.shiftIdxToAtomTypeV3K;
  }

  private molfileVersion: MOLFILE_VERSION;
  protected parseAtomAndBondCounts: () => AtomAndBondCounts;
  protected getCountsLineIdx: () => number;
  protected getAtomBlockIdx: () => number;
  protected shiftIdxToXColumn: (lineStartIdx: number) => number;
  protected shiftIdxToAtomType: (idx: number) => number;

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

  protected getBondBlockIdx(): number {
    if (this.molfileVersion === MOLFILE_VERSION.V2000)
      return this.getBondBlockIdxV2K();
    else
      return this.getBondBlockIdxV3K();
  };

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
    let idx = this.getCountsLineIdx();
    // 3 is the # of digits allocated for each of the counts
    const atomCount = parseInt(this.file.substring(idx, idx + 3));
    idx += 3;
    const bondCount = parseInt(this.file.substring(idx, idx + 3));
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
