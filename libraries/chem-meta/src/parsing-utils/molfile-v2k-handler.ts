import {AtomAndBondCounts} from './chemical-table-parser-base';
import {MolfileHandlerBase} from './molfile-handler-base';
import {V2K_CONST} from '../formats/molfile-v2k-const';

export class MolfileV2KHandler extends MolfileHandlerBase {
  constructor(molfile: string) {
    super(molfile);
  }

  public static validate(molfile: string): boolean {
    return (molfile.indexOf(V2K_CONST.HEADER) !== -1 &&
      molfile.indexOf(V2K_CONST.END) !== -1);
  }

  protected shiftIdxToAtomType(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.ATOM_TYPE_COL);
  }

  protected getCountsLineIdx(): number {
    let idx = 0;
    for (let i = 0; i < V2K_CONST.NUM_OF_HEADER_LINES; ++i)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  protected getAtomBlockIdx(): number {
    let idx = this.getCountsLineIdx();
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  protected shiftIdxToXColumn(lineStartIdx: number): number {
    return this.getNextColumnIdx(lineStartIdx);
  }

  protected shiftIdxToBondedAtomsPair(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.FIRST_BONDED_ATOM_COL);
  }

  protected shiftIdxToBondType(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V2K_CONST.BOND_TYPE_COL);
  }

  protected getBondBlockIdx(): number {
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < this.atomCount; i++)
      idx = this.getNextLineIdx(idx);
    return idx;
  }

  protected parseAtomAndBondCounts(): AtomAndBondCounts {
    let begin = this.getCountsLineIdx();
    let end = begin + V2K_CONST.NUM_OF_COUNTS_DIGITS;
    const atomCount = parseInt(this.fileContent.substring(begin, end));
    begin = end;
    end += V2K_CONST.NUM_OF_COUNTS_DIGITS;
    const bondCount = parseInt(this.fileContent.substring(begin, end));
    return {atomCount: atomCount, bondCount: bondCount};
  };

  protected queryCriterion(idx: number): boolean {
    return this.fileContent[idx] === 'R' ||
      (this.fileContent[idx] === 'L' && !this.isAlpha(this.fileContent[idx + 1]));
  }
}
