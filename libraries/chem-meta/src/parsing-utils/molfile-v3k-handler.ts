import {AtomAndBondCounts} from './chemical-table-parser-base';
import {MolfileHandlerBase} from './molfile-handler-base';
import {V3K_CONST} from '../formats/molfile-v3k-const';

export class MolfileV3KHandler extends MolfileHandlerBase {
  constructor(molfile: string) {
    super(molfile);
    this.init(molfile);
  }

  protected shiftIdxToAtomType(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.ATOM_TYPE_COL);
  }

  protected getCountsLineIdx(): number {
    return this.fileContent.indexOf(V3K_CONST.BEGIN_COUNTS_LINE);
  }

  protected getAtomBlockIdx(): number {
    let idx = this.fileContent.indexOf(V3K_CONST.BEGIN_ATOM_BLOCK);
    idx = this.getNextLineIdx(idx);
    return idx;
  }

  protected shiftIdxToXColumn(lineStartIdx: number): number {
    let idx = this.shiftIdxToAtomType(lineStartIdx);
    if (this.isQuote(idx)) {
      idx = this.getNextIdenticalChar(idx);
      idx = this.getNextColumnIdx(idx);
      return idx;
    }
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.X_COL);
  }

  protected shiftIdxToBondedAtomsPair(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.FIRST_BONDED_ATOM_COL);
  }

  protected shiftIdxToBondType(lineStartIdx: number): number {
    return this.shiftIdxToSpecifiedColumn(lineStartIdx, V3K_CONST.BOND_TYPE_COL);
  }

  protected getBondBlockIdx(): number {
    return this.getNextLineIdx(this.fileContent.indexOf(V3K_CONST.BEGIN_BOND_BLOCK));
  }

  public static validate(molfile: string): boolean {
    return (molfile.indexOf(V3K_CONST.HEADER) !== -1 &&
    molfile.indexOf(V3K_CONST.END) !== -1);
  }

  protected parseAtomAndBondCounts(): AtomAndBondCounts {
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

  protected queryCriterion(idx: number): boolean {
    return this.fileContent[idx] === 'R' || this.isQuote(idx) ||
      (this.fileContent[idx] === 'L' && !this.isAlpha(this.fileContent[idx + 1]));
  }
}
