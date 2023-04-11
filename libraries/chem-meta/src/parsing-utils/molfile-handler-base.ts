import {ChemicalTableParserBase, AtomAndBondCounts} from './chemical-table-parser-base';

export abstract class MolfileHandlerBase extends ChemicalTableParserBase {
  constructor(molfile: string) {
    super(molfile);
    this.init(molfile);
  }

  protected init(molfile: string) {
    super.init(molfile);
  }

  protected abstract parseAtomAndBondCounts(): AtomAndBondCounts;
  protected abstract getCountsLineIdx(): number;
  protected abstract getAtomBlockIdx(): number;
  protected abstract shiftIdxToXColumn(lineStartIdx: number): number;
  protected abstract shiftIdxToAtomType(lineStartIdx: number): number;
  protected abstract getBondBlockIdx(): number;
  protected abstract shiftIdxToBondedAtomsPair(lineStartIdx: number): number;
  protected abstract shiftIdxToBondType(lineStartIdx: number): number;
  protected abstract queryCriterion(idx: number): boolean;

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

  protected isQuote(idx: number): boolean {
    return this.fileContent.at(idx) === '\"' || this.fileContent.at(idx) === '\'';
  }

  protected getNextIdenticalChar(idx: number): number {
    const sym = this.fileContent.at(idx);
    return sym ? this.fileContent.indexOf(sym, idx + 1) : -1;
  }

  public isQuery(): boolean {
    const atomCount = this.atomCount;
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < atomCount; i++) {
      idx = this.shiftIdxToAtomType(idx);
      if (this.queryCriterion(idx))
        return true;
      idx = this.getNextLineIdx(idx);
    }
    return false;
  };

  protected isAlpha(sym: string): boolean {
    const charCode = sym.charCodeAt(0);
    return (charCode >= 65 && charCode <= 90) ||
      (charCode >= 97 && charCode <= 122);
  }
}
