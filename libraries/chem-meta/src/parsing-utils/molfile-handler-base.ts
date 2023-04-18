import {ChemicalTableParserBase, AtomAndBondCounts} from './chemical-table-parser-base';
import { D_QUOTE, L, R, S_QUOTE } from './const';
import { isAlpha } from './utils';

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
    const letter = this.fileContent[idx].charCodeAt(0);
    return letter === S_QUOTE || letter === D_QUOTE;
  }

  protected getNextIdenticalChar(idx: number): number {
    const sym = this.fileContent[idx];
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


  protected queryCriterion(idx: number): boolean {
    const letter = this.fileContent[idx];
    const char = letter.charCodeAt(idx);
    return char === R || !isAlpha(letter) ||
      (char === L && !isAlpha(this.fileContent[idx + 1]));
  }
}
