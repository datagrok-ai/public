import {V3K_CONST} from '../formats/molfile-const';
import {AtomAndBondCounts} from './chemical-table-parser-base';
import {MolfileHandlerBase} from './molfile-handler-base';

export class MolfileV3KHandler extends MolfileHandlerBase {
  constructor(molfile: string) {
    super(molfile);
    this.init(molfile);
  }

  getAtomLines(): string[] {
    const atomBlockIdx = this.getAtomBlockIdx();
    const end = this.getBondBlockIdx();
    return this.fileContent.substring(atomBlockIdx, end).split('\n').slice(0, this.atomCount);
  }

  getBondLines(): string[] {
    const bondBlockIdx = this.getBondBlockIdx();
    return this.fileContent.substring(bondBlockIdx).split('\n').slice(0, this.bondCount);
  }

  getRGroupIdToAtomicIdxMap(): Map<number, number> {
    const map = new Map<number, number>();
    this.getAtomLines().forEach((line, idx) => {
      const rGroupMatches = line.match(/RGROUPS=\(([\d\s]+)\)/);

      if (rGroupMatches) {
        const rGroupData = rGroupMatches[1].split(/\s+/).map((rgId) => parseInt(rgId));
        // todo: handle cases when there are more than 2 r groups
        if (rGroupData.length > 2)
          throw new Error(`R group data ${rGroupData} has more than 2 elements`);
        const rGroupId = rGroupData[1];
        if (map.has(rGroupId))
          throw new Error(`R group ${rGroupId} is already in the map`);

        map.set(rGroupId, idx);
      }
    });

    return map;
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

  static isValidMolfile(molfile: string): boolean {
    return (molfile.indexOf(V3K_CONST.TYPE) !== -1 &&
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
}
