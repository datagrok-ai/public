import {R_GROUP_ELEMENT_SYMBOL, V2K_CONST} from '../formats/molfile-const';
import {AtomAndBondCounts} from './chemical-table-parser-base';
import {MolfileHandlerBase} from './molfile-handler-base';

export class MolfileV2KHandler extends MolfileHandlerBase {
  constructor(molfile: string) {
    super(molfile);
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
    const lines = this.fileContent.split('\n');
    const rgroupLines = lines.filter((line: string) => line.startsWith(V2K_CONST.RGP_LINE_START));
    rgroupLines.forEach((line: string) => {
      const atomIdxToRgpIdxList = this.getAtomIdxToRgpIdxList(line);
      for (const [key, value] of atomIdxToRgpIdxList) {
        if (map.has(key))
          throw new Error(`R group ${key} is already in the map`);
        map.set(key, value);
      }
    });

    const atomAliasLinesIndices = lines.map((line: string, idx: number) => {
      if (line.startsWith(V2K_CONST.ATOM_ALIAS_LINE_START))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
    const atomAliasLines = atomAliasLinesIndices.map((idx) => lines[idx]);
    const atomAliasTextLines = atomAliasLinesIndices.map((idx) => lines[idx + 1]);
    atomAliasLines.forEach((line: string, idx: number) => {
      const rgpAtomIdx = parseInt(line.split(/\s+/)[1]) - 1;
      const rgpId = parseInt(atomAliasTextLines[idx].substring(1));
      if (map.has(rgpId))
        throw new Error(`R group ${rgpId} is already in the map`);
      map.set(rgpId, rgpAtomIdx);
    });

    const rGroupAtomicIndices = this.getRGroupAtomicIndices();
    const unaccounted = rGroupAtomicIndices.filter((idx) => !Array.from(map.values()).includes(idx));
    if (unaccounted.length !== 0)
      throw new Error(`Unaccounted R group indices: ${unaccounted}`);

    return map;
  }

  private getAtomIdxToRgpIdxList(rgpLine: string): [number, number][] {
    const indices = rgpLine.split(/\s+/).filter((item) => item)
      .slice(3).map((item) => parseInt(item));
    const atomIdxToRgpIdxList = new Array<[number, number]>(indices.length / 2);
    for (let i = 0; i < indices.length; i += 2)
      atomIdxToRgpIdxList[i / 2] = [indices[i + 1], indices[i] - 1];
    return atomIdxToRgpIdxList;
  }

  private getRGroupAtomicIndices(): number[] {
    return this.atomTypes.map((line: string, idx: number) => {
      if (line.includes(R_GROUP_ELEMENT_SYMBOL))
        return idx;
    }).filter((idx) => idx !== undefined) as number[];
  }

  static isValidMolfile(molfile: string): boolean {
    return (molfile.indexOf(V2K_CONST.TYPE) !== -1 &&
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
}
