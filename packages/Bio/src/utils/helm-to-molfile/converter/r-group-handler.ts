import {HYDROGEN_SYMBOL, V2K_CONST} from './const';
import {MolfileAtoms} from './mol-atoms';
import {MolfileBonds} from './mol-bonds';
import {PositionInBonds} from './types';


export class RGroupHandler {
  constructor(rGroupLines: string[], private atoms: MolfileAtoms, private bonds: MolfileBonds) {
    this.rGroupIdToAtomicIndexMap = this.getRGroupIdToAtomicIdxMap(rGroupLines);
  }

  /** Relates R group id (starting from 1) to its atomic index within the
   * molfile  */
  rGroupIdToAtomicIndexMap: Map<number, number>;

  /** Maps R group id (starting from 1) to its position in the bond block  */
  private rGroupBondPositionMap = new Map<number, PositionInBonds>();

  getAtomicIdx(rGroupId: number): number | null {
    const atomicIdx = this.rGroupIdToAtomicIndexMap.get(rGroupId);
    return atomicIdx === undefined ? null : atomicIdx;
  }

  private removeRGroupsFromAtomBlock(rGroupIds: number[]): void {
    rGroupIds.forEach((rgroupId) => {
      const atomicIdx = this.rGroupIdToAtomicIndexMap.get(rgroupId);
      if (atomicIdx === undefined)
        throw new Error(`Cannot find atomic index for R group ${rgroupId}`);
    });

    const rGroupAtomicIndices = Array.from(this.rGroupIdToAtomicIndexMap.entries()).filter(
      ([rGroupId, _]) => rGroupIds.includes(rGroupId)
    ).map(([_, atomicIdx]) => atomicIdx);
    this.atoms.deleteAtoms(rGroupAtomicIndices);
  }

  removeRGroups(rGroupIds: number[]): void {
    this.removeRGroupsFromAtomBlock(rGroupIds);

    rGroupIds.forEach((rGroupId) => {
      const dummyPosition = this.replaceRGroupInBondsByDummy(rGroupId);
      this.rGroupBondPositionMap.set(rGroupId, dummyPosition);
    });
  }

  /** Replace RGroups by -1, update associated maps, and return the position in
   * bond block */
  private replaceRGroupInBondsByDummy(rGroupId: number): PositionInBonds {
    const rGroupAtomicIdx = this.rGroupIdToAtomicIndexMap.get(rGroupId)!;

    if (this.rGroupBondPositionMap.has(rGroupId))
      throw new Error(`R group ${rGroupId} is already handled`);

    const positions = this.bonds.getPositionsInBonds(rGroupAtomicIdx + 1);
    if (positions.length === 0)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    if (positions.length > 1)
      throw new Error(`More than one position for R group ${rGroupId}`);

    const rGroupPosition = positions[0];

    this.bonds.replacePositionsInBondsByDummy([rGroupPosition]);
    this.bonds.removeAtomIdFromBonds(rGroupAtomicIdx + 1);
    this.removeRGroupFromAtomicIdxMap(rGroupId, rGroupAtomicIdx);

    return rGroupPosition;
  }

  private removeRGroupFromAtomicIdxMap(deletedId: number, deletedAtomicIdx: number): void {
    this.rGroupIdToAtomicIndexMap.delete(deletedId);
    for (const [rGroupId, rGroupAtomicIdx] of this.rGroupIdToAtomicIndexMap) {
      if (rGroupAtomicIdx > deletedAtomicIdx)
        this.rGroupIdToAtomicIndexMap.set(rGroupId, rGroupAtomicIdx - 1);
    }
  }

  private getRGroupIdToAtomicIdxMap(lines: string[]): Map<number, number> {
    function getAtomIdxToRgpIdxList(rgpLine: string): [number, number][] {
      const indices = rgpLine.split(/\s+/).filter((item) => item)
        .slice(3).map((item) => parseInt(item));
      const atomIdxToRgpIdxList = new Array<[number, number]>(indices.length / 2);
      for (let i = 0; i < indices.length; i += 2)
        atomIdxToRgpIdxList[i / 2] = [indices[i + 1], indices[i] - 1];
      return atomIdxToRgpIdxList;
    }

    const map = new Map<number, number>();

    const rgroupLines = lines.filter((line: string) => line.startsWith(V2K_CONST.RGP_LINE_START));
    rgroupLines.forEach((line: string) => {
      const atomIdxToRgpIdxList = getAtomIdxToRgpIdxList(line);
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

    const rGroupAtomicIndices = this.atoms.getRGroupAtomicIndices();
    const unaccounted = rGroupAtomicIndices.filter((idx) => !Array.from(map.values()).includes(idx));
    if (unaccounted.length !== 0)
      throw new Error(`Unaccounted R group indices: ${unaccounted}`);

    return map;
  }

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    const position = this.rGroupBondPositionMap.get(rGroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    const {bondLineIdx} = position;
    this.bonds.deleteBondLines([bondLineIdx]);
    this.rGroupBondPositionMap.delete(rGroupId);
    this.rGroupIdToAtomicIndexMap.delete(rGroupId);
    // update values of other positions
    this.rGroupBondPositionMap.forEach((position) => {
      if (position.bondLineIdx > bondLineIdx)
        position.bondLineIdx -= 1;
    });
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void {
    const position = this.rGroupBondPositionMap.get(rGroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rGroupId}`);
    const {bondLineIdx, nodeIdx} = position;
    this.bonds.bondedAtoms[bondLineIdx][nodeIdx] = externalAtom;
  }

  /** Atom id is molfile id starting from 1  */
  getAttachmentAtomIdByRGroupId(rgroupId: number): number {
    const position = this.rGroupBondPositionMap.get(rgroupId);
    if (!position)
      throw new Error(`Cannot find position for R group ${rgroupId}`);
    const {bondLineIdx, nodeIdx} = position;
    return this.bonds.bondedAtoms[bondLineIdx][(nodeIdx + 1) % 2];
  }

  /** WARNING: capping RGRoups and deletion of the bonded ones don't commute */
  capRGroups(capGroupElements: string[]): void {
    this.rGroupIdToAtomicIndexMap.forEach((atomicIdx, rGroupId) => {
      const element = capGroupElements[rGroupId - 1];
      if (element === HYDROGEN_SYMBOL) {
        this.removeRGroups([rGroupId]);
        this.deleteBondLineWithSpecifiedRGroup(rGroupId);
      } else
        this.atoms.replaceElementSymbol(atomicIdx, element);
    });
  }
}

