import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {HYDROGEN_SYMBOL} from './const';
import {MolfileAtoms} from './mol-atoms';
import {MolfileBonds} from './mol-bonds';
import {PositionInBonds} from './types';


export class RGroupHandler {
  constructor(molfileHandler: MolfileHandlerBase, private atoms: MolfileAtoms, private bonds: MolfileBonds) {
    this.rGroupIdToAtomicIndexMap = molfileHandler.getRGroupIdToAtomicIdxMap();
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

