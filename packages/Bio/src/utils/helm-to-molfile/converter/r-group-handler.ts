import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {HYDROGEN_SYMBOL} from './const';
import {MolfileAtoms} from './mol-atoms';
import {MolfileBonds} from './mol-bonds';
import {CapGroupInfo, PositionInBonds} from './types';


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
    return atomicIdx == undefined ? null : atomicIdx;
  }

  private removeRGroupsFromAtomBlock(rGroupIds: number[]): void {
    rGroupIds.forEach((rgroupId) => {
      const atomicIdx = this.rGroupIdToAtomicIndexMap.get(rgroupId);
      if (atomicIdx == undefined)
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

  /** WARNING: capping RGroups and deletion of the bonded ones don't commute */
  capRGroups(capGroupInfo: CapGroupInfo[], rdKitModule: RDModule): void {
    this.rGroupIdToAtomicIndexMap.forEach((atomicIdx, rGroupId) => {
      const info = capGroupInfo.find((info) => info.rGroupId === rGroupId) ?? capGroupInfo[rGroupId - 1];
      if (info.isSimple) {
        if (info.element === HYDROGEN_SYMBOL) {
          this.removeRGroups([rGroupId]);
          this.deleteBondLineWithSpecifiedRGroup(rGroupId);
        } else
          this.atoms.replaceRGroupSymbolByElement(atomicIdx, info.element);
      } else
        this.capWithComplexGroup(atomicIdx, info.smiles, rdKitModule);
    });
  }

  /** Cap an R-group with a multi-atom cap group by parsing the cap SMILES,
   * then inserting its atoms and bonds into the monomer molfile */
  private capWithComplexGroup(
    rGroupAtomicIdx: number, capSmiles: string, rdKitModule: RDModule
  ): void {
    // Replace [*:N] with placeholder element Xe so RDKit can parse the SMILES
    const PLACEHOLDER = 'Xe';
    const parsableSmiles = capSmiles.replace(/\[\*:\d+\]/g, `[${PLACEHOLDER}]`);
    const capMol = rdKitModule.get_mol(parsableSmiles);
    if (!capMol)
      throw new Error(`Cannot parse cap group SMILES: ${capSmiles}`);

    let capMolfile: string;
    try {
      capMolfile = capMol.get_v3Kmolblock();
    } finally {
      capMol.delete();
    }

    const capHandler = MolfileHandler.getInstance(capMolfile);
    const capAtomLines = capHandler.getAtomLines();
    const capBondPairs = capHandler.pairsOfBondedAtoms;
    const capBondLines = capHandler.getBondLines();
    const capX = capHandler.x;
    const capY = capHandler.y;
    const capAtomTypes = capHandler.atomTypes;

    // Find the placeholder atom (was the [*:N] attachment point)
    let dummyCapIdx = -1; // 0-based
    for (let i = 0; i < capAtomTypes.length; i++) {
      if (capAtomTypes[i] === PLACEHOLDER) {
        dummyCapIdx = i;
        break;
      }
    }
    if (dummyCapIdx === -1)
      throw new Error(`Cannot find placeholder atom in cap group SMILES: ${capSmiles}`);

    // Find the attachment atom (bonded to placeholder) and the bond connecting them
    let attachmentCapIdx = -1; // 0-based
    for (let i = 0; i < capBondPairs.length; i++) {
      const [a1, a2] = capBondPairs[i]; // 1-based
      if (a1 === dummyCapIdx + 1) {
        attachmentCapIdx = a2 - 1;
        break;
      }
      if (a2 === dummyCapIdx + 1) {
        attachmentCapIdx = a1 - 1;
        break;
      }
    }
    if (attachmentCapIdx === -1)
      throw new Error(`Cannot find attachment atom in cap group SMILES: ${capSmiles}`);

    // Compute coordinate translation: place cap attachment at R-group position
    const rGroupCoords = this.atoms.atomCoordinates[rGroupAtomicIdx];
    const tx = rGroupCoords.x - capX[attachmentCapIdx];
    const ty = rGroupCoords.y - capY[attachmentCapIdx];

    // Replace the R# atom symbol with the attachment atom's element
    const attachmentSymbol = capAtomTypes[attachmentCapIdx];
    this.atoms.replaceRGroupSymbolByElement(rGroupAtomicIdx, attachmentSymbol);

    // Build index mapping: cap 1-based → monomer 1-based
    const capToMonomer = new Map<number, number>();
    capToMonomer.set(attachmentCapIdx + 1, rGroupAtomicIdx + 1);

    // Append remaining cap atoms (excluding placeholder and attachment)
    let nextMonomerIdx = this.atoms.count + 1; // 1-based
    for (let i = 0; i < capAtomLines.length; i++) {
      if (i === dummyCapIdx || i === attachmentCapIdx) continue;
      const newX = capX[i] + tx;
      const newY = capY[i] + ty;
      this.atoms.appendAtomLine(capAtomLines[i], newX, newY);
      capToMonomer.set(i + 1, nextMonomerIdx);
      nextMonomerIdx++;
    }

    // Append cap bonds (excluding any bond involving the placeholder)
    for (let i = 0; i < capBondPairs.length; i++) {
      const [a1, a2] = capBondPairs[i]; // 1-based in cap
      if (a1 === dummyCapIdx + 1 || a2 === dummyCapIdx + 1) continue;
      const newA1 = capToMonomer.get(a1)!;
      const newA2 = capToMonomer.get(a2)!;
      this.bonds.appendBondLine(capBondLines[i], [newA1, newA2]);
    }
  }
}

