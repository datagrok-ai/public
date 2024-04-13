import {HELM_POLYMER_TYPE, HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibManager} from '../../monomer-lib/lib-manager';
import {MolfileWrapperBase} from './mol-wrapper-base';
import {MolfileWrapper} from './mol-wrapper-new';

export class MonomerWrapper {
  constructor(
    monomerSymbol: string,
    polymerType: HELM_POLYMER_TYPE,
  ) {
    const monomerLib = MonomerLibManager.instance.getBioLib();
    const monomer = monomerLib.getMonomer(polymerType, monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${monomerSymbol} is not found in the library`);
    this.molfileWrapper = MolfileWrapper.getInstance(monomer.molfile, monomerSymbol);
    this.capGroupElements = monomer.rgroups.map((rgroup) => {
      const smiles = rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES] ||
        // WARNING: ignore because both key variants coexist in HELM Core Library!
        // @ts-ignore
        rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE];
      // extract the element symbol
      return smiles.replace(/(\[|\]|\*|:|\d)/g, '');
    });
  }

  private molfileWrapper: MolfileWrapperBase;
  private capGroupElements: string[] = [];

  shiftCoordinates(shift: {x: number, y: number}): void {
    this.molfileWrapper.shiftCoordinates(shift);
  }

  getAtomLines(): string[] {
    return this.molfileWrapper.getAtomLines();
  }

  getBondLines(): string[] {
    return this.molfileWrapper.getBondLines();
  }

  removeBondedRGroups(rGroupIds: number[]): void {
    this.molfileWrapper.removeRGroups(rGroupIds);
  }

  capTrailingRGroups(): void {
    this.molfileWrapper.capRGroups(this.capGroupElements);
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, attachmentAtomIdx: number): void {
    this.molfileWrapper.replaceRGroupWithAttachmentAtom(rGroupId, attachmentAtomIdx);
  };

  getAttachmentAtomByRGroupId(rGroupId: number): number {
    const attachmentAtom = this.molfileWrapper.getAttachmentAtomByRGroupId(rGroupId);
    return attachmentAtom;
  }

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    this.molfileWrapper.deleteBondLineWithSpecifiedRGroup(rGroupId);
  }

  shiftBonds(shift: number): void {
    this.molfileWrapper.shiftBonds(shift);
  }
}

