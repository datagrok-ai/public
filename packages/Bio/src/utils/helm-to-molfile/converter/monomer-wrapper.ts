import {IMonomerLib, IMonomerLibBase, Monomer} from '@datagrok-libraries/bio/src/types';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

import {Helm} from './helm';
import {MolfileWrapper} from './mol-wrapper';
import {MolfileWrapperFactory} from './mol-wrapper-factory';

export class MonomerWrapper {
  private readonly molfileWrapper: MolfileWrapper;
  private capGroupElements: string[] = [];

  constructor(
    public readonly monomerSymbol: string,
    public readonly monomerIdx: number,
    private helm: Helm,
    shift: { x: number, y: number },
    rdKitModule: RDModule,
    private readonly monomerLib: IMonomerLibBase
  ) {
    const libraryMonomerObject = this.getLibraryMonomerObject();

    let molfile = libraryMonomerObject.molfile;
    if (MolfileHandler.isMolfileV2K(molfile))
      molfile = this.convertMolfileToV3KFormat(molfile, monomerSymbol, rdKitModule);

    this.molfileWrapper = MolfileWrapperFactory.getInstance(molfile, monomerSymbol);
    this.capGroupElements = this.getCapGroupElements(libraryMonomerObject);

    this.removeRGroups(helm.bondedRGroupsMap[monomerIdx]!);
    this.capRemainingRGroups();

    this.shiftCoordinates(shift);
  }

  public get atomCount() { return this.molfileWrapper.atomCount; }

  public get bondCount() { return this.molfileWrapper.bondCount; }

  private convertMolfileToV3KFormat(molfileV2K: string, monomerSymbol: string, rdKitModule: RDModule): string {
    let mol: RDMol | null = null;
    try {
      mol = rdKitModule.get_mol(molfileV2K, JSON.stringify({mergeQueryHs: true}));
      if (mol)
        return mol.get_v3Kmolblock();
      else
        throw new Error(`Cannot convert ${monomerSymbol} to molV3000`);
    } finally {
      mol?.delete();
    }
  }

  private getLibraryMonomerObject(): Monomer {
    const polymerType = this.helm.getPolymerTypeByMonomerIdx(this.monomerIdx);
    const monomer = this.monomerLib.getMonomer(polymerType, this.monomerSymbol);
    if (!monomer)
      throw new Error(`Monomer ${this.monomerSymbol} is not found in the library`);
    return monomer;
  }

  private getCapGroupElements(
    libraryMonomerObject: Monomer
  ): string[] {
    const rgroups = libraryMonomerObject.rgroups;
    const result = rgroups.map((rgroup) => {
      const smiles = rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES] ||
        // WARNING: ignore because both key variants coexist in HELM Core Library!
        // @ts-ignore
        rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE];
      // extract the element symbol
      return smiles.replace(/(\[|\]|\*|:|\d)/g, '');
    });

    return result;
  }

  private shiftCoordinates(shift: { x: number, y: number }): void {
    this.molfileWrapper.shiftCoordinates(shift);
  }

  getAtomLines(): string[] {
    return this.molfileWrapper.getAtomLines();
  }

  getBondLines(): string[] {
    return this.molfileWrapper.getBondLines();
  }

  private removeRGroups(rGroupIds: number[]): void {
    this.molfileWrapper.removeRGroups(rGroupIds);
  }

  private capRemainingRGroups(): void {
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

