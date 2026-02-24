import {IMonomerLib, IMonomerLibBase, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

import {Helm} from './helm';
import {MolfileWrapper} from './mol-wrapper';
import {MolfileWrapperFactory} from './mol-wrapper-factory';
import {CapGroupInfo} from './types';

/** Returns true if the string is a valid single element symbol (e.g. 'H', 'O', 'C', 'Cl') */
function isSimpleElement(s: string): boolean {
  return /^[A-Z][a-z]?$/.test(s);
}

export class MonomerWrapper {
  private readonly molfileWrapper: MolfileWrapper;
  private capGroupInfo: CapGroupInfo[] = [];
  private static molfileV2KToV3KCache: Map<string, string> = new Map();
  constructor(
    public readonly monomerSymbol: string,
    public readonly monomerIdx: number,
    private helm: Helm,
    shift: { x: number, y: number },
    private readonly rdKitModule: RDModule,
    private readonly monomerLib: IMonomerLibBase
  ) {
    const libraryMonomerObject = this.getLibraryMonomerObject();

    let molfile = libraryMonomerObject.molfile;
    if (MolfileHandler.isMolfileV2K(molfile))
      molfile = this.convertMolfileToV3KFormat(molfile, monomerSymbol, rdKitModule);

    this.molfileWrapper = MolfileWrapperFactory.getInstance(molfile, monomerSymbol);
    this.capGroupInfo = this.getCapGroupInfo(libraryMonomerObject);

    this.removeRGroups(helm.bondedRGroupsMap[monomerIdx]!);
    this.capRemainingRGroups();

    this.shiftCoordinates(shift);
  }

  public get atomCount() { return this.molfileWrapper.atomCount; }

  public get bondCount() { return this.molfileWrapper.bondCount; }

  private convertMolfileToV3KFormat(molfileV2K: string, monomerSymbol: string, rdKitModule: RDModule): string {
    if (MonomerWrapper.molfileV2KToV3KCache.has(molfileV2K))
      return MonomerWrapper.molfileV2KToV3KCache.get(molfileV2K)!;
    let mol: RDMol | null = null;
    try {
      mol = rdKitModule.get_mol(molfileV2K, JSON.stringify({mergeQueryHs: true}));
      if (mol) {
        const res = mol.get_v3Kmolblock();
        MonomerWrapper.molfileV2KToV3KCache.set(molfileV2K, res);
        return res;
      } else
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

  private getCapGroupInfo(
    libraryMonomerObject: Monomer
  ): CapGroupInfo[] {
    const rgroups = libraryMonomerObject.rgroups;
    return rgroups.map((rgroup, ind) => {
      const smiles = rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES] ||
        // WARNING: ignore because both key variants coexist in HELM Core Library!
        // @ts-ignore
        rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE];
      let rgroupId = rgroup[HELM_RGROUP_FIELDS.LABEL][1];
      if (!rgroupId || !parseInt(rgroupId) || isNaN(parseInt(rgroupId))) {
        // try to parse it from smiles, which can look like '[H][*:1]', 'O[*:2]', 'C=C[*:3]'
        const match = smiles?.match(/\[\*:(\d)\]/);
        if (match && match[1])
          rgroupId = match[1];
      }
      if (!rgroupId || !parseInt(rgroupId) || isNaN(parseInt(rgroupId)))
        rgroupId = `${ind + 1}`; // fallback to index-based id, starting from 1
      // extract the element symbol
      const element = smiles.replace(/(\[|\]|\*|:|\d)/g, '');
      return {element, smiles, isSimple: isSimpleElement(element), rGroupId: parseInt(rgroupId)};
    });
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
    this.molfileWrapper.capRGroups(this.capGroupInfo, this.rdKitModule);
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

