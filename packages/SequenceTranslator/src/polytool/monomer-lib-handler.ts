/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HELM_FIELDS, DUMMY_MONOMER} from '@datagrok-libraries/bio/src/utils/const';
import {
  helmFieldsToPolyToolInputFields, R_GROUP_BLOCK_DUMMY
} from './const';
import {Monomer, RGroup} from '@datagrok-libraries/bio/src/types/index';

export class PolyToolMonomerLibHandler {
  constructor(private rawLib: any[]) { }

  isValid(): boolean {
    return this.rawLib.every((entry) => {
      return typeof entry === 'object' &&
        Object.values(helmFieldsToPolyToolInputFields).every((field) => {
          return field in entry &&
            typeof entry[field] === 'string';
        });
    });
  }

  getJsonMonomerLib(): any {
    const resultLib: any[] = [];
    this.rawLib.forEach((rawMonomer) => {
      const monomer = this.prepareMonomer(rawMonomer);
      resultLib.push(monomer);
    });
    return resultLib;
  }

  private prepareMonomer(rawMonomer: any): Monomer {
    const monomer: Monomer = {...DUMMY_MONOMER};

    Object.entries(helmFieldsToPolyToolInputFields).forEach(([key, value]) => {
      const monomerSymbol = rawMonomer[value] as string;
      //@ts-ignore
      monomer[key] = monomerSymbol;
    });

    let key = HELM_FIELDS.SMILES;
    const rawSmiles = rawMonomer[helmFieldsToPolyToolInputFields[key]];
    const smilesHandler = new SmilesHandler(rawSmiles);
    const smiles = smilesHandler.getCappedSmiles();
    monomer[key] = smiles;

    key = HELM_FIELDS.RGROUPS;
    monomer[key] = RGroupHandler.getRGroups(smilesHandler.getNumberOfRGroups());

    key = HELM_FIELDS.MOLFILE;
    monomer[key] = new MolBlockHandler(smilesHandler.getSmilesWithRGroups()).getMolfile();

    return monomer;
  }
}

class SmilesHandler {
  constructor(rawSmiles: string) {
    const regex = /\[R(\d+)\]/g;
    let i = 0;
    this.smilesWithRGroups = rawSmiles.replace(regex, (_, capturedDigit) => { ++i; return `[${capturedDigit}*]`; });
    this.numberOfRGroups = i;
  }

  private numberOfRGroups: number;
  private smilesWithRGroups: string;

  getSmilesWithRGroups(): string {
    return this.smilesWithRGroups;
  }

  getCappedSmiles(): string {
    const smiles = this.capRGroups();
    return smiles;
  }

  getNumberOfRGroups(): number {
    return this.numberOfRGroups;
  }

  private capRGroups(): string {
    let result = this.smilesWithRGroups.replace('[1*]', '[H:1]');
    result = result.replace('[2*]', '[OH:2]');
    return result.replace('[3*]', '[H:3]');
  }
}

class RGroupHandler {
  private constructor() {};

  static getRGroups(numberOfRGroups: number): RGroup[] {
    return R_GROUP_BLOCK_DUMMY.slice(0, numberOfRGroups);
  }
}

class MolBlockHandler {
  constructor(private smilesWithRGroups: string) { }

  getMolfile(): string {
    let molfile = DG.chem.convert(this.smilesWithRGroups, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    molfile = this.restoreRGPLine(molfile);
    molfile = this.fixRGroupSymbols(molfile);
    return molfile;
  }

  private restoreRGPLine(rawMolfile: string): string {
    return rawMolfile.replace('M  ISO', 'M  RGP');
  }

  private fixRGroupSymbols(molfile: string): string {
    return molfile.replace(/\bR\b/g, 'R#');
  }
}
