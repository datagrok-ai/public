/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMonomerLib} from '../../package';
import {DELIMITER, SYNTHESIZERS, TECHNOLOGIES} from '../../model/const';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';

import {HELM_REQUIRED_FIELDS as REQ, HELM_OPTIONAL_FIELDS as OPT} from '@datagrok-libraries/bio/src/utils/const';

import {isValidSequence} from '../code-converter/conversion-validation-tools';

import {HardcodeTerminator} from '../hardcode-terminator';
const terminator = new HardcodeTerminator();

const TERMINAL_SMILES = 'threePrimeTerminalSmiles';

type CodesField = {
  [synthesizer: string]: {
    [technology: string]: string[]
  }
}

type MetaField = {
  [key: string]: string | CodesField,
}

export class MonomerLibWrapper {
  // todo: dependency injection of monomer lib instead of getMonomeLib
  private constructor() {
    const lib = getMonomerLib();
    if (lib === null)
      throw new Error('SequenceTranslator: monomer library is null');
    this.lib = lib!;
  }

  private lib: IMonomerLib;
  private static instance?: MonomerLibWrapper;

  static getInstance(): MonomerLibWrapper {
    if (MonomerLibWrapper.instance === undefined)
      MonomerLibWrapper.instance = new MonomerLibWrapper();
    return MonomerLibWrapper.instance!;
  }

  getMolfileByName(monomerName: string): string {
    const monomer = this.getMonomer(monomerName);
    return monomer.molfile;
  }

  getSmilesByName(monomerName: string): string {
    const monomer = this.getMonomer(monomerName);
    return monomer.smiles;
    // return terminator.getSmilesByName(monomerName);
  }

  get3PrimeTerminalSmiles(modificationName: string): string {
    if (!this.isModification(modificationName))
      throw new Error(`SequenceTranslator: ${modificationName} is not a modification`);
    const monomer = this.getMonomer(modificationName);
    return monomer[OPT.META]![TERMINAL_SMILES];
  }

  isModification(monomerName: string): boolean {
    const molfile = this.getMolfileByName(monomerName);
    return (molfile.includes('MODIFICATION')) ? true : false;
  }

  getCodeToNameMap(format: string): Map<string, string> {
    const monomers = this.getAllMonomers();
    const codeToNameMap = new Map<string, string>;

    for (const monomer of monomers) {
      const name = monomer[REQ.NAME];
      const codes = monomer[OPT.META]?.codes;
      console.log('name, codes:', name, codes);
      if (Object.keys(codes).includes(format)) {
        for (const technology of Object.keys(codes[format])) {
          for (const code of codes[format][technology])
            codeToNameMap.set(code, name);
        }
      }
    }
    return codeToNameMap;
  }

  getModificationCodes(): string[] {
    return terminator.getModificationCodes();
  }

  getTableForViewer(): DG.DataFrame {
    const monomerList = this.getAllMonomers();
    const formattedObjects = monomerList.map((monomer) => this.formatMonomerForViewer(monomer));
    const df = DG.DataFrame.fromObjects(formattedObjects)!;
    return df;
  }

  private formatMonomerForViewer(sourceObj: Monomer): {[key: string]: string} {
    const enum FIELD {
      NAME = 'name',
      MOLFILE = 'molfile',
      CODES = 'codes',
      META = 'meta',
    }

    const formattedObject: {[key: string]: string} = {};
    formattedObject[FIELD.NAME] = sourceObj[FIELD.NAME];
    formattedObject[FIELD.MOLFILE] = sourceObj[FIELD.MOLFILE];

    const meta = sourceObj[FIELD.META] as MetaField;
    const codes = meta[FIELD.CODES] as CodesField;

    for (const synthesizer of Object.values(SYNTHESIZERS)) {
      const fieldName = synthesizer;
      const valuesList = [];
      for (const technology of Object.values(TECHNOLOGIES)) {
        if (codes[synthesizer] !== undefined) {
          if (codes[synthesizer][technology] !== undefined)
            valuesList.push(codes[synthesizer][technology].toString());
        }
      }
      formattedObject[fieldName] = valuesList.toString();
    }
    return formattedObject;
  }

  private getAllMonomers(): Monomer[] {
    const monomerTypes = this.lib.getTypes();
    let result: Monomer[] = [];
    for (const monomerType of monomerTypes) {
      const monomerNames = this.lib.getMonomerNamesByType(monomerType);
      const monomersByType: Monomer[] = monomerNames
        .map((monomerName) => this.lib.getMonomer(monomerType, monomerName))
        .filter((monomer): monomer is Monomer => monomer !== null);
      result = result.concat(monomersByType);
    }
    return result;
  }

  private getMonomer(monomerName: string): Monomer {
    const monomer = this.lib.getMonomer('RNA', monomerName);
    if (monomer === undefined)
      throw new Error(`SequenceTranslator: no monomer with name ${monomerName}`);
    return monomer!;
  }
}
