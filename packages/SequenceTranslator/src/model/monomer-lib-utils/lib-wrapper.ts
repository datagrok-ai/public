/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMonomerLib} from '../../package';
import {SYNTHESIZERS, TECHNOLOGIES} from '../../model/const';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';

import {HardcodeTerminator} from '../hardcode-terminator';
const terminator = new HardcodeTerminator();

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
    const monomer = this.lib?.getMonomer('RNA', monomerName);
    return monomer?.molfile!;
  }

  getSmilesByName(monomerName: string): string {
    return terminator.getSmilesByName(monomerName);
  }

  get3PrimeTerminalSmiles(modificationName: string): string {
    return terminator.get3PrimeTerminalSmiles(modificationName);
  }

  isModification(monomerName: string): boolean {
    const monomerMolfile = this.getMolfileByName(monomerName);
    if (monomerMolfile !== undefined) {
      return (monomerMolfile.includes('MODIFICATION')) ? true : false;
    } else {
      const modifications = terminator.getModificationCodes();
      return (modifications.includes(monomerName)) ? true : false;
    }
  }

  getCodeToNameMap(sequence: string, format: string): Map<string, string> {
    return (new HardcodeTerminator().getCodeToNameMap(sequence, format));
  }

  getModificationCodes(): string[] {
    return terminator.getModificationCodes();
  }

  getTableForViewer(): DG.DataFrame {
    const monomerList = this.getAllMonomers();
    const formattedObjectsList = new Array(monomerList.length);
    for (let i = 0; i < monomerList.length; i++)
      formattedObjectsList[i] = this.formatMonomerForViewer(monomerList[i]);
    const df = DG.DataFrame.fromObjects(formattedObjectsList)!;
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
      // formattedObject['technologies'] = [...technologySet].toString();
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
}
