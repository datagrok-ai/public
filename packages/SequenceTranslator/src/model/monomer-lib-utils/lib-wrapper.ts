/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {SYNTHESIZERS, TECHNOLOGIES} from '../../model/const';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';

import {HELM_REQUIRED_FIELDS as REQ, HELM_OPTIONAL_FIELDS as OPT} from '@datagrok-libraries/bio/src/utils/const';

const TERMINAL_SMILES = 'threePrimeTerminalSmiles';

type TechnologiesObject = {
  [technology: string]: string[]
}

type Codes = {
  [synthesizer: string]: TechnologiesObject | string[]
}

type Meta = {
  [key: string]: string | Codes,
}

export class MonomerLibWrapper {
  // todo: dependency injection of monomer lib instead of getMonomeLib
  private constructor() {
    const lib = _package.monomerLib;
    if (lib === null)
      throw new Error('SequenceTranslator: monomer library is null');
    this.lib = lib!;
    this.allMonomers = this.getAllMonomers();
  }

  private lib: IMonomerLib;
  private static instance?: MonomerLibWrapper;
  private allMonomers: Monomer[];

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
    const codeToNameMap = new Map<string, string>;

    for (const monomer of this.allMonomers) {
      const name = monomer[REQ.NAME];
      const codes = this.getCodesObject(monomer);
      if (Object.keys(codes).includes(format)) {
        if (Array.isArray(codes[format])) {
          const arr = codes[format] as string[];
          arr.forEach((code) => codeToNameMap.set(code, name));
        } else {
          const obj = codes[format] as TechnologiesObject;
          for (const technology in obj) {
            for (const code of obj[technology])
              codeToNameMap.set(code, name);
          }
        }
      }
    }
    return codeToNameMap;
  }

  getModificationGCRSCodes(): string[] {
    let result: string[] = [];
    const modifications = this.getMonomersByFormat(SYNTHESIZERS.GCRS)
      .filter((monomer) => this.isModification(monomer.name));
    for (const monomer of modifications) {
      const codes = this.getCodesObject(monomer)[SYNTHESIZERS.GCRS] as TechnologiesObject;
      for (const technology in codes)
        result = result.concat(codes[technology]);
    }
    return result;
  }

  private getCodesObject(monomer: Monomer): Codes {
    const meta: Meta = monomer[OPT.META] as Meta;
    return meta['codes'] as Codes;
  };

  private getMonomersByFormat(format: string): Monomer[] {
    return this.allMonomers.filter((monomer) => Object.keys(monomer[OPT.META]?.codes).includes(format));
  }

  getCodesByFormat(format: string): string[] {
    let codes: string[] = [];
    const monomers = this.getMonomersByFormat(format);
    for (const monomer of monomers) {
      const codesObj = this.getCodesObject(monomer);
      if (Array.isArray(codesObj[format])) {
        const array = codesObj[format] as string[];
        codes = codes.concat(array);
      } else {
        for (const technology in codesObj[format]) {
          const obj = codesObj[format] as TechnologiesObject;
          codes = codes.concat(obj[technology]);
        }
      }
    }
    return codes;
  }

  getTableForViewer(): DG.DataFrame {
    const formattedObjects = this.allMonomers.map((monomer) => this.formatMonomerForViewer(monomer));
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

    const meta = sourceObj[FIELD.META] as Meta;
    const codes = meta[FIELD.CODES] as Codes;

    for (const synthesizer of Object.values(SYNTHESIZERS)) {
      const fieldName = synthesizer;
      const valuesList = [];
      if (codes[synthesizer] !== undefined) {
        if (Array.isArray(codes[synthesizer])) {
          const arr = codes[synthesizer] as string[];
          valuesList.push(arr.toString());
        } else {
          for (const technology of Object.values(TECHNOLOGIES)) {
            const obj = codes[synthesizer] as TechnologiesObject;
            if (obj[technology] !== undefined)
              valuesList.push(obj[technology].toString());
          }
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

  getCodesToWeightsMap(): Map<string, number> {
    const codesToWeightsMap = new Map<string, number>();
    const monomers = this.getAllMonomers();
    for (const monomer of monomers) {
      const codesObj = this.getCodesObject(monomer);
      const weight = monomer[OPT.META]?.['weight'];
      let codesArray: string[] = [];
      for (const synthesizer in codesObj){
        if (Array.isArray(codesObj[synthesizer])) {
          const arr = codesObj[synthesizer] as Array<string>;
          codesArray = codesArray.concat(arr);
        } else {
          for (const technology in codesObj[synthesizer]) {
            const obj = codesObj[synthesizer] as TechnologiesObject;
            codesArray = codesArray.concat(obj[technology]);
          }
        }
      }
      for (const code of codesArray) {
        codesToWeightsMap.set(code, weight);
      }
    }
    return codesToWeightsMap;
  }
}
