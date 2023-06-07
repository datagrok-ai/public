/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {FORMAT, TECHNOLOGIES} from '../const';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';

import {HELM_REQUIRED_FIELDS as REQ, HELM_OPTIONAL_FIELDS as OPT} from '@datagrok-libraries/bio/src/utils/const';
import {META_FIELDS as MET} from './const';

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

  getMolfileBySymbol(monomerSymbol: string): string {
    const monomer = this.getMonomer(monomerSymbol);
    return monomer.molfile;
  }

  getSmilesBySymbol(monomerSymbol: string): string {
    const monomer = this.getMonomer(monomerSymbol);
    return monomer.smiles;
  }

  get3PrimeTerminalSmiles(modificationSymbol: string): string {
    if (!this.isModification(modificationSymbol))
      throw new Error(`SequenceTranslator: ${modificationSymbol} is not a modification`);
    const monomer = this.getMonomer(modificationSymbol);
    return monomer[OPT.META]![MET.TERMINAL_SMILES];
  }

  // todo: a better criterion
  isModification(monomerSymbol: string): boolean {
    const molfile = this.getMolfileBySymbol(monomerSymbol);
    return (molfile.includes('MODIFICATION')) ? true : false;
  }

  getCodeToSymbolMap(format: string): Map<string, string> {
    const codeToSymbolMap = new Map<string, string>;

    for (const monomer of this.allMonomers) {
      const monomerSymbol = monomer[REQ.SYMBOL];
      const codes = this.getCodesObject(monomer);
      if (Object.keys(codes).includes(format)) {
        if (Array.isArray(codes[format])) {
          const arr = codes[format] as string[];
          arr.forEach((code) => codeToSymbolMap.set(code, monomerSymbol));
        } else {
          const obj = codes[format] as TechnologiesObject;
          for (const technology in obj) {
            for (const code of obj[technology])
              codeToSymbolMap.set(code, monomerSymbol);
          }
        }
      }
    }
    return codeToSymbolMap;
  }

  getModificationGCRSCodes(): string[] {
    let result: string[] = [];
    const modifications = this.getMonomersByFormat(FORMAT.GCRS)
      .filter((monomer) => this.isModification(monomer[REQ.SYMBOL]));
    for (const monomer of modifications) {
      const codes = this.getCodesObject(monomer)[FORMAT.GCRS] as TechnologiesObject;
      for (const technology in codes)
        result = result.concat(codes[technology]);
    }
    return result;
  }

  private getCodesObject(monomer: Monomer): Codes {
    const meta: Meta = monomer[OPT.META] as Meta;
    return meta[MET.CODES] as Codes;
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

  // todo: refactor, unify with getCodesByFormat
  getGcrsCodesWithoutLinkages(): string[] {
    let codes: string[] = [];
    const format = FORMAT.GCRS;
    const monomers = this.getMonomersByFormat(FORMAT.GCRS);
    for (const monomer of monomers) {
      if (monomer.name.includes('linkage'))
        continue;
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

  private formatMonomerForViewer(sourceObj: Monomer): {[key: string]: string} {
    const formattedObject: {[key: string]: string} = {};
    formattedObject[REQ.NAME] = sourceObj[REQ.NAME];
    formattedObject[REQ.MOLFILE] = sourceObj[REQ.MOLFILE];

    const meta = sourceObj[OPT.META] as Meta;
    const codes = meta[MET.CODES] as Codes;

    for (const synthesizer of Object.values(FORMAT)) {
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
    const polymerTypes = this.lib.getPolymerTypes();
    let result: Monomer[] = [];
    for (const polymerType of polymerTypes) {
      const monomerSymbols = this.lib.getMonomerSymbolsByType(polymerType);
      const monomersByType: Monomer[] = monomerSymbols
        .map((monomerSymbol) => this.lib.getMonomer(polymerType, monomerSymbol))
        .filter((monomer): monomer is Monomer => monomer !== null);
      result = result.concat(monomersByType);
    }
    return result;
  }

  private getMonomer(monomerSymbol: string): Monomer {
    const monomer = this.lib.getMonomer('RNA', monomerSymbol);
    if (monomer === undefined)
      throw new Error(`SequenceTranslator: no monomer with symbol ${monomerSymbol}`);
    return monomer!;
  }

  getCodesToWeightsMap(): Map<string, number> {
    const codesToWeightsMap = new Map<string, number>();
    const monomers = this.getAllMonomers();
    for (const monomer of monomers) {
      const codesObj = this.getCodesObject(monomer);
      const weight = monomer[OPT.META]?.[MET.MOLWEIGHT];
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
