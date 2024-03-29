import * as DG from 'datagrok-api/dg';

import {_package} from '../../../../package';
import {DEFAULT_FORMATS} from '../const';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';

import {HELM_REQUIRED_FIELDS as REQ, HELM_OPTIONAL_FIELDS as OPT} from '@datagrok-libraries/bio/src/utils/const';
import {META_FIELDS as MET} from './const';
import {CODES_TO_SYMBOLS_DICT} from '../data-loader/json-loader';

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

  private formatMonomerForViewer(sourceObj: Monomer): {[key: string]: string} {
    const formattedObject: {[key: string]: string} = {};
    formattedObject[REQ.NAME] = sourceObj[REQ.SYMBOL];
    formattedObject[REQ.SYMBOL] = sourceObj[REQ.SYMBOL];
    formattedObject[REQ.MOLFILE] = sourceObj[REQ.MOLFILE];
    const formats = this.getAllFormats();
    formats.forEach((format) => {
      if (format === DEFAULT_FORMATS.HELM)
        return;
      const map = CODES_TO_SYMBOLS_DICT[format];
      const codes = Object.keys(map).filter((code) => map[code] === sourceObj.symbol);
      formattedObject[format] = codes.join(', ');
    });

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

  static getInstance(): MonomerLibWrapper {
    if (MonomerLibWrapper.instance === undefined)
      MonomerLibWrapper.instance = new MonomerLibWrapper();
    return MonomerLibWrapper.instance!;
  }

  getMolfileBySymbol(monomerSymbol: string): string {
    const monomer = this.getMonomer(monomerSymbol);
    return monomer.molfile;
  }

  getNaturalAnalogBySymbol(monomerSymbol: string): string {
    const monomer = this.getMonomer(monomerSymbol);
    const naturalAnalog = monomer.naturalAnalog;
    if (!naturalAnalog)
      throw new Error(`ST: no natural analog for ${monomerSymbol}`);
    return naturalAnalog!;
  }

  // todo: a better criterion
  isModification(monomerSymbol: string): boolean {
    const molfile = this.getMolfileBySymbol(monomerSymbol);
    return (molfile.includes('MODIFICATION')) ? true : false;
  }

  getCodeToSymbolMap(format: string): Map<string, string> {
    return new Map<string, string>(Object.entries(CODES_TO_SYMBOLS_DICT[format]));
  }

  getCodesByFormat(format: string): string[] {
    return Object.keys(CODES_TO_SYMBOLS_DICT[format]);
  }

  getAllFormats(): string[] {
    return Object.keys(CODES_TO_SYMBOLS_DICT);
  }

  getTableForViewer(): DG.DataFrame {
    const formattedObjects = this.allMonomers.map((monomer) => this.formatMonomerForViewer(monomer));
    const df = DG.DataFrame.fromObjects(formattedObjects)!;
    return df;
  }

  getCodesToWeightsMap(): Map<string, number> {
    const codesToWeightsMap = new Map<string, number>();
    Object.entries(CODES_TO_SYMBOLS_DICT).forEach(([_, dict]) => {
      Object.entries(dict).forEach(([code, monomerSymbol]) => {
        const monomer = this.getMonomer(monomerSymbol);
        const weight = monomer[OPT.META]?.[MET.MOLWEIGHT];
        codesToWeightsMap.set(code, weight);
      });
    });
    return codesToWeightsMap;
  }
}
