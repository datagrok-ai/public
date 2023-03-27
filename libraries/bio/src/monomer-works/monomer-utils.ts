// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {
  HELM_FIELDS, HELM_CORE_FIELDS, RGROUP_FIELDS, jsonSdfMonomerLibDict,
  MONOMER_ENCODE_MAX, MONOMER_ENCODE_MIN, SDF_MONOMER_NAME
} from '../utils/const';
import {getSplitter, SplitterFunc, TAGS} from '../utils/macromolecule';
import {IMonomerLib, Monomer} from '../types/index';


export const expectedMonomerData = ['symbol', 'name', 'molfile', 'rgroups', 'polymerType', 'monomerType'];

export function encodeMonomers(col: DG.Column): DG.Column | null {
  let encodeSymbol = MONOMER_ENCODE_MIN;
  const monomerSymbolDict: { [key: string]: number } = {};
  const units = col.tags[DG.TAGS.UNITS];
  const sep = col.getTag(TAGS.separator);
  const splitterFunc: SplitterFunc = getSplitter(units, sep);
  const encodedStringArray = [];
  for (let i = 0; i < col.length; ++i) {
    let encodedMonomerStr = '';
    const monomers = splitterFunc(col.get(i));
    monomers.forEach((m) => {
      if (!monomerSymbolDict[m]) {
        if (encodeSymbol > MONOMER_ENCODE_MAX) {
          grok.shell.error(`Not enough symbols to encode monomers`);
          return null;
        }
        monomerSymbolDict[m] = encodeSymbol;
        encodeSymbol++;
      }
      encodedMonomerStr += String.fromCodePoint(monomerSymbolDict[m]);
    });
    encodedStringArray.push(encodedMonomerStr);
  }
  return DG.Column.fromStrings('encodedMolecules', encodedStringArray);
}

export function getMolfilesFromSeq(col: DG.Column, monomersLibObject: any[]): any[][] | null {
  const units = col.tags[DG.TAGS.UNITS];
  const sep = col.getTag('separator');
  const splitterFunc: SplitterFunc = getSplitter(units, sep);
  const monomersDict = createMomomersMolDict(monomersLibObject);
  const molFiles = [];
  for (let i = 0; i < col.length; ++i) {
    const macroMolecule = col.get(i);
    const monomers = splitterFunc(macroMolecule);
    const molFilesForSeq = [];
    for (let j = 0; j < monomers.length; ++j) {
      if (monomers[j]) {
        if (!monomersDict[monomers[j]]) {
          grok.shell.warning(`Monomer ${monomers[j]} is missing in HELM library. Structure cannot be created`);
          return null;
        }
        // what is the reason of double conversion?
        molFilesForSeq.push(JSON.parse(JSON.stringify(monomersDict[monomers[j]])));
      }
    }
    molFiles.push(molFilesForSeq);
  }
  return molFiles;
}

export function getMolfilesFromSingleSeq(cell: DG.Cell, monomersLibObject: any[]): any[][] | null {
  const units = cell.column.tags[DG.TAGS.UNITS];
  const sep = cell.column!.getTag('separator');
  const splitterFunc: SplitterFunc = getSplitter(units, sep);
  const monomersDict = createMomomersMolDict(monomersLibObject);
  const molFiles = [];
  const macroMolecule = cell.value;
  const monomers = splitterFunc(macroMolecule);
  const molFilesForSeq = [];
  for (let j = 0; j < monomers.length; ++j) {
    if (monomers[j]) {
      if (!monomersDict[monomers[j]]) {
        grok.shell.warning(`Monomer ${monomers[j]} is missing in HELM library. Structure cannot be created`);
        return null;
      }
      molFilesForSeq.push(JSON.parse(JSON.stringify(monomersDict[monomers[j]])));
    }
  }
  molFiles.push(molFilesForSeq);
  return molFiles;
}

export function createMomomersMolDict(lib: any[]): { [key: string]: string | any } {
  const dict: { [key: string]: string | any } = {};
  lib.forEach((it) => {
    if (it['polymerType'] === 'PEPTIDE') {
      const monomerObject: { [key: string]: any } = {};
      HELM_CORE_FIELDS.forEach((field) => {
        monomerObject[field] = it[field];
      });
      dict[it[HELM_FIELDS.SYMBOL]] = monomerObject;
    }
  });
  return dict;
}

export function createJsonMonomerLibFromSdf(table: DG.DataFrame): any {
  const resultLib = [];
  for (let i = 0; i < table.rowCount; i++) {
    const monomer: { [key: string]: string | any } = {};
    Object.keys(jsonSdfMonomerLibDict).forEach((key) => {
      if (key === HELM_FIELDS.SYMBOL) {
        const monomerSymbol = table.get(jsonSdfMonomerLibDict[key], i);
        monomer[key] = monomerSymbol === '.' ? table.get(SDF_MONOMER_NAME, i) : monomerSymbol;
      } else if (key === HELM_FIELDS.RGROUPS) {
        const rgroups = table.get(jsonSdfMonomerLibDict[key], i).split('\n');
        const jsonRgroups: any[] = [];
        rgroups.forEach((g: string) => {
          const rgroup: { [key: string]: string | any } = {};
          const altAtom = g.substring(g.lastIndexOf(']') + 1);
          const radicalNum = g.match(/\[R(\d+)\]/)![1];
          rgroup[RGROUP_FIELDS.CAP_GROUP_SMILES] = altAtom === 'H' ? `[*:${radicalNum}][H]` : `O[*:${radicalNum}]`;
          rgroup[RGROUP_FIELDS.ALTER_ID] = altAtom === 'H' ? `R${radicalNum}-H` : `R${radicalNum}-OH`;
          rgroup[RGROUP_FIELDS.CAP_GROUP_NAME] = altAtom === 'H' ? `H` : `OH`;
          rgroup[RGROUP_FIELDS.LABEL] = `R${radicalNum}`;
          jsonRgroups.push(rgroup);
        });
        monomer[key] = jsonRgroups;
      } else {
        if ((jsonSdfMonomerLibDict as { [key: string]: string | any })[key])
          monomer[key] = table.get((jsonSdfMonomerLibDict as { [key: string]: string | any })[key], i);
      }
    });
    resultLib.push(monomer);
  }
  return resultLib;
}

export interface IMonomerLibHelper {
  /** Singleton monomer library */
  getBioLib(): IMonomerLib;

  /** (Re)Loads libraries based on settings in user storage {@link LIB_STORAGE_NAME} to singleton.
   * @param {boolean} reload Clean {@link monomerLib} before load libraries [false]
   */
  loadLibraries(reload?: boolean): Promise<void>;

  /** Reads library from file shares, handles .json and .sdf */
  readLibrary(path: string, fileName: string): Promise<IMonomerLib>;
}

export async function getMonomerLibHelper(): Promise<IMonomerLibHelper> {
  const funcList = DG.Func.find({package: 'Bio', name: 'getMonomerLibHelper'});
  if (funcList.length === 0)
    throw new Error('Package "Bio" must be installed for MonomerLibHelper.');

  const res: IMonomerLibHelper = (await funcList[0].prepare().call()).getOutputParamValue() as IMonomerLibHelper;
  return res;
}
