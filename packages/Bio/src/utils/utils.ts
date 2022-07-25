import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '@datagrok-libraries/bio/src/viewers/web-logo';
import * as grok from 'datagrok-api/grok';
import {
  CAP_GROUP_NAME, CAP_GROUP_SMILES, jsonSdfMonomerLibDict, MONOMER_ENCODE_MAX, MONOMER_ENCODE_MIN, MONOMER_SYMBOL,
  RGROUP_ALTER_ID, RGROUP_FIELD, RGROUP_LABEL, SDF_MONOMER_NAME
} from '../const';

export const HELM_CORE_LIB_FILENAME = '/samples/HELMCoreLibrary.json';
export const HELM_CORE_LIB_MONOMER_SYMBOL = 'symbol';
export const HELM_CORE_LIB_MOLFILE = 'molfile';
export const HELM_CORE_FIELDS = ['symbol', 'molfile', 'rgroups', 'name'];


export function encodeMonomers(col: DG.Column): DG.Column | null {
  let encodeSymbol = MONOMER_ENCODE_MIN;
  const monomerSymbolDict:  { [key: string]: number }= {};
  const units = col.tags[DG.TAGS.UNITS];
  const sep = col.getTag('separator');
  const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, sep);
  const encodedStringArray = [];
  for (let i = 0; i < col.length; ++i) {
    let encodedMonomerStr = '';
    const monomers = splitterFunc(col.get(i));
    monomers.forEach(m => {
      if(!monomerSymbolDict[m]) {
        if(encodeSymbol > MONOMER_ENCODE_MAX) {
          grok.shell.error(`Not enougth symbols to encode monomers`);
          return null;
        }
        monomerSymbolDict[m] = encodeSymbol;
        encodeSymbol++;
      }
      encodedMonomerStr += String.fromCodePoint(monomerSymbolDict[m]);
    })
    encodedStringArray.push(encodedMonomerStr);
  }
  return DG.Column.fromStrings('encodedMolecules', encodedStringArray);
}

export function getMolfilesFromSeq(col: DG.Column, monomersLibObject: any[]): any[][] | null {
  const units = col.tags[DG.TAGS.UNITS];
  const sep = col.getTag('separator');
  const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, sep);
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
  const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, sep);
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
      dict[it[HELM_CORE_LIB_MONOMER_SYMBOL]] = monomerObject;
    }
  });
  return dict;
}


export function createJsonMonomerLibFromSdf(table: DG.DataFrame): any {
  const resultLib = [];
  for (let i = 0; i < table.rowCount; i++) {
    const monomer: { [key: string]: string | any } = {};
    Object.keys(jsonSdfMonomerLibDict).forEach((key) => {
      if (key === MONOMER_SYMBOL) {
        const monomerSymbol = table.get(jsonSdfMonomerLibDict[key], i);
        monomer[key] = monomerSymbol === '.' ? table.get(SDF_MONOMER_NAME, i) : monomerSymbol;
      } else if (key === RGROUP_FIELD) {
        const rgroups = table.get(jsonSdfMonomerLibDict[key], i).split('\n');
        const jsonRgroups: any[] = [];
        rgroups.forEach((g: string) => {
          const rgroup: { [key: string]: string | any } = {};
          const altAtom = g.substring(g.lastIndexOf(']') + 1);
          const radicalNum = g.match(/\[R(\d+)\]/)![1];
          rgroup[CAP_GROUP_SMILES] = altAtom === 'H' ? `[*:${radicalNum}][H]` : `O[*:${radicalNum}]`;
          rgroup[RGROUP_ALTER_ID] = altAtom === 'H' ? `R${radicalNum}-H` : `R${radicalNum}-OH`;
          rgroup[CAP_GROUP_NAME] = altAtom === 'H' ? `H` : `OH`;
          rgroup[RGROUP_LABEL] = `R${radicalNum}`;
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