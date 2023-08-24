// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {
  HELM_FIELDS, HELM_CORE_FIELDS, HELM_RGROUP_FIELDS, jsonSdfMonomerLibDict,
  MONOMER_ENCODE_MAX, MONOMER_ENCODE_MIN, SDF_MONOMER_NAME, helmFieldsToEnumeratorInputFields, rGroupsDummy
} from '../utils/const';
import {IMonomerLib} from '../types/index';
import {ISeqSplitted, SplitterFunc} from '../utils/macromolecule/types';
import {UnitsHandler} from '../utils/units-handler';
import {splitAlignedSequences} from '../utils/splitter';

export function encodeMonomers(col: DG.Column): DG.Column | null {
  let encodeSymbol = MONOMER_ENCODE_MIN;
  const monomerSymbolDict: { [key: string]: number } = {};
  const uh = UnitsHandler.getOrCreate(col);
  const splitterFunc: SplitterFunc = uh.getSplitter();
  const encodedStringArray = [];
  for (let i = 0; i < col.length; ++i) {
    let encodedMonomerStr = '';
    const monomers = splitterFunc(col.get(i));
    for (const m of monomers) {
      if (!monomerSymbolDict[m]) {
        if (encodeSymbol > MONOMER_ENCODE_MAX) {
          grok.shell.error(`Not enough symbols to encode monomers`);
          return null;
        }
        monomerSymbolDict[m] = encodeSymbol;
        encodeSymbol++;
      }
      encodedMonomerStr += String.fromCodePoint(monomerSymbolDict[m]);
    }
    encodedStringArray.push(encodedMonomerStr);
  }
  return DG.Column.fromStrings('encodedMolecules', encodedStringArray);
}

export function getMolfilesFromSeq(col: DG.Column, monomersLibObject: any[]): any[][] | null {
  const uh = UnitsHandler.getOrCreate(col);
  const splitter: SplitterFunc = uh.getSplitter();
  const monomersDict = createMomomersMolDict(monomersLibObject);
  const molFiles = [];
  for (let i = 0; i < col.length; ++i) {
    const macroMolecule = col.get(i);
    const monomers = splitter(macroMolecule);
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
  const uh = UnitsHandler.getOrCreate(cell.column);
  const splitterFunc: SplitterFunc = uh.getSplitter();
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

export function isValidEnumeratorLib(json: any[]): boolean {
  return json.every((entry) => {
    return typeof entry === 'object' &&
      Object.values(helmFieldsToEnumeratorInputFields).every((field) => {
        return field in entry &&
          typeof entry[field] === 'string';
      });
  });
}

/** Specific to Enumerator for peptides */
export function getJsonMonomerLibForEnumerator(rawLib: any[]): any {
  // todo
  function validateMonomerLibWithSmiles() {
  }

  function restoreRgroupsInSmiles(rawSmiles: string) {
    const regex = new RegExp('\\[r\\]', 'g');
    let i = 0;
    return rawSmiles.replace(regex, (match) => { ++i; return `[${i}*]`; });
  }

  function prepareOutputSmilesColValue(smilesWithRestoredRgroups: string): string {
    const result = smilesWithRestoredRgroups.replace('[1*]', '[H:1]');
    return result.replace('[2*]', '[OH:2]');
  }

  function prepareMolblock(rawMolblock: string): string {
    return rawMolblock.replace('M  ISO', 'M  RGP');
  }

  const resultLib: any[] = [];
  validateMonomerLibWithSmiles();

  // for (let i = 0; i < df.rowCount; i++) {
  rawLib.forEach((monomer) => {
    Object.keys(jsonSdfMonomerLibDict).forEach((key) => {
      if (key === HELM_FIELDS.SYMBOL) {
        const monomerSymbol = monomer[helmFieldsToEnumeratorInputFields[key]];
        monomer[key] = monomerSymbol;
      } else if (key === HELM_FIELDS.SMILES) {
        const rawSmiles = monomer[helmFieldsToEnumeratorInputFields[key]];
        const smilesWithRestoredRgroups = restoreRgroupsInSmiles(rawSmiles);
        const smiles = prepareOutputSmilesColValue(smilesWithRestoredRgroups);
        monomer[key] = smiles;
      } else if (key === HELM_FIELDS.RGROUPS)
        monomer[key] = rGroupsDummy;
      else if (key === HELM_FIELDS.MOLFILE) {
        const rawSmiles = monomer[helmFieldsToEnumeratorInputFields[HELM_FIELDS.SMILES]];
        const smiles = restoreRgroupsInSmiles(rawSmiles);
        const rawMolfile = DG.chem.convert(smiles, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
        monomer[key] = prepareMolblock(rawMolfile);
      }
    });
    resultLib.push(monomer);
  });
  return resultLib;
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
          rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_SMILES] = altAtom === 'H' ? `[*:${radicalNum}][H]` : `O[*:${radicalNum}]`;
          rgroup[HELM_RGROUP_FIELDS.ALTERNATE_ID] = altAtom === 'H' ? `R${radicalNum}-H` : `R${radicalNum}-OH`;
          rgroup[HELM_RGROUP_FIELDS.CAP_GROUP_NAME] = altAtom === 'H' ? `H` : `OH`;
          rgroup[HELM_RGROUP_FIELDS.LABEL] = `R${radicalNum}`;
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

/** Calculates chemical similarity between reference sequence and list of sequences.
 * Similarity is computed as a sum of monomer similarities on corresponding positions. Monomer similarity is calculated
 * based on Morgan fingerprints.
 * @param {DG.Column<string>[]} positionColumns List of position columns containing monomers.
 * @param {string[]} referenceSequence Reference sequence.
 * @returns {Promise<DG.Column<number>>} Column with similarity values. */
export async function sequenceChemSimilarity(sequenceCol: DG.Column<string>, referenceSequence: ISeqSplitted): Promise<DG.Column<number>>;
export async function sequenceChemSimilarity(positionColumns: DG.Column<string>[], referenceSequence: ISeqSplitted): Promise<DG.Column<number>>;
export async function sequenceChemSimilarity(positionColumns: DG.Column<string>[] | DG.Column<string>,
  referenceSequence: ISeqSplitted): Promise<DG.Column<number>> {
  if (positionColumns instanceof DG.Column)
    positionColumns = splitAlignedSequences(positionColumns).columns.toList();

  const libHelper = await getMonomerLibHelper();
  const monomerLib = libHelper.getBioLib();
  // const smilesCols: DG.Column<string>[] = new Array(monomerCols.length);
  const rawCols: {categories: string[], data: Uint32Array, emptyIndex: number}[] = new Array(positionColumns.length);
  const rowCount = positionColumns[0].length;
  const totalSimilarity = new Float32Array(rowCount);

  // Calculate base similarity
  for (let position = 0; position < positionColumns.length; ++position) {
    const referenceMonomer = referenceSequence[position];
    const referenceMol = monomerLib.getMonomer('PEPTIDE', referenceMonomer)?.smiles ?? '';

    const monomerCol = positionColumns[position];
    const monomerColData = monomerCol.getRawData() as Uint32Array;
    const monomerColCategories = monomerCol.categories;
    const emptyCategoryIdx = monomerColCategories.indexOf('');
    rawCols[position] = {categories: monomerColCategories, data: monomerColData, emptyIndex: emptyCategoryIdx};
    if (typeof referenceMonomer === 'undefined')
      continue;
    
    // Calculating similarity for 
    const molCol = DG.Column.fromStrings('smiles',
      monomerColCategories.map((cat) => monomerLib.getMonomer('PEPTIDE', cat)?.smiles ?? ''));
    const _df = DG.DataFrame.fromColumns([molCol]); // getSimilarities expects that column is in dataframe
    const similarityCol = (await grok.chem.getSimilarities(molCol, referenceMol))!;
    const similarityColData = similarityCol.getRawData();

    for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
      const monomerCategoryIdx = monomerColData[rowIdx];
      totalSimilarity[rowIdx] += referenceMonomer !== '' && monomerCategoryIdx !== emptyCategoryIdx ?
        similarityColData[monomerCategoryIdx] :
        referenceMonomer === '' && monomerCategoryIdx === emptyCategoryIdx ? 1 : 0;
    }
  }

  for (let similarityIndex = 0; similarityIndex < totalSimilarity.length; ++similarityIndex) {
    let updatedSimilarity = totalSimilarity[similarityIndex] / referenceSequence.length;
    for (let position = 0; position < positionColumns.length; ++position) {
      const currentRawCol = rawCols[position];
      if ((position >= referenceSequence.length && currentRawCol.data[similarityIndex] !== currentRawCol.emptyIndex) ||
        (currentRawCol.data[similarityIndex] === currentRawCol.emptyIndex && position < referenceSequence.length)) {
        updatedSimilarity = DG.FLOAT_NULL;
        break;
      }
    }
    totalSimilarity[similarityIndex] = updatedSimilarity;
  }

  const similarityCol = DG.Column.fromFloat32Array('Similarity', totalSimilarity);
  return similarityCol;
}
