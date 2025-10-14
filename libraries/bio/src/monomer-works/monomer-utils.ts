// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {Observable} from 'rxjs';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

import {
  HELM_FIELDS, HELM_CORE_FIELDS, HELM_RGROUP_FIELDS, jsonSdfMonomerLibDict,
  MONOMER_ENCODE_MAX, MONOMER_ENCODE_MIN, SDF_MONOMER_NAME, HELM_REQUIRED_FIELD,
} from '../utils/const';
import {IMonomerLib, IMonomerSet} from '../types/index';
import {ISeqSplitted} from '../utils/macromolecule/types';
import {splitAlignedSequences} from '../utils/splitter';
import {GAP_SYMBOL} from '../utils/macromolecule/consts';
import {ISeqHelper} from '../utils/seq-helper';

export function encodeMonomers(col: DG.Column, seqHelper: ISeqHelper): DG.Column | null {
  let encodeSymbol = MONOMER_ENCODE_MIN;
  const monomerSymbolDict: { [key: string]: number } = {};
  const sh = seqHelper.getSeqHandler(col);
  const encodedStringArray = [];
  const rowCount = col.length;
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    let encodedMonomerStr = '';
    const monomers = sh.getSplitted(rowIdx);
    for (let posIdx = 0; posIdx < rowCount; ++posIdx) {
      const cm = monomers.getCanonical(posIdx);
      if (!monomerSymbolDict[cm]) {
        if (encodeSymbol > MONOMER_ENCODE_MAX) {
          grok.shell.error(`Not enough symbols to encode monomers`);
          return null;
        }
        monomerSymbolDict[cm] = encodeSymbol;
        encodeSymbol++;
      }
      encodedMonomerStr += String.fromCodePoint(monomerSymbolDict[cm]);
    }
    encodedStringArray.push(encodedMonomerStr);
  }
  return DG.Column.fromStrings('encodedMolecules', encodedStringArray);
}

export function getMolfilesFromSeq(
  col: DG.Column, monomersLibObject: any[], seqHelper: ISeqHelper
): any[][] | null {
  const sh = seqHelper.getSeqHandler(col);
  const monomersDict = createMomomersMolDict(monomersLibObject);
  const molFiles = [];
  const rowCount = col.length;
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    const monomers = sh.getSplitted(rowIdx);
    const molFilesForSeq = [];
    for (let j = 0; j < rowCount; ++j) {
      const cm = monomers.getCanonical(j);
      if (cm) {
        if (!monomersDict[cm]) {
          grok.shell.warning(`Monomer ${cm} is missing in HELM library. Structure cannot be created`);
          return null;
        }
        // what is the reason of double conversion?
        molFilesForSeq.push(JSON.parse(JSON.stringify(monomersDict[cm])));
      }
    }
    molFiles.push(molFilesForSeq);
  }
  return molFiles;
}

export function getMolfilesFromSingleSeq(
  cell: DG.Cell, monomersLibObject: any[], seqHelper: ISeqHelper
): any[][] | null {
  const sh = seqHelper.getSeqHandler(cell.column);
  const monomersDict = createMomomersMolDict(monomersLibObject);
  const molFiles = [];
  const monomers = sh.getSplitted(cell.rowIndex);
  const molFilesForSeq = [];
  for (let j = 0; j < monomers.length; ++j) {
    const cm = monomers.getCanonical(j);
    if (cm) {
      if (!monomersDict[cm]) {
        grok.shell.warning(`Monomer ${cm} is missing in HELM library. Structure cannot be created`);
        return null;
      }
      molFilesForSeq.push(JSON.parse(JSON.stringify(monomersDict[cm])));
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
    monomer[HELM_REQUIRED_FIELD.AUTHOR] ??= '';
    monomer[HELM_REQUIRED_FIELD.ID] ??= 0;
    monomer[HELM_REQUIRED_FIELD.SMILES] ??= '';
    monomer[HELM_REQUIRED_FIELD.MONOMER_TYPE] ??= 'Backbone';
    monomer[HELM_REQUIRED_FIELD.CREATE_DATE] ??= null;

    resultLib.push(monomer);
  }
  return resultLib;
}

export interface IMonomerLibFileEventManager {
  get addLibraryFileRequested$(): Observable<void>;
  get updateUIControlsRequested$(): Observable<string[]>;
  get librarySelectionRequested$(): Observable<[string, boolean]>;

  updateLibrarySelectionStatus(libFileName: string, isSelected: boolean): void;
}

export interface IMonomerLibFileManager {
  get eventManager(): IMonomerLibFileEventManager;

  getValidLibraryPaths(): string[];
  getValidLibraryPathsAsynchronously(): Promise<string[]>;

  addLibraryFile(fileContent: string, fileName: string): Promise<void>;
  deleteLibraryFile(fileName: string): Promise<void>;

  loadLibraryFromFile(path: string, fileName: string): Promise<IMonomerLib>;
}

export interface IMonomerLibHelper {
  get eventManager(): IMonomerLibFileEventManager;

  /** Ensures files are loaded and validated, throws error after timeout */
  awaitLoaded(timeout?: number): Promise<void>;

  /** Singleton monomer library collected from various sources */
  getMonomerLib(): IMonomerLib;

  /**  @deprecated Use {@link getMonomerLib} */
  getBioLib(): IMonomerLib;

  /** Singleton monomer set collected from various sources */
  getMonomerSets(): IMonomerSet;

  getFileManager(): Promise<IMonomerLibFileManager>;

  /** (Re)Loads libraries based on settings in user storage {@link LIB_STORAGE_NAME} to singleton.
   * @param {boolean} reload Clean {@link monomerLib} before load libraries [false]
   */
  loadMonomerLib(reload?: boolean): Promise<void>;

  /** @deprecated Use {@link loadMonomerLib} */
  loadLibraries(reload?: boolean): Promise<void>;

  /** (Re)loads monomer sets based on settings in user storage {@link SETS_STORAGE_NAME} to singleton.
   * @param {boolean} reload Clean {@link monomerSets} before load sets [false]
   */
  loadMonomerSets(reload?: boolean): Promise<void>;

  /** Reads library from file shares, handles .json and .sdf */
  readLibrary(path: string, fileName: string): Promise<IMonomerLib>;

  // -- Settings --

  /** Changes user lib settings. */
  loadMonomerLibForTests(): Promise<void>;
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

/* eslint-disable max-len */
export async function sequenceChemSimilarity(
  sequenceCol: DG.Column<string>, referenceSequence: ISeqSplitted, seqHelper: ISeqHelper): Promise<DG.Column<number>>;
export async function sequenceChemSimilarity(
  positionColumns: DG.Column<string>[], referenceSequence: ISeqSplitted, seqHelper: ISeqHelper): Promise<DG.Column<number>>;
export async function sequenceChemSimilarity(
  positionColumns: DG.Column<string>[] | DG.Column<string>, referenceSequence: ISeqSplitted, seqHelper: ISeqHelper
): Promise<DG.Column<number>> {
  /* eslint-enable max-len */
  if (positionColumns instanceof DG.Column)
    positionColumns = splitAlignedSequences(positionColumns, seqHelper).columns.toList();

  const libHelper = await getMonomerLibHelper();
  const monomerLib = libHelper.getMonomerLib();
  // const smilesCols: DG.Column<string>[] = new Array(monomerCols.length);
  const rawCols: { categories: string[], data: Uint32Array, emptyIndex: number }[] = new Array(positionColumns.length);
  const rowCount = positionColumns[0].length;
  const totalSimilarity = new Float32Array(rowCount);

  // Calculate base similarity
  for (let position = 0; position < positionColumns.length; ++position) {
    const referenceMonomerCanonical = position < referenceSequence.length ?
      referenceSequence.getCanonical(position) : GAP_SYMBOL;
    const referenceMol = monomerLib.getMonomer('PEPTIDE', referenceMonomerCanonical)?.smiles ?? '';

    const monomerCol = positionColumns[position];
    const monomerColData = monomerCol.getRawData() as Uint32Array;
    const monomerColCategories = monomerCol.categories;
    const emptyCategoryIdx = monomerColCategories.indexOf('');
    rawCols[position] = {categories: monomerColCategories, data: monomerColData, emptyIndex: emptyCategoryIdx};
    if (typeof referenceMonomerCanonical === 'undefined')
      continue;

    // Calculating similarity for
    const molCol = DG.Column.fromStrings('smiles',
      monomerColCategories.map((cat) => monomerLib.getMonomer('PEPTIDE', cat)?.smiles ?? ''));
    const _df = DG.DataFrame.fromColumns([molCol]); // getSimilarities expects that column is in dataframe
    const similarityCol: DG.Column<number> | null = (await grok.chem.getSimilarities(molCol, referenceMol))!;
    const similarityColData: Float32Array | null = similarityCol ? (similarityCol.getRawData() as Float32Array) : null;

    for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
      const monomerCategoryIdx = monomerColData[rowIdx];
      if (referenceMonomerCanonical !== GAP_SYMBOL && monomerCategoryIdx !== emptyCategoryIdx) {
        totalSimilarity[rowIdx] += similarityColData![monomerCategoryIdx];
      } else if (referenceMonomerCanonical === GAP_SYMBOL && monomerCategoryIdx === emptyCategoryIdx) {
        totalSimilarity[rowIdx] += 1;
      } // Do not increase similarity on mismatch score/penalty equals 0;
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

/* eslint-disable max-len */
/** Calculates chemical similarity between each pair of monomers.
 * @param {string[]} monomerSet Set of unique monomers.
 * @returns {Promise<{scoringMatrix: number[][], alphabetIndexes: {[monomerId: string]: number}}>} Ojbect containing similarity scoring matrix and monomer to index mapping. */
export async function calculateMonomerSimilarity(monomerSet: string[],
): Promise<{ scoringMatrix: number[][], alphabetIndexes: { [monomerId: string]: number } }> {
  /* eslint-enable max-len */
  const libHelper = await getMonomerLibHelper();
  const monomerLib = libHelper.getMonomerLib();
  const scoringMatrix: number[][] = [];
  const alphabetIndexes: { [id: string]: number } = {};
  const monomerMolecules = monomerSet.map((monomer) => monomerLib.getMonomer('PEPTIDE', monomer)?.smiles ?? '');
  const monomerMoleculesCol = DG.Column.fromStrings('smiles', monomerMolecules);

  for (let monomerIndex = 0; monomerIndex < monomerMolecules.length; ++monomerIndex) {
    const monomer = monomerSet[monomerIndex];
    alphabetIndexes[monomer] = monomerIndex;
    const monomerMol = monomerMolecules[monomerIndex];
    const similarityScores = monomerMol === '' ? new Array(monomerMolecules.length).fill(0) :
      (await grok.chem.getSimilarities(monomerMoleculesCol, monomerMol))!.getRawData();
    similarityScores[monomerIndex] = 1;

    scoringMatrix[monomerIndex] = Array.from(similarityScores);
  }

  return {scoringMatrix, alphabetIndexes};
}

export async function getMonomerSubstitutionMatrix(monomerSet: string[], fingerprintType: string = 'Morgan',
): Promise<{ scoringMatrix: number[][], alphabetIndexes: { [monomerId: string]: number } }> {
  const libHelper = await getMonomerLibHelper();
  const monomerLib = libHelper.getMonomerLib();
  const scoringMatrix: number[][] =
    new Array(monomerSet.length).fill(0).map(() => new Array(monomerSet.length).fill(0));
  const alphabetIndexes: { [id: string]: number } = {};
  // note, below specifically boolean OR is used to get either molfile or smiles, because we want to skip '' molfiles
  const monomerMolecules = monomerSet.map((monomer) => monomerLib.getMonomer('PEPTIDE', monomer)?.molfile || monomerLib.getMonomer('PEPTIDE', monomer)?.smiles || '');
  const fingerprintsFunc = DG.Func.find({package: 'Chem', name: 'getFingerprints'})[0];
  if (!fingerprintsFunc) {
    console.warn('Function "Chem:getFingerprints" is not found in chem package. falling back to Morgan fingerprints');
    return await calculateMonomerSimilarity(monomerSet);
  }
  const monomerMoleculesCol = DG.Column.fromStrings('smiles', monomerMolecules);
  // needed for function to work
  const _unusedMolDf = DG.DataFrame.fromColumns([monomerMoleculesCol]);
  const fingerPrints: (BitArray | null)[] =
    (await fingerprintsFunc.apply({col: monomerMoleculesCol, fingerprintType: fingerprintType}))?.entries!;
  if (!fingerPrints) {
    console.warn(`${fingerprintType} Fingerprints could not be calculated for monomers from chem package.
      falling back to Morgan fingerprints`);
    return await calculateMonomerSimilarity(monomerSet);
  }

  for (let i = 0; i < fingerPrints.length; ++i) {
    scoringMatrix[i][i] = 1;
    alphabetIndexes[monomerSet[i]] = i;
    if (!fingerPrints[i])
      continue;
    for (let j = i + 1; j < fingerPrints.length; ++j) {
      if (!fingerPrints[j])
        continue;
      // @ts-ignore
      scoringMatrix[i][j] = scoringMatrix[j][i] = tanimotoSimilarity(fingerPrints[i]!, fingerPrints[j]!);
    }
  }
  return {scoringMatrix, alphabetIndexes};
};
