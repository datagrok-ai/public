import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {getRdKitModule, getRdKitService} from './utils/chem-common-rdkit';
import {
  chemBeginCriticalSection,
  chemEndCriticalSection,
  defaultMorganFpLength,
  defaultMorganFpRadius,
  Fingerprint,
  rdKitFingerprintToBitArray,
} from './utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {assure} from '@datagrok-libraries/utils/src/test';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {getMolSafe} from './utils/mol-creation_rdkit';
import {_package} from './package';
import {IFpResult} from './rdkit-service/rdkit-service-worker-similarity';
import {SubstructureSearchType, getSearchProgressEventName, getSearchQueryAndType, getTerminateEventName} from './constants';
import {SubstructureSearchWithFpResult} from './rdkit-service/rdkit-service';
import { RDMol } from '@datagrok-libraries/chem-meta/src/rdkit-api';

const enum FING_COL_TAGS {
  molsCreatedForVersion = '.mols.created.for.version',
}

export const enum FILTER_TYPES {
  scaffold = 'scaffold',
  substructure = 'substructure'
}

let lastColumnInvalidated: string = '';
const canonicalSmilesColName = 'canonicalSmiles';

const currentSearchSmiles: {[key: string]: {[key: string]: string}} = {
  [FILTER_TYPES.scaffold]: {},
  [FILTER_TYPES.substructure]: {}
};

function _chemFindSimilar(molStringsColumn: DG.Column, fingerprints: (BitArray | null)[],
  queryMolString: string, settings: { [name: string]: any }): DG.DataFrame {
  const len = molStringsColumn.length;
  const distances = _chemGetSimilarities(queryMolString, fingerprints);
  const limit = Math.min((settings.hasOwnProperty('limit') ? settings.limit : len), len);
  const minScore = settings.hasOwnProperty('minScore') ? settings.minScore : 0.0;
  const sortedIndices = Array.from(Array(len).keys()).sort((i1, i2) => {
    const a1 = distances[i1];
    const a2 = distances[i2];
    if (a2 < a1) return -1;
    if (a2 > a1) return +1;
    return 0;
  });
  const sortedMolStringsArr: string[] = [];
  const sortedMolIndArr: number[] = [];
  const sortedScoresArr: number[] = [];
  for (let n = 0; n < limit; n++) {
    const idx = sortedIndices[n];
    const score = distances[idx];
    if (score < minScore)
      break;

    sortedMolStringsArr.push(molStringsColumn.get(idx));
    sortedScoresArr.push(score);
    sortedMolIndArr.push(idx);
  }
  const length = sortedMolStringsArr.length;
  const sortedMolStrings = DG.Column.fromType(DG.TYPE.STRING, 'molecule', length).init((i) => sortedMolStringsArr[i]);
  const sortedMolInd = DG.Column.fromType(DG.TYPE.INT, 'index', length).init((i) => sortedMolIndArr[i]);
  sortedMolStrings.semType = DG.SEMTYPE.MOLECULE;
  const sortedScores = DG.Column.fromType(DG.TYPE.FLOAT, 'score', length).init((i) => sortedScoresArr[i]); ;
  return DG.DataFrame.fromColumns([sortedMolStrings, sortedScores, sortedMolInd]);
}

// Only this function receives {sorted} in settings
export function _chemGetSimilarities(queryMolString: string, fingerprints: (BitArray | null)[]): number[] {
  const distances = new Array(fingerprints.length).fill(0.0);
  try {
    const sample = chemGetFingerprint(queryMolString, Fingerprint.Morgan);
    for (let i = 0; i < fingerprints.length; ++i)
      distances[i] = tanimotoSimilarity(fingerprints[i] ?? BitArray.fromBytes(new Uint8Array()), sample);
  } catch {
    return distances;
  }
  return distances;
}

function _chemGetDiversities(limit: number, molStringsColumn: DG.Column, fingerprints: (BitArray | null)[]): string[] {
  limit = Math.min(limit, fingerprints.length);
  const indexes = ArrayUtils.indexesOf(fingerprints, (f) => f != null);
  const diverseIndexes = getDiverseSubset(indexes.length, limit,
    (i1: number, i2: number) => 1 - tanimotoSimilarity(fingerprints[indexes[i1]]!, fingerprints[indexes[i2]]!));

  const diversities = new Array(limit).fill('');

  for (let i = 0; i < limit; i++)
    diversities[i] = molStringsColumn.get(indexes[diverseIndexes[i]]);

  return diversities;
}


function invalidatedColumnKey(col: DG.Column): string {
  return col.dataFrame?.name + '.' + col.name;
}

function colInvalidated(col: DG.Column): Boolean {
  return (lastColumnInvalidated == invalidatedColumnKey(col) &&
    col.getTag(FING_COL_TAGS.molsCreatedForVersion) == String(col.version));
}

//updates array of RDMols in case column has been changed
async function _invalidate(molCol: DG.Column) {
  if (!colInvalidated(molCol)) {
    await (await getRdKitService()).initMoleculesStructures(molCol.toList());
    molCol.setTag(FING_COL_TAGS.molsCreatedForVersion, String(molCol.version + 1));
    lastColumnInvalidated = invalidatedColumnKey(molCol);
  }
}


function checkForSavedColumn(col: DG.Column, colTag: Fingerprint | string): DG.Column<any> | null {
  const savedColName = '~' + col.name + '.' + colTag;
  const colVerTag = '.' + colTag + '.Version';

  if (col.dataFrame && col.dataFrame?.col(savedColName) &&
      col.getTag(colVerTag) == String(col.version)) {
    const savedCol = col.dataFrame.columns.byName(savedColName);
    return savedCol;
  }
  return null;
}

function saveColumns(col: DG.Column, data: (Uint8Array | string | null)[][],
  tags: string[], colTypes: DG.COLUMN_TYPE[]): void {
  if (col.dataFrame) {
    for (let i = 0; i < data.length; i++)
      saveColumn(col, data[i], tags[i], colTypes[i]);

    const colVersion = col.version + data.length;
    for (let i = 0; i < data.length; i++)
      col.setTag('.' + tags[i] + '.Version', String(colVersion));
  }
}


function saveColumn(col: DG.Column, data: (Uint8Array | string | null)[],
  tagName: string, colType: DG.COLUMN_TYPE): void {
  if (!col.dataFrame)
    throw new Error('Column has no parent dataframe');

  const savedColumnName = '~' + col.name + '.' + tagName;
  const newCol = col.dataFrame.columns.getOrCreate(savedColumnName, colType, data.length);

  newCol.init((i) => data[i]);
  newCol.meta.includeInCsvExport = false;
  newCol.meta.includeInBinaryExport = false;
}

/**
* Creates columns with fingerprints and canonical smiles (if required) if haven't been created previously or
*  molecular column has been changed
* @async
* @param {DG.Column} molCol - column to count fps for
* @param {Fingerprint} fingerprintsType - Morgan or Pattern
* @param {boolean} returnSmiles - optional parameter, if passed column with canonical smiles is returned
additionally to fps
* */
async function invalidateAndSaveColumns(molCol: DG.Column, fingerprintsType: Fingerprint,
  returnSmiles?: boolean): Promise<IFpResult> {
  const fpRes = await (await getRdKitService()).getFingerprints(fingerprintsType, molCol.toList(), returnSmiles);
  returnSmiles ?
    saveColumns(molCol, [fpRes.fps, fpRes.smiles!], [fingerprintsType, canonicalSmilesColName],
      [DG.COLUMN_TYPE.BYTE_ARRAY, DG.COLUMN_TYPE.STRING]):
    saveColumns(molCol, [fpRes.fps], [fingerprintsType], [DG.COLUMN_TYPE.BYTE_ARRAY]);
  return fpRes;
}

async function getUint8ArrayFingerprints(
  molCol: DG.Column, fingerprintsType: Fingerprint = Fingerprint.Morgan,
  returnSmiles = false): Promise<IFpResult> {
  await chemBeginCriticalSection();
  try {
    const fgsCheck = checkForSavedColumn(molCol, fingerprintsType);
    if (returnSmiles) {
      const smilesCheck = checkForSavedColumn(molCol, canonicalSmilesColName);
      if (fgsCheck && smilesCheck)
        return {fps: fgsCheck.toList(), smiles: smilesCheck.toList()};
      else {
        const fpResult = await invalidateAndSaveColumns(molCol, fingerprintsType, returnSmiles);
        return {fps: fpResult.fps, smiles: fpResult.smiles};
      }
    } else {
      if (fgsCheck)
        return {fps: fgsCheck.toList(), smiles: null};
      else {
        const fpResult = await invalidateAndSaveColumns(molCol, fingerprintsType);
        return {fps: fpResult.fps, smiles: null};
      }
    }
  } finally {
    chemEndCriticalSection();
  }
}

export async function chemGetFingerprints(...args: [DG.Column, Fingerprint?, boolean?]):
  Promise<(BitArray | null)[]> {
  return (await getUint8ArrayFingerprints(...args)).fps.map((el) => el ? rdKitFingerprintToBitArray(el) : null);
}

export async function chemGetSimilarities(molStringsColumn: DG.Column, queryMolString = '')
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, false)!;

  return queryMolString.length != 0 ?
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      _chemGetSimilarities(queryMolString, fingerprints)) : null;
}

export async function chemGetDiversities(molStringsColumn: DG.Column, limit: number)
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(limit, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, false)!;

  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'molecule',
    _chemGetDiversities(limit, molStringsColumn, fingerprints));
}

export async function chemFindSimilar(molStringsColumn: DG.Column, queryMolString = '',
  settings: { [name: string]: any } = {}) : Promise<DG.DataFrame | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, false)!;
  return queryMolString.length != 0 ?
    _chemFindSimilar(molStringsColumn, fingerprints, queryMolString, settings) : null;
}

/**
* Performes substructure search in the given moelcular column by a given substructure
* @async
* @param {DG.Column} molStringsColumn - column search in
* @param {string} molString - smiles/molblock to filter by
* @param {string} molBlockFailover - smart to filter by (is used if creation of RDMol object from
  query parameter failed)
* @param {boolean} columnIsCanonicalSmiles - if column is canonical smiles itself, than invisible
canonical smiles column
will not be created along with pattern fp column
* @param {boolean} awaitAll - in case of true fucntion will wait for results on a whole table to be received
before returning (required for compatibility)
* */
export async function chemSubstructureSearchLibrary(
  molStringsColumn: DG.Column, molString: string, molBlockFailover: string, filterType = FILTER_TYPES.substructure,
  columnIsCanonicalSmiles = false, awaitAll = true, searchType = SubstructureSearchType.CONTAINS, similarityCutOff = 0.8,
  fp = Fingerprint.Morgan): Promise<BitArray> {
  const searchKey = `${molStringsColumn?.dataFrame?.name ?? ''}-${molStringsColumn?.name ?? ''}`;
  const currentSearch = `${molBlockFailover}_${searchType}_${similarityCutOff}_${fp}`;
  currentSearchSmiles[filterType][searchKey] = currentSearch;
  _package.logger.debug(`in chemSubstructureSearchLibrary, filterType: ${filterType}, searchkey: ${searchKey}, currentSearch: ${currentSearch}`);
  await chemBeginCriticalSection();
  _package.logger.debug(`in chemSubstructureSearchLibrary, began critical section currentSearch: ${currentSearch}`);
  const terminateEventName = getTerminateEventName(molStringsColumn.dataFrame?.name ?? '', molStringsColumn.name);
  if (currentSearchSmiles[filterType][searchKey] !== currentSearch && filterType !== FILTER_TYPES.scaffold) {
    _package.logger.debug(`in chemSubstructureSearchLibrary, ending critical section without search: ${currentSearch}`);
    grok.events.fireCustomEvent(terminateEventName, getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff));
    chemEndCriticalSection();
    _package.logger.debug(`in chemSubstructureSearchLibrary, ended critical section: ${currentSearch}`);
    return new BitArray(molStringsColumn.length);
  }

  try {
    _package.logger.debug(`in chemSubstructureSearchLibrary, search start: ${currentSearch}`);
    let invalidateCacheFlag = false;
    const rdKitService = await getRdKitService();
    await rdKitService.setTerminateFlag(false);
    const currentCol = invalidatedColumnKey(molStringsColumn);
    if (lastColumnInvalidated !== currentCol) {
      invalidateCacheFlag = true;
      lastColumnInvalidated = currentCol;
    }

    const matchesBitArray = new BitArray(molStringsColumn.length);

    const searchProgressEventName =
      getSearchProgressEventName(molStringsColumn.dataFrame?.name ?? '', molStringsColumn.name);
    const updateFilterFunc = (progress: number) => {
      grok.events.fireCustomEvent(searchProgressEventName, progress * 100);
    };

    const getSavedCol = (colName: string) => {
      const savedCol = checkForSavedColumn(molStringsColumn, colName);
      if (!savedCol)
        invalidateCacheFlag = true;
      return savedCol;
    };

    const canonicalSmilesList = getSavedCol(canonicalSmilesColName)?.toList() ??
      new Array<string | null>(molStringsColumn.length).fill(null);
    const fpType = searchType === SubstructureSearchType.IS_SIMILAR ? fp : Fingerprint.Pattern;
    const fgsList = getSavedCol(fpType)?.toList() ??
      new Array<Uint8Array | null>(molStringsColumn.length).fill(null);

    if (invalidateCacheFlag) //invalidating cache in case column has been changed
      await rdKitService.invalidateCache();

    const result: SubstructureSearchWithFpResult = {
      bitArray: matchesBitArray,
      fpsRes: {
        fps: fgsList,
        smiles: !columnIsCanonicalSmiles ? canonicalSmilesList : null,
      },
    };
    const subFuncs = await rdKitService.
      searchSubstructureWithFps(molString, molBlockFailover, result, updateFilterFunc,
        molStringsColumn.toList(), !columnIsCanonicalSmiles, searchType, similarityCutOff, fp);     
    const saveProcessedColumns = () => {
      try {
        !columnIsCanonicalSmiles ?
          saveColumns(molStringsColumn, [result.fpsRes!.fps, result.fpsRes!.smiles!],
            [fpType, canonicalSmilesColName], [DG.COLUMN_TYPE.BYTE_ARRAY, DG.COLUMN_TYPE.STRING]):
          saveColumns(molStringsColumn, [result.fpsRes!.fps], [fpType], [DG.COLUMN_TYPE.BYTE_ARRAY]);
          _package.logger.debug(`in chemSubstructureSearchLibrary, saveProcessedColumns: ${currentSearch}`);
      } catch {

      } finally {
        _package.logger.debug(`in chemSubstructureSearchLibrary, ending critical section: ${currentSearch}`);
        chemEndCriticalSection();
        _package.logger.debug(`in chemSubstructureSearchLibrary, ended critical section: ${currentSearch}`);
      }
    };
    const fireFinishEvents = () => {
      _package.logger.debug(`in chemSubstructureSearchLibrary, fireFinishEvents: ${getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff)}`);
      grok.events.fireCustomEvent(searchProgressEventName, 100);
      grok.events.fireCustomEvent(terminateEventName, getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff));
      saveProcessedColumns();
    };
    if (awaitAll) {
      await Promise.all(subFuncs.promises);
      fireFinishEvents();
    } else {
      const sub = grok.events.onCustomEvent(terminateEventName).subscribe(async (molAndSearchType: string) => {
        _package.logger.debug(`in chemSubstructureSearchLibrary, terminate event handler, ${molAndSearchType} ****** ${getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff)}`);
        if (molAndSearchType === getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff)) {
          await rdKitService.setTerminateFlag(true);
          subFuncs!.setTerminateFlag();
          await Promise.allSettled(subFuncs.promises);
          await rdKitService.setTerminateFlag(false);
          saveProcessedColumns();
          sub.unsubscribe();
        }
      });

      subFuncs?.promises && (Promise.allSettled(subFuncs?.promises).then(() => {
        _package.logger.debug(`in chemSubstructureSearchLibrary, subFuncs all settled, ${currentSearch}`);
        if (!subFuncs!.getTerminateFlag()) {
          sub.unsubscribe();
          _package.logger.debug(`in chemSubstructureSearchLibrary, subFuncs all settled firing finish events, ${currentSearch}`);
          fireFinishEvents();
        }
      }));
    }
    return result.bitArray;
  } catch (e: any) {
    grok.shell.error(e.message);
    _package.logger.debug(`in chemSubstructureSearchLibrary, ending chemEndCriticalSection in catch, ${currentSearch}`);
    chemEndCriticalSection();
    _package.logger.debug(`in chemSubstructureSearchLibrary, ended chemEndCriticalSection in catch, ${currentSearch}`);
    throw e;
  }
}

export function getRDKitFpAsUint8Array(mol: RDMol, fingerprint: Fingerprint): Uint8Array {
  if (fingerprint == Fingerprint.Morgan) {
    return mol.get_morgan_fp_as_uint8array(JSON.stringify({
      radius: defaultMorganFpRadius,
      nBits: defaultMorganFpLength,
    }));
  } else if (fingerprint == Fingerprint.Pattern)
    return mol.get_pattern_fp_as_uint8array();
  else if (fingerprint == Fingerprint.AtomPair)
    return mol.get_atom_pair_fp_as_uint8array();
  else if (fingerprint == Fingerprint.MACCS)
    return mol.get_maccs_fp_as_uint8array();
  else if (fingerprint == Fingerprint.RDKit)
    return mol.get_rdkit_fp_as_uint8array();
  else if (fingerprint == Fingerprint.TopologicalTorsion)
    return mol.get_topological_torsion_fp_as_uint8array();
  else
    throw new Error(`${fingerprint} does not match any fingerprint`);
}

export function chemGetFingerprint(molString: string, fingerprint: Fingerprint, onError?:
  (error: Error) => any, convertToBitArray = true): BitArray {
  let mol = null;
  try {
    mol = getMolSafe(molString, {}, getRdKitModule()).mol;
    if (mol) {
      let fp = getRDKitFpAsUint8Array(mol, fingerprint);
      return rdKitFingerprintToBitArray(fp) as BitArray;
    } else
      throw new Error(`Chem | Possibly a malformed molString: ${molString}`);
  } catch (e: any) {
    if (onError)
      return onError(e);
    else
      throw new Error(`Chem | Possibly a malformed molString: ${molString}`);
  } finally {
    mol?.delete();
  }
}

