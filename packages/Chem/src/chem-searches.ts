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
import {SubstructureSearchType, getSearchProgressEventName, getSearchQueryAndType, getTerminateEventName,
  getValuesChangedEventName} from './constants';
import {SubstructureSearchWithFpResult} from './rdkit-service/rdkit-service';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {filter} from 'rxjs/operators';


const enum FING_COL_TAGS {
  molsCreatedForVersion = '.mols.created.for.version',
}

const DATA_CHANGED_TAG = 'dataChanged';
const DATA_CHANGED_SUBSCRIBED_TAG = 'dataChangedSubscribed';

export const enum FILTER_TYPES {
  scaffold = 'scaffold',
  substructure = 'substructure'
}

let lastColumnInvalidated: string = '';
const canonicalSmilesColName = 'canonicalSmiles';

const currentSearchSmiles: {[key: string]: {[key: string]: string}} = {
  [FILTER_TYPES.scaffold]: {},
  [FILTER_TYPES.substructure]: {},
};

function _chemFindSimilar(molStringsColumn: DG.Column, fingerprints: (BitArray | null)[],
  queryMolString: string, settings: { [name: string]: any }): DG.DataFrame {
  const len = molStringsColumn.length;
  const distances = _chemGetSimilarities(queryMolString, fingerprints);
  const limit = Math.min(settings.limit ?? len, len);
  const minScore = settings.minScore ?? 0.0;
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
  const sortedScores = DG.Column.fromType(DG.TYPE.FLOAT, 'score', length).init((i) => sortedScoresArr[i]);
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

  // Pre-allocated to the known length; every slot is written below so no need to .fill('').
  const diversities = new Array<string>(limit);

  for (let i = 0; i < limit; i++)
    diversities[i] = molStringsColumn.get(indexes[diverseIndexes[i]]);

  return diversities;
}


function invalidatedColumnKey(col: DG.Column): string {
  return col.dataFrame?.name + '.' + col.name;
}

/**
 * Subscribes (once per column) to value changes so that the set of edited row indexes is tracked
 * in col.temp[DATA_CHANGED_TAG] and a per-column event is fired for interested consumers (e.g. the
 * substructure filter).
 */
export function subscribeToColumnChanges(col: DG.Column): void {
  if (col.temp[DATA_CHANGED_SUBSCRIBED_TAG])
    return;
  col.temp[DATA_CHANGED_SUBSCRIBED_TAG] = true;
  col.temp[DATA_CHANGED_TAG] = null; // null = nothing changed
  const valuesChangedEventName = getValuesChangedEventName(col.dataFrame?.name ?? '', col.name);
  col.dataFrame?.onValuesChanged.pipe(filter((args: any) =>
    args?.args?.column?.name === col.name && args?.args?.indexes?.length)).subscribe((args: any) => {
    // accumulate edited row indexes into a sparse Set (cheap for large columns)
    const changed: Set<number> = col.temp[DATA_CHANGED_TAG] ?? new Set<number>();
    const indexes = args.args.indexes;
    for (let i = 0; i < indexes.length; i++) {
      if (indexes[i] < col.length)
        changed.add(indexes[i]);
    }
    col.temp[DATA_CHANGED_TAG] = changed;
    grok.events.fireCustomEvent(valuesChangedEventName, {indexes: Array.from(indexes as ArrayLike<number>)});
  });
}

interface ISavedColumnsResult {
  cols: {[key: string | Fingerprint]: DG.Column<any>};
  changedRows: number[] | null; // null = nothing changed; otherwise the indexes of edited rows
}

/**
 * Returns the saved (cached) columns that exist for the given tags, together with the indexes of the
 * rows edited since the last call (or null if none). The change set is consumed (reset) on return, so
 * the caller is responsible for recomputing the affected rows. Subscribes to the column on first use.
 */
function checkForSavedColumns(col: DG.Column, colTags: (Fingerprint | string) []): ISavedColumnsResult {
  subscribeToColumnChanges(col);

  const changed: Set<number> | null = col.temp[DATA_CHANGED_TAG] ?? null;
  col.temp[DATA_CHANGED_TAG] = null;

  const cols: {[key: string | Fingerprint]: DG.Column<any>} = {};
  for (const colTag of colTags) {
    const savedColName = '~' + col.name + '.' + colTag;
    if (col.dataFrame?.col(savedColName))
      cols[colTag] = col.dataFrame.columns.byName(savedColName);
  }

  // keep only in-range indexes; treat an empty result as "nothing changed" (null)
  const inRange = changed ? Array.from(changed).filter((i) => i < col.length) : [];
  return {cols, changedRows: inRange.length > 0 ? inRange : null};
}

function saveColumns(col: DG.Column, data: (Uint8Array | string | null)[][],
  tags: string[], colTypes: DG.COLUMN_TYPE[]): void {
  if (col.dataFrame) {
    for (let i = 0; i < data.length; i++)
      saveColumn(col, data[i], tags[i], colTypes[i]);
  }
}


function saveColumn(col: DG.Column, data: (Uint8Array | string | null)[],
  tagName: string, colType: DG.COLUMN_TYPE): void {
  if (!col.dataFrame)
    throw new Error('Column has no parent dataframe');

  const savedColumnName = '~' + col.name + '.' + tagName;
  const newCol = col.dataFrame.columns.getOrCreate(savedColumnName, colType);

  newCol.init((i) => data[i]);
  newCol.meta.includeInCsvExport = false;
  newCol.meta.includeInBinaryExport = false;
}

/**
* (Re)computes and saves fingerprint (and optionally canonical smiles) columns. When saved columns exist they
* are reused, recomputing only the given changed rows; when there are no saved columns to reuse the whole
* column is computed.
* @async
* @param {DG.Column} molCol - column to count fps for
* @param {Fingerprint} fingerprintsType - Morgan or Pattern
* @param {boolean} returnSmiles - if true, canonical smiles are computed/saved alongside fps
* @param {DG.Column | null} savedFpCol - existing saved fp column to reuse (null -> nothing to reuse -> full recompute)
* @param {DG.Column | null} savedSmilesCol - existing saved canonical smiles column to reuse
* @param {number[] | null} changedRows - rows to recompute when reusing saved columns (null -> nothing changed)
* */
async function invalidateAndSaveColumns(molCol: DG.Column, fingerprintsType: Fingerprint, returnSmiles: boolean,
  savedFpCol: DG.Column | null, savedSmilesCol: DG.Column | null, changedRows: number[] | null): Promise<IFpResult> {
  const rdKitService = await getRdKitService();
  const canReuse = !!savedFpCol && (!returnSmiles || !!savedSmilesCol);

  let fps: (Uint8Array | null)[];
  let smiles: (string | null)[] | null;
  if (!canReuse) {
    // no saved columns to reuse -> compute the whole column
    const fpRes = await rdKitService.getFingerprints(fingerprintsType, molCol.toList(), returnSmiles);
    fps = fpRes.fps;
    smiles = fpRes.smiles;
  } else {
    // reuse saved fps/canonical smiles, recompute only the changed rows (none -> pure reuse)
    fps = savedFpCol!.toList();
    smiles = returnSmiles ? savedSmilesCol!.toList() : null;
    if (changedRows) {
      const changedMols = changedRows.map((i) => molCol.get(i));
      const recomputed = await rdKitService.getFingerprints(fingerprintsType, changedMols, returnSmiles);
      for (let k = 0; k < changedRows.length; k++) {
        fps[changedRows[k]] = recomputed.fps[k];
        if (returnSmiles && recomputed.smiles)
          smiles![changedRows[k]] = recomputed.smiles[k];
      }
    }
  }

  returnSmiles ?
    saveColumns(molCol, [fps, smiles!], [fingerprintsType, canonicalSmilesColName],
      [DG.COLUMN_TYPE.BYTE_ARRAY, DG.COLUMN_TYPE.STRING]):
    saveColumns(molCol, [fps], [fingerprintsType], [DG.COLUMN_TYPE.BYTE_ARRAY]);
  return {fps, smiles};
}

async function getUint8ArrayFingerprints(
  molCol: DG.Column, fingerprintsType: Fingerprint = Fingerprint.Morgan,
  returnSmiles = false): Promise<IFpResult> {
  await chemBeginCriticalSection();
  try {
    const {cols, changedRows} = checkForSavedColumns(molCol, [fingerprintsType, canonicalSmilesColName]);
    const noChanges = changedRows === null;
    const fgsCheck = cols[fingerprintsType] ?? null;
    if (returnSmiles) {
      const smilesCheck = cols[canonicalSmilesColName] ?? null;
      if (fgsCheck && smilesCheck && noChanges)
        return {fps: fgsCheck.toList(), smiles: smilesCheck.toList()};
      else {
        const fpResult = await invalidateAndSaveColumns(
          molCol, fingerprintsType, true, fgsCheck, smilesCheck, changedRows);
        return {fps: fpResult.fps, smiles: fpResult.smiles};
      }
    } else {
      if (fgsCheck && noChanges)
        return {fps: fgsCheck.toList(), smiles: null};
      else {
        const fpResult = await invalidateAndSaveColumns(molCol, fingerprintsType, false, fgsCheck, null, changedRows);
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

  return queryMolString.length !== 0 ?
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      _chemGetSimilarities(queryMolString, fingerprints)) : null;
}

export async function chemGetDiversities(molStringsColumn: DG.Column, limit: number)
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(limit, 'limit');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, false)!;

  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'molecule',
    _chemGetDiversities(limit, molStringsColumn, fingerprints));
}

export async function chemFindSimilar(molStringsColumn: DG.Column, queryMolString = '',
  settings: { [name: string]: any } = {}) : Promise<DG.DataFrame | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn, Fingerprint.Morgan, false)!;
  return queryMolString.length !== 0 ?
    _chemFindSimilar(molStringsColumn, fingerprints, queryMolString, settings) : null;
}

/**
* Performs substructure search in the given molecular column by a given substructure
* @async
* @param {DG.Column} molStringsColumn - column search in
* @param {string} molString - smiles/molblock to filter by
* @param {string} molBlockFailover - smart to filter by (is used if creation of RDMol object from
  query parameter failed)
* @param {boolean} columnIsCanonicalSmiles - if column is canonical smiles itself, than invisible
canonical smiles column
will not be created along with pattern fp column
* @param {boolean} awaitAll - in case of true function will wait for results on a whole table to be received
before returning (required for compatibility)
* */
export async function chemSubstructureSearchLibrary(
  molStringsColumn: DG.Column, molString: string, molBlockFailover: string, filterType = FILTER_TYPES.substructure,
  columnIsCanonicalSmiles = false, awaitAll = true, searchType = SubstructureSearchType.CONTAINS, similarityCutOff = 0.8,
  fp = Fingerprint.Morgan, includeMask: BitArray | null = null): Promise<BitArray> {
  const searchKey = `${molStringsColumn?.dataFrame?.name ?? ''}-${molStringsColumn?.name ?? ''}`;
  const currentSearch = `${molBlockFailover}_${searchType}_${similarityCutOff}_${fp}`;
  currentSearchSmiles[filterType][searchKey] = currentSearch;
  _package.logger.debug(`in chemSubstructureSearchLibrary, filterType: ${filterType}, searchkey: ${searchKey}, currentSearch: ${currentSearch}`);
  await chemBeginCriticalSection();
  _package.logger.debug(`in chemSubstructureSearchLibrary, began critical section currentSearch: ${currentSearch}`);
  const terminateEventName = getTerminateEventName(molStringsColumn.dataFrame?.name ?? '', molStringsColumn.name);
  if (currentSearchSmiles[filterType][searchKey] !== currentSearch && filterType !== FILTER_TYPES.scaffold) {
    _package.logger.debug(`in chemSubstructureSearchLibrary, ending critical section without search: ${currentSearch}`);
    grok.events.fireCustomEvent(terminateEventName, getSearchQueryAndType(molBlockFailover, searchType, fp,
      similarityCutOff));
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

    const fpType = searchType === SubstructureSearchType.IS_SIMILAR ? fp : Fingerprint.Pattern;
    const {cols: savedCols, changedRows} = checkForSavedColumns(molStringsColumn, [canonicalSmilesColName, fpType]);
    const savedFpCol = savedCols[fpType] ?? null;
    const savedSmilesCol = savedCols[canonicalSmilesColName] ?? null;
    const hasSaved = !!savedFpCol && (columnIsCanonicalSmiles || !!savedSmilesCol);

    // in case there are no saved columns -> recompute all
    const recomputeAll = !hasSaved;
    const rowsToRecompute: number[] = [];
    if (!recomputeAll) {
      if (includeMask) {
        for (let i = includeMask.findNext(-1); i !== -1; i = includeMask.findNext(i)) {
          if (i < molStringsColumn.length)
            rowsToRecompute.push(i);
        }
      } else if (changedRows) {
        for (const i of changedRows)
          rowsToRecompute.push(i);
      }
    }

    let canonicalSmilesList: (string | null)[];
    let fgsList: (Uint8Array | null)[];
    if (recomputeAll) {
      // no saved columns: pass empty lists; searchSubstructureWithFps computes fps lazily during the search
      canonicalSmilesList = new Array<string | null>(molStringsColumn.length).fill(null);
      fgsList = new Array<Uint8Array | null>(molStringsColumn.length).fill(null);
    } else if (rowsToRecompute.length > 0) {
      // reuse saved fps/canonical smiles, recompute only the changed rows (persisted by invalidateAndSaveColumns)
      const res = await invalidateAndSaveColumns(molStringsColumn, fpType, !columnIsCanonicalSmiles,
        savedFpCol, savedSmilesCol, rowsToRecompute);
      fgsList = res.fps;
      canonicalSmilesList = res.smiles ?? new Array<string | null>(molStringsColumn.length).fill(null);
    } else {
      // nothing changed: reuse saved fps/canonical smiles as-is
      fgsList = savedFpCol!.toList();
      canonicalSmilesList = savedSmilesCol?.toList() ??
        new Array<string | null>(molStringsColumn.length).fill(null);
    }

    if (recomputeAll || invalidateCacheFlag) //invalidating whole cache on column switch or full recompute
      await rdKitService.invalidateCache();

    const result: SubstructureSearchWithFpResult = {
      bitArray: matchesBitArray,
      fpsRes: {
        fps: fgsList,
        smiles: !columnIsCanonicalSmiles ? canonicalSmilesList : null,
      },
    };
    let numOfCalculatedFpBatches = 0;
    const updateNumOfCalculatedFpBatches = () => {
      numOfCalculatedFpBatches++;
    };
    const subFuncs = await rdKitService.
      searchSubstructureWithFps(molString, molBlockFailover, result, updateFilterFunc,
        molStringsColumn.toList(), !columnIsCanonicalSmiles, searchType, similarityCutOff, fp,
        updateNumOfCalculatedFpBatches, includeMask);
    const saveProcessedColumns = () => {
      try {
        //save procecced columns only in case at least one fp batch has been calculated.
        // (rows recomputed for individual edits are already persisted by invalidateAndSaveColumns)
        if (numOfCalculatedFpBatches) {
          !columnIsCanonicalSmiles ?
            saveColumns(molStringsColumn, [result.fpsRes!.fps, result.fpsRes!.smiles!],
              [fpType, canonicalSmilesColName], [DG.COLUMN_TYPE.BYTE_ARRAY, DG.COLUMN_TYPE.STRING]):
            saveColumns(molStringsColumn, [result.fpsRes!.fps], [fpType], [DG.COLUMN_TYPE.BYTE_ARRAY]);
          _package.logger.debug(`in chemSubstructureSearchLibrary, saveProcessedColumns: ${currentSearch}`);
        }
      } catch (e: any) {
        _package.logger.debug(e);
      } finally {
        _package.logger.debug(`in chemSubstructureSearchLibrary, ending critical section: ${currentSearch}`);
        chemEndCriticalSection();
        _package.logger.debug(`in chemSubstructureSearchLibrary, ended critical section: ${currentSearch}`);
      }
    };
    const fireFinishEvents = () => {
      _package.logger.debug(`in chemSubstructureSearchLibrary, fireFinishEvents: ${getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff)}`);
      grok.events.fireCustomEvent(searchProgressEventName, 100);
      grok.events.fireCustomEvent(terminateEventName, getSearchQueryAndType(molBlockFailover, searchType, fp,
        similarityCutOff));
      saveProcessedColumns();
    };
    if (awaitAll) {
      await Promise.all(subFuncs.promises);
      fireFinishEvents();
    } else {
      const sub = grok.events.onCustomEvent(terminateEventName).subscribe(async (molAndSearchType: string) => {
        _package.logger.debug(`in chemSubstructureSearchLibrary, terminate event handler, ${molAndSearchType} ****** ${getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff)}`);
        if (molAndSearchType ===
            getSearchQueryAndType(molBlockFailover, searchType, fp, similarityCutOff)) {
          await rdKitService.setTerminateFlag(true);
          subFuncs!.setTerminateFlag();
          await Promise.allSettled(subFuncs.promises);
          await rdKitService.setTerminateFlag(false);
          saveProcessedColumns();
          sub.unsubscribe();
        }
      });

      if (subFuncs?.promises) {
        Promise.allSettled(subFuncs.promises).then(() => {
          _package.logger.debug(`in chemSubstructureSearchLibrary, subFuncs all settled, ${currentSearch}`);
          if (!subFuncs!.getTerminateFlag()) {
            sub.unsubscribe();
            _package.logger.debug(`in chemSubstructureSearchLibrary, subFuncs all settled firing finish events, ${currentSearch}`);
            fireFinishEvents();
          }
        }).catch((err) => _package.logger.debug(`in chemSubstructureSearchLibrary, subFuncs settle failed: ${err}`));
      }
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
    mol.remove_hs_in_place(); // hydrogens can cause identical molecules to have different fingerprints
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
      const fp = getRDKitFpAsUint8Array(mol, fingerprint);
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

