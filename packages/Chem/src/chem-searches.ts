import * as DG from 'datagrok-api/dg';
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
import {tanimotoSimilarity} from '@datagrok-libraries/utils/src/similarity-metrics';
import {assure} from '@datagrok-libraries/utils/src/test';

const enum FING_COL_TAGS{
  fingerprintsColumnName = '.fingerprints.column.name',
  fingerprintsColumnVersion = '.fingerprints.column.version',
  moleculesIndexed = '.molecules.were.indexed',
  morganFingerprintsIndexed = '.morgan.fingerprints.were.indexed',
  hasMorganFingerprints = '.has.morgan.fingerprints'
}

function _chemFindSimilar(molStringsColumn: DG.Column, fingerprints: BitArray[],
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
  const sortedMolStrings = DG.Column.fromType(DG.TYPE.STRING, 'molecule', limit);
  const sortedMolInd = DG.Column.fromType(DG.TYPE.INT, 'index', limit);
  sortedMolStrings.semType = DG.SEMTYPE.MOLECULE;
  const sortedScores = DG.Column.fromType(DG.TYPE.FLOAT, 'score', limit);
  for (let n = 0; n < limit; n++) {
    const idx = sortedIndices[n];
    const score = distances[idx];
    if (score < minScore) {
      sortedMolStrings.dataFrame.rows.removeAt(n, limit - n);
      sortedScores.dataFrame.rows.removeAt(n, limit - n);
      break;
    }
    sortedMolStrings.set(n, molStringsColumn.get(idx));
    sortedScores.set(n, score);
    sortedMolInd.set(n, idx);
  }
  return DG.DataFrame.fromColumns([sortedMolStrings, sortedScores, sortedMolInd]);
}

// Only this function receives {sorted} in settings
function _chemGetSimilarities(queryMolString: string, fingerprints: BitArray[]): number[] {
  const distances = new Array(fingerprints.length).fill(0.0);
  const sample = chemGetFingerprint(queryMolString, Fingerprint.Morgan);
  for (let i = 0; i < fingerprints.length; ++i)
    distances[i] = tanimotoSimilarity(fingerprints[i], sample);
  return distances;
}

class CacheParams {
  // cachedForCol : DG.Column | null = null;
  // cachedForColVersion: number | null = null;
  // column: DG.Column | null = null;
  // query: string | null = null;
  // moleculesWereIndexed: boolean | null = false;
  // morganFingerprintsWereIndexed: boolean | null = false;
  // morganFingerprints: BitArray[] | null = null;
  morganFingerprintsData: DG.Column | null = null;
  get morganFingerprints(): BitArray[] | null {
    if (this.morganFingerprintsData == null)
      return null;
    else
      return this.morganFingerprintsData.toList().map((el) => BitArray.fromBytes(el));
  }
}

const _chemCache = new CacheParams();


async function getMorganFingerprints(
  molCol: DG.Column, queryMolString: string | null,
  includeFingerprints: boolean, endSection = true): Promise<BitArray[]> {
  await chemBeginCriticalSection();
  let morganFingerprints = null;
  try {
    const fingerprintColumnName = '~' + molCol.name + 'Fingerprints';
    const sameColumnAndVersion1 = () =>
      molCol.getTag(FING_COL_TAGS.fingerprintsColumnName) &&
      molCol.getTag(FING_COL_TAGS.fingerprintsColumnVersion) == String(molCol.version);

    if (sameColumnAndVersion1() && (molCol.getTag(FING_COL_TAGS.hasMorganFingerprints) == 'true')) {
      const fingCol = molCol.dataFrame.columns.byName(fingerprintColumnName);
      return fingCol.toList().map((el) => BitArray.fromBytes(el));
    }

    if (!sameColumnAndVersion1() || !(molCol.getTag(FING_COL_TAGS.hasMorganFingerprints) == 'true') &&
      (queryMolString === null || queryMolString.length === 0)) {
      molCol.setTag(FING_COL_TAGS.fingerprintsColumnName, fingerprintColumnName);
      molCol.setTag(FING_COL_TAGS.moleculesIndexed, 'false');
      molCol.setTag(FING_COL_TAGS.morganFingerprintsIndexed, 'false');
      molCol.setTag(FING_COL_TAGS.fingerprintsColumnVersion, String(molCol.version + 1));
    }
    if (sameColumnAndVersion1() && !(molCol.getTag(FING_COL_TAGS.moleculesIndexed) == 'true')) {
      // TODO: avoid creating an additional array here
      const {molIdxToHash, hashToMolblock} =
        await (await getRdKitService()).initMoleculesStructures(molCol.toList());
      let i = 0;
      if (molIdxToHash.length > 0) {
        let needsUpdate = false;
        for (const item of molIdxToHash) {
          const notify = (i === molIdxToHash.length - 1);
          const molStr = hashToMolblock[item];
          if (molStr) {
            molCol.setString(i, molStr, notify);
            needsUpdate = true;
          }
          ++i;
        }
        if (needsUpdate) {
          // This seems to be the only way to trigger re-calculation of categories
          molCol.compact();
        }
      }
      molCol.setTag(FING_COL_TAGS.moleculesIndexed, 'true');
      molCol.setTag(FING_COL_TAGS.morganFingerprintsIndexed, 'false');
      molCol.setTag(FING_COL_TAGS.fingerprintsColumnVersion, String(molCol.version + 1));

      console.log(`RdKit molecules for a dataset ${molCol?.name} invalidated`);
    }
    if (includeFingerprints) {
      if (sameColumnAndVersion1() && !(molCol.getTag(FING_COL_TAGS.morganFingerprintsIndexed) == 'true')) {
        await (await getRdKitService()).initMorganFingerprints();
        molCol.setTag(FING_COL_TAGS.morganFingerprintsIndexed, 'true');
        molCol.setTag(FING_COL_TAGS.fingerprintsColumnVersion, String(molCol.version + 1));
        morganFingerprints = await (await getRdKitService()).getMorganFingerprints();
        const newCol = DG.Column.fromType(DG.COLUMN_TYPE.BYTE_ARRAY,
          fingerprintColumnName,
          morganFingerprints.length);

        for (let i = 0; i < morganFingerprints.length; i++)
          newCol.set(i, new Uint8Array(morganFingerprints[i].buffer.buffer));

        const df = molCol.dataFrame;
        if (df.columns.names().includes(fingerprintColumnName))
          df.columns.replace(fingerprintColumnName, newCol);
        else
          df.columns.add(newCol);

        molCol.setTag(FING_COL_TAGS.hasMorganFingerprints, 'true');
        console.log(`Morgan fingerprints for a dataset ${molCol?.name} invalidated`);
      }
      // to delete
    }
    return morganFingerprints;
  } finally {
    if (endSection)
      chemEndCriticalSection();
  }
}

// molStringsColumn and molString can be anything  RDKit supports:
// smiles, cxsmiles, molblock, v3Kmolblock, and inchi;
// see https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h

export async function chemGetSimilarities(molStringsColumn: DG.Column, queryMolString = '')
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await getMorganFingerprints(molStringsColumn, queryMolString, true)!;

  return queryMolString.length != 0 ?
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      _chemGetSimilarities(queryMolString, fingerprints)) : null;
}

export async function chemFindSimilar(molStringsColumn: DG.Column, queryMolString = '',
  settings: { [name: string]: any } = {}) : Promise<DG.DataFrame | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await getMorganFingerprints(molStringsColumn, queryMolString, true)!;
  return queryMolString.length != 0 ?
    _chemFindSimilar(molStringsColumn, fingerprints, queryMolString, settings) : null;
}

export function chemSubstructureSearchGraph(molStringsColumn: DG.Column, molString: string): DG.BitSet {
  const len = molStringsColumn.length;
  const result = DG.BitSet.create(len);
  if (molString.length == 0)
    return result;
  let subMol = null;
  try {
    subMol = getRdKitModule().get_mol(molString);
    for (let i = 0; i < len; ++i) {
      const item = molStringsColumn.get(i);
      let mol = null;
      try {
        mol = getRdKitModule().get_mol(item);
        const match = mol.get_substruct_match(subMol);
        if (match !== '{}')
          result.set(i, true, false);
      } catch (e) {
        console.error('Possibly a malformed molString: `' + item + '`');
        // Explicitly won't rethrow
      } finally {
        mol?.delete();
      }
    }
    return result;
  } finally {
    subMol?.delete();
  }
}

export async function chemSubstructureSearchLibrary(
  molStringsColumn: DG.Column, molString: string, molStringSmarts: string): Promise<DG.BitSet> {
  // ??
  // await _invalidate(molStringsColumn, molString, false, false);
  try {
    const result = DG.BitSet.create(molStringsColumn.length);
    if (molString.length != 0) {
      const matches = await (await getRdKitService()).searchSubstructure(molString, molStringSmarts);
      for (const match of matches)
        result.set(match, true, false);
    }
    return result;
  } finally {
    chemEndCriticalSection();
  }
}

export function chemGetFingerprint(molString: string, fingerprint: Fingerprint): BitArray {
  let mol = null;
  try {
    mol = getRdKitModule().get_mol(molString);
    let fp;
    if (fingerprint == Fingerprint.Morgan)
      fp = mol.get_morgan_fp_as_uint8array(defaultMorganFpRadius, defaultMorganFpLength);
    /*else if (fingerprint == Fingerprint.RDKit)
      fp = mol.get_rdkit_fp(defaultMorganFpRadius, defaultMorganFpLength);*/
    else if (fingerprint == Fingerprint.Pattern)
      fp = mol.get_pattern_fp_as_uint8array();
    else
      throw new Error(`${fingerprint} does not match any fingerprint`);

    return rdKitFingerprintToBitArray(fp);
  } catch {
    throw new Error(`Chem | Possibly a malformed molString: ${molString}`);
  } finally {
    mol?.delete();
  }
}

export async function chemGetFingerprints(molStringsColumn: DG.Column, fingerprint: Fingerprint): Promise<BitArray[]> {
  const len = molStringsColumn.length;
  let fingerprints: BitArray[] = [];
  const fallbackCountForSyncExecution = 150;
  if (len <= fallbackCountForSyncExecution) {
    for (let i = 0; i < len; ++i) {
      try {
        fingerprints.push(chemGetFingerprint(molStringsColumn.get(i), fingerprint));
      } catch {
        fingerprints.push(new BitArray(defaultMorganFpLength));
      }
    }
  } else
    fingerprints = await getMorganFingerprints(molStringsColumn, null, true)!;
  return fingerprints;
}
