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
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {assure} from '@datagrok-libraries/utils/src/test';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

const enum FING_COL_TAGS {
  invalidatedForVersion = '.invalideted.for.version',
}

let LastColumnInvalidated: string = '';

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

function _chemGetDiversities(limit: number, molStringsColumn: DG.Column, fingerprints: BitArray[]): string[] {
  limit = Math.min(limit, fingerprints.length);
  const indexes = ArrayUtils.indexesOf(fingerprints, (f) => f != null);

  const diverseIndexes = getDiverseSubset(fingerprints, indexes.length, limit,
    (i1, i2) => 1 - tanimotoSimilarity(fingerprints[indexes[i1]], fingerprints[indexes[i2]]));

  const molIds: number[] = [];
  const diversities = new Array(limit).fill('');

  for (let i = 0; i < limit; i++)
    diversities[i] = molStringsColumn.get(indexes[diverseIndexes[i]]);

  return diversities;
}

function colInvalidated(col: DG.Column): Boolean {
  return (LastColumnInvalidated == col.name &&
    col.getTag(FING_COL_TAGS.invalidatedForVersion) == String(col.version));
}

async function _invalidate(molCol: DG.Column) {
  if (!colInvalidated(molCol)) {
    await (await getRdKitService()).initMoleculesStructures(molCol.toList());
    LastColumnInvalidated = molCol.name;
    molCol.setTag(FING_COL_TAGS.invalidatedForVersion, String(molCol.version + 1));
  }
}

function checkForFingerprintsColumn(col: DG.Column, fingerprintsType: Fingerprint): Uint8Array[] | null {
  const colNameTag = '.' + fingerprintsType + '.Column';
  const colVerTag = '.' + fingerprintsType + '.Version';

  if (col.getTag(colNameTag) &&
      col.getTag(colVerTag) == String(col.version) &&
      col.dataFrame) {
    const fingCol = col.dataFrame.columns.byName(col.getTag(colNameTag));
    if (fingCol)
      return fingCol.toList();
  }
  return null;
}

function saveFingerprintsToCol(col: DG.Column, fgs: Uint8Array[], fingerprintsType: Fingerprint): void {
  if (!col.dataFrame)
    throw new Error('Column has no parent dataframe');

  const colNameTag = '.' + fingerprintsType + '.Column';
  const colVerTag = '.' + fingerprintsType + '.Version';

  const fingerprintColumnName = col.getTag(colNameTag) ??
    '~' + col.name + fingerprintsType + 'Fingerprints';
  const df = col.dataFrame;
  const newCol: DG.Column<Uint8Array> = df.columns.getOrCreate(
    fingerprintColumnName, DG.COLUMN_TYPE.BYTE_ARRAY, fgs.length);

  newCol.init((i) => fgs[i]);

  col.setTag(colNameTag, fingerprintColumnName);
  col.setTag(colVerTag, String(col.version + 2));
  col.setTag(FING_COL_TAGS.invalidatedForVersion, String(col.version + 1));
}

async function getUint8ArrayFingerprints(
  molCol: DG.Column, fingerprintsType: Fingerprint = Fingerprint.Morgan, useSection = true): Promise<Uint8Array[]> {
  if (useSection)
    await chemBeginCriticalSection();
  try {
    const fgsCheck = checkForFingerprintsColumn(molCol, fingerprintsType);
    if (fgsCheck)
      return fgsCheck;
    else {
      await _invalidate(molCol);
      const fingerprints = await (await getRdKitService()).getFingerprints(fingerprintsType);
      saveFingerprintsToCol(molCol, fingerprints, fingerprintsType);
      chemEndCriticalSection();
      return fingerprints;
    }
  } finally {
    if (useSection)
      chemEndCriticalSection();
  }
}

function substructureSearchPatternsMatch(molString: string, querySmarts: string, fgs: Uint8Array[]): BitArray {
  const patternFpUint8Length = 256;
  const result = new BitArray(fgs.length, false);
  let queryMol = null;
  try {
    queryMol = getRdKitModule().get_mol(molString, '{"mergeQueryHs":true}');
  } catch (e2) {
    queryMol?.delete();
    queryMol = null;
    if (querySmarts !== null && querySmarts !== '')
      queryMol = getRdKitModule().get_qmol(querySmarts);
    else
      throw new Error('Chem | SMARTS not set');
  }
  const fpRdKit = queryMol.get_pattern_fp_as_uint8array();
  checkEl:
  for (let i = 0; i < fgs.length; ++i) {
    for (let j = 0; j < patternFpUint8Length; ++j) {
      if ((fgs[i][j] & fpRdKit[j]) != fpRdKit[j])
        continue checkEl;
    }
    result.setBit(i, true, false);
  }
  return result;
}

export async function chemGetFingerprints(...args: [DG.Column, Fingerprint?, boolean?]): Promise<BitArray[]> {
  return (await getUint8ArrayFingerprints(...args)).map((el) => rdKitFingerprintToBitArray(el));
}

export async function chemGetSimilarities(molStringsColumn: DG.Column, queryMolString = '')
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn)!;

  return queryMolString.length != 0 ?
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      _chemGetSimilarities(queryMolString, fingerprints)) : null;
}

export async function chemGetDiversities(molStringsColumn: DG.Column, limit: number)
      : Promise<DG.Column | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(limit, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn)!;

  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'molecule', _chemGetDiversities(limit, molStringsColumn, fingerprints));
}

export async function chemFindSimilar(molStringsColumn: DG.Column, queryMolString = '',
  settings: { [name: string]: any } = {}) : Promise<DG.DataFrame | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn)!;
  return queryMolString.length != 0 ?
    _chemFindSimilar(molStringsColumn, fingerprints, queryMolString, settings) : null;
}

export async function chemSubstructureSearchLibrary(
  molStringsColumn: DG.Column, molString: string, molBlockFailover: string, usePatternFingerprints = false)
  : Promise<DG.BitSet> {
  await chemBeginCriticalSection();
  try {
    const result = DG.BitSet.create(molStringsColumn.length);
    if (molString.length != 0) {
      let matches: number[];
      if (usePatternFingerprints) {
        const fgs: Uint8Array[] = await getUint8ArrayFingerprints(molStringsColumn, Fingerprint.Pattern, false);
        const bitset: BitArray = substructureSearchPatternsMatch(molString, molBlockFailover, fgs);
        matches = await (await getRdKitService()).searchSubstructure(molString, molBlockFailover, bitset);
      } else {
        await _invalidate(molStringsColumn);
        matches = await (await getRdKitService()).searchSubstructure(molString, molBlockFailover);
      }
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
    if (fingerprint == Fingerprint.Morgan) {
      fp = mol.get_morgan_fp_as_uint8array(JSON.stringify({
        radius: defaultMorganFpRadius,
        nBits: defaultMorganFpLength,
      }));
    } else if (fingerprint == Fingerprint.Pattern)
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
