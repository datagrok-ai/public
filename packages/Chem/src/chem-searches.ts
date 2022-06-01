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
  invalidatedForVersion = '.invalideted.for.version',
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

async function _invalidate(molCol: DG.Column) {
  if (!(molCol.getTag(FING_COL_TAGS.invalidatedForVersion) == String(molCol.version))) {
    const {molIdxToHash, hashToMolblock} =
        await (await getRdKitService()).initMoleculesStructures(molCol.toList(), false, false);
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
    molCol.setTag(FING_COL_TAGS.invalidatedForVersion, String(molCol.version + 1));
  }
}

async function getUint8ArrayFingerprints(
  molCol: DG.Column, fingerprintsType: Fingerprint = Fingerprint.Morgan, endSection = true): Promise<Uint8Array[]> {
  await chemBeginCriticalSection();
  const colNameTag = '.' + fingerprintsType + '.Column';
  const colVerTag = '.' + fingerprintsType + '.Version';

  try {
    const fingerprintColumnName = molCol.getTag(colNameTag) ?? '~' + molCol.name + fingerprintsType + 'Fingerprints';

    if (molCol.getTag(colNameTag) &&
        molCol.getTag(colVerTag) == String(molCol.version)) {
      const fingCol = molCol.dataFrame.columns.byName(fingerprintColumnName);
      return fingCol.toList();
    } else {
      await _invalidate(molCol);
      const fingerprints = await (await getRdKitService()).getFingerprints(fingerprintsType);
      if (molCol.dataFrame) {
        const df = molCol.dataFrame;
        const newCol: DG.Column<Uint8Array> = df.columns.names().includes(fingerprintColumnName) ?
          df.columns.byName(fingerprintColumnName) :
          DG.Column.fromType(DG.COLUMN_TYPE.BYTE_ARRAY, fingerprintColumnName, fingerprints.length);

        for (let i = 0; i < fingerprints.length; ++i)
          newCol.set(i, fingerprints[i]);

        if (df.columns.names().includes(fingerprintColumnName))
          df.columns.replace(fingerprintColumnName, newCol);
        else
          df.columns.add(newCol);

        molCol.setTag(colNameTag, fingerprintColumnName);
        molCol.setTag(colVerTag, String(molCol.version + 2));
        molCol.setTag(FING_COL_TAGS.invalidatedForVersion, String(molCol.version + 1));
      }
      return fingerprints;
    }
  } finally {
    if (endSection)
      chemEndCriticalSection();
  }
}

export async function chemGetFingerprints(...args: [DG.Column, Fingerprint?, boolean?]): Promise<BitArray[]> {
  return (await getUint8ArrayFingerprints(...args)).map(el => rdKitFingerprintToBitArray(el));
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

export async function chemFindSimilar(molStringsColumn: DG.Column, queryMolString = '',
  settings: { [name: string]: any } = {}) : Promise<DG.DataFrame | null> {
  assure.notNull(molStringsColumn, 'molStringsColumn');
  assure.notNull(queryMolString, 'queryMolString');

  const fingerprints = await chemGetFingerprints(molStringsColumn)!;
  return queryMolString.length != 0 ?
    _chemFindSimilar(molStringsColumn, fingerprints, queryMolString, settings) : null;
}

export async function chemSubstructureSearchLibrary(
  molStringsColumn: DG.Column, molString: string, molStringSmarts: string): Promise<DG.BitSet> {
  try {
    const result = DG.BitSet.create(molStringsColumn.length);
    if (molString.length != 0) {
      const patternFps: Uint8Array[] = await getUint8ArrayFingerprints(molStringsColumn, Fingerprint.Pattern);
      const matches = await (await getRdKitService()).searchSubstructure(molString, molStringSmarts, patternFps);
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

