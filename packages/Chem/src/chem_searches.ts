import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitService} from './chem_common_rdkit';
import BitArray from "@datagrok-libraries/utils/src/bit-array";
import {rdKitFingerprintToBitArray, tanimoto} from './chem_common';

export function _morganFP(molString: string, fp_length = 128, fp_radius = 2) {
  if (molString.length == 0) {
    console.error(
      "Possibly an empty molString: `" + molString + "`");
  } else {
    try {
      let mol = getRdKitModule().get_mol(molString);
      let mfp = mol.get_morgan_fp(fp_radius, fp_length);
      mol.delete();
      return mfp;
    } catch (e) {
      console.error(
        "Possibly a malformed molString: `" + molString + "`");
      // Won't rethrow
    }
  }
  return '0'.repeat(fp_length);
}

export function moleculesToFingerprints(molStringsColumn: DG.Column, settings: { [name: string]: number } = { }) {
  const len = molStringsColumn.length;
  const fpLength = settings.hasOwnProperty('fpLength') ? settings.fpLength : 128;
  const fpRadius = settings.hasOwnProperty('fpRadius') ? settings.fpRadius : 2;
  let fingerprints = [];
  for (let i = 0; i < molStringsColumn.length; ++i) {
    let molString = molStringsColumn.get(i);
    let morganFp = _morganFP(molString, fpLength, fpRadius);
    const fingerprint = DG.BitSet.fromString(morganFp);
    fingerprints.push(fingerprint);
  }
  return DG.Column.fromList('object', 'fingerprints', fingerprints);
}

async function _chemFindSimilar(molStringsColumn: DG.Column,
    queryMolString: string, settings: { [name: string]: any }) {
  const len = molStringsColumn.length;
  const distances =  await _chemGetSimilarities(queryMolString);
  const limit = Math.min((settings.hasOwnProperty('limit') ? settings.limit : len), len);
  const minScore = settings.hasOwnProperty('minScore') ? settings.minScore : 0.0;
  let sortedIndices = Array.from(Array(len).keys()).sort((i1, i2) => {
    const a1 = distances[i1];
    const a2 = distances[i2];
    if (a2 < a1) return -1;
    if (a2 > a1) return +1;
    return 0; // a2.compareTo(a1)
  });
  let sortedMolStrings = DG.Column.fromType(DG.TYPE.STRING, 'molecule', limit);
  let sortedMolInd = DG.Column.fromType(DG.TYPE.INT, 'index', limit);
  sortedMolStrings.semType = DG.SEMTYPE.MOLECULE;
  let sortedScores = DG.Column.fromType(DG.TYPE.FLOAT, 'score', limit);
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
async function _chemGetSimilarities(queryMolString: string) {
  const fingerprints = _chemCache.similarityFingerprints!;
  let distances = new Array(fingerprints.length).fill(0.0);
  try {
    const mol = getRdKitModule().get_mol(queryMolString);
    const fp = mol.get_morgan_fp(2, 128);
    const sample = rdKitFingerprintToBitArray(fp, 128);
    for (let i = 0; i < fingerprints.length; ++i) {
      distances[i] = tanimoto(fingerprints[i], sample);
    }
  } finally {
    return distances;
  }
}

interface CacheParams {
  cachedForCol : DG.Column | null;
  cachedForColVersion: number | null;
  column: DG.Column | null;
  query: string | null;
  moleculesWereIndexed: boolean | null;
  similarityFingerprintsWereIndexed: boolean | null;
  similarityFingerprints: BitArray[] | null;
};

const cacheParamsDefaults: CacheParams = {
  cachedForCol: null,
  cachedForColVersion: null,
  column: null,
  query: null,
  moleculesWereIndexed: false,
  similarityFingerprintsWereIndexed: false,
  similarityFingerprints: null
};

let _chemCache = { ...cacheParamsDefaults };

async function _invalidate(molStringsColumn: DG.Column, queryMolString: string, includeFingerprints: boolean) {
  const sameColumnAndVersion = () =>
    molStringsColumn === _chemCache.cachedForCol &&
    molStringsColumn.version === _chemCache.cachedForColVersion;
  if (!sameColumnAndVersion() || queryMolString === null || queryMolString.length === 0) {
    _chemCache.cachedForCol = molStringsColumn;
    _chemCache.cachedForColVersion = molStringsColumn.version;
    _chemCache.moleculesWereIndexed = false;
    _chemCache.similarityFingerprintsWereIndexed = false;
    _chemCache.similarityFingerprints = null;
  }
  if (sameColumnAndVersion() && !_chemCache.moleculesWereIndexed) {
    // TODO: avoid creating an additional array here
    const { molIdxToHash, hashToMolblock }
      = await getRdKitService().initMoleculesStructures(molStringsColumn.toList());
    let i = 0;
    let needsUpdate = false;
    for (const item of molIdxToHash) {
      const notify = (i === molIdxToHash.length - 1);
      const molStr = hashToMolblock[item];
      if (molStr) {
        molStringsColumn.setString(i, molStr, notify);
        needsUpdate = true;
      }
      ++i;
    }
    if (needsUpdate) {
      // This seems to be the only way to trigger re-calculation of categories
      molStringsColumn.compact();
    }
    _chemCache.moleculesWereIndexed = true;
    _chemCache.similarityFingerprintsWereIndexed = false;
    _chemCache.similarityFingerprints = null;
    console.log(`RdKit molecules for a dataset ${_chemCache.cachedForCol?.name} invalidated`);
  }
  if (includeFingerprints) {
    if (sameColumnAndVersion() && !_chemCache.similarityFingerprintsWereIndexed) {
      await getRdKitService().initMorganFingerprints();
      _chemCache.similarityFingerprintsWereIndexed = true;
      _chemCache.similarityFingerprints = await getRdKitService().getMorganFingerprints();
      console.log(`Morgan fingerprints for a dataset ${_chemCache.cachedForCol?.name} invalidated`);
    }
  }
}

// molStringsColumn and molString can be anything  RDKit supports:
// smiles, cxsmiles, molblock, v3Kmolblock, and inchi;
// see https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h

export async function chemGetSimilarities(
    molStringsColumn: DG.Column, queryMolString = "",
    settings: { [name: string]: any } = {}) {
  await _invalidate(molStringsColumn, queryMolString, true);
  const result = queryMolString.length != 0 ?
    await DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      await _chemGetSimilarities(queryMolString))  : null;
  return result;
}

export async function chemFindSimilar(
    molStringsColumn: DG.Column, queryMolString = "", settings: { [name: string]: any } = {}) {
  await _invalidate(molStringsColumn, queryMolString, true);
  const result = queryMolString.length != 0 ?
    _chemFindSimilar(molStringsColumn, queryMolString, settings) : null;
  return result;
}

export function chemSubstructureSearchGraph(molStringsColumn: DG.Column, molString: string) {
  const len = molStringsColumn.length;
  let result = DG.BitSet.create(len);
  if (molString.length == 0) {
    return result;
  }
  let subMol = getRdKitModule().get_mol(molString);
  for (let i = 0; i < len; ++i) {
    let item = molStringsColumn.get(i);
    try {
      let mol = getRdKitModule().get_mol(item);
      let match = mol.get_substruct_match(subMol);
      if (match !== "{}")
        result.set(i, true, false);
      mol.delete();
    } catch (e) {
      console.error(
        "Possibly a malformed molString: `" + item + "`");
      // Won't rethrow
    }
  }
  subMol.delete();
  return result;
}

export async function chemSubstructureSearchLibrary(
    molStringsColumn: DG.Column, molString: string, molStringSmarts: string) {
  await _invalidate(molStringsColumn, molString, false);
  let result = DG.BitSet.create(molStringsColumn.length);
  if (molString.length != 0) {
    const matches = await getRdKitService().searchSubstructure(molString, molStringSmarts);
    for (let match of matches)
      result.set(match, true, false);
  }
  return result;
}