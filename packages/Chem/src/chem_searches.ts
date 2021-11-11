import {RdKitParallel} from './rdkit_parallel';
import * as DG from 'datagrok-api/dg';

let searchesRdKitModule: any = null;
let rdKitParallel: RdKitParallel | null = null;

export async function setSearchesRdKitModule(module: any)
{
  searchesRdKitModule = module;
}

let _initRdKitWorkersState: {[_:string]: any} = {};
async function _initRdKitWorkers(workerWebRoot: string) {
  if (typeof _initRdKitWorkersState.initialized == 'undefined' || !_initRdKitWorkersState.initialized) {
    rdKitParallel = new RdKitParallel();
    await rdKitParallel.init(workerWebRoot);
    _initRdKitWorkersState.initialized = true;
  }
}

interface CacheParams {
  cachedForCol : DG.Column | null;
  cachedForColVersion: number | null;
  cachedStructure: any;
  column: DG.Column | null;
  query: string | null;
};

const cacheParamsDefaults: CacheParams = {
  cachedForCol: null,
  cachedForColVersion: null,
  cachedStructure: null,
  column: null,
  query: null
};

async function _cacheByAction(params: CacheParams, invalidator: (_: CacheParams) => Promise<void>) {

  let invalidateCache = false;

  if (
    params.cachedForCol === null &&
     params.cachedStructure === null) {
    invalidateCache = true;
  }

  if (params.column !== params.cachedForCol ||
      (params.column!.version !== params.cachedForColVersion) || params.query == null) {
    invalidateCache = true;
  }

  if (invalidateCache) {
    await invalidator(params);
    params.cachedForCol = params.column;
    params.cachedForColVersion = params.column!.version;
    console.log('Molecule cache was invalidated');
  }
}

function _morganFP(molString: string, fp_length = 128, fp_radius = 2) {

  if (molString.length == 0) {
    console.error(
      "Possibly an empty molString: `" + molString + "`");
  } else {
    try {
      let mol = searchesRdKitModule.get_mol(molString);
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

function _foldFingerprint(bitsetFp: DG.BitSet, newLength: number) {
  let result = DG.BitSet.create(newLength);
  for (let idx in bitsetFp.getSelectedIndexes())
    result.set(+idx % newLength, true, false);
  return result;
}

function _fingerprintSimilarity(bitsetFp1: DG.BitSet, bitsetFp2: DG.BitSet) {
  const len1 = bitsetFp1.length;
  const len2 = bitsetFp2.length;
  if (len1 < len2)
    bitsetFp2 = _foldFingerprint(bitsetFp2, len1);
  else if (len2 < len1)
    bitsetFp1 = _foldFingerprint(bitsetFp1, len2);
  return bitsetFp1.similarityTo(bitsetFp2, 'tanimoto'); // tanimotoSimilarity(fp1, fp2);
}

// Only this function receives {sorted} in settings
function _chemSimilarityScoringByFingerprints(
    fingerprintCol: DG.Column, fingerprint: DG.BitSet,
    molStringsColumn: DG.Column, settings: { [name: string]: any }) {
  const len = fingerprintCol.length;
  let distances = DG.Column.fromType(DG.TYPE.FLOAT, 'distances', len);
  for (let row = 0; row < len; ++row) {
    const fp = fingerprintCol.get(row);
    distances.set(row, fp == null ? 1.0 : _fingerprintSimilarity(fingerprint, fp));
  }
  if (settings.hasOwnProperty('sorted') && settings.sorted === true) {
    const limit = Math.min((settings.hasOwnProperty('limit') ? settings.limit : len), len);
    const minScore = settings.hasOwnProperty('minScore') ? settings.minScore : 0.0;
    let sortedIndices = Array.from(Array(len).keys()).sort((i1, i2) => {
      const a1 = distances.get(i1);
      const a2 = distances.get(i2);
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
      const score = distances.get(idx);
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
  } else {
    return distances;
  }
}

// molStringsColumn and molString can be anything  RDKit supports:
// smiles, cxsmiles, molblock, v3Kmolblock, and inchi;
// see https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h

let _chemSimilarityScoringCacheParams = { ...cacheParamsDefaults };
async function _chemSimilarityScoring(molStringsColumn: DG.Column, molString: string, settings: { [name: string]: any })
    : Promise<ReturnType<typeof _chemSimilarityScoringByFingerprints> | null> {

  // await _initRdKitWorkers();

  _chemSimilarityScoringCacheParams.column = molStringsColumn;
  _chemSimilarityScoringCacheParams.query = molString;

  await _cacheByAction(
      _chemSimilarityScoringCacheParams,
    async (params: CacheParams) => {
      params.cachedStructure = moleculesToFingerprints(molStringsColumn, settings);
    });

  if (molString.length != 0) {
    const fingerprintCol = _chemSimilarityScoringCacheParams.cachedStructure;
    const fingerprint = moleculesToFingerprints(DG.Column.fromStrings('molecules', [molString]), settings).get(0);
    return _chemSimilarityScoringByFingerprints(fingerprintCol, fingerprint, molStringsColumn, settings);
  } else {
    return null;
  }

}

export async function chemGetSimilarities(
    molStringsColumn: DG.Column, molString = "", settings: { [name: string]: any } = {}) {
  settings.sorted = false;
  return _chemSimilarityScoring(molStringsColumn, molString, settings);
}

export async function chemFindSimilar(
    molStringsColumn: DG.Column, molString = "", settings: { [name: string]: any } = {}) {
  settings.sorted = true;
  return _chemSimilarityScoring(molStringsColumn, molString, settings);
}

export function chemSubstructureSearchGraph(molStringsColumn: DG.Column, molString: string) {

  const len = molStringsColumn.length;
  let result = DG.BitSet.create(len);
  if (molString.length == 0) {
    return result;
  }
  let subMol = searchesRdKitModule.get_mol(molString);
  for (let i = 0; i < len; ++i) {
    let item = molStringsColumn.get(i);
    try {
      let mol = searchesRdKitModule.get_mol(item);
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

let _chemSubstructureSearchLibraryParams = {...cacheParamsDefaults};
export async function chemSubstructureSearchLibrary(
    molStringsColumn: DG.Column, molString: string, molStringSmarts: string, workerWebRoot: any) {

  await _initRdKitWorkers(workerWebRoot);

  _chemSubstructureSearchLibraryParams.column = molStringsColumn;
  _chemSubstructureSearchLibraryParams.query = molString + "|" + molStringSmarts;

  await _cacheByAction(
    _chemSubstructureSearchLibraryParams,
    async (params) => {
      // TODO: avoid creating an additional array here
      const { molIdxToHash, hashToMolblock }
        = await rdKitParallel!.substructInit(molStringsColumn.toList());
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
    }
  );

  let result = DG.BitSet.create(molStringsColumn.length);
  if (molString.length != 0) {
    const matches = await rdKitParallel!.substructSearch(molString, molStringSmarts);
    for (let match of matches)
      result.set(match, true, false);
  }

  return result;
}