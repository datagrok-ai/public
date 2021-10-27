var rdKitParallel = null;
async function _initRdKitWorkers() {
  let foo = _initRdKitWorkers;
  if (typeof foo.initialized == 'undefined' || !foo.initialized) {
    rdKitParallel = new RdKitParallel();
    await rdKitParallel.init(rdKitWorkerWebRoot);
    _initRdKitWorkers.initialized = true;
  }
}

async function _cacheByAction(params, invalidator) {

  let invalidateCache = false;
  let {foo, column, query} = params;

  if (
    typeof foo.cachedForCol == 'undefined' &&
    typeof foo.cachedStructure == 'undefined') {
    foo.cachedForCol = null;
    foo.cachedStructure = null;
    foo.cachedForColVersion = null;
    invalidateCache = true;
  }

  if (column !== foo.cachedForCol || (column.version !== foo.cachedForColVersion) || query == null) {
    invalidateCache = true;
  }

  if (invalidateCache) {
    await invalidator(params);
    foo.cachedForCol = column;
    foo.cachedForColVersion = column.version;
    console.log('Molecule cache was invalidated');
  }
}

function _morganFP(molString, fp_length = 128, fp_radius = 2) {

  if (molString.length == 0) {
    console.error(
      "Possibly an empty molString: `" + molString + "`");
  } else {
    try {
      let mol = rdKitModule.get_mol(molString);
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

function _moleculesToFingerprints(molStringsColumn, settings = {}) {
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

function _foldFingerprint(bitsetFp, newLength) {
  let result = DG.BitSet.create(newLength);
  for (let idx in bitsetFp.getSelectedIndexes())
    result.set(idx % newLength, true, false);
  return result;
}

function _fingerprintSimilarity(bitsetFp1, bitsetFp2) {
  const len1 = bitsetFp1.length;
  const len2 = bitsetFp2.length;
  if (len1 < len2)
    bitsetFp2 = _foldFingerprint(bitsetFp2, len1);
  else if (len2 < len1)
    bitsetFp1 = _foldFingerprint(bitsetFp1, len2);
  return bitsetFp1.similarityTo(bitsetFp2, 'tanimoto'); // tanimotoSimilarity(fp1, fp2);
}

// Only this function receives {sorted} in settings
function _chemSimilarityScoringByFingerprints(fingerprintCol, fingerprint, molStringsColumn, settings) {

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
        sortedMolStrings.removeAt(n, limit - n);
        sortedScores.removeAt(n, limit - n);
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
async function _chemSimilarityScoring(molStringsColumn, molString, settings) {

  // await _initRdKitWorkers();

  _cacheByAction(
    {foo: _chemSimilarityScoring, column: molStringsColumn, query: molString},
    (params) => {
      let {foo, column, query} = params;
      foo.cachedStructure = _moleculesToFingerprints(molStringsColumn, settings);
    });

  if (molString.length != 0) {
    const fingerprintCol = _chemSimilarityScoring.cachedStructure;
    const fingerprint = _moleculesToFingerprints(DG.Column.fromStrings('molecules', [molString]), settings).get(0);
    return _chemSimilarityScoringByFingerprints(fingerprintCol, fingerprint, molStringsColumn, settings);
  } else {
    return null;
  }

}

async function chemGetSimilarities(molStringsColumn, molString = "", settings = {}) {
  settings.sorted = false;
  return _chemSimilarityScoring(molStringsColumn, molString, settings);
}

async function chemFindSimilar(molStringsColumn, molString = "", settings = {}) {
  settings.sorted = true;
  return _chemSimilarityScoring(molStringsColumn, molString, settings);
}

function chemSubstructureSearchGraph(molStringsColumn, molString) {

  const len = molStringsColumn.length;
  let result = DG.BitSet.create(len);
  if (molString.length == 0) {
    return result;
  }
  let subMol = rdKitModule.get_mol(molString);
  for (let i = 0; i < len; ++i) {
    let item = molStringsColumn.get(i);
    try {
      let mol = rdKitModule.get_mol(item);
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

async function chemSubstructureSearchLibrary(molStringsColumn, molString, molStringSmarts) {

  await _initRdKitWorkers();

  await _cacheByAction({
      foo: chemSubstructureSearchLibrary,
      column: molStringsColumn,
      query: molString + " | " + molStringSmarts
    },
    async (params) => {
      const { molIdxToHash, hashToMolblock } = await rdKitParallel.substructInit(molStringsColumn.toList());
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
        // without the following two lines, categories are not updated
        molStringsColumn.setCategoryOrder(molStringsColumn.categories);
        molStringsColumn.setCategoryOrder(null);
      }
      // TODO: avoid creating an additional array here
    }
  );

  let result = DG.BitSet.create(molStringsColumn.length);
  if (molString.length != 0) {
    const matches = await rdKitParallel.substructSearch(molString, molStringSmarts);
    for (let match of matches)
      result.set(match, true, false);
  }

  return result;
}

onmessage = function (e) {
  var length = e.data[0], fpColumn = e.data[1], outDimensions = e.data[2], maxDistance = e.data[3], 
      maxDistanceSteps = e.data[4], radiusPercent = e.data[5], lambda0 = e.data[6], 
      lambda1 = e.data[7], steps = e.data[8], cycles = e.data[9]; 
  var a = 0, b = 0;

  const _onBitCount = Int8Array.from([
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8]);

  const wordCount = 4, byteSize = 8, wordSize = 32;

  function getIndexes() {
    b = Math.floor(Math.random() * length);
    for (a = Math.floor(Math.random() * length); b == a; b = Math.floor(Math.random() * length));
  }

  function _euclideanDistance(coordinates, row1, row2) {
    var d = 0.0;
    for (let i = 0; i < coordinates.length; ++i) {
      var diff = coordinates[i][row1] - coordinates[i][row2];
      d += diff * diff;
    }
    return Math.sqrt(d);
  }

  function calculateBitCount(mask) {
    var count = 0;
    for (let j = 0; j < wordSize / byteSize; ++j) {
      var subMask = (mask >> j * 8) & 0xff;
      count += _onBitCount[subMask];
    }
    return count;
  }

  function _tanimotoSimilarity(maskA, maskB) {
    var bothABcount = 0, onlyAcount = 0, onlyBcount = 0;
    for (let i = 0; i < wordCount; ++i) {
      var curMaskA = maskA[i], curMaskB = maskB[i];
      var bothAB = (curMaskA & curMaskB);
      var onlyA = (curMaskA ^ bothAB);
      var onlyB = (curMaskB ^ bothAB);
      bothABcount += calculateBitCount(bothAB);
      onlyAcount += calculateBitCount(onlyA);
      onlyBcount += calculateBitCount(onlyB);
    }  
    return bothABcount / (onlyAcount + onlyBcount + bothABcount);
  }

  var distance = (x, y) => 1.0 - _tanimotoSimilarity(fpColumn[x], fpColumn[y]);

  if (maxDistanceSteps == null)
    maxDistanceSteps = length * Math.floor((length - 1) / 2);
  if (maxDistance == null) {
    maxDistance = -1e9;
    for (let n = 0; n < maxDistanceSteps; ++n) {
      getIndexes();
      var d = distance(a, b);
      maxDistance = Math.max(d, maxDistance);
    }
  }

  var radius = (radiusPercent == 0.0) ? maxDistance : maxDistance * radiusPercent;
  var epsilon = 1e-9, lambda = lambda0;
  var coordinates = new Array(outDimensions);
  for (let n = 0; n < outDimensions; ++n) {
    var dim = new Float32Array(length);
    for (let m = 0; m < length; ++m)
      dim[m] = Math.random();
    coordinates[n] = dim;
  }    

  for (let i = 0; i < cycles; ++i) {
    for (let j = 0; j < steps; ++j) {
      getIndexes();
      var dab = _euclideanDistance(coordinates, a, b);
      var simab = distance(a, b);
      if (simab > radius && dab >= simab)
        continue;
      else {
        var t1 = lambda * 0.5 * (simab - dab) / (dab + epsilon);
        for (let k = 0; k < outDimensions; ++k) {
          coordinates[k][a] += t1 * (coordinates[k][a] - coordinates[k][b]);
          coordinates[k][b] += t1 * (coordinates[k][b] - coordinates[k][a]);
        }
      }
    }
    lambda -= (lambda0 - lambda1) / (cycles - 1.0);
    if (lambda < lambda1)
      break;
  }
  
  postMessage(coordinates);
};