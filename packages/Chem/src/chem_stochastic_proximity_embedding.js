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
    for (let subMask = mask; subMask != 0; subMask >>>= 8) {
      count += _onBitCount[subMask & 0xff];
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