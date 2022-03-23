onmessage = function(e) {
  const length = e.data[0]; const fpColumn = e.data[1]; const outDimensions = e.data[2]; let maxDistance = e.data[3];
  let maxDistanceSteps = e.data[4]; const radiusPercent = e.data[5]; const lambda0 = e.data[6];
  const lambda1 = e.data[7]; const steps = e.data[8]; const cycles = e.data[9];
  let a = 0; let b = 0;

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

  const wordCount = 4;

  function getIndexes() {
    b = Math.floor(Math.random() * length);
    for (a = Math.floor(Math.random() * length); b == a; b = Math.floor(Math.random() * length));
  }

  function _euclideanDistance(coordinates, row1, row2) {
    let d = 0.0;
    for (let i = 0; i < coordinates.length; ++i) {
      const diff = coordinates[i][row1] - coordinates[i][row2];
      d += diff * diff;
    }
    return Math.sqrt(d);
  }

  function calculateBitCount(mask) {
    let count = 0;
    for (let subMask = mask; subMask != 0; subMask >>>= 8) {
      count += _onBitCount[subMask & 0xff];
    }
    return count;
  }

  function _tanimotoSimilarity(maskA, maskB) {
    let bothABcount = 0; 
    let onlyAcount = 0; 
    let onlyBcount = 0;
    for (let i = 0; i < wordCount; ++i) {
      const curMaskA = maskA[i]; const curMaskB = maskB[i];
      const bothAB = (curMaskA & curMaskB);
      const onlyA = (curMaskA ^ bothAB);
      const onlyB = (curMaskB ^ bothAB);
      bothABcount += calculateBitCount(bothAB);
      onlyAcount += calculateBitCount(onlyA);
      onlyBcount += calculateBitCount(onlyB);
    }
    return bothABcount / (onlyAcount + onlyBcount + bothABcount);
  }

  const distance = (x, y) => 1.0 - _tanimotoSimilarity(fpColumn[x], fpColumn[y]);

  if (maxDistanceSteps == null) {
    maxDistanceSteps = length * Math.floor((length - 1) / 2);
  }
  if (maxDistance == null) {
    maxDistance = -1e9;
    for (let n = 0; n < maxDistanceSteps; ++n) {
      getIndexes();
      const d = distance(a, b);
      maxDistance = Math.max(d, maxDistance);
    }
  }

  const radius = (radiusPercent == 0.0) ? maxDistance : maxDistance * radiusPercent;
  const epsilon = 1e-9; let lambda = lambda0;
  const coordinates = new Array(outDimensions);
  for (let n = 0; n < outDimensions; ++n) {
    const dim = new Float32Array(length);
    for (let m = 0; m < length; ++m) {
      dim[m] = Math.random();
    }
    coordinates[n] = dim;
  }

  for (let i = 0; i < cycles; ++i) {
    for (let j = 0; j < steps; ++j) {
      getIndexes();
      const dab = _euclideanDistance(coordinates, a, b);
      const simab = distance(a, b);
      if (simab > radius && dab >= simab) {
        continue;
      } else {
        const t1 = lambda * 0.5 * (simab - dab) / (dab + epsilon);
        for (let k = 0; k < outDimensions; ++k) {
          coordinates[k][a] += t1 * (coordinates[k][a] - coordinates[k][b]);
          coordinates[k][b] += t1 * (coordinates[k][b] - coordinates[k][a]);
        }
      }
    }
    lambda -= (lambda0 - lambda1) / (cycles - 1.0);
    if (lambda < lambda1) {
      break;
    }
  }

  postMessage(coordinates);
};
