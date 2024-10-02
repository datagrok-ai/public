/* eslint-disable max-len */

// in all the functions below, the variables a and b are assumed to be arrays of uint32/f32
//values which are infered from the code this chunk is injected into
// also, we have access to computeInfo struct, which contains the following fields:
// computeInfo.entrySizes: array of arrays of u32 containing the sizes of the entries
// other fields are specific to the distance function should be injected from the main script that calls this function,
// and should be available in the supplementaryInfo struct
// like the similarity matrix for monomer chemical distance.
// the getProcessInfo function should return correct buffer allocation mechanism for the supplementaryInfo,
// for every entry list
// the maxDistance variable is also assumed to be available in the
// scope of the function, in case of knn it is the distance in the last postion of knn on this index,
// in case of sparse matrix, it can be just the threshold for the distance.

// hamming distance for sequnences of uint32 arrays of max length ${maxArraySize}
export function webGPUHamming(_maxArraySize: number, entryIndex: number) {
  return `
  let aLength: u32 = computeInfo.entrySizes[${entryIndex}][aIndex];
  let bLength: u32 = computeInfo.entrySizes[${entryIndex}][bIndex];
  let maxLength: u32 = max(aLength, bLength);
  let minLength: u32 = min(aLength, bLength);
  let sizeDiff: u32 = maxLength - minLength;
  
  let maxIntDistance = ceil(maxDistance * f32(maxLength)) - f32(sizeDiff);

  var diff: f32 = 0.0;
  for (var i = 0u; i < ${_maxArraySize}; i = i + 1u) {
      diff = diff + f32(a[i] != b[i]);
      if (diff > maxIntDistance) {
            return 1.0;
        }
  }
  diff += f32(sizeDiff);
  return diff / ${_maxArraySize};
  `;
}

export function webGPUMonomerChemicalDistance(_maxArraySize: number, entryIndex: number) {
  // it is assumet that suppInfo struct contains correct matrix called similarityMatrix${entryIndex}, (similarityMatrix0, similarityMatrix1, etc)
  // this should be guaranteed by the getProcessInfo function.
  return `
  let aLength: u32 = computeInfo.entrySizes[${entryIndex}][aIndex];
  let bLength: u32 = computeInfo.entrySizes[${entryIndex}][bIndex];
  let maxLength: u32 = max(aLength, bLength);
  let minLength: u32 = min(aLength, bLength);
  let sizeDiff: u32 = maxLength - minLength;
  
  let maxIntDistance = ceil(maxDistance * f32(maxLength)) - f32(sizeDiff);

  let simMatrix = &(suppInfo.similarityMatrix${entryIndex}); // using pointers make things faster
  var diff: f32 = 0.0;
  for (var i = 0u; i < ${_maxArraySize}; i = i + 1u) {
      diff = diff + 1.0 - (*simMatrix)[u32(a[i])][u32(b[i])];
      if (diff > maxIntDistance) {
            return 1.0;
        }
  }
  diff += f32(sizeDiff);
  return diff / ${_maxArraySize};
  `;
}

export function webGPULevenstein(maxArraySize: number, entryIndex: number) {
  return `
  let aLength: u32 = computeInfo.entrySizes[${entryIndex}][aIndex];
  let bLength: u32 = computeInfo.entrySizes[${entryIndex}][bIndex];
  let maxLength: u32 = max(aLength, bLength);
  let minLength: u32 = min(aLength, bLength);

  let maxIntDistance = ceil(maxDistance * f32(maxLength));

  // we will store two arrays as matrix and swap the working indices per pass.
  // this way we can reduce memory usage per computation to just O(aLength)
  // the grid will have aLength + 1 columns and bLength + 1 rows
  // this will be guaranteed by iteration, but the array sizes must be known at compile time, so we will use a fixed size of maxArraySize
  var dynamicPassMat: array<array<f32, ${maxArraySize + 1}u>, 2>; // initialize to 0
  
  var prevIndex: u32 = 0;
  var curIndex: u32 = 1; // we will swap these indices per pass

  // initialize the first row
  for (var i = 0u; i <= aLength; i = i + 1u) {
      dynamicPassMat[prevIndex][i] = f32(i);
  }

  // iterate over the rows
  for (var i = 1u; i <= bLength; i = i + 1u) {
      dynamicPassMat[curIndex][0] = f32(i);
      var minEntry: f32 = f32(maxLength);
      let prevRow = &dynamicPassMat[prevIndex];
      let curRow = &dynamicPassMat[curIndex];
      let bMon = u32(b[i - 1]);
      for (var j = 1u; j <= aLength; j = j + 1u) {
          var cost: f32 = f32(a[j - 1] != bMon);
          var res: f32 = min(
              min(
                (*prevRow)[j] + 1.0, // deletion
                (*curRow)[j - 1] + 1.0, // insertion
              ),
              (*prevRow)[j - 1] + cost // substitution
              );
              (*curRow)[j] = res;
          if (res < minEntry) {
              minEntry = res;
          }
      }
      // swap the indices
      let temp: u32 = prevIndex;
      prevIndex = curIndex;
      curIndex = temp;
      if (minEntry > maxIntDistance) {
          return 1.0;
      }
  }

  return dynamicPassMat[prevIndex][aLength] / f32(maxLength);
  `;
}

export function webGPUNeedlemanWunsch(maxArraySize: number, entryIndex: number) {
  // version of the levenshtain where the cost of substitution is customizable
  // it is assumet that suppInfo struct contains correct matrix called similarityMatrix${entryIndex}, (similarityMatrix0, similarityMatrix1, etc)
  // and gapOpenPenalty, gapExtensionPenalty
  // this should be guaranteed by the getProcessInfo function.
  return `
      let aLength: u32 = computeInfo.entrySizes[${entryIndex}][aIndex];
      let bLength: u32 = computeInfo.entrySizes[${entryIndex}][bIndex];
      let maxLength: u32 = max(aLength, bLength);
      let minLength: u32 = min(aLength, bLength);
  
      let maxIntDistance = ceil(maxDistance * f32(maxLength));
      // we will store two arrays as matrix and swap the working indices per pass.
      // this way we can reduce memory usage per computation to just O(aLength)
      // the grid will have aLength + 1 columns and bLength + 1 rows
      // this will be guaranteed by iteration, but the array sizes must be known at compile time, so we will use a fixed size of maxArraySize
      var dynamicPassMat: array<array<f32, ${maxArraySize + 1}u>, 2>; // initialize to 0
      
      // we need to keep track of which operation led to the current cell
      // i.e. whether we came from the left, top or diagonal to assign gap open/gap extend penalty
      var verticalGaps: array<u32, ${maxArraySize + 1}u>;
      var horizontalGaps: array<u32, ${maxArraySize + 1}u>;

      let gapOpenPenalty: f32 = suppInfo.gapOpenPenalty${entryIndex};
      let gapExtensionPenalty: f32 = suppInfo.gapExtensionPenalty${entryIndex};
      var prevIndex: u32 = 0;
      var curIndex: u32 = 1; // we will swap these indices per pass
      // initialize the first row
      for (var i = 0u; i <= aLength; i = i + 1u) {
          dynamicPassMat[prevIndex][i] = gapExtensionPenalty + f32(i - 1) * gapExtensionPenalty; // accounting for the fact that left and right gaps are less costly
          dynamicPassMat[curIndex][i] = 0.0;
      }
      dynamicPassMat[0][0] = 0.0;

      let simMatrix = &suppInfo.similarityMatrix${entryIndex}; // using pointers make things faster
      // iterate over the rows
      for (var i = 1u; i <= bLength; i = i + 1u) {
          let prevRow = &dynamicPassMat[prevIndex];
          let curRow = &dynamicPassMat[curIndex];
          (*curRow)[0] = gapExtensionPenalty + f32(i - 1) * gapExtensionPenalty;
          var minEntry: f32 = f32(maxLength);
          let monB = u32(b[i - 1]);
          for (var j = 1u; j <= aLength; j = j + 1u) {
              let monA = u32(a[j - 1]);
              
              let cost: f32 = (*prevRow)[j - 1] + 1f - (*simMatrix)[monA][monB];
              var top = (*prevRow)[j]; // deletion
              if (verticalGaps[j] > 0 || i == 1 || i == bLength) {
                  top = top + gapExtensionPenalty;
              } else {
                  top = top + gapOpenPenalty;
              }
              var left = (*curRow)[j - 1]; // insertion
              if (horizontalGaps[j - 1] > 0 || j == 1 || j == aLength) {
                  left = left + gapExtensionPenalty;
              } else {
                  left = left + gapOpenPenalty;
              }
              var res: f32 = min(
                  min(
                      top, // deletion
                      left, // insertion
                  ),
                  cost // substitution
                  );
                  (*curRow)[j] = res;
              if (res < minEntry) {
                  minEntry = res;
              }
              // update the horizontal and vertical gaps
              if (res == cost) {
                  verticalGaps[j] = 0;
                  horizontalGaps[j] = 0;
              } else if (res == left) {
                  verticalGaps[j] = 0;
                  horizontalGaps[j] = 1;
              } else {
                  verticalGaps[j] = 1;
                  horizontalGaps[j] = 0;
              }
          }
          // swap the indices
          let temp: u32 = prevIndex;
          prevIndex = curIndex;
          curIndex = temp;
          if (minEntry > maxIntDistance) {
              return 1.0;
          }
      }
      return dynamicPassMat[prevIndex][aLength] / f32(minLength);

  `;
}


export function webGPUEuclidean(maxArraySize: number, _entryIndex: number) {
  return `
  var dist: f32 = 0.0;
  for (var i = 0u; i < ${maxArraySize}; i = i + 1u) {
    dist = dist + f32(a[i] - b[i]) * f32(a[i] - b[i]);
  }
  return sqrt(dist);
  `;
}

export function webGPUManhattan(maxArraySize: number, _entryIndex: number) {
  return `
  var dist: f32 = 0.0;
  for (var i = 0u; i < ${maxArraySize}; i = i + 1u) {
    dist = dist + abs(f32(a[i] - b[i]));
  }
  return dist;
  `;
}

export function webGPUOneHotDistance(_maxArraySize: number, entryIndex: number) {
  return `
  let aLength: u32 = computeInfo.entrySizes[${entryIndex}][aIndex];
  let bLength: u32 = computeInfo.entrySizes[${entryIndex}][bIndex];
  if (aLength != bLength) {
    return 1.0;
  }
  for (var i = 0u; i < aLength; i = i + 1u) {
    if(a[i] != b[i]) {
      return 1.0;
    }
  }
  return 0.0;
  `;
}

export function webGPUNumericDistance(_maxArraySize: number, entryIndex: number) {
  // we assume that range${entryIndex} is available in the supplementaryInfo struct
  return `
    let range = suppInfo.range${entryIndex};
    return f32(abs(f32(a[0]) - f32(b[0])) / range);
  `;
}


// tanimoto distance for uint32 arrays of length ${maxArraySize}
export function webGPUTanimotoBitArray(maxArraySize: number, _entryIndex: number) {
  return `
  var onBitsa: u32 = 0u;
  var onBitsb: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      onBitsa = onBitsa + countOneBits(a[i]);
      onBitsb = onBitsb + countOneBits(b[i]);
  }

  if (onBitsa == 0u && onBitsb == 0u) {
      return 0.0;
  }

  let totalOnBits = onBitsa + onBitsb;
  var commonBits: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      commonBits = commonBits + countOneBits(a[i] & b[i]);
  }

  return 1.0 - f32(commonBits) / f32(totalOnBits - commonBits);
  `;
}

export function webGPUAsymmetricBitArray(maxArraySize: number, _entryIndex: number) {
  return `
  var onBitsa: u32 = 0u;
  var onBitsb: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      onBitsa = onBitsa + countOneBits(a[i]);
      onBitsb = onBitsb + countOneBits(b[i]);
  }
  let min = min(onBitsa, onBitsb);
  if (min == 0u) {
    return 1.0;
  }
  var commonBits: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      commonBits = commonBits + countOneBits(a[i] & b[i]);
  }
  return 1.0 - f32(commonBits) / f32(min);
  `;
}

export function webGPUCosineBitArray(maxArraySize: number, _entryIndex: number) {
  return `
  var onBitsa: u32 = 0u;
  var onBitsb: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      onBitsa = onBitsa + countOneBits(a[i]);
      onBitsb = onBitsb + countOneBits(b[i]);
  }
  let total = onBitsa * onBitsb; // p.s. here total is taken by multiplying
  if (total == 0u) {
    return 1.0;
  }
  var commonBits: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      commonBits = commonBits + countOneBits(a[i] & b[i]);
  }
  return 1.0 - f32(commonBits) / sqrt(f32(total));
  `;
}

export function webGPUSokalBitArray(maxArraySize: number, _entryIndex: number) {
  return `
  var onBitsa: u32 = 0u;
  var onBitsb: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      onBitsa = onBitsa + countOneBits(a[i]);
      onBitsb = onBitsb + countOneBits(b[i]);
  }
  let total = onBitsa + onBitsb;
  if (total == 0u) {
    return 1.0;
  }
  var commonBits: u32 = 0u;
  for (var i = 0u; i < ${maxArraySize}u; i = i + 1u) {
      commonBits = commonBits + countOneBits(a[i] & b[i]);
  }
  return 1.0 - f32(commonBits) / f32(total * 2 - commonBits * 3);
  `;
}

export enum WEBGPUDISTANCE {
  HAMMING = 'Hamming',
  EUCLIDEAN = 'Euclidean',
  MANHATTAN = 'Manhattan',
  TANIMOTO = 'Tanimoto',
  LEVENSTEIN = 'Levenshtein',
  NEEDLEMAN_WUNSCH = 'Needlemann-Wunsch',
  MONOMER_CHEMICAL_DISTANCE = 'Monomer chemical distance',
  SOKAL = 'Sokal',
  COSINE = 'Cosine',
  ASYMMETRIC = 'Asymmetric',
  Difference = 'Difference',
  OneHot = 'One-Hot',
}

export const webGPUFunctions:{ [key in WEBGPUDISTANCE]: (maxArraySize: number, entryIndex: number) => string } = {
  [WEBGPUDISTANCE.HAMMING]: webGPUHamming,
  [WEBGPUDISTANCE.EUCLIDEAN]: webGPUEuclidean,
  [WEBGPUDISTANCE.MANHATTAN]: webGPUManhattan,
  [WEBGPUDISTANCE.TANIMOTO]: webGPUTanimotoBitArray,
  [WEBGPUDISTANCE.LEVENSTEIN]: webGPULevenstein,
  [WEBGPUDISTANCE.NEEDLEMAN_WUNSCH]: webGPUNeedlemanWunsch,
  [WEBGPUDISTANCE.MONOMER_CHEMICAL_DISTANCE]: webGPUMonomerChemicalDistance,
  [WEBGPUDISTANCE.SOKAL]: webGPUSokalBitArray,
  [WEBGPUDISTANCE.COSINE]: webGPUCosineBitArray,
  [WEBGPUDISTANCE.ASYMMETRIC]: webGPUAsymmetricBitArray,
  [WEBGPUDISTANCE.Difference]: webGPUNumericDistance,
  [WEBGPUDISTANCE.OneHot]: webGPUOneHotDistance
};

export const distanceFunctionComplexity: { [key in WEBGPUDISTANCE]: (maxEntrySize: number) => number } = {//
  [WEBGPUDISTANCE.HAMMING]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 30),
  [WEBGPUDISTANCE.EUCLIDEAN]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 30),
  [WEBGPUDISTANCE.MANHATTAN]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 30),
  [WEBGPUDISTANCE.TANIMOTO]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 60),
  [WEBGPUDISTANCE.SOKAL]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 60),
  [WEBGPUDISTANCE.COSINE]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 60),
  [WEBGPUDISTANCE.ASYMMETRIC]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 60),
  [WEBGPUDISTANCE.LEVENSTEIN]: (maxEntrySize: number) => Math.ceil(maxEntrySize * maxEntrySize / 60),
  [WEBGPUDISTANCE.NEEDLEMAN_WUNSCH]: (maxEntrySize: number) => Math.ceil(maxEntrySize * maxEntrySize / 60),
  [WEBGPUDISTANCE.MONOMER_CHEMICAL_DISTANCE]: (maxEntrySize: number) => Math.ceil(maxEntrySize / 25),
  [WEBGPUDISTANCE.Difference]: (_maxEntrySize: number) => 1,
  [WEBGPUDISTANCE.OneHot]: (_maxEntrySize: number) => Math.ceil(_maxEntrySize / 40),
};

export type BioDistanceFnOptions = {
  scoringMatrix?: number[][], alphabetIndexes?: {
    [monomerId: string]: number;
}, gapOpenPenalty?: number, gapExtensionPenalty?: number
}
export type SupportedEntryTypes = string | Uint32Array | Int32Array | Float32Array | number | {_data: Uint32Array, _length: number}//last is bitarray;

export const enum WGPUENTRYTYPE {
  STRING = 'STRING',
  UINT32ARRAY = 'UINT32ARRAY',
  INT32ARRAY = 'INT32ARRAY',
  FLOAT32ARRAY = 'FLOAT32ARRAY',
  NUMBER = 'NUMBER',
  BITARRAY = 'BITARRAY'
}

export const TypeSupportedDistances: {[key in WGPUENTRYTYPE] : Set<WEBGPUDISTANCE>} = {
  [WGPUENTRYTYPE.STRING]: new Set([WEBGPUDISTANCE.HAMMING, WEBGPUDISTANCE.LEVENSTEIN, WEBGPUDISTANCE.NEEDLEMAN_WUNSCH, WEBGPUDISTANCE.MONOMER_CHEMICAL_DISTANCE, WEBGPUDISTANCE.OneHot]),
  [WGPUENTRYTYPE.UINT32ARRAY]: new Set([WEBGPUDISTANCE.HAMMING, WEBGPUDISTANCE.EUCLIDEAN, WEBGPUDISTANCE.MANHATTAN, WEBGPUDISTANCE.MONOMER_CHEMICAL_DISTANCE, WEBGPUDISTANCE.LEVENSTEIN, WEBGPUDISTANCE.NEEDLEMAN_WUNSCH, WEBGPUDISTANCE.TANIMOTO, WEBGPUDISTANCE.COSINE, WEBGPUDISTANCE.SOKAL, WEBGPUDISTANCE.ASYMMETRIC, WEBGPUDISTANCE.OneHot, WEBGPUDISTANCE.Difference]),
  [WGPUENTRYTYPE.INT32ARRAY]: new Set([WEBGPUDISTANCE.EUCLIDEAN, WEBGPUDISTANCE.MANHATTAN, WEBGPUDISTANCE.OneHot, WEBGPUDISTANCE.Difference]),
  [WGPUENTRYTYPE.FLOAT32ARRAY]: new Set([WEBGPUDISTANCE.EUCLIDEAN, WEBGPUDISTANCE.MANHATTAN, WEBGPUDISTANCE.Difference]),
  [WGPUENTRYTYPE.NUMBER]: new Set([WEBGPUDISTANCE.EUCLIDEAN, WEBGPUDISTANCE.MANHATTAN, WEBGPUDISTANCE.Difference]),
  [WGPUENTRYTYPE.BITARRAY]: new Set([WEBGPUDISTANCE.TANIMOTO, WEBGPUDISTANCE.COSINE, WEBGPUDISTANCE.SOKAL, WEBGPUDISTANCE.ASYMMETRIC])
};

