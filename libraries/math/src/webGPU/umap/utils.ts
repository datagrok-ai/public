import {levenbergMarquardt} from 'ml-levenberg-marquardt';
import {SparseKNNInfo} from './types';

export function linear(x1: number, x2: number, n: number) {
  const step = (x2 - x1) / (n - 1);
  return Array.from({length: n}, (_, i) => x1 + step * i);
}

/**
 * Fit a, b params for the differentiable curve used in lower
 * dimensional fuzzy simplicial complex construction. We want the
 * smooth curve (from a pre-defined family with simple gradient) that
 * best matches an offset exponential decay.
 */
export function findABParams(spread: number, minDist: number) {
  const curve = ([a, b]: number[]) => (x: number) => {
    return 1.0 / (1.0 + a * x ** (2 * b));
  };

  const xv = linear(0, spread * 3, 300)
    .map((val) => (val < minDist ? 1.0 : val));

  const yv = new Array<number>(xv.length).fill(0).map((val, index) => {
    const gte = xv[index] >= minDist;
    return gte ? Math.exp(-(xv[index] - minDist) / spread) : val;
  });

  const initialValues = [0.5, 0.5];
  const data = {x: xv, y: yv};

  // Default options for the algorithm (from github example)
  const options = {
    damping: 1.5,
    initialValues,
    gradientDifference: 10e-2,
    maxIterations: 100,
    errorTolerance: 10e-3,
  };

  // eslint-disable-next-line new-cap
  const {parameterValues} = levenbergMarquardt(data, curve, options);
  const [a, b] = parameterValues as number[];
  return {a, b};
}

/**
 * heuristic for number of epochs of umap, differs from cpu version, because of the way gpu passes data
 * @param entryLen number of entries
 * @returns
 */
export function getNEpochs(entryLen: number) {
  //
  const length = entryLen;
  if (length <= 2500)
    return 1000;
  else if (length <= 5000)
    return 800;
  else if (length <= 7500)
    return 650;
  else
    return 400;
}

/**  helper function that transforms array of sizes to offset form */
export function toOffsetForm(ar: Uint32Array | number[] | Int32Array) {
  const res = new Uint32Array(ar.length + 1);
  let offset = 0;
  for (let i = 0; i < ar.length; i++) {
    res[i] = offset;
    offset += ar[i];
  }
  res[ar.length] = offset;
  return res;
}

/** multiplies flat representation of sparse knn info by scalar.
 * also no point of doing this on gpu, can be done on cpu pretty quickly.
 * if there will be need for it, it can be done on gpu as well very very easily.
 * @param values values array
 * @param scalar
 */
export function multiplyScalar(values: Int32Array | number[] | Uint32Array | Float32Array, scalar: number) {
  for (let i = 0; i < values.length; i++)
    values[i] = values[i] * scalar;
}
/**
 * generate transpose of knn from the transposition info taken from sparse knn info
 * using this information avoids billions of array.push calls, is much cleaner, allows typed array use and is faster.
 * no point of doing this on gpu, can be done on cpu as we already have all the info needed.
 * @param knnIndexes
 * @param knnDistances
 * @param knnSparseInfo
 * @param knnSize
 * @returns
 */
export async function transposeKNN(knnIndexes: number[][] | Int32Array[] | Uint32Array[],
  knnDistances: number[][] | Float32Array[] | Uint32Array[],
  knnSparseInfo: SparseKNNInfo, knnSize = 15) {
  const numOfEntries = knnIndexes.length;
  const fullKnnEntrySize = knnSize * numOfEntries;
  const transposeSizes = knnSparseInfo.resultTransposeSizesArray;
  const transposeOffsets = toOffsetForm(transposeSizes);
  const offsetsCopy = new Uint32Array(transposeOffsets.length);
  offsetsCopy.set(transposeOffsets); // will be used for dynamic indexing
  const transposeKNNIndexes = new Int32Array(fullKnnEntrySize).fill(0);
  const transposeKNNDistances = new Float32Array(fullKnnEntrySize).fill(0);

  for (let i = 0; i < numOfEntries; i++) {
    for (let j = 0; j < knnSize; j++) {
      const otherIndex = knnIndexes[i][j];
      const otherIndexOffset = offsetsCopy[otherIndex];
      transposeKNNIndexes[otherIndexOffset] = i;
      transposeKNNDistances[otherIndexOffset] = knnDistances[i][j];
      offsetsCopy[otherIndex] += 1;
    }
  }
  return {transposeKNNIndexes, transposeKNNDistances, transposeOffsets};
}

