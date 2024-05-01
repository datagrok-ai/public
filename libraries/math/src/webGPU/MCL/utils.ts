import {SparseMatrix} from '../types';
import {toOffsetForm} from '../umap/utils';
import {MCLOpReturnType} from './types';

export function countOccurancesAddLoops(sparseMatrix: SparseMatrix, nRows: number) {
  const res = new Uint32Array(nRows);
  for (let i = 0; i < sparseMatrix.i.length; i++) {
    res[sparseMatrix.i[i]]++;
    res[sparseMatrix.j[i]]++;
  }
  // add loops
  for (let i = 0; i < nRows; i++)
    res[i]++;

  return res;
}

/**
 * BE WARE, this function is ment for half sparse matrices only!!!
 * @param sparseMatrix
 * @param nRows
 * @returns
 */
export function toSparseKNNSimilarityForm(sparseMatrix: SparseMatrix, nRows: number): MCLOpReturnType {
  const occurenceCount = countOccurancesAddLoops(sparseMatrix, nRows);
  const offsetForm = toOffsetForm(occurenceCount);
  const insertCounter = new Uint32Array(nRows).fill(1); // for counting how many things we inserted into the knn arrays.
  // we initialize these counters to one because first we insert the self loops
  // dont be confused by the name knn, k here is variable, its just the thing that in this form, we know that at index i
  // we have some k(variable) ammount of sparse neigbours
  // we know for a fact that sparse matrix will not contain any duplicates, we add nRows to account for self loops
  const KNNIndexes = new Uint32Array(sparseMatrix.i.length * 2 + nRows);
  const KNNSimilarities = new Float32Array(sparseMatrix.i.length * 2 + nRows);

  // insert self loops at the start of each offset
  for (let i = 0; i < nRows; i++) {
    KNNIndexes[offsetForm[i]] = i;
    KNNSimilarities[offsetForm[i]] = 1;
  }
  for (let i = 0; i < sparseMatrix.i.length; i++) {
    const row = sparseMatrix.i[i];
    const col = sparseMatrix.j[i];
    const similarity = 1 - sparseMatrix.distance[i];
    KNNIndexes[offsetForm[row] + insertCounter[row]] = col;
    KNNSimilarities[offsetForm[row] + insertCounter[row]] = similarity;
    insertCounter[row]++;
    KNNIndexes[offsetForm[col] + insertCounter[col]] = row;
    KNNSimilarities[offsetForm[col] + insertCounter[col]] = similarity;
    insertCounter[col]++;
  }
  return {KNNIndexes, KNNSimilarities, indexOffsets: offsetForm};
}


export function knnToSparseForm(knnIndexes: number[][], knnSims: number[][]): MCLOpReturnType {
  const lengths = knnIndexes.map((x) => x.length);
  const totalLength = lengths.reduce((a, b) => a + b, 0);
  const offsets = toOffsetForm(lengths);
  const outIndexes = new Uint32Array(totalLength);
  const outSims = new Float32Array(totalLength);
  for (let i = 0; i < knnIndexes.length; i++) {
    outIndexes.set(knnIndexes[i], offsets[i]);
    outSims.set(knnSims[i], offsets[i]);
  }
  return {KNNIndexes: outIndexes, KNNSimilarities: outSims, indexOffsets: offsets};
}
