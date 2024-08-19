import {SparseMatrix} from '@datagrok-libraries/math';
import {multSparseMatrix} from './mcl-sparse-matrix-mult';

export function prune(sparseMatrix: SparseMatrix, pruneValue: number): SparseMatrix {
  const is = [];
  const js = [];
  const ds = [];
  for (let i = 0; i < sparseMatrix.i.length; i++) {
    if (sparseMatrix.distance[i] > pruneValue) {
      is.push(sparseMatrix.i[i]);
      js.push(sparseMatrix.j[i]);
      ds.push(sparseMatrix.distance[i]);
    }
  }
  return {i: new Uint32Array(is), j: new Uint32Array(js), distance: new Float32Array(ds)};
}

export function mclInflate(sparseMatrix: SparseMatrix, pow: number) {
  const d = sparseMatrix.distance;
  for (let i = 0; i < sparseMatrix.i.length; i++)
    d[i] = Math.pow(d[i], pow);

  return sparseMatrix;
}

export function colwiseNormilize(
  sparseMatrix: SparseMatrix, nRows: number, pruneValue: number, doPrune?: boolean
): SparseMatrix {
  if (doPrune)
    sparseMatrix = prune(sparseMatrix, pruneValue);

  const colSums = new Float32Array(nRows);
  const d = sparseMatrix.distance;
  const js = sparseMatrix.j;
  for (let i = 0; i < sparseMatrix.i.length; i++)
    colSums[js[i]] += d[i];

  for (let i = 0; i < sparseMatrix.i.length; i++) {
    const sum = colSums[js[i]];
    if (sum !== 0)
      d[i] /= sum;
  }
  return sparseMatrix;
}

export function initMCLMatrix(sparseMatrix: SparseMatrix, nRows: number, pruneValue: number): SparseMatrix {
  // here we get sparse matrix that is only a half of the matrix
  // also, the matrix contains distances, need to convert those into similarities
  // also, we need to add self loops
  const is = new Uint32Array(sparseMatrix.i.length * 2 + nRows);
  const js = new Uint32Array(sparseMatrix.j.length * 2 + nRows);
  const ds = new Float32Array(sparseMatrix.distance.length * 2 + nRows);
  for (let i = 0; i < nRows; i++) {
    is[i] = i;
    js[i] = i;
    ds[i] = 1;
  }
  for (let i = 0; i < sparseMatrix.i.length; i++) {
    const start = i * 2 + nRows;
    is[start] = sparseMatrix.i[i];
    js[start] = sparseMatrix.j[i];
    ds[start] = 1 - sparseMatrix.distance[i];

    is[start + 1] = sparseMatrix.j[i];
    js[start + 1] = sparseMatrix.i[i];
    ds[start + 1] = 1 - sparseMatrix.distance[i];
  }

  // normalize
  colwiseNormilize({i: is, j: js, distance: ds}, nRows, pruneValue, false);

  return {i: is, j: js, distance: ds};
}

export async function runMarkovClustering(
  sparseMatrix: SparseMatrix, nRows: number, inflate: number, maxIter: number, pruneValue: number
) {
  let matrix = initMCLMatrix(sparseMatrix, nRows, pruneValue);
  for (let i = 0; i < maxIter; i++) {
    matrix = await multSparseMatrix(matrix, nRows, pruneValue);
    // console.time('inflate');
    //@ts-ignore
    matrix = colwiseNormilize(matrix, nRows, pruneValue, false);
    matrix = mclInflate(matrix, inflate);

    //@ts-ignore
    matrix = colwiseNormilize(matrix, nRows, pruneValue, true);
    console.log('MCL iteration', i);
  }

  return matrix;
}
