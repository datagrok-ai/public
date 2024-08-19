import {getKnnSparseOpInfo} from './knn-sparse-info';
import {pairwiseOpSparse} from './pairwise-sparse-ops';
import {multiplyScalar, toOffsetForm, transposeKNN} from './utils';

export async function performMatrixOps(
  knnIndexes: number[][] | Int32Array[] | Uint32Array[],
  knnDistances: number[][] | Float32Array[] | Uint32Array[],
  setOpMixRatio = 1.0) {
  const flatKnnIndexes = new Int32Array(knnIndexes.length * knnIndexes[0].length);
  const flatKnnDistances = new Float32Array(knnDistances.length * knnDistances[0].length);
  for (let i = 0; i < knnIndexes.length; i++) {
    flatKnnIndexes.set(knnIndexes[i], i * knnIndexes[i].length);
    flatKnnDistances.set(knnDistances[i], i * knnDistances[i].length);
  }
  const unionInfo = await getKnnSparseOpInfo(knnIndexes);
  const numOfEntries = knnIndexes.length;
  const transposed = await transposeKNN(knnIndexes, knnDistances, unionInfo!, knnDistances[0].length);
  const sourceOffsets = toOffsetForm(new Uint32Array(numOfEntries).fill(knnIndexes[0].length));
  const unionMatrixOffsets = toOffsetForm(unionInfo!.resultUnionSizesArray);

  const transposeMultRes = await pairwiseOpSparse(flatKnnIndexes, flatKnnDistances, sourceOffsets,
    transposed.transposeKNNIndexes, transposed.transposeKNNDistances,
    transposed.transposeOffsets, numOfEntries, unionMatrixOffsets, '*');
  const addRes = await pairwiseOpSparse(flatKnnIndexes, flatKnnDistances, sourceOffsets,
    transposed.transposeKNNIndexes, transposed.transposeKNNDistances,
    transposed.transposeOffsets, numOfEntries, unionMatrixOffsets, '+');
    // these offsets will be same as union matrix offsets
  const subtRes = await pairwiseOpSparse(addRes!.resultKnnIndexes!, addRes!.resultKnnDistances, unionMatrixOffsets,
        transposeMultRes!.resultKnnIndexes, transposeMultRes!.resultKnnDistances,
        unionMatrixOffsets, numOfEntries, unionMatrixOffsets, '-');

  if (setOpMixRatio !== 1.0) {
    multiplyScalar(subtRes!.resultKnnDistances!, setOpMixRatio);
    multiplyScalar(transposeMultRes!.resultKnnDistances, 1.0 - setOpMixRatio);
    const res = await pairwiseOpSparse(subtRes!.resultKnnIndexes!, subtRes!.resultKnnDistances, unionMatrixOffsets,
            transposeMultRes!.resultKnnIndexes, transposeMultRes!.resultKnnDistances, unionMatrixOffsets,
            numOfEntries, unionMatrixOffsets, '+');
    return {res, unionMatrixOffsets};
  }
  return {res: subtRes, unionMatrixOffsets, unionSizes: unionInfo!.resultUnionSizesArray};
}
