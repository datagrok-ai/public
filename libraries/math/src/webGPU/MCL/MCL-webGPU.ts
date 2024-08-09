/* eslint-disable max-len */
import {getGPUDevice} from '../getGPUDevice';
import {SparseMatrix} from '../types';
import {sparseKNNNoralizeColwise} from './colwise-normalize';
import {expandNoRevive, inflate} from './inflate-expand';
import {MCLOpReturnType} from './types';
import {toSparseKNNSimilarityForm} from './utils';

export async function markovClusterWebGPU(
  sparseMatrix: SparseMatrix, nRows: number, maxIterations: number = 5, inflateFactor: number = 2
): Promise<MCLOpReturnType> {
  const device = await getGPUDevice();
  if (!device)
    throw new Error('no gpu device found');

  // self loops are already added here
  const sparseKNNForm = toSparseKNNSimilarityForm(sparseMatrix, nRows);
  // first, we normilize the similarities

  await sparseKNNNoralizeColwise(device, sparseKNNForm.KNNSimilarities, sparseKNNForm.indexOffsets, nRows);

  let res: MCLOpReturnType = sparseKNNForm;
  // then we expand the similarities
  for (let i = 0; i < maxIterations; i++) {
    //console.log(checkSorted(res.KNNIndexes, res.indexOffsets));

    const expandRes = await expandNoRevive(device, res.KNNSimilarities, res.KNNIndexes, res.indexOffsets, nRows);

    await sparseKNNNoralizeColwise(device, expandRes.KNNSimilarities, expandRes.indexOffsets, nRows);
    // then we inflate the similarities
    inflate(expandRes.KNNSimilarities, inflateFactor);

    // then we normilize the similarities again
    await sparseKNNNoralizeColwise(device, expandRes.KNNSimilarities, expandRes.indexOffsets, nRows);
    res = expandRes;
  }
  return res;
}

