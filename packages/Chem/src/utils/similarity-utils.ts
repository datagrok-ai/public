import * as DG from 'datagrok-api/dg';

import {getSimilarityFromDistance} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {_chemGetSimilarities, chemGetFingerprints} from '../chem-searches';
import {Fingerprint} from './chem-common';

export async function getSimilaritiesMarix(dim: number, smiles: DG.Column, dfSmiles: DG.DataFrame,
  colName: string, simArr: (DG.Column | null)[]): Promise<(DG.Column | null)[]> {
  const fingerprints = await chemGetFingerprints(dfSmiles.col(colName)!, Fingerprint.Morgan, true, false);
  for (let i = 0; i != dim - 1; ++i) {
    const fp = fingerprints.shift();
    simArr[i] = fp ?
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
        _chemGetSimilarities(smiles.get(i)!, fingerprints)) : null;
  }
  return simArr;
}

export function getSimilaritiesMarixFromDistances(dim: number, distances: Matrix, simArr: DG.Column[])
  : DG.Column[] {
  for (let i = 0; i < dim - 1; ++i) {
    const similarityArr = new Float32Array(dim - i - 1).fill(0);
    for (let j = i + 1; j < dim; ++j)
      similarityArr[j - i - 1] = getSimilarityFromDistance(distances[i][j]);
    simArr[i] = DG.Column.fromFloat32Array('similarity', similarityArr);
  }
  return simArr;
}
