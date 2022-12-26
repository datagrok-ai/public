import * as DG from 'datagrok-api/dg';

import {getSimilarityFromDistance} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import * as chemSearches from '../chem-searches';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';

export async function getSimilaritiesMarix(dim: number, smiles: DG.Column, dfSmiles: DG.DataFrame, colName: string, simArr: DG.Column[])
  : Promise<DG.Column[]> {
  for (let i = 0; i != dim - 1; ++i) {
    const mol = smiles.get(i);
    dfSmiles.rows.removeAt(0, 1, false);
    simArr[i] = (await chemSearches.chemGetSimilarities(dfSmiles.col(colName)!, mol))!;
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
