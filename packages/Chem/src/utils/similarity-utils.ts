import * as DG from 'datagrok-api/dg';

import {getSimilarityFromDistance} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import * as chemSearches from '../chem-searches';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import { _chemGetSimilarities, chemGetFingerprints } from '../chem-searches';

export async function getSimilaritiesMarix(dim: number, smiles: DG.Column, dfSmiles: DG.DataFrame,
  colName: string, simArr: (DG.Column | null)[]): Promise<(DG.Column | null)[]> {
  let fingerprints = await chemGetFingerprints(dfSmiles.col(colName)!);
  for (let i = 0; i != dim - 1; ++i) {
    fingerprints.shift();
    const queryMolString = smiles.get(i)!;
    simArr[i] = queryMolString.length != 0 ?
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
      _chemGetSimilarities(queryMolString, fingerprints)) : null;
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
