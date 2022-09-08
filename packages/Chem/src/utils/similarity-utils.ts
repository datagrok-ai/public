import {getSimilarityFromDistance} from '@datagrok-libraries/utils/src/similarity-metrics';
import * as chemSearches from '../chem-searches';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import * as DG from 'datagrok-api/dg';

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
    const similarityArr = [];
    for (let j = i + 1; j < dim; ++j)
      similarityArr.push(getSimilarityFromDistance(distances[i][j]));
    simArr[i] = DG.Column.fromFloat32Array('similarity', Float32Array.from(similarityArr));
  }
  return simArr;
}
