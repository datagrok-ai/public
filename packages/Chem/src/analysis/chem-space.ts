import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
import {createDimensinalityReducingWorker, IReduceDimensionalityResult} from
  '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {normalize} from '@datagrok-libraries/utils/src/vector-operations';
import {Fingerprint} from '../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { Matrix } from '@datagrok-libraries/utils/src/type-declarations';

// export interface IChemSpaceResults {
//
// }
export interface IChemSpaceResult {
  distance: Matrix;
  coordinates: DG.ColumnList;
}

export async function chemSpace(molColumn: DG.Column, methodName: string, similarityMetric: string,
        axes: string[], options?: any): Promise<IChemSpaceResult> {

  const fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
  const dimensionalityReduceRes: IReduceDimensionalityResult =
    await createDimensinalityReducingWorker(
      // @ts-ignore
      {data: fpColumn, metric: similarityMetric as BitArrayMetrics},
      methodName, options);

  const cols: DG.Column[] = [];

  dimensionalityReduceRes.embedding[0] = normalize(dimensionalityReduceRes.embedding[0]);
  //coordinates[1] = normalize(coordinates[1]);
  for (let i = 0; i < axes.length; ++i) {
    const name = axes[i];
    cols[i] = (DG.Column.fromFloat32Array(name, dimensionalityReduceRes.embedding[i]));
  }
  return {distance: dimensionalityReduceRes.distance, coordinates: new DG.ColumnList(cols)};
}

export function getEmbeddingColsNames(df: DG.DataFrame){
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
