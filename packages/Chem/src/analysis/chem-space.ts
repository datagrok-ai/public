import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
import {createDimensinalityReducingWorker} from
  '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {normalize} from '@datagrok-libraries/utils/src/vector-operations';
import {Fingerprint} from '../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

// export interface IChemSpaceResults {
//
// }

export async function chemSpace(molColumn: DG.Column, methodName: string, similarityMetric: string,
        axes: string[], options?: any, returnDistances?: boolean): Promise<any> {

  const fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
  let distance;
  let coordinates;
  const dimensionalityReduceRes: any =
    await createDimensinalityReducingWorker(
      // @ts-ignore
      {data: fpColumn, metric: similarityMetric as BitArrayMetrics},
      methodName, options, returnDistances);
  if (returnDistances) {
    distance = dimensionalityReduceRes.distance;
    coordinates = dimensionalityReduceRes.embedding;
  } else
    coordinates = dimensionalityReduceRes;

  const cols: DG.Column[] = [];

  coordinates[0] = normalize(coordinates[0]);
  //coordinates[1] = normalize(coordinates[1]);
  for (let i = 0; i < axes.length; ++i) {
    const name = axes[i];
    cols[i] = (DG.Column.fromFloat32Array(name, coordinates[i]));
  }
  return returnDistances ? {distance: distance, coordinates: new DG.ColumnList(cols)} : new DG.ColumnList(cols);
}
