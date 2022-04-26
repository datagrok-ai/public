import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
//import {getRdKitWebRoot} from '../chem-common-rdkit';
import {Coordinates} from '@datagrok-libraries/utils/src/type-declarations';
import {createDimensinalityReducingWorker} from
  '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {normalize} from '@datagrok-libraries/utils/src/vector-operations';
import {Fingerprint} from '../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export async function chemSpace(
  molColumn: DG.Column,
  methodName: string,
  similarityMetric: string,
  axes: string[],
  returnDistances?: boolean): Promise<any> {
  const fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
  const {distance: distance, embedding: coordinates}: any =
    await createDimensinalityReducingWorker(
      {data: fpColumn as BitArray[], metric: similarityMetric as BitArrayMetrics},
      methodName, undefined, returnDistances);
  const cols: DG.Column[] = [];

  coordinates[0] = normalize(coordinates[0]);
  //coordinates[1] = normalize(coordinates[1]);
  for (let i = 0; i < axes.length; ++i) {
    const name = axes[i];
    cols[i] = (DG.Column.fromFloat32Array(name, coordinates[i]));
  }
  return returnDistances ? {distance: distance, coordinates: new DG.ColumnList(cols)} : new DG.ColumnList(cols);
}
