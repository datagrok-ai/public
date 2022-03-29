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
  similarityMetric: string): Promise<DG.ColumnList> {
  const fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
  const coordinates: Coordinates =
    await createDimensinalityReducingWorker(
      {data: fpColumn as BitArray[], metric: similarityMetric as BitArrayMetrics}, methodName) as Coordinates;
  const axes = ['Embed_X', 'Embed_Y'];
  const cols: DG.Column[] = [];

  coordinates[0] = normalize(coordinates[0]);
  for (let i = 0; i < axes.length; ++i) {
    const name = axes[i];
    cols[i] = (DG.Column.fromFloat32Array(name, coordinates[i]));
  }
  return new DG.ColumnList(cols);
}
