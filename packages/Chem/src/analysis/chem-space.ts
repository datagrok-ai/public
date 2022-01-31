import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
//import {getRdKitWebRoot} from '../chem-common-rdkit';
import {Coordinates, Vector} from '@datagrok-libraries/utils/src/type-declarations';
import {createDimensinalityReducingWorker,} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {BitArrayMetrics, isBitArrayMetric, KnownMetrics, StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Fingerprint} from "../utils/chem-common";
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export async function chemSpace(molColumn: DG.Column, methodName: string, similarityMetric: string): Promise<DG.ColumnList> {

  function normalize(data: Vector) {
    let mean = 0;
    let std = 0;

    data.forEach((val) => {
      mean += val;
    });
    mean /= data.length;

    data.forEach((val) => {
      std += (val - mean) * (val - mean);
    });
    std = Math.sqrt(std / data.length);

    data.forEach((val, idx) => {
      data[idx] = (val - mean) / std;
    });
  }
  
  const fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
  const coordinates: Coordinates = await createDimensinalityReducingWorker({data: fpColumn as BitArray[], metric: similarityMetric as BitArrayMetrics}, methodName) as Coordinates;
  const axes = ['Embed_X', 'Embed_Y'];
  const cols: DG.Column[] = [];
  
  normalize(coordinates[0]);
  for (let i = 0; i < axes.length; ++i) {
    const name = axes[i];
    cols[i] = (DG.Column.fromFloat32Array(name, coordinates[i]));
  }
  return new DG.ColumnList(cols);
}
