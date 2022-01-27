import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
//import {getRdKitWebRoot} from '../chem-common-rdkit';
import {Coordinates} from '@datagrok-libraries/utils/src/type-declarations';
import {createDimensinalityReducingWorker,} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {BitArrayMetrics, isBitArrayMetric, KnownMetrics, StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Fingerprint} from "../utils/chem-common";
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export async function chemSpace(table: DG.DataFrame, molColumn: DG.Column, methodName: string, similarityMetric: string) {
  let fpColumn = molColumn.toList();
  let coordinates: Coordinates;
  if (isBitArrayMetric(similarityMetric as KnownMetrics)) {
    fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
    coordinates = await createDimensinalityReducingWorker({data: fpColumn as BitArray[], metric: similarityMetric as BitArrayMetrics}, methodName) as Coordinates;
  } else {
    coordinates = await createDimensinalityReducingWorker({data: fpColumn as string[], metric: similarityMetric as StringMetrics}, methodName) as Coordinates;
  }
  const axes = ['Embed_X', 'Embed_Y'];

  for (let i = 0; i < axes.length; ++i) {
    const name = axes[i];
    const cols = (table.columns as DG.ColumnList);

    if (table.col(name))
      cols.remove(name);

    cols.insert(DG.Column.fromFloat32Array(name, coordinates[i]));
  }

  const view = grok.shell.addTableView(table);
  view.scatterPlot({x: axes[0], y: axes[1]});
}
