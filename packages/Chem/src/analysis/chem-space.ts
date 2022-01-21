import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetMorganFingerprints} from '../chem-searches';
//import {getRdKitWebRoot} from '../chem-common-rdkit';
import {Coordinates} from '@datagrok-libraries/utils/src/type-declarations';
import {
  createDimensinalityReducingWorker,
} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {MetricDataTypes} from '@datagrok-libraries/ml/src/string-measure'

export async function chemSpace(table: DG.DataFrame, molColumn: DG.Column, methodName: string, similarityMetric: string) {
  let fpColumn = molColumn.toList();
  if (MetricDataTypes['BitArray'].includes(similarityMetric))
    fpColumn = await chemGetMorganFingerprints(molColumn);
  const coordinates = await createDimensinalityReducingWorker(fpColumn, methodName, similarityMetric) as Coordinates;
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

  /*if (window.Worker) {
    const myWorker = new Worker(getRdKitWebRoot() + 'src/analysis/chem-spe.js');
    const fpBuffers = new Array(fpColumn.length);
    for (let i = 0; i < fpColumn.length; ++i) {
      //@ts-ignore
      const buffer = fpColumn[i].getRawData();
      fpBuffers[i] = buffer;
    }
    myWorker.postMessage([fpColumn.length, fpBuffers,
      2, null, null, 1.0, 2.0, 0.01, fpColumn.length * 100, 100]);
    return new Promise<void>((resolve, reject) => {
      myWorker.onmessage = function(event) {
        const coordinates = event.data;
        const coords = [
          DG.Column.fromFloat32Array('SPE_X', coordinates[0]),
          DG.Column.fromFloat32Array('SPE_Y', coordinates[1]),
        ];
        table = DG.DataFrame.fromColumns(table.columns.toList().concat(coords));
        const view = grok.shell.addTableView(table);
        view.scatterPlot({
          x: 'SPE_X',
          y: 'SPE_Y',
        });
        resolve();
      };
      myWorker.onerror = function(error) {
        reject(error.message);
      };
    });
  } else
    throw new Error('Your browser doesn\'t support web workers.');*/
}
