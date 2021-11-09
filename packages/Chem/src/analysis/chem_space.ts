import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getMorganFingerprints, rdKitWorkerWebRoot} from '../package';

export async function chemSpace(table: DG.DataFrame, molColumn: DG.Column) {
    const fpColumn = getMorganFingerprints(molColumn);
    if (fpColumn.stats.missingValueCount > 0) {
      grok.shell.info('Molecule column has null entries which would be droped out.');
      grok.functions.call(`ClinicalCase:filterNulls`, {'df': table, 'colName': molColumn.name});
    }
  
    if (window.Worker) {
      const myWorker = new Worker(rdKitWorkerWebRoot + 'src/analysis/chem_stochastic_proximity_embedding.js');
      const fpBuffers = new Array(fpColumn.length);
  
      for (let i = 0; i < fpColumn.length; ++i) {
        const buffer = fpColumn.get(i).getBuffer();
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
    } else {
      throw new Error('Your browser doesn\'t support web workers.');
    }
  }