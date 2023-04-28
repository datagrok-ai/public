import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {delay} from '@datagrok-libraries/utils/src/test';
import {step} from './utils';

const dataFn = 'data/sample_FASTA_DNA.csv';

export async function demoBio01UI(funcPath: string) {
  let view: DG.TableView;
  let df: DG.DataFrame;

  try {
    await step(`Loading DNA notation 'fasta'.`, async () => {
      df = await _package.files.readCsv(dataFn);
      view = grok.shell.addTableView(df);
      view.path = view.basePath = funcPath;
    })();

    await step('Sequence similarity search.', async () => {
      const simViewer = await df.plot.fromType('Sequence Similarity Search') as DG.Viewer;
      view.dockManager.dock(simViewer, DG.DOCK_TYPE.RIGHT, null, 'Similarity search', 0.35);
    })();

    await step('Sequence diversity search.', async () => {
      const divViewer = await df.plot.fromType('Sequence Diversity Search') as DG.Viewer;
      view.dockManager.dock(divViewer, DG.DOCK_TYPE.DOWN, null, 'Diversity search', 0.27);
    })();

    await step('Current row 3.', async () => {
      df.currentRowIdx = 3;
    })();

    await step('Current row  7', async () => {
      df.currentRowIdx = 7;
    });
  } catch (err: any) {
    if (err instanceof Error)
      _package.logger.error(err.message, undefined, err.stack);
    else
      _package.logger.error(err.toString());
  }
}
