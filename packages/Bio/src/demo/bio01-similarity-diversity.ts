import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {delay} from '@datagrok-libraries/utils/src/test';
import {handleError} from './utils';

const dataFn = 'data/sample_FASTA_DNA.csv';

export async function demoBio01UI() {
  let view: DG.TableView;
  let df: DG.DataFrame;

  try {
    const demoScript = new DemoScript('Demo', 'Sequence similarity / diversity search');
    await demoScript
      .step(`Loading DNA notation 'fasta'`, async () => {
        df = await _package.files.readCsv(dataFn);
        view = grok.shell.addTableView(df);
      }, {
        description: `Load dataset with macromolecules of 'fasta' notation, 'DNA' alphabet.`,
        delay: 1600
      })
      .step('Sequence similarity search', async () => {
        const simViewer = await df.plot.fromType('Sequence Similarity Search') as DG.Viewer;
        view.dockManager.dock(simViewer, DG.DOCK_TYPE.RIGHT, null, 'Similarity search', 0.35);
      }, {
        description: `Add 'Sequence Similarity Search' viewer.`,
        delay: 1600
      })
      .step('Sequence diversity search', async () => {
        const divViewer = await df.plot.fromType('Sequence Diversity Search') as DG.Viewer;
        view.dockManager.dock(divViewer, DG.DOCK_TYPE.DOWN, null, 'Diversity search', 0.27);
      }, {
        description: `Add 'Sequence Deversity Search' viewer.`,
        delay: 1600
      })
      .step('Set current row 3', async () => {
        df.currentRowIdx = 3;
      }, {
        description: 'Handling current row changed of data frame showing update of similar sequences.',
        delay: 1600,
      })
      .step('Set current row 7', async () => {
        df.currentRowIdx = 7;
      }, {
        description: 'Changing current row to another.',
        delay: 1600,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
