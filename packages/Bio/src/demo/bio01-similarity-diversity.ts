import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {handleError} from './utils';
import {SequenceDiversityViewer} from '../analysis/sequence-diversity-viewer';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';

const dataFn: string = 'samples/FASTA_PT_activity.csv';

export async function demoBio01UI() {
  let view: DG.TableView;
  let df: DG.DataFrame;

  try {
    const demoScript = new DemoScript('Similarity, Diversity',
      'Sequence similarity tracking and evaluation dataset diversity', false, {autoStartFirstStep: true});
    await demoScript
      .step(`Load DNA sequences`, async () => {
        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;

        df = await _package.files.readCsv(dataFn);
        view = grok.shell.addTableView(df);

        view.grid.columns.byName('cluster')!.visible = false;
        view.grid.columns.byName('sequence_id')!.visible = false;
        view.grid.columns.byName('sequence')!.width = 300;
        view.grid.columns.byName('activity')!.visible = false;
        view.grid.columns.byName('is_cliff')!.visible = false;

        // TODO: Fix column width
      }, {
        description: `Load dataset with macromolecules of 'fasta' notation, 'DNA' alphabet.`,
        delay: 2000,
      })
      .step('Find the most similar sequences to the current one', async () => {
        const simViewer = new SequenceSimilarityViewer(true);
        view.addViewer(simViewer, {
          moleculeColumnName: 'sequence',
          similarColumnLabel: 'Similar to current',
        });
        view.dockManager.dock(simViewer, DG.DOCK_TYPE.RIGHT, null, 'Similarity search', 0.35);
      }, {
        description: `Add 'Sequence Similarity Search' viewer.`,
        delay: 2000,
      })
      .step('Explore most diverse sequences in a dataset', async () => {
        const divViewer = await df.plot.fromType('Sequence Diversity Search', {
          moleculeColumnName: 'sequence',
          diverseColumnLabel: 'Top diverse sequences of all data',
        }) as SequenceDiversityViewer;
        view.dockManager.dock(divViewer, DG.DOCK_TYPE.DOWN, null, 'Diversity search', 0.27);
      }, {
        description: `Add 'Sequence Deversity Search' viewer.`,
        delay: 2000,
      })
      .step('Choose another sequence for similarity search', async () => {
        df.currentRowIdx = 3;
      }, {
        description: 'Handling current row changed of data frame showing update of similar sequences.',
        delay: 2000,
      })
      .step('One more sequence for similarity search', async () => {
        df.currentRowIdx = 7;
      }, {
        description: 'Just one more sequence to search similar ones.',
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
