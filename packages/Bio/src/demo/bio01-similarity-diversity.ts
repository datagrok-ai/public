import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {handleError} from './utils';
import {SequenceDiversityViewer} from '../analysis/sequence-diversity-viewer';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {adjustGridcolAfterRender} from '../utils/ui-utils';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';

const dataFn: string = 'samples/FASTA_PT_activity.csv';

export async function demoBioSimDiv() {
  const t = await _package.files.readCsv('samples/peptides-non-natural.csv');
  t.name = 'Similarity and Diversity Demo';
  t.col('activity')!.setTag('format', '3 significant digits');
  t.col('sequence')!.temp[MmcrTemps.maxMonomerLength] = 4;
  const tv = grok.shell.addTableView(t);
  await t.meta.detectSemanticTypes();
  await grok.data.detectSemanticTypes(t);
  const simV = tv.addViewer('Sequence Similarity Search', {limit: 20});
  const dn = tv.dockManager.dock(simV, DG.DOCK_TYPE.RIGHT, null, 'Similarity search', 0.45);
  adjustGridcolAfterRender(tv.grid, 'sequence', 500, 30);
  const divV = tv.addViewer('Sequence Diversity Search', {limit: 20});
  tv.dockManager.dock(divV, DG.DOCK_TYPE.DOWN, dn, 'Diversity search', 0.4);
  grok.functions.call('Dendrogram:HierarchicalClustering',
    {df: grok.shell.t, colNameList: ['sequence'], distance: 'euclidian', linkage: 'complete'});
}

export async function demoBio01UISteps() {
  const seqHelper: ISeqHelper = await getSeqHelper();

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
        const simViewer = new SequenceSimilarityViewer(seqHelper, true);
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
