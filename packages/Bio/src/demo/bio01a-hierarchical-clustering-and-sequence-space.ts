import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

import * as lev from 'fastest-levenshtein';
import {DistanceMatrix} from '@datagrok-libraries/ml/src/distance-matrix';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {getDendrogramService, IDendrogramService} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {demoSequenceSpace, handleError} from './utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

const dataFn = 'data/sample_FASTA_PT_activity.csv';
const seqColName = 'sequence';

export async function demoBio01aUI() {
  let treeHelper: ITreeHelper;
  let dendrogramSvc: IDendrogramService;
  let view: DG.TableView;
  let df: DG.DataFrame;
  let _spViewer: DG.ScatterPlotViewer;

  const dimRedMethod: string = 'UMAP';
  const activityColName = 'activity';

  try {
    const demoScript = new DemoScript(
      'Sequence Space',
      'Exploring sequence space of Macromolecules, comparison with hierarchical clustering results');
    await demoScript
      .step(`Load DNA sequences`, async () => {
        [df, treeHelper, dendrogramSvc] = await Promise.all([
          _package.files.readCsv(dataFn),
          getTreeHelper(),
          getDendrogramService(),
        ]);
        view = grok.shell.addTableView(df);
        view.grid.props.rowHeight = 22;
        view.grid.columns.byName('cluster')!.visible = false;
        view.grid.columns.byName('sequence')!.width = 200;
        view.grid.columns.byName('is_cliff')!.visible = false;

        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;
      }, {
        description: `Load dataset with macromolecules of 'fasta' notation, 'DNA' alphabet.`,
        delay: 2000,
      })
      .step('Build sequence space', async () => {
        _spViewer = await demoSequenceSpace(view, df, seqColName, dimRedMethod);
        _spViewer.setOptions({color: activityColName});
      }, {
        description: `Reduce sequence space dimensionality to display on 2D representation.`,
        delay: 2000,
      })
      .step('Cluster sequences', async () => {
        const seqCol: DG.Column<string> = df.getCol(seqColName);
        const seqList = seqCol.toList();
        const distance: DistanceMatrix = DistanceMatrix.calc(seqList, (aSeq: string, bSeq: string) => {
          const levDistance = lev.distance(aSeq, bSeq);
          return levDistance / ((aSeq.length + bSeq.length) / 2);
        });
        const treeRoot = await treeHelper.hierarchicalClusteringByDistance(distance, 'ward');
        dendrogramSvc.injectTreeForGrid(view.grid, treeRoot, undefined, 150, undefined);
      }, {
        description: `Perform hierarchical clustering to reveal relationships between sequences.`,
        delay: 2000,
      })
      .step('Select a sequence', async () => {
        df.selection.set(65, true);
      }, {
        description: `Handling selection of data frame row reflecting on linked viewers.`,
        delay: 2000,
      })
      .step('Select a bunch of sequences', async () => {
        [67, 72, 77, 82, 83, 84, 85, 91, 93].forEach((idx) => df.selection.set(idx, true));
        df.currentRowIdx = 27;
      }, {
        // eslint-disable-next-line max-len
        description: 'Selecting a group of rows from a data frame to show their similarity and proximity to each other on a viewer..',
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
