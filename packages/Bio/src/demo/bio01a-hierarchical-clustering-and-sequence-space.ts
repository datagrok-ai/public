import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

import * as lev from 'fastest-levenshtein';
import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {getDendrogramService, IDendrogramService} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {demoSequenceSpace, handleError} from './utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

const dataFn = 'data/sample_FASTA_DNA.csv';
const seqColName = 'sequence';

export async function demoBio01aUI() {
  let treeHelper: ITreeHelper;
  let dendrogramSvc: IDendrogramService;
  let view: DG.TableView;
  let df: DG.DataFrame;
  let spViewer: DG.ScatterPlotViewer;

  const method: string = 'UMAP';
  const idRows: { [id: number]: number } = {};
  const embedCols: { [colName: string]: DG.Column<number> } = {};

  try {
    const demoScript = new DemoScript('Demo', 'Exploring sequence space');
    await demoScript
      .step(`Loading DNA notation 'fasta'`, async () => {
        [df, treeHelper, dendrogramSvc] = await Promise.all([
          _package.files.readCsv(dataFn),
          getTreeHelper(),
          getDendrogramService()
        ]);
        view = grok.shell.addTableView(df);
        view.grid.props.rowHeight = 22;
      }, {
        description: `Load dataset with macromolecules of 'fasta' notation, 'DNA' alphabet.`,
        delay: 1600,
      })
      .step('Building sequence space', async () => {
        spViewer = await demoSequenceSpace(view, df, seqColName, method);
      }, {
        description: `Reduce sequence space dimensionality to display on 2D representation.`,
        delay: 1600
      })
      .step('Hierarchical clustering', async () => {
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
        delay: 1600,
      })
      .step('Selection', async () => {
        df.selection.init((idx: number) => [15].includes(idx));
      }, {
        description: `Handling selection of data frame row reflecting on linked viewers.`,
        delay: 1600,
      })
      .step('Select a bunch of sequences', async () => {
        df.selection.init((idx: number) => [21, 9, 58].includes(idx));
        df.currentRowIdx = 27;
      }, {
        description: 'Selecting a group of rows from a data frame to show their similarity and proximity to each other on a viewer..',
        delay: 1600,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
