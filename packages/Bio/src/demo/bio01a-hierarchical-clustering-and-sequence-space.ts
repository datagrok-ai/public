import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package, sequenceSpaceTopMenu} from '../package';
import {delay} from '@datagrok-libraries/utils/src/test';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import * as lev from 'fastest-levenshtein';
import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {getDendrogramService, IDendrogramService} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {demoSequenceSpace, step} from './utils';

const dataFn = 'data/sample_FASTA_DNA.csv';


export async function demoBio01aUI(funcPath: string) {
  let treeHelper: ITreeHelper;
  let dendrogramSvc: IDendrogramService;
  let view: DG.TableView;
  let df: DG.DataFrame;
  let spViewer: DG.ScatterPlotViewer;

  const method: string = 'UMAP';
  const idRows: { [id: number]: number } = {};
  const embedCols: { [colName: string]: DG.Column<number> } = {};

  try {
    await step(`Loading DNA notation 'fasta'.`, async () => {
      [df, treeHelper, dendrogramSvc] = await Promise.all([
        _package.files.readCsv(dataFn),
        getTreeHelper(),
        getDendrogramService()
      ]);
      view = grok.shell.addTableView(df);
      view.grid.props.rowHeight = 22;
      view.path = view.basePath = funcPath;
    })();

    await step('Building sequence space.', async () => {
      spViewer = await demoSequenceSpace(view, df, method);
    })();

    await step('Hierarchical clustering.', async () => {
      const seqCol: DG.Column<string> = df.getCol('sequence');
      const seqList = seqCol.toList();
      const distance: DistanceMatrix = DistanceMatrix.calc(seqList, (aSeq: string, bSeq: string) => {
        const levDistance = lev.distance(aSeq, bSeq);
        return levDistance / ((aSeq.length + bSeq.length) / 2);
      });
      const treeRoot = await treeHelper.hierarchicalClusteringByDistance(distance, 'ward');
      dendrogramSvc.injectTreeForGrid(view.grid, treeRoot, undefined, 150, undefined);
    })();

    await step('Selection.', async () => {
      df.selection.init((idx: number) => [15].includes(idx));
    })();

    await step('Select bunch of sequences.', async () => {
      df.selection.init((idx: number) => [21, 9, 58].includes(idx));
      df.currentRowIdx = 27;
    })();
  } catch (err: any) {
    if (err instanceof Error)
      _package.logger.error(err.message, undefined, err.stack);
    else
      _package.logger.error(err.toString());
  }
}
