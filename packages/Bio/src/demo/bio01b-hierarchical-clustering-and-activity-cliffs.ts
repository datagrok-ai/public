import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package, activityCliffs, sequenceSpaceTopMenu} from '../package';
import wu from 'wu';
import $ from 'cash-dom';

import {TAGS as acTAGS, TEMPS as acTEMPS} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {delay} from '@datagrok-libraries/utils/src/test';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import * as lev from 'fastest-levenshtein';
import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {getDendrogramService, IDendrogramService} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {step} from './utils';

const dataFn = 'samples/sample_FASTA.csv';

export async function demoBio01bUI(funcPath: string) {
  let treeHelper: ITreeHelper;
  let dendrogramSvc: IDendrogramService;
  let view: DG.TableView;
  let df: DG.DataFrame;
  let activityCliffsViewer: DG.ScatterPlotViewer;

  const method: string = 'UMAP';
  const idRows: { [id: number]: number } = {};

  try {
    await step('Loading DNA notation \'fasta\'.', async () => {
      [df, treeHelper, dendrogramSvc] = await Promise.all([
        _package.files.readCsv(dataFn),
        getTreeHelper(),
        getDendrogramService()
      ]);

      view = grok.shell.addTableView(df);
      view.path = view.basePath = funcPath;
      view.grid.props.rowHeight = 22;
      const uniProtKbGCol = view.grid.columns.byName('UniProtKB')!;
      uniProtKbGCol.width = 75;
      const lengthGCol = view.grid.columns.byName('Length')!;
      lengthGCol.width = 0;
    })();

    await step('Analyze for activity cliffs.', async () => {
      activityCliffsViewer = (await activityCliffs(
        df, df.getCol('Sequence'), df.getCol('Activity'),
        80, method)) as DG.ScatterPlotViewer;
      view.dockManager.dock(activityCliffsViewer, DG.DOCK_TYPE.RIGHT, null, 'Activity Cliffs', 0.35);

      // Show grid viewer with the cliffs
      const cliffsLink: HTMLButtonElement = $(activityCliffsViewer.root)
        .find('button.scatter_plot_link,cliffs_grid').get()[0] as HTMLButtonElement;
      cliffsLink.click();
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

      // adjust for visual
      const activityGCol = view.grid.columns.byName('Activity')!;
      activityGCol.scrollIntoView();
    })();

    await step('Browse the cliff.', async () => {
      const currentCliffIdx = 1;
      const cliffsDfGrid: DG.Grid = activityCliffsViewer.dataFrame.temp[acTEMPS.cliffsDfGrid];
      cliffsDfGrid.dataFrame.selection.init((i) => i == currentCliffIdx);
      cliffsDfGrid.dataFrame.currentRowIdx = currentCliffIdx;

      /* workaround to select rows of the cliff */
      const entryCol: DG.Column = df.getCol('Entry');
      df.selection.init((rowIdx) => ['UPI00000BFE1D', 'UPI00000BFE17'].includes(entryCol.get(rowIdx)));

      const selectionIdxList: Int32Array = df.selection.getSelectedIndexes();
      if (selectionIdxList.length > 0) {
        df.currentRowIdx = selectionIdxList[0];
        view.grid.scrollToCell('UniProtKB', view.grid.tableRowToGrid(selectionIdxList[0]));
      }
    })();
  } catch (err: any) {
    if (err instanceof Error)
      _package.logger.error(err.message, undefined, err.stack);
    else
      _package.logger.error(err.toString());
  }
}
