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

const dataFn = 'samples/sample_FASTA.csv';

export async function demoBio01bUI(funcPath: string) {
  let treeHelper: ITreeHelper;
  let dendrogramSvc: IDendrogramService;
  let view: DG.TableView;
  let df: DG.DataFrame;
  let activityCliffsViewer: DG.ScatterPlotViewer;

  const method: string = 'UMAP';
  const idRows: { [id: number]: number } = {};

  /** Load data and get helpers */
  async function step01() {
    grok.shell.info(`Loading DNA sequences notation 'fasta'.`);
    const pi = DG.TaskBarProgressIndicator.create('Loading DNA notation \'fasta\'.');
    try {
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
    } finally { pi.close(); }
  }

  async function step02() {
    grok.shell.info(`Activity Cliffs.`);
    const pi = DG.TaskBarProgressIndicator.create('Analyze for activity cliffs...');
    try {
      activityCliffsViewer = (await activityCliffs(
        df, df.getCol('Sequence'), df.getCol('Activity'),
        80, method)) as DG.ScatterPlotViewer;
      view.dockManager.dock(activityCliffsViewer, DG.DOCK_TYPE.RIGHT, null, 'Activity Cliffs', 0.35);

      // Show grid viewer with the cliffs
      const cliffsLink: HTMLButtonElement = $(activityCliffsViewer.root)
        .find('button.scatter_plot_link,cliffs_grid').get()[0] as HTMLButtonElement;
      cliffsLink.click();
    } finally { pi.close(); }
  }

  async function step03() {
    grok.shell.info('Hierarchical clustering.');
    const pi = DG.TaskBarProgressIndicator.create('Hierarchical clustering...');
    try {
      const seqCol: DG.Column<string> = df.getCol('sequence');
      const seqList = seqCol.toList();

      const t1: number = Date.now();
      const distance: DistanceMatrix = DistanceMatrix.calc(seqList, (aSeq: string, bSeq: string) => {
        const levDistance = lev.distance(aSeq, bSeq);
        return levDistance / ((aSeq.length + bSeq.length) / 2);
      });
      const t2: number = Date.now();
      console.debug('Bio: demoBio01aUI(), ' +
        `, calc distance matrix ET: ${((t2 - t1) / 1000)} s`);

      const treeRoot = await treeHelper.hierarchicalClusteringByDistance(distance, 'ward');
      dendrogramSvc.injectTreeForGrid(view.grid, treeRoot, undefined, 150, undefined);

      // adjust for visual
      const activityGCol = view.grid.columns.byName('Activity')!;
      activityGCol.scrollIntoView();
    } finally {
      pi.close();
    }
  }

  /** Switch selection */
  async function step04a() {
    grok.shell.info(`Selection ...`);
    df.selection.init((idx: number) => [2].includes(idx));
  }

  /** Switch selection */
  async function step04b() {
    grok.shell.info(`Selection ...`);
    df.selection.init((idx: number) => [7, 8, 9].includes(idx));
  }

  /** Switch selection */
  async function step04c() {
    grok.shell.info(`Selection ...`);
    //df.selection.init((idx: number) => [15].includes(idx));
    // TODO: Select some rows
  }

  /** Select some bunch of similar sequences */
  async function step05() {
    grok.shell.info(`Select bunch of sequences.`);
    df.selection.init((idx: number) => [21, 9, 58].includes(idx));
    df.currentRowIdx = 27;
  }

  /** To make current the second cliff in the cliff grid */
  async function step06() {
    grok.shell.info(`Browse the cliff`);

    const cliffsDfGrid: DG.Grid = activityCliffsViewer.dataFrame.temp[acTEMPS.cliffsDfGrid];
    // Selecting the cliff leads to select objects of the cliff
    cliffsDfGrid.dataFrame.selection.init((i) => i == 1);

    /* workaround to select rows of the cliff */
    const entryCol: DG.Column = df.getCol('Entry');
    df.selection.init((rowIdx) => ['UPI00000BFE1D', 'UPI00000BFE17'].includes(entryCol.get(rowIdx)));

    const selectionIdxList: Int32Array = df.selection.getSelectedIndexes();
    if (selectionIdxList.length > 0)
      view.grid.scrollToCell('UniProtKB', selectionIdxList[0]);
  }

  Promise.resolve()
    .then(async () => { await step01(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step02(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step03(); })
    .then(async () => { await delay(1600); })
    // .then(async () => { await step04a(); })
    // .then(async () => { await delay(1600); })
    // .then(async () => { await step04b(); })
    // .then(async () => { await delay(1600); })
    // .then(async () => { await step04c(); })
    // .then(async () => { await delay(1600); })
    // .then(async () => { await step05(); })
    // .then(async () => { await delay(1600); })
    .then(async () => { await step06(); })
    .catch((err) => _package.logger.error(err));
}
