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

const dataFn = 'data/sample_FASTA_DNA.csv';

enum EMBED_COL_NAMES {
  X = 'Embed_X',
  Y = 'Embed_Y'
}

export async function demoBio01aUI(funcPath: string) {
  let treeHelper: ITreeHelper;
  let dendrogramSvc: IDendrogramService;
  let view: DG.TableView;
  let df: DG.DataFrame;
  let spViewer: DG.ScatterPlotViewer;

  const method: string = 'UMAP';
  const idRows: { [id: number]: number } = {};
  const embedCols: { [colName: string]: DG.Column<number> } = {};

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
      view.grid.props.rowHeight = 22;
      view.path = view.basePath = funcPath;
    } finally { pi.close(); }
  }

  async function step02() {
    grok.shell.info(`Sequence space.`);
    const pi = DG.TaskBarProgressIndicator.create('Analyze for sequence space...');
    try {
      if (true) {
        // Custom sequence space implementation for closer resembling of hierarchical clustering results.
        const embedColNameList = Object.values(EMBED_COL_NAMES);
        // ensure embed columns exist
        for (let embedI: number = 0; embedI < embedColNameList.length; embedI++) {
          const embedColName: string = embedColNameList[embedI];
          const embedCol: DG.Column | null = df.col(embedColName);
          if (!embedCol) {
            // Notification is required to reflect added data frame Embed_<X> columns to grid columns
            // MolecularLiabilityBrowser.setView() corrects grid columns' names with .replace('_', ' ');
            const notify: boolean = embedI == embedColNameList.length - 1; // notify on adding last Embed_<X> column
            df.columns.add(DG.Column.float(embedColName, df.rowCount), notify);
          }
        }

        if (df.rowCount >= 1) {
          const seqCol: DG.Column<string> = df.getCol('sequence');
          const seqList = seqCol.toList();

          const t1: number = Date.now();
          _package.logger.debug('Bio: demoBio01aUI(), calc reduceDimensionality start...');
          const redDimRes = await reduceDimensinalityWithNormalization( // TODO: Rename method typo
            seqList, method, StringMetricsNames.Levenshtein, {});
          const t2: number = Date.now();
          _package.logger.debug('Bio: demoBio01aUI(), calc reduceDimensionality ' +
            `ET: ${((t2 - t1) / 1000)} s`);

          for (let embedI: number = 0; embedI < embedColNameList.length; embedI++) {
            const embedColName: string = embedColNameList[embedI];
            const embedCol: DG.Column = df.getCol(embedColName);
            const embedColData: Float32Array = redDimRes.embedding[embedI];
            // TODO: User DG.Column.setRawData()
            // embedCol.setRawData(embedColData);
            embedCol.init((rowI) => { return embedColData[rowI]; });
          }

          const rowCount: number = df.rowCount;
          const idCol: DG.Column = df.getCol('id');
          for (let idRowI = 0; idRowI < rowCount; idRowI++) {
            const id = idCol.get(idRowI);
            idRows[id] = idRowI;
          }

          for (const embedColName of Object.values(EMBED_COL_NAMES)) {
            const embedCol: DG.Column<number> = df.getCol(embedColName);
            embedCols[embedColName] = embedCol;
          }

          const t3: number = Date.now();
          _package.logger.debug('MLB: MlbVrSpaceBrowser.buildView(), postprocess reduceDimensionality ' +
            `ET: ${((t3 - t2) / 1000)} s`);

          spViewer = (await df.plot.fromType(DG.VIEWER.SCATTER_PLOT, {
            'xColumnName': EMBED_COL_NAMES.X,
            'yColumnName': EMBED_COL_NAMES.Y,
            'lassoTool': true,
          })) as DG.ScatterPlotViewer;
          view.dockManager.dock(spViewer, DG.DOCK_TYPE.RIGHT, null, 'Sequence space', 0.3);
        }
      } else {
        const spaceViewer = (await sequenceSpaceTopMenu(df, df.getCol('sequence'),
          'UMAP', StringMetricsNames.Levenshtein, true)) as DG.Viewer;
        view.dockManager.dock(spaceViewer, DG.DOCK_TYPE.RIGHT, null, 'Sequence Space', 0.35);
      }
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
    df.selection.init((idx: number) => [15].includes(idx));
  }

  /** Select some bunch of similar sequences */
  async function step05() {
    grok.shell.info(`Select bunch of sequences.`);
    df.selection.init((idx: number) => [21, 9, 58].includes(idx));
    df.currentRowIdx = 27;
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
    .then(async () => { await step04c(); })
    .then(async () => { await delay(1600); })
    .then(async () => { await step05(); })
    .catch((err) => _package.logger.error(err));
}
