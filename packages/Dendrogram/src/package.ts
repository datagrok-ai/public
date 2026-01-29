/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Dendrogram} from './viewers/dendrogram';
import {TreeHelper} from './utils/tree-helper';
import {DendrogramApp} from './apps/dendrogram-app';
import {HierarchicalClusteringApp} from './apps/hierarchical-clustering-app';
import {hierarchicalClusteringDialog, hierarchicalClusteringUI} from './utils/hierarchical-clustering';
import {TreeForGridApp} from './apps/tree-for-grid-app';
import {TreeForGridFilterApp} from './apps/tree-for-grid-filter-app';
import {TreeForGridCutApp} from './apps/tree-for-grid-cut-app';
import {DendrogramService} from './utils/dendrogram-service';
import {DistanceMetric, NodeType} from '@datagrok-libraries/bio/src/trees';
import {IDendrogramService} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {HierarchicalClusteringSequencesApp} from './apps/hierarchical-clustering-sequences-app';
import {heatmapDemo} from './demos/heatmapDemo';
import {DendrogramPackage} from './package-types';

export * from './package.g';
export const _package = new DendrogramPackage(/*{debug: true}/**/);

/*
Scripting parameter types
https://datagrok.ai/help/compute/scripting
 */

type DendrogramWindowType = Window & { $dendrogramService?: IDendrogramService }
declare const window: DendrogramWindowType;

export class PackageFunctions {
  @grok.decorators.func({
    'name': 'info'
    })
  static info() {
    grok.shell.info(_package.webRoot);
  }


  @grok.decorators.func({
    'meta': {'icon': 'files/icons/dendrogram-viewer.svg', role: 'viewer'},
    'outputs': [
    {
    'name': 'result',
    'type': 'viewer'
    }
    ],
    'name': 'Dendrogram',
    'description': 'Dendrogram tree visualization'
    })
  static dendrogram(): DG.JsViewer {
    return new Dendrogram();
  }


  @grok.decorators.func({
    'meta': {},
    'outputs': [
    {
    'name': 'result',
    'type': 'object'
    }
    ],
    'name': 'getTreeHelper'
    })
  static getTreeHelper(): ITreeHelper {
    return new TreeHelper();
  }


  @grok.decorators.func({
    'outputs': [
    {
    'name': 'result',
    'type': 'object'
    }
    ],
    'name': 'getDendrogramService'
    })
  static getDendrogramService(): IDendrogramService {
    if (!(window.$dendrogramService)) {
      const svc: IDendrogramService = new DendrogramService();
      window.$dendrogramService = svc;
    }
    return window.$dendrogramService;
  }


  @grok.decorators.func({
    'name': 'generateTreeDialog'
    })
  static generateTreeDialog() {
    const sizeInput = ui.input.int('Tree size (node count)', {value: 10000});
    const filenameInput = ui.input.string('File name', {value: 'tree-gen-10000'});

    return ui.dialog('Generate tree')
      .add(ui.divV([sizeInput, filenameInput]))
      .onOK(async () => {
        const th: ITreeHelper = new TreeHelper();
        const treeRoot: NodeType = th.generateTree(sizeInput.value!);
        const treeNwk = th.toNewick(treeRoot);

        const leafList = th.getLeafList(treeRoot);
        const leafCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Leaf',
          leafList.map((n) => n.name));
        const activityCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Activity',
          leafList.map((_n) => Math.random()));

        const df = DG.DataFrame.fromColumns([leafCol, activityCol]);
        await _package.files.writeAsText(filenameInput.value + '.nwk', treeNwk);
        await _package.files.writeAsText(filenameInput.value + '.csv', df.toCsv());
      })
      .show();
  }


  @grok.decorators.func({
    'meta': {},
    'name': 'dendrogramApp',
    'description': 'Test/demo app for Dendrogram'
    })
  static async dendrogramApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open Dendrogram app');
    try {
      const app = new DendrogramApp();
      await app.init();
    } catch (err: unknown) {
      const msg: string = 'Dendrogram: dendrogramApp() error: ' +
        `${err instanceof Error ? err.message : (err as Object).toString()}`;
      grok.shell.error(msg);
      console.error(msg);
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'dendrogramLargeApp',
    'description': 'Test/demo app for Dendrogram Large'
    })
  static async dendrogramLargeApp(): Promise<void> {
    const largeDataSize: number = 100000;
    const largeDataFn: string = 'data/tree-gen-100000.nwk';

    const pi = DG.TaskBarProgressIndicator.create('open Dendrogram Large app');
    try {
      const th = new TreeHelper();
      let largeNewickStr: string;
      if (await _package.files.exists(largeDataFn)) {
        largeNewickStr = await _package.files.readAsText('data/tree-gen-100000.nwk');
      } else {
        grok.shell.warning(`File '${largeDataFn}' does not exist, generating data...`);
        largeNewickStr = th.toNewick(th.generateTree(largeDataSize));
      }

      const largeTreeDf: DG.DataFrame = th.newickToDf(largeNewickStr, 'large');
      const app = new DendrogramApp();
      await app.init(largeTreeDf, 'dendrogramLargeApp');
    } catch (err: unknown) {
      const msg: string = 'Dendrogram: dendrogramLargeApp() error: ' +
        `${err instanceof Error ? err.message : (err as Object).toString()}`;
      grok.shell.error(msg);
      console.error(msg);
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'treeForGridApp',
    'description': 'Test/demo app for TreeForGrid (custom renderer)'
    })
  static async treeForGridApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open treeForGridFilter app');
    try {
      const app = new TreeForGridApp();
      await app.init();
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'treeForGridFilterApp',
    'description': 'Test/demo app for TreeForGridFilter (custom renderer)'
    })
  static async treeForGridFilterApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open treeForGrid large app');
    try {
      const app = new TreeForGridFilterApp();
      await app.init();
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'treeForGridCutApp',
    'description': 'Test/demo app for TreeForGridCutApp (custom renderer, cutting slider)'
    })
  static async treeForGridCutApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open treeForGridCut large app');
    try {
      const app = new TreeForGridCutApp();
      await app.init();
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'hierarchicalClusteringApp',
    'description': 'Test/demo app for hierarchical clustering (inject tree to grid)'
    })
  static async hierarchicalClusteringApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('opem Hierarchical Clustering app');
    try {
      const app = new HierarchicalClusteringApp();
      await app.init();
    } catch (err: unknown) {
      const msg: string = 'Dendrogram: hierarchicalClusteringApp() error: ' +
        `${err instanceof Error ? err.message : (err as Object).toString()}`;
      grok.shell.error(msg);
      console.error(msg);
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'hierarchicalClusteringSequencesApp',
    'description': 'Test/demo app for hierarchical clustering (inject tree to grid)'
    })
  static async hierarchicalClusteringSequencesApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('opem Hierarchical Clustering app for sequences');
    try {
      const app = new HierarchicalClusteringSequencesApp();
      await app.init();
    } catch (err: unknown) {
      const msg: string = 'Dendrogram: hierarchicalClusteringApp() error: ' +
        `${err instanceof Error ? err.message : (err as Object).toString()}`;
      grok.shell.error(msg);
      console.error(msg);
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    'name': 'Hierarchical Clustering',
    'description': 'Calculates hierarchical clustering on features and injects tree to grid'
    })
  static async hierarchicalClustering(
    df: DG.DataFrame,
    colNameList: DG.ColumnList,
    @grok.decorators.param({type: 'string', options: {initialValue: 'euclidean', choices:['euclidean', 'manhattan']}}) distance: DistanceMetric = DistanceMetric.Euclidean,
    @grok.decorators.param({options: {initialValue: 'ward', choices:['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']}})linkage: string,
  ): Promise<void> {
    await hierarchicalClusteringUI(df, colNameList.names(), distance, linkage);
  }


  @grok.decorators.fileHandler({
    'ext': 'nwk, newick',
    'outputs': [
    {
    'name': 'tables',
    'type': 'list'
    }
    ],
    'name': 'importNwk',
    'description': 'Opens Newick file'
    })
  static async importNewick(
    @grok.decorators.param({'name':'fileContent', 'type':'string'}) fileContent: string): Promise<DG.DataFrame[]> {
    const th: ITreeHelper = new TreeHelper();
    const df: DG.DataFrame = th.newickToDf(fileContent, '');

    const app = new DendrogramApp();
    await app.init(df);

    return [];
  }


  @grok.decorators.fileViewer({
    'fileViewer': 'nwk,newick',
    'outputs': [
    {
    'name': 'preview',
    'type': 'view'
    }
    ]
    })
  static async previewNewick(
    @grok.decorators.param({'name':'file', 'type':'file'}) file: DG.FileInfo) {
    const newickString = await file.readAsString();
    const treeHelper = await PackageFunctions.getTreeHelper();
    const df = treeHelper.newickToDf(newickString, file.fileName.slice(0, -4));

    const viewerRoot = ((await df.plot.fromType('PhylocanvasGL', {})) as DG.JsViewer).root;
    viewerRoot.style.setProperty('width', '100%', 'important');
    viewerRoot.style.setProperty('height', '100%', 'important');

    return DG.View.fromRoot(viewerRoot);
  }


  @grok.decorators.func({
    'meta': {},
    'top-menu': 'Bio | Analyze | Hierarchical Clustering...',
    'name': 'Hierarchical Clustering (Sequences)',
    'description': 'Calculates hierarchical clustering on features and injects tree to grid'
    })
  static async hierarchicalClusteringSequences(): Promise<void> {
    hierarchicalClusteringDialog((t) => t.columns.bySemType(DG.SEMTYPE.MACROMOLECULE));
  }


  @grok.decorators.func({
    'top-menu': 'Chem | Analyze | Hierarchical Clustering...',
    'name': 'Hierarchical Clustering (Molecules)',
    'description': 'Calculates hierarchical clustering on features and injects tree to grid'
    })
  static async hierarchicalClusteringMolecules(): Promise<void> {
    hierarchicalClusteringDialog((t) => t.columns.bySemType(DG.SEMTYPE.MOLECULE));
  }


  @grok.decorators.func({
    'top-menu': 'ML | Cluster | Hierarchical Clustering...',
    'name': 'Hierarchical Clustering (All)',
    'description': 'Calculates hierarchical clustering on features and injects tree to grid'
    })
  static async hierarchicalClustering2(): Promise<void> {
    hierarchicalClusteringDialog((t) => t.columns.bySemType(DG.SEMTYPE.MOLECULE) ?? t.columns.bySemType(DG.SEMTYPE.MACROMOLECULE));
  }


  @grok.decorators.func({
    'meta': {
    'demoPath': 'Visualization | General | Heatmap'
    },
    'name': 'heatMapDemo',
    'description': 'Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.'
    })
  static async _heatMapDemo() {
    await heatmapDemo();
  }
}
