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

export const _package = new DendrogramPackage(/*{debug: true}/**/);

/*
Scripting parameter types
https://datagrok.ai/help/compute/scripting
 */

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


// -- Viewers --

//name: Dendrogram
//description: Dendrogram tree visualization
//meta.icon: files/icons/dendrogram-viewer.svg
//tags: viewer
//output: viewer result
export function dendrogram(): DG.JsViewer {
  return new Dendrogram();
}


// -- Utils --

//name: getTreeHelper
//output: object result
export function getTreeHelper(): ITreeHelper {
  return new TreeHelper();
}

type DendrogramWindowType = Window & { $dendrogramService?: IDendrogramService }
declare const window: DendrogramWindowType;

//name: getDendrogramService
//output: object result
export function getDendrogramService(): IDendrogramService {
  if (!(window.$dendrogramService)) {
    const svc: IDendrogramService = new DendrogramService();
    window.$dendrogramService = svc;
  }
  return window.$dendrogramService;
}

//name: generateTreeDialog
export function generateTreeDialog() {
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


// -- Apps for tests --

//name: dendrogramApp
//description: Test/demo app for Dendrogram
export async function dendrogramApp(): Promise<void> {
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

//name: dendrogramLargeApp
//description: Test/demo app for Dendrogram Large
export async function dendrogramLargeApp(): Promise<void> {
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

//name: treeForGridApp
//description: Test/demo app for TreeForGrid (custom renderer)
export async function treeForGridApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open treeForGridFilter app');
  try {
    const app = new TreeForGridApp();
    await app.init();
  } finally {
    pi.close();
  }
}

//name: treeForGridFilterApp
//description: Test/demo app for TreeForGridFilter (custom renderer)
export async function treeForGridFilterApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open treeForGrid large app');
  try {
    const app = new TreeForGridFilterApp();
    await app.init();
  } finally {
    pi.close();
  }
}

//name: treeForGridCutApp
//description: Test/demo app for TreeForGridCutApp (custom renderer, cutting slider)
export async function treeForGridCutApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open treeForGridCut large app');
  try {
    const app = new TreeForGridCutApp();
    await app.init();
  } finally {
    pi.close();
  }
}

//name:hierarchicalClusteringApp
//description: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringApp(): Promise<void> {
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

//name:hierarchicalClusteringSequencesApp
//description: Test/demo app for hierarchical clustering (inject tree to grid)
export async function hierarchicalClusteringSequencesApp(): Promise<void> {
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

//name: Hierarchical Clustering
//description: Calculates hierarchical clustering on features and injects tree to grid
//input: dataframe df
//input: column_list colNameList
//input: string distance = 'euclidean' {choices: ['euclidean', 'manhattan']}
//input: string linkage = 'ward' {choices: ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']}
export async function hierarchicalClustering(
  df: DG.DataFrame,
  colNameList: DG.ColumnList,
  distance: DistanceMetric = DistanceMetric.Euclidean,
  linkage: string,
): Promise<void> {
  await hierarchicalClusteringUI(df, colNameList.names(), distance, linkage);
}


// -- File handlers --

//name: importNwk
//description: Opens Newick file
//tags: file-handler
//meta.ext: nwk, newick
//input: string fileContent
//output: list tables
export async function importNewick(fileContent: string): Promise<DG.DataFrame[]> {
  const th: ITreeHelper = new TreeHelper();
  const df: DG.DataFrame = th.newickToDf(fileContent, '');

  const app = new DendrogramApp();
  await app.init(df);

  return [];
}

// -- File preview --

//tags: fileViewer
//meta.fileViewer: nwk,newick
//input: file file
//output: view preview
export async function previewNewick(file: DG.FileInfo) {
  const newickString = await file.readAsString();
  const treeHelper = await getTreeHelper();
  const df = treeHelper.newickToDf(newickString, file.fileName.slice(0, -4));

  const preview = DG.View.create();
  const host = ui.divH([
    ui.button('Load dataframe', () => {
      const view = grok.shell.addTableView(df);
      const viewer = DG.Viewer.fromType('PhylocanvasGL', df);
      view.addViewer(viewer);
      view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT);
      return view;
    }, 'View in a dataframe'),
    DG.Viewer.fromType('Dendrogram', df).root,
  ], 'd4-ngl-viewer');

  preview.append(host);
  return preview;
}

// -- Top menu --

//top-menu: Bio | Analyze | Hierarchical Clustering...
//name: Hierarchical Clustering (Sequences)
//description: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClusteringSequences(): Promise<void> {
  hierarchicalClusteringDialog();
}

//top-menu: Chem | Analyze | Hierarchical Clustering...
//name: Hierarchical Clustering (Molecules)
//description: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClusteringMolecules(): Promise<void> {
  hierarchicalClusteringDialog();
}

//top-menu: ML | Cluster | Hierarchical Clustering...
//name: Hierarchical Clustering (All)
//description: Calculates hierarchical clustering on features and injects tree to grid
export async function hierarchicalClustering2(): Promise<void> {
  hierarchicalClusteringDialog();
}

// -- Demo --

//name: heatMapDemo
// eslint-disable-next-line max-len
//description: Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.
//meta.demoPath: Visualization | General | Heatmap
//test: _heatMapDemo() //wait: 2000
export async function _heatMapDemo() {
  await heatmapDemo();
};
