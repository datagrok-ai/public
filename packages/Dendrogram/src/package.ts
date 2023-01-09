/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Dendrogram, MyViewer} from './viewers/dendrogram';
import {TreeHelper} from './utils/tree-helper';
import {DendrogramApp} from './apps/dendrogram-app';
import {HierarchicalClusteringApp} from './apps/hierarchical-clustering-app';
import {IDendrogramService, ITreeHelper, NodeType, TreeCutOptions} from '@datagrok-libraries/bio';
import {hierarchicalClusteringUI} from './utils/hierarchical-clustering';
import {TreeForGridFilterApp} from './apps/tree-for-grid-filter-app';
import {TreeForGridApp} from './apps/tree-for-grid-app';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {injectTreeForGridUI2} from './viewers/inject-tree-for-grid2';
import {DendrogramService} from './utils/dendrogram-service';

export const _package = new DG.Package();

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
  const sizeInput = ui.intInput('Tree size (node count)', 10000);
  const filenameInput = ui.stringInput('File name', 'tree-gen-10000');

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
        leafList.map((n) => Math.random()));

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
      largeNewickStr = th.toNewick(th.generateTree(100000));
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

// -- Top menu --

//top-menu: ML | Hierarchical Clustering ...
//name: hierarchicalClustering
//description: Calculates hierarchical clustering on features and injects tree to grid
//input: dataframe table
//input: column_list features {type: numerical}
//input: string distance = 'euclidean' {choices: ['euclidean', 'manhattan']}
//input: string linkage = 'ward' {choices: ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']}
export async function hierarchicalClustering(
  table: DG.DataFrame, features: DG.ColumnList, distance: string, linkage: string
): Promise<void> {
  const colNameList: string[] = features.names();
  await hierarchicalClusteringUI(table, colNameList, distance, linkage);
}

