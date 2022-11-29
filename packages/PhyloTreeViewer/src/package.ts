/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {newickToDf} from './utils';
import {PhyloTreeViewer} from './tree-viewer';
import {PhylocanvasGlViewer} from './viewers/phylocanvas-gl-viewer';
import {PhylocanvasGlViewerApp} from './apps/phylocanvas-gl-viewer-app';
import {GridWithTreeViewer} from './viewers/grid-with-tree-viewer';
import {TreeToGridApp} from './apps/tree-to-grid-app';
import {injectTreeToGridUI} from './viewers/inject-tree-to-grid';
import {NewickHelper} from './utils/newick-helper';
import {TreeInGridCellApp} from './apps/tree-in-grid-cell-app';
import {PhylocanvasGlService} from './utils/phylocanvas-gl-service';
import {TreeHelper} from './utils/tree-helper';
import {generateTree} from './utils/tree-generator';
import {TreeForGridApp} from './apps/tree-for-grid-app';
import {TreeCutAsTreeApp} from './apps/tree-cut-as-tree-app';
import {TreeForGridFilterApp} from './apps/tree-for-grid-filter-app';
import {Dendrogram, MyViewer} from './viewers/dendrogram';
import {DendrogramApp} from './apps/dendrogram-app';


export const _package = new DG.Package();


//name: newickToDf
//input: string newick
//input: string name
//output: dataframe df
export function _newickToDf(newick: string, name: string): DG.DataFrame {
  return newickToDf(newick, name);
}

let _newickHelper: NewickHelper | null = null;

//name: getNewickHelper
//description: Get object of interface bio.NewickHelper
//output: object result
export function getNewickHelper() {
  if (!_newickHelper)
    _newickHelper = new NewickHelper();
  return _newickHelper;
}


//tags: fileViewer, fileViewer-nwk, fileViewer-newick
//input: file file
//output: view preview
export async function nwkTreeViewer(file: DG.FileInfo) {
  const newickString = await file.readAsString();
  const df = newickToDf(newickString, file.fileName.slice(0, -4));

  const preview = DG.View.create();
  const host = ui.divH([
    ui.button('Load dataframe', () => {
      const view = grok.shell.addTableView(df);
      const viewer = DG.Viewer.fromType('PhyloTree', df);
      view.addViewer(viewer);
      view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT);
      return view;
    }, 'View in a dataframe'),
    DG.Viewer.fromType('Tree', df).root,
  ], 'd4-ngl-viewer');

  preview.append(host);
  return preview;
}

// -- Viewers --

//name: Dendrogram
//description: Dendrogram tree visualization
//tags: viewer
//output: viewer result
export function dendrogram(): DG.JsViewer {
  return new Dendrogram();
}

//name: PhylocanvasGL
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function phylocanvasGlViewer(): PhylocanvasGlViewer {
  return new PhylocanvasGlViewer();
}

//name: PhyloTree
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function phyloTreeViewer(): PhyloTreeViewer {
  return new PhyloTreeViewer();
}

// //name: GridWithTree
// //tags: viewer
// //output: viewer result
// export function gridWithTreeViewer(): GridWithTreeViewer {
//   return new GridWithTreeViewer();
// }


// -- apps for tests --

//name: dendrogramApp
//description: Test/demo app for Dendrogram
export async function dendrogramApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open Dendrogram app');
  try {
    const app = new DendrogramApp();
    await app.init();
  } catch (err: unknown) {
    const msg: string = 'PhyloTreeViewer dendrogramApp() error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(msg);
  } finally {
    pi.close();
  }
}

//name: phylocanvasGlViewerApp
//description: Test/demo app for PhylocanvasGlViewer
export async function phylocanvasGlViewerApp() {
  const pi = DG.TaskBarProgressIndicator.create('open PhylocanvasGlViewer app');
  try {
    const app = new PhylocanvasGlViewerApp();
    await app.init();
  } catch (err: unknown) {
    const msg: string = 'PhyloTreeViewer phylocanvasGlViewerApp() error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(msg);
  } finally {
    pi.close();
  }
}

//name: TreeToGrid
//description: Test/demo app for TreeToGrid (PhylocanvasGL based)
export async function treeToGridApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open treeInGrid app');
  try {
    const app = new TreeToGridApp();
    await app.init();
  } catch (err: unknown) {
    const msg: string = 'PhyloTreeViewer gridWithTreeViewerApp() error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    //@ts-ignore
    console.error(err);
    // if ('stack' in err)
    //   console.error(err['stack']);
  } finally {
    pi.close();
  }
}

//name: TreeForGrid
//description: Test/demo app for TreeForGrid (custom renderer)
export async function treeForGridApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open treeForGrid app');
  try {
    const app = new TreeForGridApp();
    await app.init();
  } finally {
    pi.close();
  }
}

//name: TreeForGridFilter
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

//name: TreeCutAsTree
//description: Test/demo app for TreeCutAsTree
export async function treeCutAsTreeApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open treeCutAsTree app');
  try {
    const app = new TreeCutAsTreeApp();
    await app.init();
  } finally {
    pi.close();
  }
}


//name: TreeInGridCell
//description: Test/demo app for TreeInGridCell
export async function treeInGridCellApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open TreeInGridCell app');
  try {
    const app = new TreeInGridCellApp();
    await app.init();
  } catch (err: unknown) {
    const msg: string = 'PhyloTreeViewer treeInGridCellApp() error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    //@ts-ignore
    console.error(err);
    // if ('stack' in err)
    //   console.error(err['stack']);
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
  const df: DG.DataFrame = newickToDf(fileContent, '');

  // const app = new PhylocanvasGlViewerApp();
  // await app.init(df);

  const app = new DendrogramApp();
  await app.init(df);

  return [];
}

// -- Custom helpers --

//name: injectTree
//description: Opens Newick file
//input: viewer grid
//input: string newickText
export async function injectTreeToGrid(grid: DG.Grid, newickText: string, leafColName?: string) {
  const colNameList: string[] = grid.dataFrame.columns.names();
  leafColName = leafColName ??
    grid.dataFrame.getTag('.newickLeafColumn') ??
    colNameList.find((colName) => colName.toLowerCase() == 'node') ??
    colNameList.find((colName) => colName.toLowerCase() == 'leaf') ??
    colNameList.find((colName) => colName.toLowerCase() == 'id');
  if (!leafColName)
    throw new Error('The leaf column name can not be inferred. Specify it as an argument.');

  injectTreeToGridUI(grid, newickText, leafColName!);
}

type PtvWindowType = Window & { $phylocanvasGlService?: PhylocanvasGlService };
declare var window: PtvWindowType;

//name: getPhylocanvasGlService
//output: object result
export function getPhylocanvasGlService(): bio.PhylocanvasGlServiceBase {
  if (!(window.$phylocanvasGlService)) {
    const svc: PhylocanvasGlService = new PhylocanvasGlService();
    window.$phylocanvasGlService = svc;
  }

  return window.$phylocanvasGlService;
}

//name: getTreeHelper
//output: object result
export function getTreeHelper(): bio.ITreeHelper {
  return new TreeHelper();
}

//name: generateTreeDialog
export function generateTreeDialog() {
  const sizeInput = ui.intInput('Tree size (node count)', 10000);
  const filenameInput = ui.stringInput('File name', 'tree-gen-10000');

  return ui.dialog('Generate tree')
    .add(ui.divV([sizeInput, filenameInput]))
    .onOK(async () => {
      const treeRoot: bio.NodeType = generateTree(sizeInput.value!);
      const th = new TreeHelper();
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

//name: MyViewer
//description: MyViewer test
//tags: viewer
//output: viewer result
export function myViewer(): DG.JsViewer {
  return new MyViewer();
}