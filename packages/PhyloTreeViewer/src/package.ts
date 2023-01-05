/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {newickToDf} from './utils';
import {PhyloTreeViewer} from './tree-viewer';
import {PhylocanvasGlViewer} from './viewers/phylocanvas-gl-viewer';
import {PhylocanvasGlViewerApp} from './apps/phylocanvas-gl-viewer-app';
import {TreeToGridApp} from './apps/tree-to-grid-app';
import {injectTreeToGridUI} from './viewers/inject-tree-to-grid';
import {TreeInGridCellApp} from './apps/tree-in-grid-cell-app';
import {PhylocanvasGlService} from './utils/phylocanvas-gl-service';
import {TreeCutAsTreeApp} from './apps/tree-cut-as-tree-app';
import {findNewick} from './scripts-api';
import {ITreeHelper, NodeType, PhylocanvasGlServiceBase} from '@datagrok-libraries/bio';


export const _package = new DG.Package();


//name: newickToDf
//input: string newick
//input: string name
//output: dataframe df
export function _newickToDf(newick: string, name: string): DG.DataFrame {
  return newickToDf(newick, name);
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
declare const window: PtvWindowType;

//name: getPhylocanvasGlService
//output: object result
export function getPhylocanvasGlService(): PhylocanvasGlServiceBase {
  if (!(window.$phylocanvasGlService)) {
    const svc: PhylocanvasGlService = new PhylocanvasGlService();
    window.$phylocanvasGlService = svc;
  }

  return window.$phylocanvasGlService;
}

//name: newickRepresentation
//input: dataframe data
//input: column col
//output: string newick
export async function newickRepresentation(data: DG.DataFrame) {
  const columns = Array.from(data.columns.numerical);
  for (let i = 0; i < columns.length; ++i) {
    if (columns[i].type === DG.TYPE.DATE_TIME)
      columns.splice(i, 1);
  }
  return await findNewick(DG.DataFrame.fromColumns(columns));
}

