/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PhylocanvasGlViewer} from './viewers/phylocanvas-gl-viewer';
import {PhylocanvasGlViewerApp} from './apps/phylocanvas-gl-viewer-app';
import {TreeToGridApp} from './apps/tree-to-grid-app';
import {injectTreeToGridUI} from './viewers/inject-tree-to-grid';
import {TreeInGridCellApp} from './apps/tree-in-grid-cell-app';
import {PhylocanvasGlService} from './utils/phylocanvas-gl-service';
import {TreeCutAsTreeApp} from './apps/tree-cut-as-tree-app';
import {findNewick} from './scripts-api';
import {PhylocanvasGlServiceBase} from '@datagrok-libraries/bio/src/viewers/phylocanvas-gl-viewer';
import {getTreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';


export const _package = new DG.Package();


// -- Viewers --

//name: PhylocanvasGL
//description: Phylogenetic tree visualization
//tags: viewer
//meta.icon: files/icons/phylocanvasgl-viewer.svg
//output: viewer result
export function phylocanvasGlViewer(): PhylocanvasGlViewer {
  return new PhylocanvasGlViewer();
}

// -- Apps for tests --

//name: phylocanvasGlViewerApp
//description: Test/demo app for PhylocanvasGlViewer
export async function phylocanvasGlViewerApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open PhylocanvasGlViewer app');
  try {
    const app = new PhylocanvasGlViewerApp('phylocanvasGlViewerApp');
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
//output: object result
export async function treeToGridApp(): Promise<TreeToGridApp | undefined> {
  const pi = DG.TaskBarProgressIndicator.create('open treeInGrid app');
  try {
    const app = new TreeToGridApp();
    await app.init();
    return app;
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
//output: object result
export async function injectTreeToGrid(grid: DG.Grid, newickText: string, leafColName?: string): Promise<GridNeighbor> {
  const colNameList: string[] = grid.dataFrame.columns.names();
  leafColName = leafColName ??
    grid.dataFrame.getTag('.newickLeafColumn') ??
    colNameList.find((colName) => colName.toLowerCase() == 'node') ??
    colNameList.find((colName) => colName.toLowerCase() == 'leaf') ??
    colNameList.find((colName) => colName.toLowerCase() == 'id');
  if (!leafColName)
    throw new Error('The leaf column name can not be inferred. Specify it as an argument.');

  const neighbor: GridNeighbor = await injectTreeToGridUI(grid, newickText, leafColName!);
  return neighbor;
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
