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
import {GridWithTreeApp} from './apps/grid-with-tree-app';
import {injectTreeToGridUI} from './viewers/inject-tree-to-grid';
import {NewickHelper} from './utils/newick-helper';
import {TreeInGridCellApp} from './apps/tree-in-grid-cell-app';
import {PhylocanvasGlService} from './utils/phylocanvas-gl-service';


export const _package = new DG.Package();


//name: newickToDf
//input: string newick
//input: string name
//output: dataframe df
export function _newickToDf(newick: string, name: string): DG.DataFrame {
  return newickToDf(newick, name);
};

let _newickHelper: NewickHelper | null = null;

//name: getNewickHelper
//description: Get object of interface bio.NewickHelper
//output: object result
export function getNewickHelper() {
  if (!_newickHelper)
    _newickHelper = new NewickHelper();
  return _newickHelper;
}

//name: PhyloTree
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function phyloTreeViewer(): PhyloTreeViewer {
  return new PhyloTreeViewer();
}


//name: PhylocanvasGL
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function phylocanvasGlViewer(): PhylocanvasGlViewer {
  return new PhylocanvasGlViewer();
}


//name: GridWithTree
//tags: viewer
//output: viewer result
export function gridWithTreeViewer(): GridWithTreeViewer {
  return new GridWithTreeViewer();
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


//name: PhylocanvasGlViewer
//tags: app
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

//name: GridWithTree
//tags: app
export async function gridWithTreeApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open GridWithTreeViewer app');
  try {
    const app = new GridWithTreeApp();
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

//name: TreeInGridCell
//tags: app
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

//name: importNwk
//description: Opens Newick file
//tags: file-handler
//meta.ext: nwk, newick
//input: string fileContent
//output: list tables
export async function importNewick(fileContent: string): Promise<DG.DataFrame[]> {
  const df: DG.DataFrame = newickToDf(fileContent, '');

  const app = new PhylocanvasGlViewerApp();
  await app.init(df);

  return [];
}

//name: injectTree
//description: Opens Newick file
//input: viewer grid
//input: string newickText
export async function injectTreeToGrid(grid: DG.Grid, newickText: string) {
  injectTreeToGridUI(grid, newickText);
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