/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {newickToDf} from './utils';
import {PhyloTreeViewer} from './tree-viewer';
import {PhylocanvasGlViewer} from './viewers/phylocanvas-gl-viewer';
import {PhylocanvasGlViewerApp} from './apps/phylocanvas-gl-viewer-app';
import {GridWithTreeViewer} from './viewers/grid-with-tree-viewer';
import {GridWithTreeViewerApp} from './apps/grid-with-tree-viewer-app';


export const _package = new DG.Package();

//name: newickToDf
//input: string newick
//input: string name
//output: dataframe df
export function _newickToDf(newick: string, name: string): DG.DataFrame {
  return newickToDf(newick, name);
};

//name: PhyloTree
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function phyloTreeViewer(): PhyloTreeViewer {
  return new PhyloTreeViewer();
}


//name: PhylocanvasGl
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


//tags: fileViewer, fileViewer-nwk
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

//name: GridWithTreeViewer
//tags: app
export async function gridWithTreeViewerApp() {
  const pi = DG.TaskBarProgressIndicator.create('open GridWithTreeViewer app');
  try {
    const app = new GridWithTreeViewerApp();
    await app.init();
  } catch (err: unknown) {
    const msg: string = 'PhyloTreeViewer gridWithTreeViewerApp() error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(msg);
  } finally {
    pi.close();
  }
}