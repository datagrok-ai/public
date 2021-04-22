/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { newickToDf } from './utils.js';
import { PhyloTreeViewer } from './tree-viewer.js';


export const _package = new DG.Package();

//name: PhyloTree
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function tree() {
  return new PhyloTreeViewer();
}

//tags: fileViewer, fileViewer-nwk
//input: file file
//output: view preview
export async function nwkTreeViewer(file) {
  const newickString = await file.readAsString();
  const df = newickToDf(newickString, file.fileName);

  const preview = DG.View.create();
  const host = ui.divH([
    ui.button('Load dataframe', () => {
      const view = grok.shell.addTableView(df);
      const viewer = DG.Viewer.fromType('PhyloTree', df);
      view.addViewer(viewer);
      view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT);
      return view;
    }, 'View in a dataframe'),
    DG.Viewer.fromType('PhyloTree', df).root
  ], 'd4-ngl-viewer');

  preview.append(host);
  return preview;
}
