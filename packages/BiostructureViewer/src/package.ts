/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DockingApp} from './apps/docking-app';
import {byId, byData} from './viewers/molstar-viewer';
import {PdbRenderer} from './utils/cell-renderer';
import {NglGlService} from './utils/ngl-gl-service';
import {NglForGridTestApp} from './apps/ngl-for-grid-test-app';
import {nglViewerGen as _nglViewerGen} from './utils/ngl-viewer-gen';

import {IPdbHelper, NglGlServiceBase} from '@datagrok-libraries/bio';
import {TwinPviewer} from './viewers/twin-p-viewer';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {PROPS as nglPROPS, NglViewer} from './viewers/ngl-viewer';
import {NglViewerApp} from './apps/ngl-viewer-app';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {PdbHelper} from './utils/pdb-helper';
import {PdbApp} from './apps/pdb-app';

export const _package = new DG.Package();

//name: pdbCellRenderer
//tags: cellRenderer
//meta.cellType: xray
//meta.columnTags: quality=xray
//output: grid_cell_renderer result
export function xRayCellRenderer(): PdbRenderer {
  return new PdbRenderer();
}

//name: Docking
//tags: app
export async function dockingApp() {
  const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure Viewer');
  try {
    const app = new DockingApp();
    await app.init();
  } finally {
    pi.close();
  }
}

//name: Mol* BioStructure Viewer Id.
//tags: viewer
//input: string pdbID
export async function molstarViewId(pdbID: string) {
  const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure* Viewer Id');
  await byId(pdbID);
  pi.close();
}

//name: Mol* BioStructure Viewer Data.
//tags: viewer
export async function molstarViewData() {
  const pdbData = await _package.files.readAsText('1bdq.pdb');
  const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure* Viewer Data');
  await byData(pdbData);
  pi.close();
}


type BsvWindowType = Window & {
  $pdbHelper?: IPdbHelper,
  $nglGlService?: NglGlServiceBase,
};
declare const window: BsvWindowType;

//name: getNglGlService
//output: object result
export function getNglGlService(): NglGlServiceBase {
  if (!(window.$nglGlService)) {
    const svc: NglGlService = new NglGlService();
    window.$nglGlService = svc;
  }

  return window.$nglGlService;
}

// -- File handlers --

//name: importPdb
//description: Opens PDB file
//tags: file-handler
//meta.ext: pdb
//input: string fileContent
//output: list tables
export async function importPdb(fileContent: string): Promise<DG.DataFrame[]> {
  const ph: IPdbHelper = new PdbHelper();
  const df: DG.DataFrame = await ph.pdbToDf(fileContent, '');

  const app = new PdbApp();
  await app.init(df);

  return [];
}

// -- Test apps --

//name: nglForGridTestApp
//description: Example app for NGL drawing in grid cells
export async function nglForGridTestApp() {
  const pi = DG.TaskBarProgressIndicator.create('open nglForGridTest app');
  try {
    const app = new NglForGridTestApp();
    await app.init();
  } finally {
    pi.close();
  }
}

//name: nglViewerApp
//description: Test app for NglViewer
export async function nglViewerApp() {
  const pi = DG.TaskBarProgressIndicator.create('open nglViewer app');
  try {
    const app = new NglViewerApp('nglViewerApp');
    await app.init();
  } finally {
    pi.close();
  }
}

// -- Viewers --

//name: NglViewer
//description: 3D structure viewer for large biological molecules (proteins, DNA, and RNA)
//tags: viewer, panel
//output: viewer result
export function nglViewer(): DG.JsViewer {
  return new NglViewer();
}

// -- Top menu --

// //name: PDB Viewer
// //description: 3D structure data for large biological molecules (proteins, DNA, and RNA)
// //top-menu: Bio | PDB ...
// //output: viewer result
// export async function pdbViewer(): Promise<void> {
//   const view: DG.TableView = grok.shell.tv;
//   const pdbTag = view.dataFrame.getTag(pdbTAGS.PDB);
//   if (pdbTag) {
//     const viewer = (await view.dataFrame.plot.fromType('NglViewer', {})) as DG.JsViewer;
//     view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL viewer', 0.4);
//   } else {
//     const bsView: DG.TableView = grok.shell.tv;
//
//     const ligandSelection = {};
//     let ligands = ['R', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'S', 'T', 'U', 'V', 'W'];
//     for (let i = 0; i < ligands.length; i++)
//       ligandSelection[ligands[i]] = [false, 400 + i];
//
//     return new Promise<void>((resolve, reject) => {
//       const fileBrowser = ui.fileBrowser({path: `System:AppData/${_package.name}/samples`});
//       const dlg: DG.Dialog = ui.dialog({title: 'Open PDB file'})
//         .add(fileBrowser.root)
//         .addButton('OK', () => {
//           setTimeout(async () => {
//             const filePath: string = fileBrowser.props.file;
//             const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('PDB Viewer');
//             try {
//               const pdbStr: string = await grok.dapi.files.readAsText(filePath);
//
//               const viewer: DG.JsViewer = (await view.dataFrame.plot.fromType('NglViewer',
//                 {[nglPROPS.pdb]: pdbStr})) as DG.JsViewer;
//               view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL viewer', 0.4);
//
//               resolve();
//             } catch (err: any) {
//               const errMsg: string = errorToConsole(err);
//               console.error(errMsg);
//               reject(err.toString());
//             } finally {
//               pi.close();
//               dlg.close();
//             }
//           }, 0 /* next event cycle */);
//         })
//         .show();
//     });
//   }
// }

// -- Utils --


//name: getPdbHelper
//output: object result
export function getPdbHelper(): IPdbHelper {
  if (!(window.$pdbHelper)) {
    const ph = new PdbHelper();
    window.$pdbHelper = ph;
  }
  return window.$pdbHelper;
}

export async function nglViewerGen(): Promise<void> {
  _nglViewerGen();
}