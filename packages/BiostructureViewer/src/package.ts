/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DockingApp} from './apps/docking-app';
import {byId, byData, MolstarViewer} from './viewers/molstar-viewer';
import {PdbGridCellRenderer} from './utils/pdb-grid-cell-renderer';
import {NglGlService} from './utils/ngl-gl-service';
import {NglForGridTestApp} from './apps/ngl-for-grid-test-app';
import {nglViewerGen as _nglViewerGen} from './utils/ngl-viewer-gen';

import {TwinPviewer} from './viewers/twin-p-viewer';
import {PROPS as nglPROPS, NglViewer} from './viewers/ngl-viewer';
import {NglViewerApp} from './apps/ngl-viewer-app';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {PdbHelper} from './utils/pdb-helper';
import {PdbApp} from './apps/pdb-app';
import {nglViewUI, nglWidgetUI} from './viewers/ngl-ui';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {NglGlServiceBase} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';
import {MolstarViewerApp} from './apps/molstar-viewer-app';
import {TaskBarProgressIndicator} from 'datagrok-api/dg';
import {dockingDemoApp} from './demo/docking';
import {biostructureInGridApp} from './demo/biostructure-in-grid';

class Package extends DG.Package {
  private _pLogger: DG.PackageLogger;

  get logger(): DG.PackageLogger {
    if (!this._pLogger) {
      this._pLogger = new class extends DG.PackageLogger {
        private logPrefix: string = 'BsV: ';

        constructor(_package: DG.Package) { super(_package); }

        debug(message: string, params?: object): void { super.debug(this.logPrefix + message, params); }
      }(this);
    }

    return this._pLogger;
  }
}

export const _package = new Package();

//name: pdbCellRenderer
//tags: cellRenderer
//meta.cellType: Molecule3D
//meta.columnTags: quality=Molecule3D
//output: grid_cell_renderer result
export function Molecule3dCellRenderer(): PdbGridCellRenderer {
  return new PdbGridCellRenderer();
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

//name: viewPdbById
//input: string pdbId
export async function viewPdbById(pdbId: string) {
  const pi = DG.TaskBarProgressIndicator.create(`Opening PDB ${pdbId}...`);
  try {
    await byId(pdbId);
  } finally {
    pi.close();
  }
}

//name: viewPdbByData
//input: string pdbData
//input: string name
export async function viewPdbByData(pdbData: string, name: string): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create(`Opening PDB '${name ?? '[data]'}'...`);
  try {
    await byData(pdbData, name);
  } finally {
    pi.close();
  }
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
  const ph: IPdbHelper = await getPdbHelper();
  const df: DG.DataFrame = await ph.pdbToDf(fileContent, '');

  const app = new PdbApp();
  await app.init(df);

  return [];
}

// -- File (pre)viewers --

// eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-mol, fileViewer-cif, fileViewer-mcif, fileViewer-mmcif, fileViewer-gro, fileViewer-pdb, fileViewer-ent, fileViewer-pqr, fileViewer-mmtf, fileViewer-mtl, fileViewer-sd
//input: file file
//output: view v
export function molecule3dNglView1(file: any): DG.View {
  return nglViewUI(file);
}

//tags: fileViewer, fileViewer-ply, fileViewer-obj
//input: file file
//output: view v
export function molecule3dNglView2(file: any): DG.View {
  return nglViewUI(file);
}

//tags: fileViewer, fileViewer-prmtop, fileViewer-parm7, fileViewer-psf, fileViewer-top
//input: file file
//output: view v
export function molecule3dNglView3(file: any): DG.View {
  return nglViewUI(file);
}

// eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-dsn6, fileViewer-brix, fileViewer-cube, fileViewer-cub, fileViewer-dx, fileViewer-dxbin, fileViewer-xplor, fileViewer-cns, fileViewer-mrc, fileViewer-map, fileViewer-ccp4
//input: file file
//output: view v
export function molecule3dNglView4(file: any): DG.View {
  return nglViewUI(file);
}

// -- Panel widgets --

//name: PDB id viewer
//tags: panel
//input: string pdbId {semType: PDB_ID}
//output: widget result
export function pdbIdNglPanelWidget(pdbId: string): DG.Widget {
  return nglWidgetUI(pdbId);
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

//name: molstarViewerApp
//description: Test app for MolstarViewer
export async function molstarViewerApp() {
  const pi = DG.TaskBarProgressIndicator.create('open molstarViewer app');
  try {
    const app = new MolstarViewerApp('molstarViewerApp');
    await app.init();
  } finally {
    pi.close();
  }
}

// -- Viewers --

//name: Ngl
// eslint-disable-next-line max-len
//description: 3D structure viewer for large biological molecules (proteins, DNA, and RNA)
//meta.keywords: PDB, Biostructure
//tags: viewer, panel
//output: viewer result
export function nglViewer(): DG.JsViewer {
  return new NglViewer();
}

//name: Biostructure
// eslint-disable-next-line max-len
//description: 3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)
//meta.keywords: Molstar, PDB
//tags: viewer, panel
//output: viewer result
export function molstarViewer(): DG.JsViewer {
  return new MolstarViewer();
}

// -- Top menu --

// code pdbViewer moved to utils/pdb-viewer.ts because of annotations in comments

// -- Utils --


//name: getPdbHelper
//output: object result
export async function getPdbHelper(): Promise<IPdbHelper> {
  if (!(window.$pdbHelper)) {
    const ph = await PdbHelper.getInstance(); // getPdbHelper
    window.$pdbHelper = ph;
  }
  return window.$pdbHelper;
}

export async function nglViewerGen(): Promise<void> {
  _nglViewerGen();
}

//name: dockingDemo
//meta.demoPath: Viewers | Docking
export async function dockingDemo() {
  const piMsg: string = 'Opening docking demo app ...';
  const pi: TaskBarProgressIndicator = TaskBarProgressIndicator.create(piMsg);
  try {
    await dockingDemoApp('dockingDemo', pi);
  } finally {
    pi.close();
  }
}

//name: inGridDemo
//meta.demoPath: Viewers | In Grid
export async function inGridDemo() {
  const piMsg: string = 'Opening biostructure in grid demo app ...';
  const pi: TaskBarProgressIndicator = TaskBarProgressIndicator.create(piMsg);
  try {
    await biostructureInGridApp('inGridDemo', pi);
  } finally {
    pi.close();
  }
}
