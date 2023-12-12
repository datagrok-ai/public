/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {INglViewer} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';
import {NglGlServiceBase} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';
import {IBiostructureViewer} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {IBiotrackViewer} from '@datagrok-libraries/bio/src/viewers/biotrack';
import {AutoDockRunResult, IAutoDockService} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {delay} from '@datagrok-libraries/utils/src/test';

import {byData, byId, MolstarViewer} from './viewers/molstar-viewer';
import {SaguaroViewer} from './viewers/saguaro-viewer';
import {PdbGridCellRenderer, PdbGridCellRendererBack} from './utils/pdb-grid-cell-renderer';
import {NglForGridTestApp} from './apps/ngl-for-grid-test-app';
import {nglViewerGen as _nglViewerGen} from './utils/ngl-viewer-gen';
import {NglViewer} from './viewers/ngl-viewer';
import {NglViewerApp} from './apps/ngl-viewer-app';
import {PdbResDataFrame} from './utils/pdb-helper';
import {nglWidgetUI} from './viewers/ngl-ui';
import {dockingDemoApp} from './demo/docking';
import {biostructureInGridApp} from './demo/biostructure-in-grid';
import {BiotrackViewerApp} from './apps/biotrack-viewer-app';
import {BiostructureAndTrackViewerApp} from './apps/biostructure-and-track-viewer-app';
import {previewBiostructure, previewNgl, viewBiostructure, viewNgl} from './viewers/view-preview';
import {BiostructureViewerApp} from './apps/biostructure-viewer-app';
import {LigandsWithBiostructureApp, LigandsWithNglApp} from './apps/ligands-with-base-app';
import {addContextMenuUI} from './utils/context-menu';
import {importPdbqtUI} from './utils/pdbqt/import-pdbqt';
import {IPdbGridCellRenderer} from './utils/types';
import {AutoDockApp} from './apps/auto-dock-app';
import {_getNglGlService, _getPdbHelper, Package} from './package-utils';
import {AutoDockService, _runAutodock, _runAutodock2} from './utils/auto-dock-service';
import {demoBio07NoScript} from './demo/bio07-molecule3d-in-grid';
import {demoBio06NoScript} from './demo/bio06-docking-ngl';
import {Pdbqt} from './utils/pdbqt/pdbqt-parser';

export const _package: Package = new Package();

//name: init
//tags: init
//description:
export async function init() {
  _package.logger.debug('BiostructureViewer.initBiostructureViewer() init package start');
}

//name: pdbCellRenderer
//tags: cellRenderer
//meta.cellType: Molecule3D
//meta.columnTags: quality=Molecule3D
//output: grid_cell_renderer result
export function Molecule3dCellRenderer(): PdbGridCellRenderer {
  return new PdbGridCellRenderer();
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

//name: getNglGlService
//output: object result
export function getNglGlService(): NglGlServiceBase {
  return _getNglGlService();
}

// -- File handlers --

/* The Chem package is opening formats 'mol2', 'sdf', 'mol' for small molecules */
//name: importPdb
//description: Opens PDB file
//tags: file-handler
//meta.ext: mmcif, cifCore, pdb, gro
//input: string fileContent
//output: list tables
export async function importPdb(fileContent: string): Promise<DG.DataFrame[]> {
  // Do not build up data frame from PDB file, allows to open various formats
  // const ph: IPdbHelper = await _getPdbHelper();
  // const df: DG.DataFrame = await ph.pdbToDf(fileContent, '');
  // const app = new BiostructureApp();
  // await app.init(df);

  await viewBiostructure(fileContent);
  return [];
}

/* as file is handled as string we don't know its extension, thus we need a separate handler **/
//name: importXYZ
//description: Opens XYZ file
//tags: file-handler
//meta.ext: xyz
//input: string fileContent
//output: list tables
export async function importXYZ(fileContent: string): Promise<DG.DataFrame[]> {
  await viewBiostructure(fileContent, 'xyz');
  return [];
}

//name: importWithNgl
//description: Opens biostructure files supported with NGL
//tags: file-handler
//meta.ext: mmtf, cns, top, prmtop, ply, obj, ccp4
//input: string fileContent
//output: list tables
export async function importWithNgl(fileContent: string): Promise<DG.DataFrame[]> {
  await viewNgl(fileContent);
  return [];
}

//name: importPdbqt
//description: Opens .pdbqt file with docking result ligand poses
//tags: file-handler
//meta.ext: pdbqt
//input: string fileContent
//input: bool test =false {optional: true}
//output: list tables
export async function importPdbqt(fileContent: string, test: boolean): Promise<DG.DataFrame[]> {
  return importPdbqtUI(fileContent, test);
}

// -- File (pre)viewers --

/** Structure formats not supported with Molstar are handled with NGL viewer
 * fileViewer-mtl is from NglViewer package, but NGL can not open it
 * //TODO: Support preview .mol2 with Molstar (hangs on sp-after.mol2)
 * //TODO: Fix preview .pqr
 */
//tags: fileViewer, fileViewer-mmtf, fileViewer-cns, fileViewer-top, fileViewer-prmtop
//input: file file
//output: view result
export function previewNglStructure(file: any): DG.View {
  return previewNgl(file);
}

/** Shape / surface formats not supported with Molstar are handled with NGL viewer
 * TODO: Support preview .ply with Molstar
 */
//tags: fileViewer, fileViewer-ply, fileViewer-obj
//input: file file
//output: view v
export function previewNglSurface(file: any) {
  return previewNgl(file);
}

// TODO: Support preview .ccp4 with Molstar
//eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-ccp4
//input: file file
//output: view v
export function previewNglDensity(file: any) {
  return previewNgl(file);
}


// eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-mol, fileViewer-mol2, fileViewer-cif, fileViewer-mcif, fileViewer-mmcif, fileViewer-gro, fileViewer-pdb, fileViewer-ent, fileViewer-sd, fileViewer-xyz
//input: file file
//output: view v
export function previewBiostructureStructure(file: DG.FileInfo): DG.View {
  return previewBiostructure(file);
}

//tags: fileViewer, fileViewer-parm7, fileViewer-psf
//input: file file
//output: view v
export function previewBiostructureTopology(file: DG.FileInfo): DG.View {
  return previewBiostructure(file);
}

// eslint-disable-next-line max-len
//tags: fileViewer, fileViewer-dsn6, fileViewer-brix, fileViewer-cube, fileViewer-cub, fileViewer-dx, fileViewer-dxbin, fileViewer-xplor, fileViewer-mrc, fileViewer-map
//input: file file
//output: view v
export function previewBiostructureDensity(file: DG.FileInfo): DG.View {
  return previewBiostructure(file);
}

//name: openPdbResidues
//input: file fi
export async function openPdbResidues(fi: DG.FileInfo): Promise<void> {
  const ph = await _getPdbHelper();
  const pdbStr: string = await fi.readAsString();
  const pdbDf: PdbResDataFrame = await ph.pdbToDf(pdbStr, fi.fileName);
  const view = grok.shell.addTableView(pdbDf);
  const viewer = await pdbDf.plot.fromType('NGL', {});
  view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.40);
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

//name: biostructureViewerApp
//description: Test app for BiostructureViewer (molstar)
export async function biostructureViewerApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open biostructureViewer app');
  try {
    const app = new BiostructureViewerApp('biostructureViewerApp');
    await app.init();
  } finally {
    pi.close();
  }
}

//name: biotrackViewerApp
//description: Test app for BiotrackViewer (saguaro)
export async function biotrackViewerApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open biotrackViewer app');
  try {
    const app = new BiotrackViewerApp('biotrackViewerApp');
    await app.init();
  } finally {
    pi.close();
  }
}

//name: biostructureAndTrackViewerApp
//description: Test app for twin BiostructureViewer (molstar) and BiotrackViewer (saguaro)
export async function biostructureAndTrackViewerApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open biostructureAndTrackViewer app');
  try {
    const app = new BiostructureAndTrackViewerApp('biostructureAndTrackViewerApp');
    await app.init();
  } finally {
    pi.close();
  }
}

//name: ligandsWithNglApp
export async function ligandsWithNglApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Ligands with NGL app');
  try {
    const app = new LigandsWithNglApp('ligandsWithNglApp');
    await app.init();
  } finally {
    pi.close();
  }
}

//name: ligandsWithBiostructureApp
export async function ligandsWithBiostructureApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Ligands with Biostructure app');
  try {
    const app = new LigandsWithBiostructureApp('ligandsWithBiostructureApp');
    await app.init();
  } finally {
    pi.close();
  }
}


//name: autoDockApp
export async function autoDockApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('AutoDock app...');
  try {
    const app = new AutoDockApp('autoDockApp');
    await app.init();
  } finally {
    pi.close();
  }
}


// -- Viewers --

//name: NGL
// eslint-disable-next-line max-len
//description: 3D structure viewer for large biological molecules (proteins, DNA, and RNA)
//meta.keywords: PDB, Biostructure
//meta.icon: files/icons/ngl-viewer.svg
//tags: viewer, panel
//output: viewer result
export function nglViewer(): DG.JsViewer & INglViewer {
  return new NglViewer();
}

//name: Biostructure
// eslint-disable-next-line max-len
//description: 3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)
//meta.keywords: Molstar, PDB
//meta.icon: files/icons/biostructure-viewer.svg
//tags: viewer, panel
//output: viewer result
export function molstarViewer(): DG.JsViewer & IBiostructureViewer {
  return new MolstarViewer();
}

//name: Biotrack
//description: structure polymer annotation tracks
//meta.keywords: PDB, track
//tags: viewer, panel
//output: viewer result
export function saguaroViewer(): DG.JsViewer & IBiotrackViewer {
  return new SaguaroViewer();
}

// -- Top menu --

// code pdbViewer moved to utils/pdb-viewer.ts because of annotations in comments

// -- Utils --


//name: getPdbHelper
//output: object result
export async function getPdbHelper(): Promise<IPdbHelper> {
  return _getPdbHelper();
}

export async function nglViewerGen(): Promise<void> {
  _nglViewerGen();
}

//name: getPdbGridCellRenderer
//input: object gridCol
//output: object result
export function getPdbGridCellRenderer(gridCol: DG.GridColumn): IPdbGridCellRenderer {
  return PdbGridCellRendererBack.getOrCreate(gridCol);
}

//name: dockingDemo
//description:
export async function dockingDemo() {
  const piMsg: string = 'Opening docking demo app ...';
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(piMsg);
  try {
    await dockingDemoApp('dockingDemo', pi);
  } finally {
    pi.close();
  }
}

//name: inGridDemo
//description:
export async function inGridDemo() {
  const piMsg: string = 'Opening biostructure in grid demo app ...';
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(piMsg);
  try {
    await biostructureInGridApp('inGridDemo', pi);
  } finally {
    pi.close();
  }
}

// -- Handle context menu --

//name: addContextMenu
//input: object event
export function addContextMenu(event: DG.EventData): void {
  addContextMenuUI(event);
}

// -- Demo --

//  // demoBio06
//  //name: demoBioDockingConformationsOld
//  //meta.demoPath: Bioinformatics | Docking Conformations Old
//  //description: Docking ligands along the structure
//  //meta.path: /apps/Tutorials/Demo/Bioinformatics/Docking%20Conformations%20Old
//  //test: demoBioDockingConformationsOld() //wait: 3000
//  export async function demoBioDockingConformationsOld(): Promise<void> {
//    // Do not use any script for this demo (askalkin, 2023-05-17)
//    //await demoBio06UI();
//    await demoBio06NoScript();
//  }

// demoBio06b
//name: demoBioDockingConformations
//meta.demoPath: Bioinformatics | Docking Conformations
//descriptions: Display ligand poses along the structure
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Docking%20Conformations
//test: demoBioDockingConformations() //wait: 1000
export async function demoBioDockingConformations(): Promise<void> {
  await demoBio06NoScript();
}

// demoBio07
//name: demoBioProteins
//meta.demoPath: Bioinformatics | Proteins
//description: View structures PDB in grids
//meta.path: /apps/Tutorials/Demo/Bioinformatics/Proteins
//test: demoBioProteins() //wait: 3000
export async function demoBioProteins(): Promise<void> {
  await demoBio07NoScript();
}


//name: getAutoDockService
//output: object result
export async function getAutoDockService(): Promise<IAutoDockService> {
  const resSvc: IAutoDockService = await AutoDockService.getSvc();
  return resSvc;
}

//name: runAutodock
//input: file receptor
//input: file ligand
//input: int x
//input: int y
//input: int z
//output: object dockingResults
export async function runAutodock(
  receptor: DG.FileInfo, ligand: DG.FileInfo, x: number, y: number, z: number,
): Promise<AutoDockRunResult | null> {
  return await _runAutodock(receptor, ligand, x, y, z);
}

//name: runAutodock2
// //input: dataframe df
// //input: column molCol { semType: Molecule }
// //input: file receptorFi { optional: true }
export async function runAutodock2(): Promise<void> {
  const [csv, receptorPdb] = await Promise.all([
    grok.functions.call(`${_package.name}:readAsText`,
      {file: 'CHEMBL2366517/ic50.mol.csv'}),
    grok.functions.call(`${_package.name}:readAsText`, {file: 'CHEMBL2366517/1bdq.pdb'}),
  ]);
  if (!csv || !receptorPdb)
    throw new Error('Empty input data');

  const df = DG.DataFrame.fromCsv(csv);
  const molCol: DG.Column<string> = df.getCol('molecule');

  return await _runAutodock2(molCol, receptorPdb);
}

//name: runAutodock3
//input: dataframe name
//input: column ligandCol
export async function runAutodock3(df: DG.DataFrame, ligandCol: DG.Column): Promise<void> {

}

const dataDir = 'Admin:Data/PDB/CHEMBL2366517/';

//name: demoFix1
export async function demoFix1(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('demoFix1 loading ...');
  try {
    const [ligandCsv, poseCsv] = await Promise.all([
      grok.dapi.files.readAsText(dataDir + 'ic50.mol.csv'),
      grok.dapi.files.readAsText(dataDir + 'ic50.pose.csv'),
      getPdbHelper(),
    ]);
    const ligandDf = DG.DataFrame.fromCsv(ligandCsv);
    const poseDf = DG.DataFrame.fromCsv(poseCsv);

    if (ligandDf.rowCount !== poseDf.rowCount / 30)
      throw new Error('Unsupported input data');

    // region Cols
    const Cols = {
      id: 'Molecule ChEMBL ID',
      molName: 'Molecule Name',
      mw: 'Molecular Weight',
      ic50: 'IC50',
      ic50units: 'IC50 Units',

      affinity: 'affinity',
      intermolecular: 'intermolecular (1)',
      electrostatic: 'electrostatic',
      ligandFixed: 'ligand fixed',
      ligandMoving: 'ligand moving',
      totalInternal: 'total internal (2)',
      torsionalFree: 'torsional free (3)',
      unboundSystems: 'unbound systems (4)',

      mol: 'mol',
    };

    const lIdCol = ligandDf.getCol(Cols.id);
    const lMolNameCol = ligandDf.getCol(Cols.molName);
    const lMwCol = ligandDf.getCol(Cols.mw);

    const lIc50ValCol = ligandDf.getCol('Standard Value');
    const lIc50UnitsCol = ligandDf.getCol('Standard Units');

    const pModelCol = poseDf.getCol('pdbqt_model');
    const pIdCol = poseDf.col(Cols.id) ?? poseDf.columns.addNewString(Cols.id);
    const pMolNameCol = poseDf.col(Cols.molName) ?? poseDf.columns.addNewString(Cols.molName);
    const pMwCol = poseDf.col(Cols.mw) ?? poseDf.columns.addNewFloat(Cols.mw);
    const pIc50ValCol = poseDf.col(Cols.ic50) ?? poseDf.columns.addNewFloat(Cols.ic50);
    const pIc50UnitsCol = poseDf.col(Cols.ic50units) ?? poseDf.columns.addNewString(Cols.ic50units);

    const ppAffinityCol = poseDf.col(Cols.affinity) ?? poseDf.columns.addNewFloat(Cols.affinity);
    const ppIntermolecularCol = poseDf.col(Cols.intermolecular) ?? poseDf.columns.addNewFloat(Cols.intermolecular);
    const ppElectrostaticCol = poseDf.col(Cols.electrostatic) ?? poseDf.columns.addNewFloat(Cols.electrostatic);
    const ppLigandFixedCol = poseDf.col(Cols.ligandFixed) ?? poseDf.columns.addNewFloat(Cols.ligandFixed);
    const ppLigandMovingCol = poseDf.col(Cols.ligandMoving) ?? poseDf.columns.addNewFloat(Cols.ligandMoving);
    const ppTotalInternalCol = poseDf.col(Cols.totalInternal) ?? poseDf.columns.addNewFloat(Cols.totalInternal);
    const ppTorsionalFreeCol = poseDf.col(Cols.torsionalFree) ?? poseDf.columns.addNewFloat(Cols.torsionalFree);
    const ppUnboundSystemsCol = poseDf.col(Cols.unboundSystems) ?? poseDf.columns.addNewFloat(Cols.unboundSystems);

    const pMolCol = poseDf.col(Cols.mol) ?? poseDf.columns.addNewString(Cols.mol);
    // endregion Cols
    let lastProgress: number = 0;
    const poseCount = poseDf.rowCount;
    for (let poseI = 0; poseI < poseCount; ++poseI) {
      const ligandI = Math.trunc(poseI / 30);

      pIdCol.set(poseI, lIdCol.get(ligandI));
      pMolNameCol.set(poseI, lMolNameCol.get(ligandI));
      pMwCol.set(poseI, lMwCol.get(ligandI));
      pIc50ValCol.set(poseI, lIc50ValCol.get(ligandI));
      pIc50UnitsCol.set(poseI, lIc50UnitsCol.get(ligandI));

      const modelPdbqtStr = pModelCol.get(poseI);
      const pm = Pdbqt.parse(modelPdbqtStr).models[0];
      // region Docking props
      ppAffinityCol.set(poseI, pm.affinity);
      ppIntermolecularCol.set(poseI, pm.intermolecular);
      ppElectrostaticCol.set(poseI, pm.electrostatic);
      ppLigandFixedCol.set(poseI, pm.ligandFixed);
      ppLigandMovingCol.set(poseI, pm.ligandMoving);
      ppTotalInternalCol.set(poseI, pm.totalInternal);
      ppTorsionalFreeCol.set(poseI, pm.torsionalFree);
      ppUnboundSystemsCol.set(poseI, pm.unboundSystems);
      // endregion Docking props

      const progress = poseI / poseCount;
      if (progress - lastProgress >= 0.01) {
        pi.update(100 * progress, 'demoFix1 parsing ...');
        lastProgress = progress;
        await delay(0);
      }
    }

    pi.update(1, 'demoFix1 downloading ...');
    DG.Utils.download('ic50.pose-props.csv', poseDf.toCsv());
  } finally {
    pi.close();
  }
}

//name: demoFix2
export async function demoFix2(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('demoFix2 loading ...');
  try {
    const [poseCsv] = await Promise.all([
      grok.dapi.files.readAsText(dataDir + 'ic50.pose-props.flt2.csv'),
    ]);
    if (!poseCsv)
      throw new Error('Empty input data');

    const poseDf = DG.DataFrame.fromCsv(poseCsv);
    const pdbqtModelCol = poseDf.getCol('pdbqt_model');
    const pdbqtVal = pdbqtModelCol.toList().join('\n');

    pi.update(1, 'demoFix2 downloading .pdbqt ...');
    DG.Utils.download('ic50.pose.flt2.pdbqt', pdbqtVal);
  } finally {
    pi.close();
  }
}

//name: demoFix3
export async function demoFix3(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('demoFix3 loading ...');
  try {
    const [poseCsv, ph] = await Promise.all([
      grok.dapi.files.readAsText(dataDir + 'ic50.pose-props.flt2.csv'),
      getPdbHelper(),
    ]);
    if (!poseCsv)
      throw new Error('Empty input data');

    const poseDf = DG.DataFrame.fromCsv(poseCsv);
    const pdbqtModelCol = poseDf.getCol('pdbqt_model');
    const pMolCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'mol', poseDf.rowCount);
    poseDf.columns.insert(pMolCol, 0);

    let lastProgress: number = 0;
    const poseCount: number = poseDf.rowCount;
    for (let poseI = 0; poseI < poseCount; ++poseI) {
      const pdbqtModelStr = pdbqtModelCol.get(poseI);
      const pMol: string = await ph.pdbqtToMol(pdbqtModelStr);
      pMolCol.set(poseI, pMol);

      const progress = poseI / poseCount;
      if (progress - lastProgress >= 0.01) {
        pi.update(100 * progress, 'demoFix3 pdbqt -> mol ...');
        lastProgress = progress;
        await delay(0);
      }
    }

    pi.update(1, 'demoFix3 downloading .pdbqt ...');
    DG.Utils.download('ic50.pose.flt2.mol.csv', poseDf.toCsv());
  } finally {
    pi.close();
  }
}


// -- cache --
//name: readAsText
//meta.cache: client
//meta.invalidateOn: 0 0 1 * * ?
//input: string file
//output: string result
export async function readAsText(file: string): Promise<string> {
  return await _package.files.readAsText(file);
}

//name: readAsTextDapi
//meta.cache: client
//meta.invalidateOn: 0 0 1 * * ?
//input: string file
//output: string result
export async function readAsTextDapi(file: string): Promise<string> {
  const [resStr, exists] = await Promise.all([
    grok.dapi.files.readAsText(file),
    grok.dapi.files.exists(file)
  ]);
  if (!exists)
    throw new Error(`File not found '${file}'.`);
  return resStr;
}
