/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RcsbGraphQLAdapter} from './utils/rcsb-gql-adapter';
import '@datagrok-libraries/bio/src/types/ngl'; // To enable import from the NGL module declared in bio lib

import {BuiltInTrajectoryFormat} from 'molstar/lib/mol-plugin-state/formats/trajectory';

import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {INglViewer} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';
import {NglGlServiceBase} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';
import {IBiostructureViewer} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {IBiotrackViewer} from '@datagrok-libraries/bio/src/viewers/biotrack';
import {BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';
import {delay} from '@datagrok-libraries/utils/src/test';

import {byData, byId, MolstarViewer} from './viewers/molstar-viewer';
import {SaguaroViewer} from './viewers/saguaro-viewer';
import {PdbGridCellRenderer, PdbGridCellRendererBack, PdbIdGridCellRenderer} from './utils/pdb-grid-cell-renderer';
import {NglForGridTestApp} from './apps/ngl-for-grid-test-app';
import {nglViewerGen as _nglViewerGen} from './utils/ngl-viewer-gen';
import {NglViewer} from './viewers/ngl-viewer';
import {NglViewerApp} from './apps/ngl-viewer-app';
import {PdbHelper, PdbResDataFrame} from './utils/pdb-helper';
import {nglWidgetUI} from './viewers/ngl-ui';
import {dockingDemoApp} from './demo/docking';
import {biostructureInGridApp} from './demo/biostructure-in-grid';
import {BiotrackViewerApp} from './apps/biotrack-viewer-app';
import {BiostructureAndTrackViewerApp} from './apps/biostructure-and-track-viewer-app';
import {previewBiostructure, previewNgl, viewNgl} from './viewers/view-preview';
import {BiostructureViewerApp} from './apps/biostructure-viewer-app';
import {LigandsWithBiostructureApp, LigandsWithNglApp} from './apps/ligands-with-base-app';
import {importPdbqtUI} from './utils/import-pdbqt';
import {_getNglGlService, BsvPackage} from './package-utils';
import {demoBio07NoScript} from './demo/bio07-molecule3d-in-grid';
import {demoBio06NoScript} from './demo/bio06-docking-ngl';
import {Pdbqt} from './utils/pdbqt-parser';
import {viewMolstarUI} from './viewers/molstar-viewer/utils';
import {BiostructureDataProviderApp} from './apps/biostructure-data-provider-app';
import {copyRawValue, downloadRawValue, showBiostructureViewer, showNglViewer} from './utils/context-menu';
import {defaultErrorHandler} from './utils/err-info';
import {extractSequenceColumns} from './utils/sequence-handler';

export * from './package.g';
export const _package: BsvPackage = new BsvPackage();


const dataDir = 'Admin:Data/PDB/CHEMBL2366517/';

export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    _package.logger.debug('BiostructureViewer.initBiostructureViewer() init package start');
  }

  @grok.decorators.func({
    name: 'pdbCellRenderer',
    meta: {cellType: 'Molecule3D', columnTags: 'quality=Molecule3D', role: 'cellRenderer'},
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static Molecule3dCellRenderer(): PdbGridCellRenderer {
    return new PdbGridCellRenderer();
  }

  @grok.decorators.func({
    name: 'chemCellRenderer',
    meta: {
      cellType: 'PDB_ID',
      role: 'cellRenderer'
    },
    outputs: [{type: 'grid_cell_renderer', name: 'result'}]
  })
  static pdbIdCellRenderer(): PdbIdGridCellRenderer {
    return new PdbIdGridCellRenderer();
  }

  @grok.decorators.func()
  static async viewPdbById(pdbId: string) {
    const pi = DG.TaskBarProgressIndicator.create(`Opening PDB ${pdbId}...`);
    try {
      await byId(pdbId);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async viewPdbByData(pdbData: string, name: string): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create(`Opening PDB '${name ?? '[data]'}'...`);
    try {
      await byData(pdbData, name);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static getNglGlService(): NglGlServiceBase {
    return _getNglGlService();
  }

  @grok.decorators.func()
  static async viewBiostructure(content: string, format?: string,
    @grok.decorators.param({options: {optional: true}}) name?: string): Promise<void> {
    await viewMolstarUI(content, name, format as BuiltInTrajectoryFormat);
  }

  // -- File handlers --

  /* The Chem package is opening formats 'mol2', 'sdf', 'mol' for small molecules */
  @grok.decorators.fileHandler({
    description: 'Opens PDB file',
    ext: 'mmcif, cifCore, pdb, gro'
  })
  static async importPdb(fileContent: string): Promise<DG.DataFrame[]> {
    // Do not build up data frame from PDB file, allows to open various formats
    // const ph: IPdbHelper = await _getPdbHelper();
    // const df: DG.DataFrame = await ph.pdbToDf(fileContent, '');
    // const app = new BiostructureApp();
    // await app.init(df);

    await PackageFunctions.viewBiostructure(fileContent, 'pdb');
    return [];
  }

  /* as file is handled as string we don't know its extension, thus we need a separate handler **/

  @grok.decorators.fileHandler({
    description: 'Opens XYZ file',
    ext: 'xyz'
  })
  static async importXYZ(fileContent: string): Promise<DG.DataFrame[]> {
    await PackageFunctions.viewBiostructure(fileContent, 'xyz');
    return [];
  }

  @grok.decorators.fileHandler({
    description: 'Opens biostructure files supported with NGL',
    ext: 'mmtf, cns, top, prmtop, ply, obj, ccp4'
  })
  static async importWithNgl(fileContent: string): Promise<DG.DataFrame[]> {
    await viewNgl(fileContent);
    return [];
  }

  @grok.decorators.fileHandler({
    description: 'Opens .pdbqt file with docking result ligand poses',
    ext: 'pdbqt'
  })
  static async importPdbqt(fileContent: string,
    @grok.decorators.param({options: {optional: true, initialValue: 'false'}}) test: boolean): Promise<DG.DataFrame[]> {
    return importPdbqtUI(fileContent, test);
  }

  // -- File (pre)viewers --

  /** Structure formats not supported with Molstar are handled with NGL viewer
   * fileViewer-mtl is from NglViewer package, but NGL can not open it
   * //TODO Support preview .mol2 with Molstar (hangs on sp-after.mol2)
   * //TODO Fix preview .pqr
   */
  @grok.decorators.fileViewer({fileViewer: 'mmtf,cns,top,prmtop,pqr'})
  static previewNglStructure( @grok.decorators.param({type: 'file'}) file: any): DG.View {
    return previewNgl(file);
  }

  /** Shape / surface formats not supported with Molstar are handled with NGL viewer
   * TODO Support preview .ply with Molstar
   */
  @grok.decorators.fileViewer({fileViewer: 'ply,obj'})
  static previewNglSurface( @grok.decorators.param({type: 'file'}) file: any): DG.View {
    return previewNgl(file);
  }

  //TODO Support preview .ccp4 with Molstar
  //eslint-disable-next-line max-len
  @grok.decorators.fileViewer({fileViewer: 'ccp4'})
  static previewNglDensity( @grok.decorators.param({type: 'file'}) file: any): DG.View {
    return previewNgl(file);
  }


  // eslint-disable-next-line max-len
  @grok.decorators.fileViewer({fileViewer: 'mol,mol2,cif,mcif,mmcif,gro,pdb,pdbqt,ent,sd,xyz'})
  static previewBiostructureStructure(file: DG.FileInfo): DG.View {
    return previewBiostructure(file);
  }

  @grok.decorators.fileViewer({fileViewer: 'parm7,psf'})
  static previewBiostructureTopology(file: DG.FileInfo): DG.View {
    return previewBiostructure(file);
  }

  @grok.decorators.fileViewer({fileViewer: 'dsn6,brix,cube,cub,dx,dxbin,xplor,mrc,map'})
  static previewBiostructureDensity(file: DG.FileInfo): DG.View {
    return previewBiostructure(file);
  }

  @grok.decorators.func()
  static async openPdbResidues(fi: DG.FileInfo): Promise<void> {
    const ph = await PdbHelper.getInstance();
    const pdbStr: string = await fi.readAsString();
    const pdbDf: PdbResDataFrame = await ph.pdbToDf(pdbStr, fi.fileName);
    const view = grok.shell.addTableView(pdbDf);
    const viewer = await pdbDf.plot.fromType('NGL', {});
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.40);
  }

  // -- Panel widgets --
  @grok.decorators.panel({
    name: 'PDB id viewer'
  })
  static pdbIdNglPanelWidget(
    @grok.decorators.param({options: {semType: 'PDB_ID'}}) pdbId: string
  ): DG.Widget {
    return nglWidgetUI(pdbId);
  }

  // -- Test apps --

  @grok.decorators.func({description: 'Example app for NGL drawing in grid cells'})
  static async nglForGridTestApp() {
    const pi = DG.TaskBarProgressIndicator.create('open nglForGridTest app');
    try {
      const app = new NglForGridTestApp();
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({description: 'Test app for NglViewer'})
  static async nglViewerApp() {
    const pi = DG.TaskBarProgressIndicator.create('open nglViewer app');
    try {
      const app = new NglViewerApp('nglViewerApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({description: 'Test app for BiostructureViewer (molstar)'})
  static async biostructureViewerApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open biostructureViewer app');
    try {
      const app = new BiostructureViewerApp('biostructureViewerApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({description: 'Test app for BiotrackViewer (saguaro)'})
  static async biotrackViewerApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open biotrackViewer app');
    try {
      const app = new BiotrackViewerApp('biotrackViewerApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({description: 'Test app for twin BiostructureViewer (molstar) and BiotrackViewer (saguaro)'})
  static async biostructureAndTrackViewerApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('open biostructureAndTrackViewer app');
    try {
      const app = new BiostructureAndTrackViewerApp('biostructureAndTrackViewerApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async ligandsWithNglApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Ligands with NGL app');
    try {
      const app = new LigandsWithNglApp('ligandsWithNglApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async ligandsWithBiostructureApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Ligands with Biostructure app');
    try {
      const app = new LigandsWithBiostructureApp('ligandsWithBiostructureApp');
      await app.init();
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async biostructureDataProviderApp(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Biostructure data provider app');
    try {
      const app = new BiostructureDataProviderApp();
      await app.init();
    } finally {
      pi.close();
    }
  }

  // -- Viewers --

  // eslint-disable-next-line max-len
  @grok.decorators.panel({
    name: 'NGL',
    description: '3D structure viewer for large biological molecules (proteins, DNA, and RNA)',
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {
      keywords: 'PDB, Biostructure',
      icon: 'files/icons/ngl-viewer.svg',
      role: 'viewer'
    }
  })
  static nglViewer(): DG.JsViewer & INglViewer {
    return new NglViewer();
  }

  @grok.decorators.panel({
    name: 'Biostructure',
    description: '3D structure molstar RCSB viewer for large biological molecules (proteins, DNA, and RNA)',
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {
      keywords: 'Molstar, PDB',
      icon: 'files/icons/biostructure-viewer.svg',
      role: 'viewer'
    }
  })
  static molstarViewer(): DG.JsViewer & IBiostructureViewer {
    return new MolstarViewer();
  }

  @grok.decorators.panel({
    name: 'Biotrack',
    description: 'structure polymer annotation tracks',
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {
      keywords: 'PDB, track',
      showInGallery: 'false',
      role: 'viewer'
    }
  })
  static saguaroViewer(): DG.JsViewer & IBiotrackViewer {
    return new SaguaroViewer();
  }

  // -- Top menu --
  // code pdbViewer moved to utils/pdb-viewer.ts because of annotations in comments

  // -- Utils --

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static async getPdbHelper(): Promise<IPdbHelper> {
    return PdbHelper.getInstance();
  }

  @grok.decorators.func()
  static async dockingDemo() {
    const piMsg: string = 'Opening docking demo app ...';
    const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(piMsg);
    try {
      await dockingDemoApp('dockingDemo', pi);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async inGridDemo() {
    const piMsg: string = 'Opening biostructure in grid demo app ...';
    const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(piMsg);
    try {
      await biostructureInGridApp('inGridDemo', pi);
    } finally {
      pi.close();
    }
  }

  // -- Handle context menu --

  @grok.decorators.func({name: 'Copy Biostructure raw value'})
  static async copyRawBiostructureValue(
    @grok.decorators.param({type: 'object'}) gridCell: DG.GridCell): Promise<void> {
    copyRawValue(gridCell);
  }

  @grok.decorators.func({name: 'Download Biostructure raw value'})
  static async downloadRawBiostructureValue(
    @grok.decorators.param({type: 'object'}) gridCell: DG.GridCell): Promise<void> {
    downloadRawValue(gridCell);
  }

  @grok.decorators.func({name: 'Show Biostructure Viewer menu item'})
  static async showBiostructureViewerMenuItem(
    @grok.decorators.param({type: 'object'}) gridCell: DG.GridCell) {
    await showBiostructureViewer(gridCell);
  }

  @grok.decorators.func({name: 'Show NGL Viewer menu item'})
  static async showNglViewerMenuItem(
    @grok.decorators.param({type: 'object'}) gridCell: DG.GridCell) {
    await showNglViewer(gridCell);
  }

  @grok.decorators.func({name: 'Open PDB residues table menu item'})
  static async openTableResiduesMenuItem(
    @grok.decorators.param({type: 'object'})fi: DG.FileInfo) {
    try {
      await PackageFunctions.openPdbResidues(fi);
    } catch (err: any) {
      defaultErrorHandler(err);
    }
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

  //name:
  //meta.demoPath:
  //meta.path:
  //meta.demoWait: 3000
  //meta.demoSkip: GROK-15250
  //test: demoBioDockingConformations() //wait: 3000, timeout: 60000, skip: GROK-15250
  @grok.decorators.demo({
    name: 'demoBioDockingConformations',
    description: 'Display ligand poses along the structure',
    demoPath: 'Bioinformatics | Docking Conformations',
    path: '/apps/Tutorials/Demo/Bioinformatics/Docking%20Conformations',
    demoWait: '3000',
    test: {test: 'demoBioDockingConformations()', wait: '3000', timeout: '60000'}
  })
  static async demoBioDockingConformations(): Promise<void> {
    await demoBio06NoScript();
  }

  @grok.decorators.demo({
    name: 'demoBioProteins',
    description: 'View structures PDB in grids',
    demoPath: 'Bioinformatics | Proteins',
    path: '/apps/Tutorials/Demo/Bioinformatics/Proteins',
    demoWait: '3000',
    test: {
      test: 'demoBioProteins()',
      wait: '3000',
      timeout: '60000',
    }
  })
  static async demoBioProteins(): Promise<void> {
    const t1: number = window.performance.now();
    await demoBio07NoScript();
    const t2: number = window.performance.now();
    _package.logger.debug(`demoBioProteins(), end ET: ${t2 - t1}`);
  }

  @grok.decorators.func()
  static async demoFix1(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('demoFix1 loading ...');
    try {
      const [ligandCsv, poseCsv] = await Promise.all([
        grok.dapi.files.readAsText(dataDir + 'ic50.mol.csv'),
        grok.dapi.files.readAsText(dataDir + 'ic50.pose.csv'),
        PackageFunctions.getPdbHelper(),
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

      const _pMolCol = poseDf.col(Cols.mol) ?? poseDf.columns.addNewString(Cols.mol);
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

  @grok.decorators.func()
  static async demoFix2(): Promise<void> {
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

  @grok.decorators.func()
  static async demoFix3(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('demoFix3 loading ...');
    try {
      const [poseCsv, ph] = await Promise.all([
        grok.dapi.files.readAsText(dataDir + 'ic50.pose-props.flt2.csv'),
        PackageFunctions.getPdbHelper(),
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

  @grok.decorators.func({
    meta: {
      cache: 'client',
      cacheInvalidateOn: '0 * * * *'
    }
  })
  static async readAsText(file: string): Promise<string> {
    return await _package.files.readAsText(file);
  }

  @grok.decorators.func({
    meta: {
      cache: 'client',
      cacheInvalidateOn: '0 * * * *'
    }
  })
  static async readAsTextDapi(file: string): Promise<string> {
    const [resStr, exists] = await Promise.all([
      grok.dapi.files.readAsText(file),
      grok.dapi.files.exists(file)
    ]);
    if (!exists)
      throw new Error(`File not found '${file}'.`);
    return resStr;
  }

  // -- Data provider --

  @grok.decorators.func({
    name: 'RCSB PDB',
    description: 'Get biostructure by id as PDB',
    meta: {
      dataProvider: 'Molecule3D'
    },
  })
  static async getBiostructureRcsbPdb(
    id: string
  ): Promise<string> {
    const url = `https://files.rcsb.org/download/${id}.pdb`;
    const response = await fetch(url);
    if (!response.ok)
      throw new Error(response.statusText);
    const data: string = await response.text();
    return BiostructureDataJson.fromData({binary: false, data: data, ext: 'pdb', options: {name: id}});
  }

  @grok.decorators.func({
    name: 'RCSB mmCIF',
    description: 'Get biostructure by id as mmCIF',
    meta: {
      dataProvider: 'Molecule3D',
      cache: 'client',
      cacheInvalidateOn: '0 * * * *',
    },
  })
  static async getBiostructureRcsbMmcif(id: string): Promise<string> {
    const url = `https://files.rcsb.org/download/${id}.cif`;
    const response = await fetch(url);
    if (!response.ok)
      throw new Error(response.statusText);
    const data: string = await response.text();
    return BiostructureDataJson.fromData({binary: false, data: data, ext: 'cif', options: {name: id}});
  }

  @grok.decorators.func({
    name: 'RCSB bCIF',
    description: 'Get biostructure by id as BinaryCIF',
    meta: {
      dataProvider: 'Molecule3D',
      cache: 'client',
      cacheInvalidateOn: '0 * * * *',
    },
  })
  static async getBiostructureRcsbBcif(
    @grok.decorators.param({options: {initialValue: '\'1QBS\''}}) id: string
  ): Promise<string> {
    const url = `https://models.rcsb.org/${id}.bcif`;
    const response = await fetch(url);
    if (!response.ok)
      throw new Error(response.statusText);
    const data: ArrayBuffer = await response.arrayBuffer();
    const dataA = new Uint8Array(data, 0, data.byteLength);
    return BiostructureDataJson.fromData({binary: true, data: dataA, ext: 'bcif', options: {name: id}});
  }

  @grok.decorators.func({
    description: 'Packs BiostructureData value into JSON string',
  })
  static biostructureDataToJson(
    binary: boolean,
    @grok.decorators.param({type: 'object'}) data: string | Uint8Array,
    ext: string,
    @grok.decorators.param({type: 'map', options: {optional: true}}) options?: object
  ): string {
    return BiostructureDataJson.fromData({binary, data, ext, options});
  }

  @grok.decorators.panel({
    name: '3D Structure',
    meta: {role: 'widgets', domain: 'bio'},
  })
  static structure3D(
    @grok.decorators.param({options: {semType: 'Molecule3D'}}) molecule: DG.SemanticValue
  ): DG.Widget {
    const widget = new DG.Widget(ui.div([]));
    widget.root.append(ui.loader());
    const {dataFrame, column, rowIndex} = molecule.cell;
    const tableView = grok.shell.getTableView(dataFrame.name);
    const {grid} = tableView;
    const gridCell = grid.cell(column.name, rowIndex);
    const renderer = new PdbGridCellRendererBack(null, column);

    renderer.createViewer(gridCell, tableView).then(async ({tview, viewer}) => {
      if (tview && viewer) {
        viewer.root.classList.add('bsv-container-info-panel');
        ui.empty(widget.root);
        widget.root.appendChild(viewer.root);
      }
    });
    return widget;
  }

  @grok.decorators.func({
    'name': 'Fetch PDB Sequences',
    //eslint-disable-next-line max-len
    'description': 'For a user-selected table and PDB ID column, fetches protein sequences and adds them as new columns.',
    'top-menu': 'Bio | Transform | Fetch PDB Sequences...'
  })
  static async fetchSequencesFromPdb(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'PDB_ID'}}) pdbColumn: DG.Column
  ): Promise<void> {
    await extractSequenceColumns(pdbColumn);
  }
}


