import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BiostructureDataJson, BiostructureData,} from '@datagrok-libraries/bio/src/pdb/types';
import {BiostructureProps, IBiostructureViewer,} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {
  IAutoDockService, getAutoDockService, GridSize,
} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {delay} from '@datagrok-libraries/utils/src/test';

import {PROPS as msvPROPS} from '../viewers/molstar-viewer/molstar-viewer';

import {_package, getPdbHelper} from '../package';


export type AutoDockDataType = {
  ligandDf: DG.DataFrame,
  ligandMolColName: string,
  receptor: BiostructureData,
};
type TargetViewerType = DG.Viewer<BiostructureProps> & IBiostructureViewer;

export class AutoDockApp {
  private readonly appFuncName: string;
  private readonly poseColName: string = 'pose';

  constructor(appFuncName: string = 'autoDockApp') {
    this.appFuncName = appFuncName;
  }

  async init(data?: AutoDockDataType): Promise<void> {
    let v = data;
    if (!v)
      v = await AutoDockApp.loadData();

    await this.setData(v);
  }

  static async loadData(): Promise<AutoDockDataType> {
    let ligandFi!: DG.FileInfo;
    let targetFi!: DG.FileInfo;
    await Promise.all([
      (async () => { ligandFi = (await _package.files.list('samples', false, '1bdq.sdf'))[0]; })(),
      (async () => { targetFi = (await _package.files.list('samples', false, '1bdq-wo-ligands.pdb'))[0]; })(),
    ]);

    const ligandSdfA = await _package.files.readAsBytes(ligandFi);
    const ligandDf: DG.DataFrame = (await grok.functions.call('Chem:importSdf', {bytes: ligandSdfA}))[0];
    const ligandMolColName = 'molecule'; // importSdf generates

    const receptorData: BiostructureData = {
      binary: false, ext: targetFi.extension, data: await _package.files.readAsText(targetFi),
      options: {name: targetFi.name}
    };

    return {ligandDf, ligandMolColName, receptor: receptorData};
  }

  // -- Data --

  private data: AutoDockDataType;

  async setData(data: AutoDockDataType): Promise<void> {
    this.data = data;

    await this.buildView();
  }

  // -- View --

  private view: DG.TableView;
  private targetViewer: TargetViewerType;
  private posesGrid: DG.Grid;

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.data.ligandDf);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.appFuncName}`;

    const poseDf: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('pose', []),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'ligand0', []) /* ligand row index */,
    ]);

    this.posesGrid = await poseDf.plot.fromType(DG.VIEWER.GRID, {}) as DG.Grid;
    this.targetViewer = await poseDf.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: BiostructureDataJson.fromData(this.data.receptor),
      [msvPROPS.ligandColumnName]: this.poseColName,
    }) as unknown as TargetViewerType;
    const k = 42;

    const posesNode = this.view.dockManager.dock(this.posesGrid, DG.DOCK_TYPE.RIGHT, null, 'Poses', 0.5);
    this.view.dockManager.dock(this.targetViewer, DG.DOCK_TYPE.LEFT, posesNode, 'Target', 0.35);

    this.setRibbonPanels();

    const adSvc: IAutoDockService = await getAutoDockService();
    if (!adSvc.ready) {
      this.runBtn.disabled = true;
      ui.tooltip.bind(this.runBtn, 'AutoDock docker container is not ready');

      await adSvc.startDockerContainer();
      if (adSvc.ready) {
        grok.shell.info('AutoDock docker container is ready.');
        this.runBtn.disabled = false;
        ui.tooltip.bind(this.runBtn, 'Run AutoDock');
      }
    }
  }

  private runBtn: HTMLButtonElement;
  private downloadPosesBtn: HTMLButtonElement;

  setRibbonPanels(): void {
    const runIcon = ui.iconFA('play');
    runIcon.classList.add('fas');
    this.runBtn = ui.button(runIcon, this.runBtnOnClick.bind(this));
    ui.tooltip.bind(this.runBtn, 'Run AutoDock');

    const downloadPosesIcon = ui.iconFA('download');
    this.downloadPosesBtn = ui.button(downloadPosesIcon, this.downloadPosesBtnOnClick.bind(this));
    this.downloadPosesBtn.disabled = true;
    ui.tooltip.bind(this.downloadPosesBtn, 'Download poses (*.csv)');

    this.view.setRibbonPanels([
      [this.runBtn],
      [this.downloadPosesBtn],
    ]);
  }

  // -- Handle controls' events --

  /** Handles {@link runBtn} click */
  async runBtnOnClick(): Promise<void> {
    this.runBtn.disabled = true;
    const pi = DG.TaskBarProgressIndicator.create('AutoDock running...');
    try {
      // TODO: Convert receptor data to PDB string
      if (this.data.receptor.ext !== 'pdb' || typeof this.data.receptor.data !== 'string')
        throw new Error('Unsupported receptor data');
      const receptorPdb: string = this.data.receptor.data as string;

      const ligandCol = this.data.ligandDf.getCol(this.data.ligandMolColName);

      const posesAllDf = await runAutoDock(receptorPdb, ligandCol, this.poseColName, pi);
      if (posesAllDf !== undefined) {
        posesAllDf.getCol(this.poseColName).setTag('cell.renderer', 'xray');

        this.posesGrid.dataFrame = posesAllDf!;
        this.posesGrid.props.rowHeight = 150;
        this.posesGrid.columns.byName(this.poseColName)!.width = 150;

        this.targetViewer.dataFrame = posesAllDf!;
        // @ts-ignore
        this.targetViewer.setData('AutoDock app crutch');

        this.downloadPosesBtn.disabled = false;
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      grok.shell.error(errMsg);
      _package.logger.error(errMsg, undefined, errStack);
    } finally {
      pi.close();
      this.runBtn.disabled = false;
    }
  }

  /** Handles {@link downloadPosesBtn } click */
  async downloadPosesBtnOnClick(): Promise<void> {
    const posesCsv = this.posesGrid.dataFrame.toCsv();
    const posesCsvFn = `${this.data.receptor.options?.name ?? 'receptor'}-poses.csv`;
    DG.Utils.download(posesCsvFn, posesCsv);
  }
}

// -- Routines --
async function runAutoDock(
  receptorPdb: string, ligandCol: DG.Column<string>, poseColName: string, pi: DG.ProgressIndicator
): Promise<DG.DataFrame | undefined> {
  let adSvc!: IAutoDockService;
  try {
    adSvc = await getAutoDockService();
  } catch (err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.warning(errMsg);
    return;
  }

  const pdbHelper: IPdbHelper = await getPdbHelper();
  let posesAllDf: DG.DataFrame | undefined = undefined;

  const npts: GridSize = {x: 40, y: 40, z: 40};
  const ligandRowCount = ligandCol.length;
  for (let lRowI = 0; lRowI < ligandRowCount && lRowI < 3; ++lRowI) {
    const ligandMol = ligandCol.get(lRowI);
    // const ligandData: BiostructureData = {binary: false, ext: 'mol', data: ligandMol};
    const ligandPdb = await pdbHelper.molToPdb(ligandMol!);
    const ligandData: BiostructureData = {binary: false, ext: 'pdb', data: ligandPdb};

    const adRunRes = await adSvc.run(
      receptorPdb, ligandData, npts, 10, poseColName);

    if (posesAllDf === undefined) {
      posesAllDf = adRunRes.posesDf.clone(
        DG.BitSet.create(adRunRes.posesDf.rowCount, (_i) => false));
      for (let colI = 0; colI < posesAllDf!.columns.length; ++colI) {
        const srcCol = adRunRes.posesDf.columns.byIndex(colI);
        const tgtCol = posesAllDf!.columns.byIndex(colI);
        // for (const [tagName, tagValue] of Object.entries(srcCol.tags)) {
        //   tgtCol.setTag(tagName, tagValue);
        // }
      }
    }
    posesAllDf!.append(adRunRes.posesDf, true);
    pi.update(100 * lRowI / ligandRowCount, 'AutoDock running...');
  }

  return posesAllDf;
}
