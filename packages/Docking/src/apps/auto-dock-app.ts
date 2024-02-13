import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BiostructureDataJson, BiostructureData,} from '@datagrok-libraries/bio/src/pdb/types';
import {BiostructureProps, IBiostructureViewer,} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {
  IAutoDockService, getAutoDockService, GridSize,
} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {IPdbHelper, getPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {delay} from '@datagrok-libraries/utils/src/test';

import {_package, CACHED_DOCKING} from '../utils/constants';
import {buildDefaultAutodockGpf} from '../utils/auto-dock-service';

export type AutoDockDataType = {
  ligandDf: DG.DataFrame,
  ligandMolColName: string,
  receptor: BiostructureData,
  gpfFile?: string,
  confirmationNum?: number,
};
type TargetViewerType = DG.Viewer<BiostructureProps> & IBiostructureViewer;

export class AutoDockApp {
  private readonly appFuncName: string;
  private readonly poseColName: string = 'pose';

  constructor(appFuncName: string = 'autoDockApp') {
    this.appFuncName = appFuncName;
  }

  async init(data?: AutoDockDataType): Promise<DG.DataFrame | undefined> {
    let v = data;
    if (!v)
      v = await AutoDockApp.loadData();

    await this.setData(v);
    if (!!data)
      return await this.getAutodockResults();
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

  private data!: AutoDockDataType;

  async setData(data: AutoDockDataType): Promise<void> {
    this.data = data;

    await this.buildView();
  }

  // -- View --

  private view!: DG.TableView;
  private posesGrid!: DG.Grid;

  async buildView(): Promise<void> {
    this.view = grok.shell.tv;
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

  private runBtn!: HTMLButtonElement;
  private downloadPosesBtn!: HTMLButtonElement;

  setRibbonPanels(): void {
    const runIcon = ui.iconFA('play');
    runIcon.classList.add('fas');
    this.runBtn = ui.button(runIcon, this.getAutodockResults.bind(this));
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
  async getAutodockResults(): Promise<DG.DataFrame | undefined> {
    this.runBtn.disabled = true;
    const pi = DG.TaskBarProgressIndicator.create('AutoDock running...');
    try {
      const ligandCol = this.data.ligandDf.getCol(this.data.ligandMolColName);
      const posesAllDf = await runAutoDock(this.data.receptor, ligandCol, this.data.gpfFile!, this.data.confirmationNum!, this.poseColName, pi);
      if (posesAllDf !== undefined) {
        this.downloadPosesBtn.disabled = false;
        CACHED_DOCKING.set(this.data, posesAllDf);
        return posesAllDf;
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
  receptorData: BiostructureData, ligandCol: DG.Column<string>, gpfFile: string, confirmationNum: number, poseColName: string, pi: DG.ProgressIndicator
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

  const ligandRowCount = ligandCol.length;
  for (let lRowI = 0; lRowI < ligandRowCount; ++lRowI) {
    const ligandMol = ligandCol.semType === DG.SEMTYPE.MOLECULE 
      ? await grok.functions.call('Chem:convertMolNotation',
    {molecule: ligandCol.get(lRowI), sourceNotation: 'unknown', targetNotation: 'v3Kmolblock'})
      : ligandCol.get(lRowI);
    // const ligandData: BiostructureData = {binary: false, ext: 'mol', data: ligandMol};
    const ligandPdb = await pdbHelper.molToPdb(ligandMol!);
    const ligandData: BiostructureData = {binary: false, ext: 'pdb', data: ligandPdb};

    const npts: GridSize = {x: 40, y: 40, z: 40};
    const autodockGpf: string = buildDefaultAutodockGpf(receptorData.options!.name!, npts);
    const posesDf = await adSvc.dockLigand(
      receptorData, ligandData, gpfFile ?? autodockGpf, confirmationNum ?? 10, poseColName);
    
    posesDf.rows.removeWhere((row) => row.get('affinity') !== posesDf.col('affinity')?.min);
    // region: add extra columns to AutoDock output

    const pdbqtCol = posesDf.getCol(poseColName);
    pdbqtCol.name = poseColName;
    pdbqtCol.semType = DG.SEMTYPE.MOLECULE3D;
    pdbqtCol.setTag(DG.TAGS.UNITS, 'pdbqt');
    // endregion: add extra columns to AutoDock output

    if (posesAllDf === undefined) {
      posesAllDf = posesDf.clone(
        DG.BitSet.create(posesDf.rowCount, (_i) => false));
      for (let colI = 0; colI < posesAllDf!.columns.length; ++colI) {
        const srcCol = posesDf.columns.byIndex(colI);
        const tgtCol = posesAllDf!.columns.byIndex(colI);
        // for (const [tagName, tagValue] of Object.entries(srcCol.tags)) {
        //   tgtCol.setTag(tagName, tagValue);
        // }
      }
    }
    posesAllDf!.append(posesDf, true);
    pi.update(100 * lRowI / ligandRowCount, 'AutoDock running...');
  }

  return posesAllDf;
}
