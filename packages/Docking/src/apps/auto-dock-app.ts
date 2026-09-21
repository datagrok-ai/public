import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BiostructureData} from '@datagrok-libraries/bio/src/pdb/types';
import {
  IAutoDockService, getAutoDockService, GridSize,
} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {IPdbHelper, getPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

import {_package, BINDING_ENERGY_COL, POSE_COL, ERROR_COL_NAME} from '../utils/constants';
import {buildDefaultAutodockGpf} from '../utils/auto-dock-service';

export type AutoDockDataType = {
  ligandDf: DG.DataFrame,
  ligandMolColName: string,
  receptor: BiostructureData,
  gpfFile?: string,
  posesNum?: number,
  ligandDfString?: string;
};

export class AutoDockApp {
  private readonly poseColName: string = 'pose';
  private data!: AutoDockDataType;

  async init(data: AutoDockDataType): Promise<DG.DataFrame | undefined> {
    this.data = data;
    return await this.getAutodockResults();
  }

  async getAutodockResults(): Promise<DG.DataFrame | undefined> {
    const pi = DG.TaskBarProgressIndicator.create('AutoDock running...');
    try {
      const ligandCol = this.data.ligandDf.getCol(this.data.ligandMolColName);
      const result = await runAutoDock(this.data.receptor, ligandCol, this.data.gpfFile!, this.data.posesNum!, this.poseColName, pi);
      let posesAllDf = result?.posesAllDf;
      const errorValues = result?.errorValues;
      if (!posesAllDf) {
        posesAllDf = DG.DataFrame.create();
        posesAllDf.columns.addNewString(POSE_COL);
      }

      errorValues?.forEach(({index, value}) => {
        posesAllDf.rows.insertAt(index, 1);
        posesAllDf.set(POSE_COL, index, value);
      });
      return posesAllDf;
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      grok.shell.error(errMsg);
      _package.logger.error(errMsg, undefined, errStack);
    } finally {
      pi.close();
    }
  }
}

// -- Routines --
async function runAutoDock(
  receptorData: BiostructureData, ligandCol: DG.Column<string>, gpfFile: string, posesNum: number, poseColName: string, pi: DG.ProgressIndicator
) {
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
  const errorValues: { index: number, value: string }[] = []

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
      receptorData, ligandData, gpfFile ?? autodockGpf, posesNum ?? 10, poseColName);

    if (posesDf.col(ERROR_COL_NAME)) {
      errorValues[errorValues.length] = {index: lRowI, value: posesDf.get(ERROR_COL_NAME, 0)};
      continue;
    }

    posesDf.rows.removeWhere((row) => row.get(BINDING_ENERGY_COL) !== posesDf.col(BINDING_ENERGY_COL)?.min);
    if (posesDf!.rowCount > 1) {
      posesDf!.rows.removeAt(1, posesDf!.rowCount - 1);
    }
    // region: add extra columns to AutoDock output

    const colNames = posesDf.columns.names();
    const nameWithExtension = receptorData.options?.name;
    const nameWithoutExtension = nameWithExtension?.replace(/\.[^/.]+$/, '');
    let remarkString = `REMARK   1 receptor.    ${nameWithoutExtension} J.\n`;
    for (let i = 1; i < colNames.length; ++i) {
      if (colNames[i] === poseColName) continue;
      const rowValue = posesDf.get(colNames[i], 0);
      remarkString += `REMARK   ${i + 1} ${colNames[i]}.    ${rowValue} J.`;
      if (i !== colNames.length - 1) remarkString += '\n';
    }

    const currentPosesValue = posesDf.get(poseColName, 0);
    const newPosesValue = currentPosesValue.replace(
      /^COMPND.*\n/,
      (match: string) => match + `${remarkString}\n`
    );
    posesDf.set(poseColName, 0, newPosesValue);

    const pdbqtCol = posesDf.getCol(poseColName);
    pdbqtCol.name = poseColName;
    pdbqtCol.semType = DG.SEMTYPE.MOLECULE3D;

    // endregion: add extra columns to AutoDock output

    if (posesAllDf === undefined)
      posesAllDf = posesDf;
    else
      posesAllDf!.append(posesDf, true);

    pi.update(100 * lRowI / ligandRowCount, 'AutoDock running...');
  }

  return {posesAllDf, errorValues};
}
