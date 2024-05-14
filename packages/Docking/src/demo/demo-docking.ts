import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BINDING_ENERGY_COL, CACHED_DOCKING, POSE_COL, TARGET_PATH, _package, setPose} from '../utils/constants';
import { BiostructureData } from '@datagrok-libraries/bio/src/pdb/types';
import { AutoDockDataType } from '../apps/auto-dock-app';
import { addColorCoding } from '../package';
import { delay } from '@datagrok-libraries/utils/src/test';

async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export async function _demoDocking(): Promise<void> {
  const tv = await openMoleculeDataset('System:AppData/Docking/demo_files/demo_dataset_small.csv');
  const grid = tv.grid;
  const table = tv.dataFrame;
  const desirableHeight = 100;
  const desirableWidth = 100;
  addColorCoding(grid.columns.byName(BINDING_ENERGY_COL)!);
  grid.sort([BINDING_ENERGY_COL]);
  await grok.data.detectSemanticTypes(table);
  grid.onCellRender.subscribe((args) => {
    grid.setOptions({ 'rowHeight': desirableHeight });
    grid.col(POSE_COL)!.width = desirableWidth;
    grid.col(BINDING_ENERGY_COL)!.width = desirableWidth + 50;
  });

  await delay(100);
  const layoutString = await grok.dapi.files.readAsText('System:AppData/Docking/demo_files/demo.layout');
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);

  table.col(POSE_COL)!.setTag('docking.role', 'ligand');
  grid.invalidate();
  setPose(POSE_COL);
  const autodockResults = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText('System:AppData/Docking/demo_files/autodock_results.csv'));
  const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(TARGET_PATH, true);
  const receptor = targetsFiles.filter(file => file.isFile).find(file => file.name === 'BACE1.pdbqt');
  const receptorData: BiostructureData = {
    binary: false,
    data: (await grok.dapi.files.readAsText('System:AppData/Docking/targets/BACE1/BACE1.pdbqt')),
    ext: receptor!.extension,
    options: {name: receptor!.name,},
  };
  const data: AutoDockDataType = {
    ligandDf: table,
    ligandMolColName: 'SMILES',
    receptor: receptorData,
    gpfFile: (await grok.dapi.files.readAsText('System:AppData/Docking/targets/BACE1/BACE1.gpf')),
    confirmationNum: 10,
    ligandDfString: table.columns.byName('SMILES').toString(),
  };
  //@ts-ignore
  CACHED_DOCKING.K.push(data);
  //@ts-ignore
  CACHED_DOCKING.V.push(autodockResults);

  table.currentCell = table.cell(0, POSE_COL);
}