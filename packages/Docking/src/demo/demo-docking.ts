import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BINDING_ENERGY_COL, CACHED_DOCKING, POSE_COL, _package, setPose} from '../utils/constants';
import { addColorCoding, prepareAutoDockData } from '../package';
import { delay } from '@datagrok-libraries/utils/src/test';

export async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export async function prepareDockingData(
  table: DG.DataFrame, resultsPath: string, target: string, smilesCol: string, poses: number
): Promise<void> {
  const autodockResults = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(resultsPath));
  const dockingData = await prepareAutoDockData(target, table, smilesCol, poses);
  //@ts-ignore
  CACHED_DOCKING.K.push(dockingData);
  //@ts-ignore
  CACHED_DOCKING.V.push(autodockResults);
}

export async function _demoDocking(): Promise<void> {
  let tv = await openMoleculeDataset('System:AppData/Docking/demo_files/demo_dataset_small.csv');
  tv = (grok.shell.view('Browse')! as DG.BrowseView)!.preview! as DG.TableView;
  grok.shell.tv.dataFrame = tv.dataFrame;
  const grid = tv.grid;
  const table = tv.dataFrame;
  const desirableHeight = 100;
  const desirableWidth = 100;
  addColorCoding(grid.columns.byName(BINDING_ENERGY_COL)!);
  grid.sort([BINDING_ENERGY_COL]);
  await grok.data.detectSemanticTypes(table);
  grid.onCellRender.subscribe((args: any) => {
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

  await prepareDockingData(table, 'System:AppData/Docking/demo_files/autodock_results.csv', 'BACE1', 'SMILES', 10);
  table.currentCell = table.cell(0, POSE_COL);
}