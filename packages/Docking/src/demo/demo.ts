import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import { delay } from '@datagrok-libraries/utils/src/test';

import {BINDING_ENERGY_COL, POSE_COL, _package, setPose} from '../utils/constants';
import { addColorCoding } from '../utils/utils';

export async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export function showHelpPanel(info: string) {
  grok.shell.windows.help.visible = true;
  grok.shell.windows.help.showHelp(ui.markdown(info));
    
  grok.shell.windows.context.visible = true;  
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showProperties = true;
  grok.shell.windows.help.visible = true;
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
  table.currentCell = table.cell(0, POSE_COL);
}

export async function _demoBoltzFolding(): Promise<void> {
  await boltzDemo('folding');
  showHelpPanel(BOLTZ_HELP_FOLDING);
}

export async function _demoBoltzDocking(): Promise<void> {
  await boltzDemo('docking');
  showHelpPanel(BOLTZ_HELP_DOCKING);
}

async function boltzDemo(type: string): Promise<void> {
  const tv: DG.TableView = await openMoleculeDataset(`System:AppData/Docking/demo_files/boltz_${type}_demo.csv`);
  const layoutString = await grok.dapi.files.readAsText(`System:AppData/Docking/demo_files/boltz_${type}_demo.layout`);
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);

  const {dataFrame} = tv;
  dataFrame.getCol('pdb').semType = DG.SEMTYPE.MOLECULE3D;
  dataFrame.currentCell = dataFrame.cell(0, 'pdb');
}

// Constants
export const BOLTZ_HELP_FOLDING = `# Folding Demo

This demo showcases **Boltz-1**, a state-of-the-art open-source model to predict biomolecular structures containing combinations of proteins, RNA, DNA, and other molecules.

To run folding for the entire column:
  - Navigate to the top menu and select **Bio | Folding | Boltz-1**.
  - This will apply **Boltz-1** to all sequences in the selected column.
  
Click on any predicted structure to view the details, including confidence score, ptm, iptm and other relevant characteristics.`

export const BOLTZ_HELP_DOCKING = `# Docking Demo

This demo showcases **Boltz-1**, a state-of-the-art open-source model to predict biomolecular structures containing combinations of proteins, RNA, DNA, and other molecules.

To run docking for the entire column:
  - Navigate to the top menu and select **Chem | Docking | Boltz-1**.
  - This will apply **Boltz-1** to all molecules in the selected column.
  
Click on any predicted structure to view the details, including confidence score, ptm, iptm and other relevant characteristics.`