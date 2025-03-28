import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {_package} from '../utils/constants';

export async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export function showHelpPanel() {
  grok.shell.windows.showContextPanel = true;
   grok.shell.windows.showHelp = true;
  grok.shell.windows.help.showHelp('/help/develop/domains/chem/docking');
  
  grok.shell.windows.context.visible = true;  
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showProperties = true;
}

export async function _demoDocking(): Promise<void> {
  await demo('docking', ['Autodock poses']);
  showHelpPanel();
}

async function demo(type: "docking", columnNames: string[]): Promise<void> {
  const datasetPath = `System:AppData/Docking/demo_files/${type}_demo.csv`;
  const layoutPath = `System:AppData/Docking/demo_files/${type}_demo.layout`;

  const tv: DG.TableView = await openMoleculeDataset(datasetPath);
  const layout = DG.ViewLayout.fromJson(await grok.dapi.files.readAsText(layoutPath));
  tv.loadLayout(layout);

  const { dataFrame } = tv;
  await grok.data.detectSemanticTypes(dataFrame);

  for (let i = 0; i < columnNames.length; i++) {
    const name = columnNames[i];
    const column = dataFrame.getCol(name);
    column.semType = DG.SEMTYPE.MOLECULE3D;
    column.meta.units = 'pdb';
    column.setTag('docking.role', 'ligand');
  }

  dataFrame.currentCell = dataFrame.cell(0, columnNames[0]);
} 
