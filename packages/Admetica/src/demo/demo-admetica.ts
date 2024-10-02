import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export async function _demoAdmetica(): Promise<void> {
  let tv = await openMoleculeDataset('System:AppData/Admetica/demo_files/demo_dataset.csv');
  tv.dataFrame.name = 'Table';

  const layoutString = await grok.dapi.files.readAsText('System:AppData/Admetica/demo_files/demo.layout');
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);
  grok.shell.windows.showContextPanel = false;
}