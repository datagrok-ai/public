import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../utils/constants';

export async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export function showHelpPanel(): void {
  grok.shell.windows.showHelp = true;
  grok.shell.windows.help.showHelp('/help/develop/domains/chem/docking');
  grok.shell.windows.context.visible = true;
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showProperties = true;
}

export async function _demoDocking(): Promise<void> {
  await demo('docking', ['AutoDock poses']);
  showHelpPanel();
}

async function demo(type: 'docking', columnNames: string[]): Promise<void> {
  const datasetPath = `System:AppData/Docking/demo_files/${type}_demo.csv`;
  const layoutPath = `System:AppData/Docking/demo_files/${type}_demo.layout`;

  // semType isn't stored in layout JSON — pin Molecule (SMILES) and rawPng
  // (PL Diagram, BSV PL object handler) at import time. Ligand columns are
  // pinned below in the same loop that adds their docking-role tag.
  const df = DG.DataFrame.fromCsv(
    await grok.dapi.files.readAsText(datasetPath),
    {columnImportOptions: [
      {name: 'SMILES', semType: DG.SEMTYPE.MOLECULE},
      {name: 'PL Diagram', semType: 'rawPng'},
    ]},
  );

  const layout = DG.ViewLayout.fromJson(await grok.dapi.files.readAsText(layoutPath));
  layout.columns.forEach((c) => {
    const col = df.col(c.name);
    if (col)
      Object.entries(c.tags).forEach(([k, v]) => col.setTag(k, v));
  });
  await df.meta.detectSemanticTypes();

  for (const name of columnNames) {
    const column = df.getCol(name);
    column.semType = DG.SEMTYPE.MOLECULE3D;
    column.meta.units = 'pdb';
    column.setTag('docking.role', 'ligand');
  }

  const tv = grok.shell.addTableView(df);
  await DG.delay(100);
  tv.loadLayout(layout, true);
  df.currentCell = df.cell(0, columnNames[0]);
}
