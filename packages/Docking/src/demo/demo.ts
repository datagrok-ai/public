import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../utils/constants';

export function showHelpPanel(): void {
  grok.shell.windows.showHelp = true;
  grok.shell.windows.help.showHelp('/help/develop/domains/chem/docking');
  grok.shell.windows.context.visible = true;
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showProperties = true;
}

export async function _demoDocking(): Promise<void> {
  // semType isn't stored in layout JSON — pin it at import time so PowerGrid's
  // rawPng renderer and BSV's PL Diagram object handler wire up on first render.
  const df = DG.DataFrame.fromCsv(
    await _package.files.readAsText('demo_files/docking_demo.csv'),
    {columnImportOptions: [
      {name: 'AutoDock poses', semType: DG.SEMTYPE.MOLECULE3D},
      {name: 'PL Diagram', semType: 'rawPng'},
    ]},
  );

  const layout = DG.ViewLayout.fromJson(
    await _package.files.readAsText('demo_files/docking_demo.layout'));
  layout.columns.forEach((c) => {
    const col = df.col(c.name);
    if (col)
      Object.entries(c.tags).forEach(([k, v]) => col.setTag(k, v));
  });
  await df.meta.detectSemanticTypes();

  const tv = grok.shell.addTableView(df);
  await DG.delay(100);
  tv.loadLayout(layout, true);
  df.currentCell = df.cell(0, 'AutoDock poses');

  showHelpPanel();
}
