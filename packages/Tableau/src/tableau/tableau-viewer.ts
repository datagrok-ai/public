import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TwbFile, TwbDatasource, TwbWorksheet} from './tableau-types';
import {twbDatasourceToDataFrame} from './tableau-to-dataframe';


function findDatasource(twbFile: TwbFile, ws: TwbWorksheet): TwbDatasource | null {
  return twbFile.datasources.find((ds) => ds.name === ws.datasourceName) || null;
}


function buildWorksheetPane(twbFile: TwbFile, ws: TwbWorksheet): HTMLElement {
  const ds = findDatasource(twbFile, ws);
  if (!ds) {
    return ui.divText(`Datasource '${ws.datasourceName}' not found.`);
  }

  const df = twbDatasourceToDataFrame(ds);
  const grid = DG.Viewer.grid(df);
  grid.root.style.width = '100%';
  grid.root.style.height = '100%';
  return grid.root;
}


export function buildTableauView(twbFile: TwbFile): HTMLElement {
  const headerItems: HTMLElement[] = [];
  if (twbFile.version)
    headerItems.push(ui.divText(`Tableau Workbook v${twbFile.version}`));
  headerItems.push(ui.divText(
    `${twbFile.datasources.length} datasource(s), ${twbFile.worksheets.length} worksheet(s)`
  ));
  const header = ui.div(headerItems, {style: {padding: '8px', borderBottom: '1px solid var(--grey-2)'}});

  const tabControl = ui.tabControl();
  for (const ws of twbFile.worksheets)
    tabControl.addPane(ws.name, () => buildWorksheetPane(twbFile, ws));
  tabControl.root.style.width = '100%';
  tabControl.root.style.flex = '1';

  const container = ui.divV([header, tabControl.root]);
  container.style.width = '100%';
  container.style.height = '100%';
  container.style.display = 'flex';
  container.style.flexDirection = 'column';
  return container;
}


export function openAllDatasources(twbFile: TwbFile): void {
  for (const ds of twbFile.datasources) {
    const df = twbDatasourceToDataFrame(ds);
    grok.shell.addTableView(df);
  }
}
