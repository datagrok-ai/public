import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TwbFile, TwbWorksheet} from './tableau-types';
import {twbDatasourceToDataFrame} from './tableau-to-dataframe';


function buildWorksheetPane(ws: TwbWorksheet): HTMLElement {
  const items: HTMLElement[] = [];

  if (ws.datasourceName)
    items.push(ui.divText(`Datasource: ${ws.datasourceName}`));
  if (ws.markClass)
    items.push(ui.divText(`Mark type: ${ws.markClass}`));
  if (ws.rows)
    items.push(ui.divText(`Rows: ${ws.rows}`));
  if (ws.cols)
    items.push(ui.divText(`Cols: ${ws.cols}`));

  if (ws.usedColumns.length > 0) {
    items.push(ui.h3('Used Columns'));
    items.push(ui.list(ws.usedColumns));
  }

  const container = ui.divV(items, {style: {padding: '8px', overflow: 'auto'}});
  container.style.width = '100%';
  container.style.height = '100%';
  return container;
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
    tabControl.addPane(ws.name, () => buildWorksheetPane(ws));
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
