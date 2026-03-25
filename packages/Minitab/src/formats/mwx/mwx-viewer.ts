import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MwxWorksheet} from './mwx-types';
import {mwxWorksheetToDataFrame} from './mwx-to-dataframe';


/** Builds a preview for a parsed MWX worksheet. */
export function buildMwxView(ws: MwxWorksheet): HTMLElement {
  const df = mwxWorksheetToDataFrame(ws);

  const header = ui.divText(
    `${ws.name} \u2014 ${ws.rowCount} rows, ${ws.columns.length} columns`,
    {style: {padding: '8px', fontWeight: 'bold', borderBottom: '1px solid var(--grey-2)'}},
  );

  const grid = DG.Viewer.grid(df);
  grid.root.style.width = '100%';
  grid.root.style.flex = '1';

  const container = ui.divV([header, grid.root]);
  container.style.width = '100%';
  container.style.height = '100%';
  container.style.display = 'flex';
  container.style.flexDirection = 'column';
  return container;
}
