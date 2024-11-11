import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {delay, expect, testEvent} from '@datagrok-libraries/utils/src/test';
import {asRenderer, IRenderer, isRenderer} from '@datagrok-libraries/bio/src/types/renderer';

import {_package} from '../package-test';
import {CellRendererBackBase, getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';

export async function loadFileAsText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}

export async function readDataframe(tableName: string): Promise<DG.DataFrame> {
  const file = await loadFileAsText(tableName);
  const df = DG.DataFrame.fromCsv(file);
  df.name = tableName.replace('.csv', '');
  return df;
}

export async function createTableView(tableName: string): Promise<DG.TableView> {
  const df = await readDataframe(tableName);
  df.name = tableName.replace('.csv', '');
  const view = grok.shell.addTableView(df);
  return view;
}


/**
 * Tests if a table has non-zero rows and columns.
 *
 * @param {DG.DataFrame} table Target table. */
export function _testTableIsNotEmpty(table: DG.DataFrame): void {
  expect(table.columns.length > 0 && table.rowCount > 0, true);
}

export async function awaitGrid(grid: DG.Grid, timeout: number = 5000): Promise<void> {
  await delay(0);
  await testEvent(grid.onAfterDrawContent, () => {},
    () => { grid.invalidate(); }, timeout);

  const colCount = grid.columns.length;
  for (let colI = 0; colI < colCount; ++colI) {
    const gridCol = grid.columns.byIndex(colI);
    if (gridCol) {
      const gridCell = grid.cell(gridCol.name, 0);
      const [_gridCol, _tableCol, temp] =
        getGridCellColTemp<void, CellRendererBackBase<void>>(gridCell);

      const renderer = asRenderer(temp.rendererBack);
      if (renderer) await renderer.awaitRendered();
    }
  }
}
