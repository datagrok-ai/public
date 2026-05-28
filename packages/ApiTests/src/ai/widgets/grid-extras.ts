// DG.Grid internals — core/client/d4/lib/src/viewers/grid/*.dart, grid/grid_column*.dart (scenario: grid-extras)
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, withTableView, wait, expectNoThrow} from '../helpers';

category('AI: Widgets: Grid Extras', () => {
  test('grid cell access + identity', async () => {
    const df = demog();
    const colName = 'age';
    await withTableView(df, async (tv) => {
      const grid = tv.grid;
      await wait();
      const gc = grid.cell(colName, 0);
      expect(gc instanceof DG.GridCell, true);
      expect(gc.gridRow, 0);
      expect(gc.isTableCell, true);
      expect(gc.isColHeader, false);
      expect(gc.isRowHeader, false);
      expect(gc.tableRowIndex, 0);
      expect(gc.gridColumn.name, colName);
      expect(gc.tableColumn != null, true);
      expect(gc.tableColumn!.name, colName);
      expect(gc.cell instanceof DG.Cell, true);
      expect(gc.value, df.get(colName, 0));
      expect(gc.cell.value, df.get(colName, 0));
    });
  });

  test('table<->grid row index mapping', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const grid = tv.grid;
      await wait();
      for (const i of [0, 1, 5, 10]) {
        const g = grid.tableRowToGrid(i);
        expect(typeof g, 'number');
        expect(grid.gridRowToTable(g), i);
      }
    });
  });

  test('GridColumnList byName/byIndex/length + GridColumn.idx alignment', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const grid = tv.grid;
      await wait();
      const cols = grid.columns;
      expect(cols instanceof DG.GridColumnList, true);
      expect(cols.length > df.columns.length, true);
      const colName = 'age';
      const byName = cols.byName(colName);
      expect(byName != null, true);
      expect(byName!.name, colName);
      expect(byName!.column != null, true);
      expect(byName!.column!.name, colName);
      const idx = byName!.idx;
      expect(typeof idx, 'number');
      expect(cols.byIndex(idx)!.name, colName);
      const viaCol = grid.col(colName);
      expect(viaCol != null, true);
      expect(viaCol!.name, colName);
      const rowHeader = cols.rowHeader;
      expect(rowHeader != null, true);
      expect(rowHeader!.idx, 0);
    });
  });

  test('column move reorders idx', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const grid = tv.grid;
      await wait();
      const a = grid.col('age')!;
      const target = grid.col('race')!.idx;
      expect(a.idx !== target, true);
      a.move(target);
      await wait();
      const moved = grid.col('age')!;
      expect(moved.idx, target);
      expect(grid.columns.byIndex(target)!.name, 'age');
    });
  });

  test('GridColumn.getDataWidth + checkEditable', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const grid = tv.grid;
      await wait();
      const col = grid.col('age')!;
      const dataWidth = col.getDataWidth();
      expect(typeof dataWidth, 'number');
      expect(dataWidth >= 0, true);
      const reason = col.checkEditable();
      expect(reason == null || typeof reason === 'string', true);
      if (col.editable)
        expect(reason == null, true);
    });
  });

  test('visible-cell enumeration + invalidate/scrollToCell smoke', async () => {
    const df = demog(1);
    await withTableView(df, async (tv) => {
      const grid = tv.grid;
      await wait();
      let gridCellCount = 0;
      for (const gc of grid.getVisibleCells()) {
        expect(gc instanceof DG.GridCell, true);
        expect(gc.gridColumn != null, true);
        gridCellCount++;
      }
      expect(gridCellCount >= 0, true);
      const col = grid.col('age')!;
      let colCellCount = 0;
      for (const gc of col.getVisibleCells()) {
        expect(gc instanceof DG.GridCell, true);
        expect(gc.gridColumn != null, true);
        colCellCount++;
      }
      expect(colCellCount >= 0, true);
      expectNoThrow(() => grid.invalidate());
      expectNoThrow(() => grid.scrollToCell('age', 0));
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
