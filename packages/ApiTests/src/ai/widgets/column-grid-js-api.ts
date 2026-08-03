import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectBoolGetSet} from '../helpers';

// DG.ColumnGrid — columnSelector/columnManager factories, passesFilter, checkAll, getChecked*,
// get/setSelectedColumns, showSearch.
category('AI: Widgets: ColumnGrid', () => {
  function typedDf(): DG.DataFrame {
    return df([
      ['a', DG.COLUMN_TYPE.INT, [1, 2, 3]],
      ['b', DG.COLUMN_TYPE.STRING, ['x', 'y', 'z']],
      ['c', DG.COLUMN_TYPE.FLOAT, [1.5, 2.5, 3.5]],
    ]);
  }

  // Build a ColumnGrid, run body, always detach+close.
  function withGrid(make: () => DG.ColumnGrid, body: (cg: DG.ColumnGrid) => void): void {
    const cg = make();
    try {
      body(cg);
    } finally {
      cg.detach();
      cg.close();
    }
  }

  test('columnSelector construction yields a ColumnGrid', async () => {
    withGrid(() => DG.ColumnGrid.columnSelector(demog()), (cg) => {
      expect(cg instanceof DG.ColumnGrid, true);
      expect(cg.dfSource instanceof DG.DataFrame, true);
    });
  });

  test('passesFilter respects a type filter', async () => {
    const source = typedDf();
    // The filter predicate is invoked from Dart with a raw Dart column handle (not a
    // wrapped DG.Column), so DG.toJs() is required for the .type getter to resolve.
    withGrid(() => DG.ColumnGrid.columnSelector(source, {filter: (c) => DG.toJs(c).type === DG.COLUMN_TYPE.STRING}),
      (cg) => {
        expect(cg.passesFilter(source.col('b')!), true);
        expect(cg.passesFilter(source.col('a')!), false);
      });
  });

  test('checkAll(true) checks every column', async () => {
    const source = typedDf();
    withGrid(() => DG.ColumnGrid.columnSelector(source), (cg) => {
      cg.checkAll(true);
      const names = cg.getCheckedColumnNames();
      for (const c of source.columns.names())
        expect(names.indexOf(c) >= 0, true);
      expect(names.length, source.columns.length);
    });
  });

  test('checkAll(false) clears all checks', async () => {
    withGrid(() => DG.ColumnGrid.columnSelector(typedDf(), {checkAll: true}), (cg) => {
      cg.checkAll(false);
      expect(cg.getCheckedColumnNames().length, 0);
    });
  });

  test('setSelectedColumns scopes the checked set to a column subset', async () => {
    // setSelectedColumns scopes the checked set by name; getSelectedColumns reads grid row selection
    // (a different state), so read back via the checked set instead.
    const source = typedDf();
    withGrid(() => DG.ColumnGrid.columnSelector(source), (cg) => {
      cg.setSelectedColumns(['a', 'c']);
      const checked = cg.getCheckedColumns();
      for (const c of checked)
        expect(c instanceof DG.Column, true);
      const names = cg.getCheckedColumnNames();
      expect(names.length, 2);
      expect(names.indexOf('a') >= 0, true);
      expect(names.indexOf('c') >= 0, true);
      expect(names.indexOf('b') >= 0, false);
    });
  });

  test('columnManager exposes dfSource', async () => {
    const source = demog();
    withGrid(() => DG.ColumnGrid.columnManager(source), (cg) => {
      expect(cg instanceof DG.ColumnGrid, true);
      expect(cg.dfSource instanceof DG.DataFrame, true);
      expect(cg.dfSource.rowCount, source.rowCount);
    });
  });

  test('showSearch get/set round-trips', async () => {
    withGrid(() => DG.ColumnGrid.columnSelector(demog()), (cg) => expectBoolGetSet(cg, 'showSearch'));
  });
}, {owner: 'agolovko@datagrok.ai'});
