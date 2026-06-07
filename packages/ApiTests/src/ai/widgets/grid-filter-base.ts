import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {addAndGetFilter, demog, expectNoThrow, withFilterGroup} from '../helpers';

// DG.GridFilterBase: saveState/applyState, filterColumnName/filterType, isActive/isFiltering,
// setSorting, caption. Built-in filters surface as GridFilterBase instances in FilterGroup.filters.
category('AI: Widgets: GridFilterBase', () => {
  test('filterColumnName + filterType reflect the added filter', async () => {
    await withFilterGroup(demog(), async (fg) => {
      const f = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      expect(f.filterColumnName, 'sex');
      expect(typeof f.filterType, 'string');
      expect(f.filterType.length > 0, true);
    });
  });

  test('saveState / applyState round-trip', async () => {
    await withFilterGroup(demog(), async (fg) => {
      const f = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      const state = f.saveState();
      expect(state != null, true);
      expect(typeof state, 'object');
      expectNoThrow(() => f.applyState(state));
    });
  });

  test('isActive get/set toggles', async () => {
    await withFilterGroup(demog(), async (fg) => {
      const f = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      const original = f.isActive;
      f.isActive = false;
      expect(f.isActive, false);
      f.isActive = true;
      expect(f.isActive, true);
      f.isActive = original;
    });
  });

  test('isFiltering is false on a fresh filter with no excluded categories', async () => {
    await withFilterGroup(demog(), async (fg) => {
      const f = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      // No category deselected yet, so the filter passes every row.
      expect(f.isFiltering, false);
    });
  });

  test('caption is a readable string and the setter does not throw', async () => {
    await withFilterGroup(demog(), async (fg) => {
      const f = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      // A built-in categorical filter derives its caption from the column, so the setter does not
      // override the readback; assert it's a non-empty string and assigning is a no-throw operation.
      expect(typeof f.caption, 'string');
      expect(f.caption.length > 0, true);
      expectNoThrow(() => {
        f.caption = 'Gender';
      });
    });
  });

  test('setSorting by name and by count does not throw', async () => {
    await withFilterGroup(demog(), async (fg) => {
      const f = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      expectNoThrow(() => {
        f.setSorting('name', true);
        f.setSorting('count', false);
      });
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
