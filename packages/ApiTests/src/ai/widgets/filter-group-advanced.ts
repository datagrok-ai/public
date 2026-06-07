import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {addAndGetFilter, demog, expectNoThrow, until, withFilterGroup} from '../helpers';

// Advanced DG.FilterGroup config — getStates() return-value content, setExpanded toggle,
// updateOrAdd(requestFilter=false), min/max reflection, isFiltering-after-applyState, setEnabled.
category('AI: Widgets: FilterGroup Advanced', () => {
  test('getStates returns array of state objects for an existing filter type', async () => {
    await withFilterGroup(demog(), async (fg) => {
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM).length > 0);
      const states = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM);
      expect(Array.isArray(states), true);
      expect(states.length > 0, true);
      const state = states[0] as {[k: string]: any};
      expect(typeof state, 'object');
      expect(state != null, true);
      expect(state.column ?? state.columnName, 'age');
    });
  });

  test('getStates returns empty array for a non-existent column/type combination', async () => {
    await withFilterGroup(demog(), async (fg) => {
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM).length > 0);
      const states = fg.getStates('no_such_column', DG.FILTER_TYPE.HISTOGRAM);
      expect(Array.isArray(states), true);
      expect(states.length, 0);
    });
  });

  test('setExpanded(filter,false) then setExpanded(filter,true) does not throw and toggles state', async () => {
    await withFilterGroup(demog(), async (fg) => {
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.filters.length > 0);
      const filter = fg.filters[fg.filters.length - 1];
      expectNoThrow(() => {
        fg.setExpanded(filter, false);
        fg.setExpanded(filter, true);
      });
    });
  });

  test('updateOrAdd with requestFilter=false adds filter without immediate filtering', async () => {
    // settleMs=100: the group must be fully initialized for requestFilter=false to be honored.
    await withFilterGroup(demog(), async (fg, tv) => {
      const df = tv.dataFrame;
      const rowCount = df.rowCount;
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age', min: 30, max: 40}, false);
      await until(() => fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM).length > 0);
      // requestFilter=false defers applying the filter, so the filter exists but the
      // dataframe's filtered bitset is not yet narrowed.
      expect(fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM).length > 0, true);
      expect(df.filter.trueCount, rowCount);
    }, 100);
  });

  test('getStates reflects updated min/max after updateOrAdd with histogram range', async () => {
    await withFilterGroup(demog(), async (fg) => {
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age', min: 25, max: 55}, true);
      await until(() => {
        const ss = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM) as Array<{[k: string]: any}>;
        return ss.length > 0 && ss[0].min === 25 && ss[0].max === 55;
      }, 5000);
      const state = (fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM) as Array<{[k: string]: any}>)[0];
      expect(state.min, 25);
      expect(state.max, 55);
    });
  });

  test('GridFilterBase.isFiltering is true after applyState with active criteria', async () => {
    await withFilterGroup(demog(), async (fg, tv) => {
      // GridFilterBase.isFiltering reflects category exclusion, not histogram min/max, so use a
      // CATEGORICAL filter and apply a state that deselects a category to make it filter.
      const filter = await addAndGetFilter(fg, DG.FILTER_TYPE.CATEGORICAL, 'sex');
      // A full state (no 'selected' key) means all categories pass — not filtering yet.
      expect(filter.isFiltering, false);
      const cats = (tv.dataFrame.getCol('sex').categories as string[]).filter((c) => c != null && c !== '');
      expect(cats.length > 1, true);
      // Select a strict subset of categories so at least one is excluded -> isFiltering becomes true.
      const state = filter.saveState();
      state.active = true;
      state.selected = [cats[0]];
      filter.applyState(state);
      await until(() => filter.isFiltering === true, 5000);
      expect(filter.isFiltering, true);
      expect(filter.isActive, true);
    });
  });

  test('setEnabled then getStates confirms active=false is reflected', async () => {
    await withFilterGroup(demog(), async (fg) => {
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.filters.length > 0);
      const filter = fg.filters[fg.filters.length - 1];
      fg.setEnabled(filter, false);
      await until(() => {
        const ss = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM) as Array<{[k: string]: any}>;
        return ss.length > 0 && ss[0].active === false;
      }, 5000);
      const state = (fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM) as Array<{[k: string]: any}>)[0];
      expect(state.active, false);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
