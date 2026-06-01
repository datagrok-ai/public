import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, subscribeAll, until, wait, withTableView} from '../helpers';

// The six event observables on the built-in filter panel (DG.FilterGroup).
category('AI: Viewers: Filters Events', () => {
  test('onFilterAdded fires when a filter is added', async () => {
    await withTableView(demog(), async (tv) => {
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      await wait(100);
      await expectFiresWithin(fg.onFilterAdded,
        () => fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'}), 3000);
    });
  });

  test('onFilterRemoved fires when a filter is removed', async () => {
    await withTableView(demog(), async (tv) => {
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      await wait(100);
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.filters.length > 0);
      expect(fg.filters.length > 0, true);
      const filter = fg.filters[fg.filters.length - 1];
      await expectFiresWithin(fg.onFilterRemoved, () => fg.remove(filter), 3000);
    });
  });

  test('onFilterEnabledChanged fires when a filter is disabled', async () => {
    await withTableView(demog(), async (tv) => {
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      await wait(100);
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.filters.length > 0);
      expect(fg.filters.length > 0, true);
      const filter = fg.filters[fg.filters.length - 1];
      await expectFiresWithin(fg.onFilterEnabledChanged, () => fg.setEnabled(filter, false), 3000);
    });
  });

  test('onFilterCriteriaChanged fires when filter criteria change', async () => {
    await withTableView(demog(), async (tv) => {
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      await wait(100);
      fg.add({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
      await until(() => fg.filters.length > 0);
      await expectFiresWithin(fg.onFilterCriteriaChanged,
        () => fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age', min: 20, max: 60}, true), 5000);
    });
  });

  test('onFilterSync and onFormulaFilterChanged are well-formed rxjs Observables', async () => {
    await withTableView(demog(), async (tv) => {
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      await wait(100);
      subscribeAll([fg.onFilterSync, fg.onFormulaFilterChanged])();
    });
  });

  test('all six event getters are well-formed (subscribe is a function)', async () => {
    await withTableView(demog(), async (tv) => {
      const fg = tv.getFiltersGroup({createDefaultFilters: false});
      await wait(100);
      subscribeAll([fg.onFilterAdded, fg.onFilterRemoved, fg.onFilterCriteriaChanged,
        fg.onFilterEnabledChanged, fg.onFilterSync, fg.onFormulaFilterChanged])();
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
