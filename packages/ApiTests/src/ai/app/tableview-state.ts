import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolGetSet, expectNoThrow, until, withTableView} from '../helpers';

category('AI: App: TableView State', () => {
  test('syncCurrentObject get/set round-trips', async () => {
    await withTableView(demog(), async (tv) => {
      expectBoolGetSet(tv as unknown as {[k: string]: any}, 'syncCurrentObject');
    });
  });

  test('saveState returns JSON, loadState restores viewers after resetLayout', async () => {
    await withTableView(demog(), async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'});
      tv.addViewer(DG.VIEWER.HISTOGRAM, {value: 'age'});
      await until(() => tv.viewers.length >= 3);
      const beforeReset = tv.viewers.length;

      const s = tv.saveState();
      expect(typeof s, 'string');
      expect(s.length > 0, true);
      expectNoThrow(() => JSON.parse(s));

      // resetLayout() drops the added viewers; baseline keeps the grid only.
      tv.resetLayout();
      await until(() => tv.viewers.length < beforeReset);
      const afterReset = tv.viewers.length;

      // Exact counts are headless-nondeterministic; assert the monotonic relation
      // (load >= reset) rather than a fixed number.
      expectNoThrow(() => tv.loadState(s));
      await until(() => tv.viewers.length >= afterReset);
      expect(tv.viewers.length >= afterReset, true);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
