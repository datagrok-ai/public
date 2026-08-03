import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolGetSet, expectNoThrow, until, withTableView} from '../helpers';

function hasViewer(tv: DG.TableView, type: string): boolean {
  for (const v of tv.viewers)
    if (v.type === type)
      return true;
  return false;
}

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
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && hasViewer(tv, DG.VIEWER.HISTOGRAM));
      const beforeReset = tv.viewers.length;

      const s = tv.saveState();
      expect(typeof s, 'string');
      expect(s.length > 0, true);
      expectNoThrow(() => JSON.parse(s));

      // resetLayout() drops the added viewers; baseline keeps the grid only.
      tv.resetLayout();
      await until(() => !hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && !hasViewer(tv, DG.VIEWER.HISTOGRAM));

      // loadState must rebuild the exact arrangement: both viewer types return and the total
      // viewer count is restored to what it was before the reset.
      tv.loadState(s);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && hasViewer(tv, DG.VIEWER.HISTOGRAM));
      expect(hasViewer(tv, DG.VIEWER.SCATTER_PLOT), true);
      expect(hasViewer(tv, DG.VIEWER.HISTOGRAM), true);
      await until(() => tv.viewers.length === beforeReset);
      expect(tv.viewers.length, beforeReset);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
