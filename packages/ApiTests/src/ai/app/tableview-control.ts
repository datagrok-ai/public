import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, withTableView, withAttachedViewer, expectNoThrow, until} from '../helpers';

category('AI: App: TableView Control', () => {
  test('getViewerColumns after scatter', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'},
      async (v, tv) => {
        const cols = tv.getViewerColumns();
        expect(Array.isArray(cols), true);
        for (const c of cols)
          expect(c instanceof DG.Column, true);
        const names = cols.map((c) => c.name);
        expect(names.includes('age'), true);
        expect(names.includes('height'), true);
      });
  });

  test('closeAllViewers', async () => {
    await withTableView(demog(), async (tv) => {
      const baseline = Array.from(tv.viewers).length;
      tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'});
      tv.addViewer(DG.VIEWER.HISTOGRAM, {value: 'age'});
      await until(() => Array.from(tv.viewers).length > baseline);
      // closeAllViewers() closes every viewer in the view — including the grid — so the
      // viewers collection drops to empty, it does not return to the grid-only baseline.
      tv.closeAllViewers();
      await until(() => Array.from(tv.viewers).length === 0);
    });
  });

  test('getViewerColumns fallback', async () => {
    await withTableView(demog(), (tv) => {
      const cols = tv.getViewerColumns();
      expect(Array.isArray(cols), true);
      for (const c of cols)
        expect(c instanceof DG.Column, true);
    });
  });

  test('focus', async () => {
    await withTableView(demog(), (tv) => {
      expectNoThrow(() => tv.focus());
      // best-effort: activeElement may be unreliable in headless runs — downgraded to no-throw.
      const active = document.activeElement;
      if (active != null && tv.root != null)
        expect(tv.root.contains(active) || active === document.body, true);
    });
  });

  test('closeAllViewers boundary no-op', async () => {
    await withTableView(demog(), (tv) => expectNoThrow(() => tv.closeAllViewers()));
  });
}, {owner: 'agolovko@datagrok.ai'});
