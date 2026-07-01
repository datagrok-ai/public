import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, until, withTableView} from '../helpers';

// Note: View.getInfo() intentionally untested — TS labels its return ViewLayout but the Dart reg
// returns ViewInfo (type-label mismatch, flagged for human review).

function countViewers(tv: DG.TableView, type: string): number {
  let n = 0;
  for (const v of tv.viewers) {
    if (v.type === type)
      n++;
  }
  return n;
}

function hasViewer(tv: DG.TableView, type: string): boolean {
  return countViewers(tv, type) > 0;
}

category('AI: App: View Layout Roundtrip', () => {
  test('saveLayout returns a ViewLayout with TableView type and viewer types in json', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT);
      tv.addViewer(DG.VIEWER.HISTOGRAM);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && hasViewer(tv, DG.VIEWER.HISTOGRAM));
      const layout = tv.saveLayout();
      expect(layout instanceof DG.ViewLayout, true);
      const json = layout.toJson();
      expect(typeof json, 'string');
      expect(json.length > 0, true);
      const parsed = JSON.parse(json);
      // ViewLayout has no `type` getter on the TS wrapper; assert TableView via the serialized json.
      expect(parsed['type'], 'TableView');
      expect(json.indexOf(DG.VIEWER.SCATTER_PLOT) >= 0, true);
      expect(json.indexOf(DG.VIEWER.HISTOGRAM) >= 0, true);
    });
  });

  test('toJson / fromJson string round-trip re-creates the viewer arrangement on a fresh view', async () => {
    const df = demog();
    let json = '';
    await withTableView(df, async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT);
      tv.addViewer(DG.VIEWER.HISTOGRAM);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && hasViewer(tv, DG.VIEWER.HISTOGRAM));
      json = tv.saveLayout().toJson();
    });
    await withTableView(demog(), async (tv2) => {
      const restored = DG.ViewLayout.fromJson(json);
      expect(restored instanceof DG.ViewLayout, true);
      tv2.loadLayout(restored);
      await until(() => hasViewer(tv2, DG.VIEWER.SCATTER_PLOT) && hasViewer(tv2, DG.VIEWER.HISTOGRAM));
      expect(hasViewer(tv2, DG.VIEWER.SCATTER_PLOT), true);
      expect(hasViewer(tv2, DG.VIEWER.HISTOGRAM), true);
    });
  });

  test('loadLayout re-creates viewers on the same view after resetLayout', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT));
      const saved = tv.saveLayout();
      tv.resetLayout();
      await until(() => !hasViewer(tv, DG.VIEWER.SCATTER_PLOT));
      tv.loadLayout(saved);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT));
      expect(hasViewer(tv, DG.VIEWER.SCATTER_PLOT), true);
    });
  });

  test('resetLayout leaves only the grid', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT);
      tv.addViewer(DG.VIEWER.HISTOGRAM);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && hasViewer(tv, DG.VIEWER.HISTOGRAM));
      tv.resetLayout();
      await until(() => !hasViewer(tv, DG.VIEWER.SCATTER_PLOT) && !hasViewer(tv, DG.VIEWER.HISTOGRAM));
      expect(countViewers(tv, DG.VIEWER.SCATTER_PLOT), 0);
      expect(countViewers(tv, DG.VIEWER.HISTOGRAM), 0);
      expect(countViewers(tv, DG.VIEWER.GRID) >= 1, true);
    });
  });

  test('ViewLayout.viewState get/set round-trip and fromViewState', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT));
      const layout = tv.saveLayout();
      const state = layout.viewState;
      expect(typeof state, 'string');
      expect(state.length > 0, true);
      const fromState = DG.ViewLayout.fromViewState(state);
      expect(fromState instanceof DG.ViewLayout, true);
      const newState = fromState.viewState;
      layout.viewState = newState;
      expect(typeof layout.viewState, 'string');
      expect(layout.viewState.length > 0, true);
    });
  });

  test('ViewLayout user data get/set', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const layout = tv.saveLayout();
      layout.setUserDataValue('roundtrip-key', 'roundtrip-value');
      expect(layout.getUserDataValue('roundtrip-key'), 'roundtrip-value');
      const missing = layout.getUserDataValue('no-such-key');
      expect(missing == null || missing === '', true);
    });
  });

  test('ViewLayout.columns returns a ColumnInfo list', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const layout = tv.saveLayout();
      const columns = layout.columns;
      // columns is server-persistence metadata, empty for an in-memory saveLayout() never saved to
      // the server. Assert the wrap returns an array of ColumnInfo (don't require a specific column).
      expect(Array.isArray(columns), true);
      for (const c of columns)
        expect(typeof c.name, 'string');
    });
  });

  test('onViewLayoutGenerated fires when saveLayout is called', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      tv.addViewer(DG.VIEWER.SCATTER_PLOT);
      await until(() => hasViewer(tv, DG.VIEWER.SCATTER_PLOT));
      await expectFiresWithin(grok.events.onViewLayoutGenerated, () => {tv.saveLayout();});
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
