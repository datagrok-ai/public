import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, wait, withTableView} from '../helpers';

// Regression coverage for github.com/datagrok-ai/public/issues/3469
// (commit a91c4539af, 2025-09-16): Box plot — when the user picks a column for
// "Markers Color" and then changes the X-axis (Category 1 / Category 2)
// column, the platform used to silently overwrite the user's color column
// with the new category column. The fix: once `markerColorColumnName` (or
// `markerColorMap`) is changed *manually*, the internal
// `_allowColorSynchronization` flag flips to false and the auto-sync stops.
//
// Auto-sync still runs on every category change while the user has not
// touched the color column manually — that's the documented behavior.
// `BoxPlotLook.auto()` pre-seeds both `category1ColumnName` and
// `markerColorColumnName` to the same auto-picked categorical column, so a
// brand-new viewer has both set on attach. Tests drive *changes* on top of
// that initial state. The auto-set runs in a `Timer.run`, so every category
// change is followed by a `DG.delay` before reading state.
category('AI: gh-3469: Box plot color column auto-sync vs. manual override', () => {
  test('auto-sync: changing Category 1 propagates to markerColorColumnName when untouched', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await wait();
      expectLook(v, {category1ColumnName: 'race'});
      v.setOptions({category1ColumnName: 'sex'});
      await wait();
      expectLook(v, {category1ColumnName: 'sex', markerColorColumnName: 'sex'});
    });
  });

  test('manual override sticks: user-picked color column survives a Category 1 change', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await wait();
      // User explicitly picks a *different* column for the marker color.
      // This flips _allowColorSynchronization to false inside the core.
      v.setOptions({markerColorColumnName: 'sex'});
      await wait();
      expectLook(v, {markerColorColumnName: 'sex'});
      v.setOptions({category1ColumnName: 'started'});
      await wait();
      // The bug from the issue: this used to flip to 'started'.
      expectLook(v, {category1ColumnName: 'started', markerColorColumnName: 'sex'});
    });
  });

  test('manual override sticks across Category 2 changes too', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race', category2: 'sex'});
      tv.addViewer(v);
      await wait();
      v.setOptions({markerColorColumnName: 'started'});
      await wait();
      expectLook(v, {markerColorColumnName: 'started'});
      v.setOptions({category2ColumnName: 'disease'});
      await wait();
      expectLook(v, {category2ColumnName: 'disease', markerColorColumnName: 'started'});
    });
  });

  test('global escape hatch: allowColorSynchronization=false freezes the color column on category change', async () => {
    await withTableView(demog(), async (tv) => {
      // The settable `allowColorSynchronization` look property (default true)
      // gates the same auto-sync block. Setting it to false up front means
      // a category change must NOT propagate to markerColorColumnName.
      // BoxPlotLook.auto() pre-seeds the color column at attach, so we pin
      // the *delta* across a category swap, not the absolute value.
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await wait();
      expectLook(v, {markerColorColumnName: 'race'});
      v.setOptions({allowColorSynchronization: false});
      await wait(100);
      v.setOptions({category1ColumnName: 'sex'});
      await wait();
      expectLook(v, {category1ColumnName: 'sex', markerColorColumnName: 'race'});
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
