import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for github.com/datagrok-ai/public/issues/3469
// (commit a91c4539af, 2025-09-16): Box plot — when the user picks a column for
// "Markers Color" and then changes the X-axis (Category 1 / Category 2)
// column, the platform used to silently overwrite the user's color column
// with the new category column. The fix: once `markerColorColumnName` (or
// `markerColorMap`) is changed *manually*, the internal
// `_allowColorSynchronization` flag flips to false and the auto-sync stops.
//
// Auto-sync still runs on every category change while the user has not
// touched the color column manually — that's the documented behavior, called
// out in the inline help for `markerColorColumnName`: "Changing *Category 1*
// or *Category 2* sets the color scheme to categorical (same as selected
// category column)". Note that `BoxPlotLook.auto()` already pre-seeds both
// `category1ColumnName` and `markerColorColumnName` to the same auto-picked
// categorical column, so a brand-new viewer has both set on attach. Tests
// drive *changes* on top of that initial state.
//
// Implementation: core/client/d4/lib/src/viewers/box_plot/box_plot_core.dart
// — `_allowColorSynchronization` (line ~414) and the Timer-driven
// `markerColorColumnName.set(...)` block (line ~475). Because the auto-set
// runs in a `Timer.run`, every category change is followed by a `DG.delay`
// before reading state.
category('AI: gh-3469: Box plot color column auto-sync vs. manual override', () => {
  test('auto-sync: changing Category 1 propagates to markerColorColumnName when untouched', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await DG.delay(300);

      // Pin the starting state — auto() sets both to the same column.
      const initial = v.getOptions(true).look;
      expect(initial['category1ColumnName'], 'race');

      // Change Category 1 to a different column. With no manual color
      // override yet, markerColorColumnName must follow via Timer.run.
      v.setOptions({category1ColumnName: 'sex'});
      await DG.delay(300);

      const after = v.getOptions(true).look;
      expect(after['category1ColumnName'], 'sex');
      expect(after['markerColorColumnName'], 'sex');
    }
    finally {
      tv.close();
    }
  });

  test('manual override sticks: user-picked color column survives a Category 1 change', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await DG.delay(300);

      // User explicitly picks a *different* column for the marker color.
      // This flips _allowColorSynchronization to false inside the core.
      v.setOptions({markerColorColumnName: 'sex'});
      await DG.delay(300);
      expect(v.getOptions(true).look['markerColorColumnName'], 'sex');

      // User changes Category 1 to a third column. The color must NOT follow.
      v.setOptions({category1ColumnName: 'started'});
      await DG.delay(300);

      const afterCategoryChange = v.getOptions(true).look;
      expect(afterCategoryChange['category1ColumnName'], 'started');
      // The bug from the issue: this used to flip to 'started'.
      expect(afterCategoryChange['markerColorColumnName'], 'sex');
    }
    finally {
      tv.close();
    }
  });

  test('manual override sticks across Category 2 changes too', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race', category2: 'sex'});
      tv.addViewer(v);
      await DG.delay(300);

      v.setOptions({markerColorColumnName: 'started'});
      await DG.delay(300);
      expect(v.getOptions(true).look['markerColorColumnName'], 'started');

      // Swap Category 2 — color must stay on 'started'.
      v.setOptions({category2ColumnName: 'disease'});
      await DG.delay(300);

      const after = v.getOptions(true).look;
      expect(after['category2ColumnName'], 'disease');
      expect(after['markerColorColumnName'], 'started');
    }
    finally {
      tv.close();
    }
  });

  test('global escape hatch: allowColorSynchronization=false freezes the color column on category change', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      // The settable `allowColorSynchronization` look property (default true)
      // gates the same auto-sync block. Setting it to false up front means
      // a category change must NOT propagate to markerColorColumnName.
      // BoxPlotLook.auto() pre-seeds the color column at attach, so we pin
      // the *delta* across a category swap, not the absolute value.
      const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await DG.delay(300);

      const initialColor = v.getOptions(true).look['markerColorColumnName'];
      // Sanity: auto() seeded the color to match category1.
      expect(initialColor, 'race');

      v.setOptions({allowColorSynchronization: false});
      await DG.delay(100);

      // Now flip Category 1 — the color column must stay on 'race'.
      v.setOptions({category1ColumnName: 'sex'});
      await DG.delay(300);

      const after = v.getOptions(true).look;
      expect(after['category1ColumnName'], 'sex');
      expect(after['markerColorColumnName'], initialColor);
    }
    finally {
      tv.close();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
