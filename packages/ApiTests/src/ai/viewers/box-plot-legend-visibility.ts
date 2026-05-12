import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Source: core/client/d4/lib/src/viewers/box_plot/box_plot_core.dart:94
// (canShowLegend) and :28 (legendCol => markerColorCol). When a single
// category1 column is selected and markerColorColumnName resolves to the
// SAME column, the box plot's auto legend is suppressed because the
// categories would already be rendered along the X axis (duplicate info).
// Logic: canShowLegend() returns false iff
//   categoryColumnNames.length == 1 && sameLegendColumnDisplayed(category1)
// where sameLegendColumnDisplayed(col) checks col == legendCol == markerColorCol.
// IMPORTANT: box_plot_core.dart:475-487 auto-syncs markerColorColumnName to
// the last category column whenever category1/2 changes (Timer.run),
// regardless of the order setOptions applies properties — the look-prop
// allowColorSynchronization is read each tick. We force-resync in two
// stages: first attach with the same column (auto-sync no-op), then flip
// allowColorSynchronization off and assign markerColorColumnName for the
// "different" case so the deferred sync sees the override and respects it
// (the private _allowColorSynchronization is also cleared on user-driven
// markerColorColumnName change, see :425). legendVisibility=Always pins
// super.canShowLegend() to true so a small headless viewport doesn't
// suppress the legend independently — the same-column branch is gated by
// the FIRST half of canShowLegend (line :94) which doesn't read
// legendVisibility, so the suppression we're testing still applies.
// We use a real TableView and a 300ms delay so the deferred legend refresh
// (Timer.run in onLookChanged) lands. Same-column case: .d4-legend should
// be removed from the DOM (legend.remove() at legend_mixin.dart:530 →
// legend.dart:168 'root?.remove()'). Different-column case: .d4-legend
// should be present in the DOM.
category('AI: Viewers: BoxPlot legend visibility', () => {
  test('legend hidden when single category equals markerColor column', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BOX_PLOT,
        {valueColumnName: 'age', category1ColumnName: 'race',
          markerColorColumnName: 'race', legendVisibility: 'Always'}) as DG.BoxPlot;
      expect(v instanceof DG.BoxPlot, true);
      await DG.delay(300);
      const look = v.getOptions(true).look;
      expect(look['category1ColumnName'], 'race');
      expect(look['markerColorColumnName'], 'race');
      const legendEl = v.root.querySelector('.d4-legend');
      expect(legendEl == null, true);
    }
    finally {
      tv.close();
    }
  });

  test('legend visible when category and markerColor columns differ', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BOX_PLOT,
        {valueColumnName: 'age', category1ColumnName: 'race',
          markerColorColumnName: 'race', legendVisibility: 'Always'}) as DG.BoxPlot;
      expect(v instanceof DG.BoxPlot, true);
      await DG.delay(200);
      v.setOptions({allowColorSynchronization: false, markerColorColumnName: 'sex'});
      await DG.delay(300);
      const look = v.getOptions(true).look;
      expect(look['category1ColumnName'], 'race');
      expect(look['markerColorColumnName'], 'sex');
      const legendEl = v.root.querySelector('.d4-legend');
      expect(legendEl != null, true);
    }
    finally {
      tv.close();
    }
  });

  test('legend toggles when markerColor switches to/from category column', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BOX_PLOT,
        {valueColumnName: 'age', category1ColumnName: 'race',
          markerColorColumnName: 'race', legendVisibility: 'Always'}) as DG.BoxPlot;
      await DG.delay(200);
      v.setOptions({allowColorSynchronization: false, markerColorColumnName: 'sex'});
      await DG.delay(300);
      expect(v.root.querySelector('.d4-legend') != null, true);

      v.setOptions({markerColorColumnName: 'race'});
      await DG.delay(300);
      expect(v.root.querySelector('.d4-legend') == null, true);

      v.setOptions({markerColorColumnName: 'sex'});
      await DG.delay(300);
      expect(v.root.querySelector('.d4-legend') != null, true);
    }
    finally {
      tv.close();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
