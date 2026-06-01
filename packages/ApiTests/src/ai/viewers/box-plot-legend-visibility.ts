import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, look, until, withAttachedViewer} from '../helpers';

// Source: core/client/d4/lib/src/viewers/box_plot/box_plot_core.dart:94
// (canShowLegend) and :28 (legendCol => markerColorCol). When a single
// category1 column is selected and markerColorColumnName resolves to the
// SAME column, the box plot's auto legend is suppressed because the
// categories would already be rendered along the X axis (duplicate info).
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
// suppress the legend independently.
category('AI: Viewers: BoxPlot legend visibility', () => {
  const opts = {valueColumnName: 'age', category1ColumnName: 'race',
    markerColorColumnName: 'race', legendVisibility: 'Always'};
  const legendPresent = (v: DG.Viewer): boolean => v.root.querySelector('.d4-legend') != null;

  test('legend hidden when single category equals markerColor column', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, opts, async (v) => {
      await until(() => look(v)['markerColorColumnName'] === 'race');
      expectLook(v, {category1ColumnName: 'race', markerColorColumnName: 'race'});
      expect(legendPresent(v), false);
    });
  });

  test('legend visible when category and markerColor columns differ', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, opts, async (v) => {
      await until(() => look(v)['category1ColumnName'] === 'race');
      v.setOptions({allowColorSynchronization: false, markerColorColumnName: 'sex'});
      await until(() => look(v)['markerColorColumnName'] === 'sex' && legendPresent(v));
      expectLook(v, {category1ColumnName: 'race', markerColorColumnName: 'sex'});
      expect(legendPresent(v), true);
    });
  });

  test('legend toggles when markerColor switches to/from category column', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, opts, async (v) => {
      await until(() => look(v)['category1ColumnName'] === 'race');
      v.setOptions({allowColorSynchronization: false, markerColorColumnName: 'sex'});
      await until(() => legendPresent(v) === true);
      expect(legendPresent(v), true);
      v.setOptions({markerColorColumnName: 'race'});
      await until(() => legendPresent(v) === false);
      expect(legendPresent(v), false);
      v.setOptions({markerColorColumnName: 'sex'});
      await until(() => legendPresent(v) === true);
      expect(legendPresent(v), true);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
