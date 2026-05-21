import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectPropAndLook, expectRoundTrip, subscribeAll, withAttachedViewer} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:231 (DG.Viewer.boxPlot),
// public/js-api/src/viewer.ts:713 (DG.BoxPlot — bare class name, no Viewer
// suffix), public/js-api/src/interfaces/d4.ts:999 (IBoxPlotSettings),
// public/js-api/src/dataframe/data-frame.ts:568 (df.plot.box).
// Box-only JS surface: typed factory returning DG.BoxPlot vs df.plot.box
// returning a base Viewer, friendly-key aliasing (value/category1/category2),
// the combined category1ColumnName + category2ColumnName JSON envelope,
// statisticsFormat choices via getProperties() (it's a numeric-format list:
// 'int', 'one/two/three/four digits after comma', 'scientific', etc.),
// the four event Observables (onResetView / onAfterDrawScene /
// onBeforeDrawScene / onPointClicked), and view.addViewer(VIEWER.BOX_PLOT)
// returning a typed instance.
category('AI: Viewers: BoxPlot JS API', () => {
  test('factory DG.Viewer.boxPlot returns typed DG.BoxPlot; df.plot.box returns base Viewer', async () => {
    const df = demog();
    const v = DG.Viewer.boxPlot(df, {valueColumnName: 'age', category1ColumnName: 'race'});
    expect(v instanceof DG.BoxPlot, true);
    expect(v.dataFrame === df, true);
    expectPropAndLook(v, {valueColumnName: 'age', category1ColumnName: 'race'});
    expect(df.plot.box({valueColumnName: 'age', category1ColumnName: 'race'}).type, DG.VIEWER.BOX_PLOT);
  });

  test('friendly-key aliases value/category1/category2 collapse to canonical names', async () => {
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race', category2: 'sex'});
    expectPropAndLook(v, {valueColumnName: 'age', category1ColumnName: 'race', category2ColumnName: 'sex'});
  });

  test('combined category1ColumnName + category2ColumnName JSON envelope round-trip', async () => {
    const v = DG.Viewer.boxPlot(demog(), {valueColumnName: 'age', category1ColumnName: 'race'});
    expectRoundTrip(v, {category1ColumnName: 'race', category2ColumnName: 'sex'});
  });

  test('showStatistics + statisticsFormat combined round-trip with getProperties choices', async () => {
    const v = DG.Viewer.boxPlot(demog(), {valueColumnName: 'age', category1ColumnName: 'race'});
    expectRoundTrip(v, {showStatistics: true, statisticsFormat: 'two digits after comma'});
    expectChoices(v, 'statisticsFormat', ['two digits after comma']);
  });

  test('onResetView, onAfterDrawScene, onBeforeDrawScene, onPointClicked are rxjs Observables', async () => {
    const v = DG.Viewer.boxPlot(demog(20), {valueColumnName: 'age', category1ColumnName: 'sex'});
    subscribeAll([v.onResetView, v.onAfterDrawScene, v.onBeforeDrawScene, v.onPointClicked])();
  });

  test('view.addViewer attaches a typed DG.BoxPlot', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT,
      {valueColumnName: 'age', category1ColumnName: 'race'}, (v, tv) => {
        expect(v instanceof DG.BoxPlot, true);
        var found: DG.Viewer | undefined;
        for (var x of tv.viewers)
          if (x.type === DG.VIEWER.BOX_PLOT) { found = x; break; }
        expect(found instanceof DG.BoxPlot, true);
      });
  });
});
