import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, expectPropAndLook, subscribeAll, withAttachedViewer} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:223 (DG.Viewer.barChart),
// public/js-api/src/viewer.ts:682 (BarChartViewer — resetView,
// onCategoryClicked, onCategoryHovered, onResetView),
// public/js-api/src/dataframe/data-frame.ts:566 (df.plot.bar),
// public/js-api/src/interfaces/d4.ts:786 (IBarChartSettings).
// Bar-only JS surface: typed factory + df.plot.bar shorthand, friendly-key
// aliasing (value/split/stack/color → *ColumnName), resetView() wrapper +
// onResetView round-trip, and rxjs Observable shape of onCategoryClicked /
// onCategoryHovered (firing requires real canvas mouse events).
// Detached viewers are never close()-d (close() on a never-attached viewer
// throws inside Dart interop).
category('AI: Viewers: BarChart JS API', () => {
  test('factory typed via DG.Viewer.barChart', async () => {
    const v = DG.Viewer.barChart(demog(100), {value: 'age', valueAggrType: 'avg', split: 'race'});
    expect(v instanceof DG.BarChartViewer, true);
    expectPropAndLook(v, {valueColumnName: 'age', valueAggrType: 'avg', splitColumnName: 'race'});
  });

  test('factory via df.plot.bar shorthand', async () => {
    const df = demog();
    const v = df.plot.bar({value: 'height', split: 'race', stack: 'sex'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.BAR_CHART);
    expect(v.dataFrame === df, true);
    expectPropAndLook(v, {valueColumnName: 'height', splitColumnName: 'race', stackColumnName: 'sex'});
  });

  test('resetView does not throw on attached BarChart', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART, {value: 'age', split: 'race'},
      (v) => { expect(v.resetView() === undefined, true); });
  });

  test('onCategoryClicked and onCategoryHovered are rxjs Observables', async () => {
    const v = DG.Viewer.barChart(demog(20), {value: 'age', split: 'race'}) as DG.BarChartViewer;
    subscribeAll([v.onCategoryClicked, v.onCategoryHovered])();
  });

  test('onResetView round-trip via resetView', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART, {value: 'age', split: 'race'},
      (v) => expectFiresWithin(v.onResetView, () => v.resetView()));
  });
});
