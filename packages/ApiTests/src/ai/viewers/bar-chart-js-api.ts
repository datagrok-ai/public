import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, expectPropAndLook, subscribeAll, withAttachedViewer} from '../helpers';

// BarChart JS surface: typed factory, df.plot.bar shorthand, resetView, and category/reset Observables.
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
      (v) => {expect(v.resetView() === undefined, true);});
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
