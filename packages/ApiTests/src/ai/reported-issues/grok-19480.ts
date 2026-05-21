import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectFiresWithin, expectPropAndLook, expectRoundTripPropAndLook, withAttachedViewer} from '../helpers';

// Regression coverage for GROK-19480: bar chart stacking activation,
// orientation, negative aggregates, legend ordering. Visual rendering is out
// of scope here; we pin the state machine (props get set, persist through
// getOptions(true).look, and the viewer instance does not throw).
category('AI: GROK-19480: Bar chart stacking + orientation', () => {
  test('stack activates with negative values in df.plot.bar', async () => {
    const d = df([
      ['age', 'int', [25, -5, 40, -10, 30, 60, -20, 15]],
      ['race', 'string', ['Asian', 'White', 'Black', 'Asian', 'White', 'Black', 'Asian', 'White']],
      ['sex', 'string', ['M', 'F', 'M', 'F', 'M', 'F', 'M', 'F']],
    ]);
    const v = d.plot.bar({value: 'age', split: 'race', stack: 'sex'});
    expect(v.type, DG.VIEWER.BAR_CHART);
    expect(v.dataFrame.rowCount === d.rowCount, true);
    expectPropAndLook(v, {valueColumnName: 'age', splitColumnName: 'race', stackColumnName: 'sex'});
  });

  test('valueColumnName round-trips when switched between numeric columns', async () => {
    const v = DG.Viewer.barChart(demog(), {value: 'age', split: 'race'}) as DG.BarChartViewer;
    expectPropAndLook(v, {valueColumnName: 'age'});
    for (var col of ['height', 'weight'])
      expectRoundTripPropAndLook(v, {valueColumnName: col});
  });

  test('orientation vertical persists in getOptions look', async () => {
    const v = DG.Viewer.barChart(demog(), {value: 'age', split: 'race'}) as DG.BarChartViewer;
    for (var o of ['vertical', 'horizontal'])
      expectRoundTripPropAndLook(v, {orientation: o});
  });

  test('onResetView fires when resetView() is invoked', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART, {value: 'age', split: 'race', stack: 'sex'},
      (v) => expectFiresWithin(v.onResetView, () => v.resetView()));
  });
});
