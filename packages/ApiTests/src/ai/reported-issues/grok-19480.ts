import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {take} from 'rxjs/operators';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19480: bar chart stacking activation,
// orientation, negative aggregates, legend ordering. Visual rendering is out
// of scope here; we pin the state machine (props get set, persist through
// getOptions(true).look, and the viewer instance does not throw).
category('AI: GROK-19480: Bar chart stacking + orientation', () => {
  test('stack activates with negative values in df.plot.bar', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'age', [25, -5, 40, -10, 30, 60, -20, 15]),
      DG.Column.fromList('string', 'race', ['Asian', 'White', 'Black', 'Asian', 'White', 'Black', 'Asian', 'White']),
      DG.Column.fromList('string', 'sex', ['M', 'F', 'M', 'F', 'M', 'F', 'M', 'F']),
    ]);
    const v = df.plot.bar({value: 'age', split: 'race', stack: 'sex'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.BAR_CHART);
    // df.plot.bar returns a generic Viewer wrapper at the TS type level; runtime
    // instance is BarChartViewer via toJs (round-7 finding). dataFrame identity
    // via === is unreliable since the getter wraps a fresh JS shell each call.
    expect(v.dataFrame.rowCount === df.rowCount, true);
    expect(v.props['valueColumnName'], 'age');
    expect(v.props['splitColumnName'], 'race');
    expect(v.props['stackColumnName'], 'sex');
    const look = v.getOptions(true).look;
    expect(look['valueColumnName'], 'age');
    expect(look['splitColumnName'], 'race');
    expect(look['stackColumnName'], 'sex');
  });

  test('valueColumnName round-trips when switched between numeric columns', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'}) as DG.BarChartViewer;
    expect(v.props['valueColumnName'], 'age');
    expect(v.getOptions(true).look['valueColumnName'], 'age');
    v.setOptions({valueColumnName: 'height'});
    expect(v.props['valueColumnName'], 'height');
    expect(v.getOptions(true).look['valueColumnName'], 'height');
    v.setOptions({valueColumnName: 'weight'});
    expect(v.props['valueColumnName'], 'weight');
    expect(v.getOptions(true).look['valueColumnName'], 'weight');
  });

  test('orientation vertical persists in getOptions look', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'}) as DG.BarChartViewer;
    v.setOptions({orientation: 'vertical'});
    expect(v.props['orientation'], 'vertical');
    expect(v.getOptions(true).look['orientation'], 'vertical');
    v.setOptions({orientation: 'horizontal'});
    expect(v.props['orientation'], 'horizontal');
    expect(v.getOptions(true).look['orientation'], 'horizontal');
  });

  test('onResetView fires when resetView() is invoked', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    let sub: Subscription | undefined;
    try {
      const v = tv.addViewer(DG.VIEWER.BAR_CHART,
        {value: 'age', split: 'race', stack: 'sex'}) as DG.BarChartViewer;
      expect(v instanceof DG.BarChartViewer, true);
      const fired: Promise<boolean> = new Promise((resolve) => {
        sub = v.onResetView.pipe(take(1)).subscribe(() => resolve(true));
      });
      const timeout: Promise<boolean> = (async () => {
        await DG.delay(1500);
        return false;
      })();
      v.resetView();
      const ok = await Promise.race([fired, timeout]);
      expect(ok, true);
    }
    finally {
      if (sub != null)
        sub.unsubscribe();
      tv.close();
    }
  });
});
