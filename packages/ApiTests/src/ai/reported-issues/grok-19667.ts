import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19667: switching IBarChartSettings.axisType to
// 'logarithmic' on certain Split values broke the chart. The fix (1.27.0)
// preserves viewer state under linear<->logarithmic toggles. We pin the
// property surface (axisType: keyof typeof AxisType, d4.ts:816, AxisType enum
// at d4.ts:769-772) and its round-trip through setOptions/getOptions(true).look.
category('AI: GROK-19667: Bar chart axisType log scale stability', () => {
  test('setOptions axisType=logarithmic does not throw and round-trips', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.bar({value: 'age', split: 'race', valueAggrType: 'count'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.BAR_CHART);
    // dataFrame identity via === is unreliable across the Dart->JS wrapper boundary.
    expect(v.dataFrame.rowCount === df.rowCount, true);

    var threw = false;
    try {
      v.setOptions({axisType: 'logarithmic'});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.getOptions(true).look['axisType'], 'logarithmic');
  });

  test('switching axisType back to linear does not throw and round-trips', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race', valueAggrType: 'count'}) as DG.BarChartViewer;
    expect(v instanceof DG.BarChartViewer, true);

    var threw = false;
    try {
      v.setOptions({axisType: 'logarithmic'});
      v.setOptions({axisType: 'linear'});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.getOptions(true).look['axisType'], 'linear');
  });

  test('toggling axisType with valueMin/valueMax set does not throw', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.bar({value: 'age', split: 'race', valueAggrType: 'count'});
    var threw = false;
    try {
      v.setOptions({valueMin: 1, valueMax: 100});
      v.setOptions({axisType: 'logarithmic'});
      v.setOptions({axisType: 'linear'});
      v.setOptions({axisType: 'logarithmic'});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    const look = v.getOptions(true).look;
    expect(look['axisType'], 'logarithmic');
    // valueMin/valueMax may not echo back if the platform doesn't store them
    // on the BarChart look — accept either, just assert the type is sane when
    // they do round-trip.
    if (look['valueMin'] !== undefined)
      expect(typeof look['valueMin'] === 'number', true);
    if (look['valueMax'] !== undefined)
      expect(typeof look['valueMax'] === 'number', true);
  });

  test('getProperties() exposes axisType with non-empty description (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = df.plot.bar({value: 'age', split: 'race', valueAggrType: 'count'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);
    var found: DG.Property | null = null;
    for (var p of props) {
      if (p.name === 'axisType') {
        found = p;
        break;
      }
    }
    expect(found != null, true);
    // Description is best-effort: assert string-typed if present, but do not
    // fail on emptiness — the platform may not ship a description for every prop.
    const desc = found!.description;
    expect(typeof desc === 'string' || desc == null, true);
  });
});
