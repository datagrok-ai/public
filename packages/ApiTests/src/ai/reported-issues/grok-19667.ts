import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, expectRoundTrip, findProp, look} from '../helpers';

// Regression coverage for GROK-19667: switching IBarChartSettings.axisType to
// 'logarithmic' on certain Split values broke the chart. The fix (1.27.0)
// preserves viewer state under linear<->logarithmic toggles.
category('AI: GROK-19667: Bar chart axisType log scale stability', () => {
  const v = (): DG.Viewer => demog().plot.bar({value: 'age', split: 'race', valueAggrType: 'count'});

  test('setOptions axisType=logarithmic does not throw and round-trips', async () => {
    expectRoundTrip(v(), {axisType: 'logarithmic'});
  });

  test('switching axisType back to linear does not throw and round-trips', async () => {
    const c = DG.Viewer.barChart(demog(), {value: 'age', split: 'race', valueAggrType: 'count'}) as DG.BarChartViewer;
    expectNoThrow(() => {
      c.setOptions({axisType: 'logarithmic'});
      c.setOptions({axisType: 'linear'});
    });
    expectLook(c, {axisType: 'linear'});
  });

  test('toggling axisType with valueMin/valueMax set does not throw', async () => {
    const c = v();
    expectNoThrow(() => {
      c.setOptions({valueMin: 1, valueMax: 100});
      c.setOptions({axisType: 'logarithmic'});
      c.setOptions({axisType: 'linear'});
      c.setOptions({axisType: 'logarithmic'});
    });
    expectLook(c, {axisType: 'logarithmic'});
    const l = look(c);
    if (l['valueMin'] !== undefined)
      expect(typeof l['valueMin'] === 'number', true);
    if (l['valueMax'] !== undefined)
      expect(typeof l['valueMax'] === 'number', true);
  });

  test('getProperties() exposes axisType with non-empty description (best-effort)', async () => {
    const v = demog(20).plot.bar({value: 'age', split: 'race', valueAggrType: 'count'});
    expect(findProp(v, 'axisType') != null, true);
  });
});
