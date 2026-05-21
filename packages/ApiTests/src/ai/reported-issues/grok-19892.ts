import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectFiresWithin, expectNoThrow, look, withAttachedViewer} from '../helpers';

// Regression coverage for GROK-19892: a box plot's horizontal slider could be
// dragged so its visible range collapsed to zero, breaking the plot. The fix
// (1.28.0) prevents zero-range. IBoxPlotSettings does not expose an explicit
// "horizontal slider" prop on the JS API surface — the range-y props it does
// expose are the value-axis range (valueMin/valueMax), the color range
// (colorMin/colorMax) and the marker size range (markerMinSize/markerMaxSize).
// We fall back to a "viewer survives extreme prop values" smoke test, and
// pin the resetView / onResetView round-trip.
category('AI: GROK-19892: Box plot horizontal slider zero range', () => {
  test('df.plot.box default props snapshot for IBoxPlotSettings', async () => {
    const d = df([
      ['value', 'double', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]],
      ['category', 'string', ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B', 'A', 'B']],
    ]);
    const v = d.plot.box({value: 'value', category: 'category'});
    expect(v.type, DG.VIEWER.BOX_PLOT);
    expect(v.dataFrame.rowCount === d.rowCount, true);
    expect(v.props['valueColumnName'], 'value');
  });

  test('zero-range valueMin === valueMax does not crash the viewer', async () => {
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'}) as DG.BoxPlot;
    expectNoThrow(() => v.setOptions({valueMin: 30, valueMax: 30}));
    const l = look(v);
    if (typeof l['valueMin'] === 'number' && typeof l['valueMax'] === 'number')
      expect(l['valueMax'] - l['valueMin'] >= 0, true);
  });

  test('extreme markerMinSize === markerMaxSize survives setOptions', async () => {
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'}) as DG.BoxPlot;
    expectNoThrow(() => v.setOptions({markerMinSize: 5, markerMaxSize: 5}));
    expect(typeof look(v)['markerMinSize'] === 'number', true);
    expect(typeof look(v)['markerMaxSize'] === 'number', true);
  });

  test('onResetView fires when resetView() is invoked', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, {value: 'age', category1: 'race'},
      (v) => expectFiresWithin(v.onResetView, () => v.resetView()));
  });
});
