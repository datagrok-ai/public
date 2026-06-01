import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectChoices, expectRoundTrip, until} from '../helpers';

// ScatterPlot prop JSON-shape round-trips not covered by scatter-plot-js-api.ts.
category('AI: Viewers: ScatterPlot extras', () => {
  const v = (): DG.ScatterPlotViewer => demog().plot.scatter({x: 'age', y: 'height'}) as DG.ScatterPlotViewer;

  test('xAxisType + yAxisType (AxisType) round-trip + getProperties choices on yAxisType', async () => {
    const c = v();
    expectRoundTrip(c, {xAxisType: 'logarithmic', yAxisType: 'logarithmic'});
    expectRoundTrip(c, {xAxisType: 'linear', yAxisType: 'linear'});
    expectChoices(c, 'yAxisType', ['linear', 'logarithmic']);
  });

  test('zoomAndFilter four-value choices round-trip + getProperties choices', async () => {
    const c = v();
    const values = ['no action', 'filter by zoom', 'zoom by filter', 'pack and zoom by filter'];
    for (const val of values)
      expectRoundTrip(c, {zoomAndFilter: val});
    expectChoices(c, 'zoomAndFilter', values);
  });

  test('xAxisLabelOrientation four-value choices including 45 degrees round-trip', async () => {
    const c = v();
    const values = ['Auto', 'Horz', 'Vert', '45 degrees'];
    for (const val of values)
      expectRoundTrip(c, {xAxisLabelOrientation: val});
    expectChoices(c, 'xAxisLabelOrientation', values);
  });

  test('showXHistogram + showYHistogram + histogramBins (5..100) round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {showXHistogram: true, showYHistogram: true, histogramBins: 25});
    for (const n of [5, 100])
      expectRoundTrip(c, {histogramBins: n});
    expectRoundTrip(c, {showXHistogram: false, showYHistogram: false});
  });

  test('regression line family combined bool round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {
      showRegressionLine: true, showRegressionLineEquation: false,
      showSpearmanCorrelation: true, showPearsonCorrelation: true,
      showMeanAbsoluteError: true, showRootMeanSquareError: true, regressionPerCategory: false,
    });
    expectRoundTrip(c, {
      showRegressionLine: false, showRegressionLineEquation: true,
      showSpearmanCorrelation: false, showPearsonCorrelation: false,
      showMeanAbsoluteError: false, showRootMeanSquareError: false, regressionPerCategory: true,
    });
  });

  test('markerOpacity/jitterSize/jitterSizeY/markerBorderWidth boundary round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {markerOpacity: 0, jitterSize: 0, jitterSizeY: 0, markerBorderWidth: 1});
    expectRoundTrip(c, {markerOpacity: 100, jitterSize: 50, jitterSizeY: 50, markerBorderWidth: 10});
  });

  test('selection visibility bools combined round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {showCurrentPoint: false, showMouseOverPoint: false,
      showMouseOverRowGroup: false, showSelectedRows: false, resetSelectionOnBackgroundClick: false});
    expectRoundTrip(c, {showCurrentPoint: true, showMouseOverPoint: true,
      showMouseOverRowGroup: true, showSelectedRows: true, resetSelectionOnBackgroundClick: true});
  });

  // Changing xColumnName after auto-detect must re-run the per-axis whisker detector.
  test('whiskers re-detect when xColumnName changes after auto-detect', async () => {
    const d = df([
      ['a', 'double', [1, 2, 3, 4, 5, 6, 7, 8]],
      ['a min', 'double', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]],
      ['a max', 'double', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]],
      ['b', 'double', [10, 20, 30, 40, 50, 60, 70, 80]],
      ['b min', 'double', [9, 18, 28, 38, 48, 58, 68, 78]],
      ['b max', 'double', [11, 22, 32, 42, 52, 62, 72, 82]],
    ]);
    const c = d.plot.scatter() as DG.ScatterPlotViewer;
    const initialX = c.props['xColumnName'] as string;
    expect(initialX === 'a' || initialX === 'b', true);
    expect(c.props['xWhiskerMinColumnName'], initialX + ' min');
    expect(c.props['xWhiskerMaxColumnName'], initialX + ' max');

    const otherBase = initialX === 'a' ? 'b' : 'a';
    c.setOptions({xColumnName: otherBase});
    await until(() => c.props['xColumnName'] === otherBase);

    expect(c.props['xColumnName'], otherBase);
    expect(c.props['xWhiskerMinColumnName'], otherBase + ' min');
    expect(c.props['xWhiskerMaxColumnName'], otherBase + ' max');
  });
}, {owner: 'agolovko@datagrok.ai'});
