import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectChoices, expectRoundTrip, until} from '../helpers';

// JS API source: public/js-api/src/viewer.ts (DG.ScatterPlotViewer / df.plot.scatter),
// public/js-api/src/interfaces/d4.d.ts (IScatterPlotSettings),
// core/client/d4/lib/src/viewers/scatterplot/scatterplot_look.dart (ScatterPlotLook).
// Scatter-plot-only JSON-shape coverage for props the existing
// src/ai/viewers/scatter-plot-js-api.ts (factory aliases, zoom keys, events,
// formula/annotation helpers, addViewer + viewBox/getInfo) and reported-
// issues regressions (GROK-16994 whiskers, GROK-19511 linesBy) don't pin:
// xAxisType/yAxisType (AxisType), zoomAndFilter (4-value enum),
// xAxisLabelOrientation (4-value, includes '45 degrees'), the histogram-
// overlay envelope (showXHistogram/showYHistogram/histogramBins 5..100),
// the regression line family bools, the marker numeric boundary envelope
// (markerOpacity/jitterSize/jitterSizeY/markerBorderWidth), and the
// selection visibility bool envelope. All assertions read state via
// getOptions(true).look — no canvas-geometry getters (round 9's
// xAxisBox/yAxisBox first-paint pitfall stays clear).
category('AI: Viewers: ScatterPlot extras', () => {
  // no negative case: ScatterPlot setOptions has no defined failure mode for
  // these look props — invalid values are coerced or ignored, not thrown.

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

  // Property-dependency coverage (commit 8dd87cb743, "Scatter plot whiskers:
  // Whiskers auto detection optimized for 10k columns, column detection on
  // axis property change in property panel added"). When the dataframe was
  // attached without explicit X/Y and the platform auto-picked whisker
  // columns by suffix (`_whiskersAutoDetected = true` in scatterplot_look.dart),
  // a later xColumnName / yColumnName change must re-run the per-axis detector
  // (scatterplot_core.dart `onLookChanged`, lines ~370-395). This is the
  // companion to the existing GROK-16994 test, which only pins the initial
  // auto-detect — here we pin the *update on axis change* path.
  test('whiskers re-detect when xColumnName changes after auto-detect', async () => {
    // Two independent <base, base min, base max> triples. Auto-detect picks
    // one for X and one for Y at attach time; swapping xColumnName to the
    // other triple must move the X whiskers to the matching min/max.
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
