import * as DG from 'datagrok-api/dg';
import {category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectChoices, expectRoundTrip, look, reloadLayoutWithLook, until,
  withTableView} from '../helpers';

// ScatterPlot prop JSON-shape round-trips not covered by scatter-plot-js-api.ts.
category('AI: Viewers: ScatterPlot extras', () => {
  const markerColors = ['filteredRowsColor', 'selectedRowsColor', 'filteredOutRowsColor', 'missingValueColor',
    'whiskerColor'] as const;

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

  test('jitterSize/jitterSizeY/markerBorderWidth boundary round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {jitterSize: 0, jitterSizeY: 0, markerBorderWidth: 1});
    expectRoundTrip(c, {jitterSize: 50, jitterSizeY: 50, markerBorderWidth: 10});
  });

  test('removed markerOpacity / line transparencies move into the color alpha', async () => {
    const c = v();
    c.setOptions({filteredRowsColor: '#ff0000', markerOpacity: 50, movingAverageLineTransparency: 0.75});
    expect((c.props.filteredRowsColor >>> 0).toString(16), '7fff0000');
    expect(DG.Color.a(c.props.movingAverageLineColor), 64);
    for (const removed of ['markerOpacity', 'movingAverageLineTransparency'])
      expect(removed in look(c), false);

    await withTableView(demog(), (tv) => {
      tv.scatterPlot({x: 'height', y: 'weight'});
      const sp = reloadLayoutWithLook(tv, DG.VIEWER.SCATTER_PLOT, {markerOpacity: 50, filteredRowsColor: 0xFFFF0000});
      expect((sp.props.filteredRowsColor >>> 0).toString(16), '7fff0000');
      expect(DG.Color.a(sp.props.selectedRowsColor), 127);
    });
  });

  test('setMarkerOpacity sets the alpha of all marker colors', async () => {
    const c = v();
    c.setMarkerOpacity(50);
    for (const color of markerColors)
      expect(DG.Color.a(c.props[color]), 127);
    expect(DG.Color.a(c.props.regressionLineColor), 255);
  });

  test('transparencies are not viewer properties; markerOpacity is a hidden marker color alpha', async () => {
    const c = v();
    for (const removed of ['regressionLineTransparency', 'movingAverageLineTransparency'])
      expect(removed in c.props, false);
    await expectExceptionAsync(async () => {
      Reflect.set(c.props, 'regressionLineTransparency', 0.5);
    });
    Reflect.set(c.props, 'markerOpacity', 40);
    for (const color of markerColors)
      expect(DG.Color.a(c.props[color]), 102);
    expect(Reflect.get(c.props, 'markerOpacity'), 40);
    expect(c.getProperties().find((p) => p.name === 'markerOpacity')?.userEditable, false);
    expect('markerOpacity' in look(c), false);
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
