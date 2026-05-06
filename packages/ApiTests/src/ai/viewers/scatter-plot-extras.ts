import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

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

  test('xAxisType + yAxisType (AxisType) round-trip + getProperties choices on yAxisType', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    v.setOptions({xAxisType: 'logarithmic', yAxisType: 'logarithmic'});
    const log = v.getOptions(true).look;
    expect(log['xAxisType'], 'logarithmic');
    expect(log['yAxisType'], 'logarithmic');
    v.setOptions({xAxisType: 'linear', yAxisType: 'linear'});
    const lin = v.getOptions(true).look;
    expect(lin['xAxisType'], 'linear');
    expect(lin['yAxisType'], 'linear');
    const props = v.getProperties();
    var yProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'yAxisType') { yProp = p; break; }
    expect(yProp != null, true);
    expect(Array.isArray(yProp!.choices), true);
    expect(yProp!.choices.indexOf('linear') >= 0, true);
    expect(yProp!.choices.indexOf('logarithmic') >= 0, true);
  });

  test('zoomAndFilter four-value choices round-trip + getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    const values = ['no action', 'filter by zoom', 'zoom by filter', 'pack and zoom by filter'];
    for (var val of values) {
      v.setOptions({zoomAndFilter: val});
      expect(v.getOptions(true).look['zoomAndFilter'], val);
    }
    const props = v.getProperties();
    var zfProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'zoomAndFilter') { zfProp = p; break; }
    expect(zfProp != null, true);
    expect(Array.isArray(zfProp!.choices), true);
    for (var val of values)
      expect(zfProp!.choices.indexOf(val) >= 0, true);
  });

  test('xAxisLabelOrientation four-value choices including 45 degrees round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    const values = ['Auto', 'Horz', 'Vert', '45 degrees'];
    for (var val of values) {
      v.setOptions({xAxisLabelOrientation: val});
      expect(v.getOptions(true).look['xAxisLabelOrientation'], val);
    }
    const props = v.getProperties();
    var orientProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'xAxisLabelOrientation') { orientProp = p; break; }
    expect(orientProp != null, true);
    expect(Array.isArray(orientProp!.choices), true);
    for (var val of values)
      expect(orientProp!.choices.indexOf(val) >= 0, true);
  });

  test('showXHistogram + showYHistogram + histogramBins (5..100) round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    v.setOptions({showXHistogram: true, showYHistogram: true, histogramBins: 25});
    const on = v.getOptions(true).look;
    expect(on['showXHistogram'], true);
    expect(on['showYHistogram'], true);
    expect(on['histogramBins'], 25);
    v.setOptions({histogramBins: 5});
    expect(v.getOptions(true).look['histogramBins'], 5);
    v.setOptions({histogramBins: 100});
    expect(v.getOptions(true).look['histogramBins'], 100);
    v.setOptions({showXHistogram: false, showYHistogram: false});
    const off = v.getOptions(true).look;
    expect(off['showXHistogram'], false);
    expect(off['showYHistogram'], false);
  });

  test('regression line family combined bool round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    v.setOptions({
      showRegressionLine: true,
      showRegressionLineEquation: false,
      showSpearmanCorrelation: true,
      showPearsonCorrelation: true,
      showMeanAbsoluteError: true,
      showRootMeanSquareError: true,
      regressionPerCategory: false,
    });
    const set = v.getOptions(true).look;
    expect(set['showRegressionLine'], true);
    expect(set['showRegressionLineEquation'], false);
    expect(set['showSpearmanCorrelation'], true);
    expect(set['showPearsonCorrelation'], true);
    expect(set['showMeanAbsoluteError'], true);
    expect(set['showRootMeanSquareError'], true);
    expect(set['regressionPerCategory'], false);
    v.setOptions({
      showRegressionLine: false,
      showRegressionLineEquation: true,
      showSpearmanCorrelation: false,
      showPearsonCorrelation: false,
      showMeanAbsoluteError: false,
      showRootMeanSquareError: false,
      regressionPerCategory: true,
    });
    const cleared = v.getOptions(true).look;
    expect(cleared['showRegressionLine'], false);
    expect(cleared['showRegressionLineEquation'], true);
    expect(cleared['showSpearmanCorrelation'], false);
    expect(cleared['showPearsonCorrelation'], false);
    expect(cleared['showMeanAbsoluteError'], false);
    expect(cleared['showRootMeanSquareError'], false);
    expect(cleared['regressionPerCategory'], true);
  });

  test('markerOpacity/jitterSize/jitterSizeY/markerBorderWidth boundary round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    v.setOptions({markerOpacity: 0, jitterSize: 0, jitterSizeY: 0, markerBorderWidth: 1});
    const lo = v.getOptions(true).look;
    expect(lo['markerOpacity'], 0);
    expect(lo['jitterSize'], 0);
    expect(lo['jitterSizeY'], 0);
    expect(lo['markerBorderWidth'], 1);
    v.setOptions({markerOpacity: 100, jitterSize: 50, jitterSizeY: 50, markerBorderWidth: 10});
    const hi = v.getOptions(true).look;
    expect(hi['markerOpacity'], 100);
    expect(hi['jitterSize'], 50);
    expect(hi['jitterSizeY'], 50);
    expect(hi['markerBorderWidth'], 10);
  });

  test('selection visibility bools combined round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'age', y: 'height'});
    v.setOptions({
      showCurrentPoint: false, showMouseOverPoint: false,
      showMouseOverRowGroup: false, showSelectedRows: false,
      resetSelectionOnBackgroundClick: false,
    });
    const off = v.getOptions(true).look;
    expect(off['showCurrentPoint'], false);
    expect(off['showMouseOverPoint'], false);
    expect(off['showMouseOverRowGroup'], false);
    expect(off['showSelectedRows'], false);
    expect(off['resetSelectionOnBackgroundClick'], false);
    v.setOptions({
      showCurrentPoint: true, showMouseOverPoint: true,
      showMouseOverRowGroup: true, showSelectedRows: true,
      resetSelectionOnBackgroundClick: true,
    });
    const on = v.getOptions(true).look;
    expect(on['showCurrentPoint'], true);
    expect(on['showMouseOverPoint'], true);
    expect(on['showMouseOverRowGroup'], true);
    expect(on['showSelectedRows'], true);
    expect(on['resetSelectionOnBackgroundClick'], true);
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
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'a max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'b max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const initialX = v.props['xColumnName'];
    expect(initialX === 'a' || initialX === 'b', true);
    const initialMin = v.props['xWhiskerMinColumnName'];
    const initialMax = v.props['xWhiskerMaxColumnName'];
    expect(initialMin, initialX + ' min');
    expect(initialMax, initialX + ' max');

    // Flip xColumnName to the *other* base. Whiskers must follow.
    const otherBase = initialX === 'a' ? 'b' : 'a';
    v.setOptions({xColumnName: otherBase});
    await DG.delay(200);

    expect(v.props['xColumnName'], otherBase);
    expect(v.props['xWhiskerMinColumnName'], otherBase + ' min');
    expect(v.props['xWhiskerMaxColumnName'], otherBase + ' max');
  });
}, {owner: 'agolovko@datagrok.ai'});
