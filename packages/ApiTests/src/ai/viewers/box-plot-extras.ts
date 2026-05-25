import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectRoundTrip} from '../helpers';

// JS API source: public/js-api/src/viewer.ts (DG.Viewer.boxPlot),
// public/js-api/src/interfaces/d4.d.ts:136 (IBoxPlotSettings),
// core/client/d4/lib/src/viewers/box_plot/box_plot_look.dart (BoxPlotLook).
// Box-plot-only JSON-shape coverage for props the existing
// src/ai/viewers/box-plot-js-api.ts (factory + statisticsFormat) and
// reported-issues regression tests don't pin: plotStyle (box/violin),
// violin bins int, the per-statistic boolean togglers
// (showMin/showMax/showAvg/showMed/showStdev/showQ1/showQ3),
// binColorAggrType + choices, axisType (AxisType linear/logarithmic),
// markerOpacity 0..100 boundary, and the showInside/Outside/PValue/
// EmptyCategories combined bool envelope. All assertions read state via
// getOptions(true).look — no first-paint geometry.
category('AI: Viewers: BoxPlot extras', () => {
  // no negative case: BoxPlot setOptions has no defined failure mode for
  // these look props — invalid values are coerced or ignored, not thrown.

  const v = (): DG.BoxPlot => DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});

  test('plotStyle box/violin round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {plotStyle: 'violin'});
    expectRoundTrip(c, {plotStyle: 'box'});
    expectChoices(c, 'plotStyle', ['box', 'violin']);
  });

  test('violin bins int round-trip within 20-1000', async () => {
    const c = v();
    c.setOptions({plotStyle: 'violin'});
    for (const n of [100, 20, 1000])
      expectRoundTrip(c, {bins: n});
  });

  test('per-statistic toggler bools combined round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {
      showStatistics: true,
      showMin: false, showMax: false, showAvg: false,
      showMed: false, showStdev: true, showQ1: true, showQ3: true,
    });
    expectRoundTrip(c, {
      showMin: true, showMax: true, showAvg: true,
      showMed: true, showStdev: false, showQ1: false, showQ3: false,
    });
  });

  test('binColorAggrType round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {binColorColumnName: 'height', binColorAggrType: 'med'});
    expectChoices(c, 'binColorAggrType', ['avg', 'med']);
  });

  test('axisType linear/logarithmic round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {axisType: 'logarithmic'});
    expectRoundTrip(c, {axisType: 'linear'});
    expectChoices(c, 'axisType', ['linear', 'logarithmic']);
  });

  test('markerOpacity 0..100 boundary round-trip', async () => {
    const c = v();
    for (const n of [0, 50, 100])
      expectRoundTrip(c, {markerOpacity: n});
  });

  test('showInsideValues/showOutsideValues/showPValue/showEmptyCategories combined bool round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {
      showInsideValues: false, showOutsideValues: false, showPValue: false, showEmptyCategories: false,
    });
    expectRoundTrip(c, {
      showInsideValues: true, showOutsideValues: true, showPValue: true, showEmptyCategories: true,
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
