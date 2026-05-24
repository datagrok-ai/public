import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectLook, expectRoundTrip, look} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:259 (DG.Viewer.densityPlot),
// public/js-api/src/interfaces/d4.d.ts:357 (IDensityPlotSettings),
// core/client/d4/lib/src/viewers/density_plot/density_plot_look.dart (DensityPlotLook).
// Comprehensive look round-trip for DensityPlot props not pinned by prior AI
// regression tests (grok-13206 covers backColor init; grok-17118 covers
// binShape rectangle/hex + x/yColumnName). Targeted gaps: bins (1-200 slider),
// linearColorScheme List<int>, invertColorScheme, colorTransformType
// (AxisType), showColorScale, showXAxis/showYAxis, x/yMin x/yMax, and
// getProperties choices introspection. All assertions read state via
// getOptions(true).look — no first-paint geometry.
category('AI: Viewers: DensityPlot JS API', () => {
  // no negative case: DensityPlot setOptions has no defined failure mode for
  // the look props under test — invalid values are coerced or ignored, not
  // thrown. Negative-mode coverage for binShape vs string columns is already
  // pinned in src/ai/reported-issues/grok-17118.ts.

  const v = (): DG.Viewer => DG.Viewer.densityPlot(demog(), {x: 'age', y: 'height'});

  test('factory returns typed Viewer with DENSITY_PLOT type', async () => {
    const c = v();
    expect(c.type, DG.VIEWER.DENSITY_PLOT);
    expectLook(c, {xColumnName: 'age', yColumnName: 'height'});
  });

  test('bins integer round-trip within 1-200 range', async () => {
    const c = v();
    for (const n of [120, 1, 200])
      expectRoundTrip(c, {bins: n});
  });

  test('linearColorScheme number-array round-trip', async () => {
    const c = v();
    const scheme = [0xFF112233, 0xFF445566, 0xFFAABBCC];
    c.setOptions({linearColorScheme: scheme});
    expectArray(look(c)['linearColorScheme'] as number[], scheme);
  });

  test('showXAxis/showYAxis/showColorScale combined bool round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {showXAxis: false, showYAxis: false, showColorScale: false});
    expectRoundTrip(c, {showXAxis: true, showYAxis: true, showColorScale: true});
  });

  test('x/yMin x/yMax numeric round-trip', async () => {
    expectRoundTrip(v(), {xMin: 10, xMax: 80, yMin: 120, yMax: 200});
  });

  test('invertColorScheme + colorTransformType logarithmic round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {invertColorScheme: true, colorTransformType: 'logarithmic'});
    expectRoundTrip(c, {invertColorScheme: false, colorTransformType: 'linear'});
  });

  test('getProperties choices introspection on colorTransformType', async () => {
    expectChoices(v(), 'colorTransformType', ['linear', 'logarithmic']);
  });
}, {owner: 'agolovko@datagrok.ai'});
