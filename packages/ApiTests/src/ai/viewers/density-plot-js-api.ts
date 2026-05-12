import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

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

  test('factory returns typed Viewer with DENSITY_PLOT type', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.DENSITY_PLOT);
    const look = v.getOptions(true).look;
    expect(look['xColumnName'], 'age');
    expect(look['yColumnName'], 'height');
  });

  test('bins integer round-trip within 1-200 range', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    v.setOptions({bins: 120});
    expect(v.props['bins'], 120);
    expect(v.getOptions(true).look['bins'], 120);
    v.setOptions({bins: 1});
    expect(v.getOptions(true).look['bins'], 1);
    v.setOptions({bins: 200});
    expect(v.getOptions(true).look['bins'], 200);
  });

  test('linearColorScheme number-array round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    const scheme = [0xFF112233, 0xFF445566, 0xFFAABBCC];
    v.setOptions({linearColorScheme: scheme});
    const look = v.getOptions(true).look;
    expect(Array.isArray(look['linearColorScheme']), true);
    expect((look['linearColorScheme'] as number[]).length, 3);
    expectArray(look['linearColorScheme'] as number[], scheme);
  });

  test('showXAxis/showYAxis/showColorScale combined bool round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    v.setOptions({showXAxis: false, showYAxis: false, showColorScale: false});
    const off = v.getOptions(true).look;
    expect(off['showXAxis'], false);
    expect(off['showYAxis'], false);
    expect(off['showColorScale'], false);
    v.setOptions({showXAxis: true, showYAxis: true, showColorScale: true});
    const on = v.getOptions(true).look;
    expect(on['showXAxis'], true);
    expect(on['showYAxis'], true);
    expect(on['showColorScale'], true);
  });

  test('x/yMin x/yMax numeric round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    v.setOptions({xMin: 10, xMax: 80, yMin: 120, yMax: 200});
    const look = v.getOptions(true).look;
    expect(look['xMin'], 10);
    expect(look['xMax'], 80);
    expect(look['yMin'], 120);
    expect(look['yMax'], 200);
  });

  test('invertColorScheme + colorTransformType logarithmic round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    v.setOptions({invertColorScheme: true, colorTransformType: 'logarithmic'});
    const look = v.getOptions(true).look;
    expect(look['invertColorScheme'], true);
    expect(look['colorTransformType'], 'logarithmic');
    v.setOptions({invertColorScheme: false, colorTransformType: 'linear'});
    const back = v.getOptions(true).look;
    expect(back['invertColorScheme'], false);
    expect(back['colorTransformType'], 'linear');
  });

  test('getProperties choices introspection on colorTransformType', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    const props = v.getProperties();
    var colorProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'colorTransformType') { colorProp = p; break; }
    expect(colorProp != null, true);
    expect(Array.isArray(colorProp!.choices), true);
    expect(colorProp!.choices.length >= 2, true);
    expect(colorProp!.choices.indexOf('linear') >= 0, true);
    expect(colorProp!.choices.indexOf('logarithmic') >= 0, true);
  });
}, {owner: 'agolovko@datagrok.ai'});
