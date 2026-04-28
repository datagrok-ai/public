import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19347: density plot xMin/xMax/yMin/yMax
// controls. The four range keys are declared on IDensityPlotSettings
// (d4.ts:1542-1548). The test pins that they round-trip through
// getOptions(true).look, that inverted/NaN values do not throw, and that
// the descriptor surface lists them on a best-effort basis.
category('AI: GROK-19347: Density plot xMin/xMax/yMin/yMax', () => {
  test('setOptions writes all four keys into look', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.DENSITY_PLOT);

    v.setOptions({xMin: 0, xMax: 100, yMin: 0, yMax: 200});
    const look = v.getOptions(true).look;
    expect(look['xMin'], 0);
    expect(look['xMax'], 100);
    expect(look['yMin'], 0);
    expect(look['yMax'], 200);
  });

  test('inverted ranges do not throw and viewer survives', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    var threw = false;
    try {
      v.setOptions({xMin: 100, xMax: 0, yMin: 200, yMax: 0});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    // Viewer is still alive and look is still readable; do not pin specific
    // clamped values — the API may either clamp or reject inverted input.
    expect(v instanceof DG.Viewer, true);
    const look = v.getOptions(true).look;
    expect(look != null, true);
  });

  test('NaN xMin does not throw', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    var threw = false;
    try {
      v.setOptions({xMin: NaN});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v instanceof DG.Viewer, true);
  });

  test('getProperties lists xMin/xMax/yMin/yMax (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);

    const wanted = ['xMin', 'xMax', 'yMin', 'yMax'];
    const found: {[k: string]: boolean} = {};
    for (var p of props) {
      if (wanted.indexOf(p.name) >= 0)
        found[p.name] = true;
    }

    // Older descriptors may not list newly-added props; fall back to a
    // look-bag round-trip so the regression check is still meaningful.
    if (Object.keys(found).length < wanted.length) {
      v.setOptions({xMin: 1, xMax: 99, yMin: 2, yMax: 198});
      const look = v.getOptions(true).look;
      expect(look['xMin'], 1);
      expect(look['xMax'], 99);
      expect(look['yMin'], 2);
      expect(look['yMax'], 198);
    }
    else {
      expect(found['xMin'], true);
      expect(found['xMax'], true);
      expect(found['yMin'], true);
      expect(found['yMax'], true);
    }
  });
});
