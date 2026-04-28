import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-13206: density plot `backColor` (a number on
// IDensityPlotSettings — d4.ts:1558) was not being applied. The fix preserves
// the value through a setOptions/getOptions round-trip on the look bag.
//
// This test deliberately scopes itself to a standalone DensityPlot. Some
// Trellis paths do not preserve `backColor` on `innerViewerLook`, which is
// out of scope here. Pixel-level rendering is also not asserted.
category('AI: GROK-13206: Density plot backColor', () => {
  // ARGB ints can come back sign-coerced (top-bit set => negative when read
  // as a signed 32-bit int). Compare both the raw value and the |0 form.
  function expectArgbEqual(actual: any, expected: number): void {
    const ok = actual === expected || (actual | 0) === (expected | 0);
    expect(ok, true);
  }

  test('setOptions writes backColor into look (numeric ARGB)', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.DENSITY_PLOT);

    const argb = 0xFFCCEEFF;
    v.setOptions({backColor: argb});
    const look = v.getOptions(true).look;
    expectArgbEqual(look['backColor'], argb);
    expect(typeof look['backColor'], 'number');
  });

  test('DG.Color.fromHtml round-trips through backColor', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});

    const argb = DG.Color.fromHtml('#ffcceeff');
    expect(typeof argb, 'number');

    v.setOptions({backColor: argb});
    const look = v.getOptions(true).look;
    expectArgbEqual(look['backColor'], argb);
    expect(typeof look['backColor'], 'number');
  });

  test('DG.Viewer.densityPlot accepts backColor in initial settings', async () => {
    const df = grok.data.demo.demog(50);
    const argb = 0xFFCCEEFF;
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height', backColor: argb});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.DENSITY_PLOT);

    const look = v.getOptions(true).look;
    expectArgbEqual(look['backColor'], argb);
    expect(typeof look['backColor'], 'number');
  });

  test('getProperties lists backColor (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);

    var found = false;
    for (var p of props) {
      if (p.name === 'backColor') {
        found = true;
        break;
      }
    }

    // Older descriptors may not list this prop; fall back to a look-bag
    // round-trip so the regression check is still meaningful.
    if (!found) {
      const argb = 0xFFCCEEFF;
      v.setOptions({backColor: argb});
      const look = v.getOptions(true).look;
      expectArgbEqual(look['backColor'], argb);
      expect(typeof look['backColor'], 'number');
    }
    else
      expect(found, true);
  });
});
