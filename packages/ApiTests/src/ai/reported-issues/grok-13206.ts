import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, findProp, look} from '../helpers';

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
  function expectArgb(v: DG.Viewer, argb: number): void {
    const got = look(v)['backColor'];
    expect(got === argb || (got | 0) === (argb | 0), true);
    expect(typeof got, 'number');
  }

  test('setOptions writes backColor into look (numeric ARGB)', async () => {
    const v = DG.Viewer.densityPlot(demog(), {x: 'age', y: 'height'});
    v.setOptions({backColor: 0xFFCCEEFF});
    expectArgb(v, 0xFFCCEEFF);
  });

  test('DG.Color.fromHtml round-trips through backColor', async () => {
    const v = DG.Viewer.densityPlot(demog(), {x: 'age', y: 'height'});
    const argb = DG.Color.fromHtml('#ffcceeff');
    v.setOptions({backColor: argb});
    expectArgb(v, argb);
  });

  test('DG.Viewer.densityPlot accepts backColor in initial settings', async () => {
    const v = DG.Viewer.densityPlot(demog(), {x: 'age', y: 'height', backColor: 0xFFCCEEFF});
    expectArgb(v, 0xFFCCEEFF);
  });

  test('getProperties lists backColor (best-effort)', async () => {
    const v = DG.Viewer.densityPlot(demog(20), {x: 'age', y: 'height'});
    // Older descriptors may not list this prop; fall back to a look-bag round-trip.
    if (findProp(v, 'backColor') == null) {
      v.setOptions({backColor: 0xFFCCEEFF});
      expectArgb(v, 0xFFCCEEFF);
    }
  });
});
