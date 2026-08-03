import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, findProp, look} from '../helpers';

// Regression coverage for GROK-13206: density plot backColor round-trips on the look bag.
category('AI: GROK-13206: Density plot backColor', () => {
  // ARGB ints can come back sign-coerced, so compare both raw and |0 forms.
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

  test('getProperties lists backColor', async () => {
    const v = DG.Viewer.densityPlot(demog(20), {x: 'age', y: 'height'});
    expect(findProp(v, 'backColor') != null, true);
  });
});
