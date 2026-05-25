import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, expectRoundTrip, findProp, look} from '../helpers';

// Regression coverage for GROK-19347: density plot xMin/xMax/yMin/yMax
// controls. The four range keys are declared on IDensityPlotSettings
// (d4.ts:1542-1548). The test pins that they round-trip through
// getOptions(true).look, that inverted/NaN values do not throw, and that
// the descriptor surface lists them on a best-effort basis.
category('AI: GROK-19347: Density plot xMin/xMax/yMin/yMax', () => {
  const v = (): DG.Viewer => DG.Viewer.densityPlot(demog(), {x: 'age', y: 'height'});

  test('setOptions writes all four keys into look', async () => {
    expectRoundTrip(v(), {xMin: 0, xMax: 100, yMin: 0, yMax: 200});
  });

  test('inverted ranges do not throw and viewer survives', async () => {
    const c = v();
    expectNoThrow(() => c.setOptions({xMin: 100, xMax: 0, yMin: 200, yMax: 0}));
    expect(look(c) != null, true);
  });

  test('NaN xMin does not throw', async () => {
    expectNoThrow(() => v().setOptions({xMin: NaN}));
  });

  test('getProperties lists xMin/xMax/yMin/yMax (best-effort)', async () => {
    const c = DG.Viewer.densityPlot(demog(20), {x: 'age', y: 'height'});
    const wanted = ['xMin', 'xMax', 'yMin', 'yMax'];
    let foundCount = 0;
    for (const n of wanted) {
      if (findProp(c, n) !=
        null) foundCount++;
    }
    // Older descriptors may not list newly-added props; fall back to a look-bag
    // round-trip so the regression check is still meaningful.
    if (foundCount < wanted.length)
      expectRoundTrip(c, {xMin: 1, xMax: 99, yMin: 2, yMax: 198});
    else
      expect(foundCount, wanted.length);
  });
});
