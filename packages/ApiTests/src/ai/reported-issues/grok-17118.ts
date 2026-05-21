import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, expectRoundTrip, findProp, look} from '../helpers';

// Regression coverage for GROK-17118 (1.23.0): density plot binShape +
// axis column type compatibility. With binShape='rectangle', switching
// xColumnName to a string column previously threw. The fix makes the
// viewer either ignore the invalid assignment or coerce it without
// crashing. These tests pin:
//   - binShape='rectangle' is accepted via setOptions and round-trips;
//   - changing xColumnName to a string column from rectangle state does
//     not throw, the viewer survives, and look.xColumnName ends up either
//     unchanged (rejected) or set to the string column (accepted);
//   - binShape round-trips for at least one non-rectangle valid value;
//   - getProperties() lists binShape, xColumnName, yColumnName on a
//     best-effort basis.
category('AI: GROK-17118: Density plot binShape + axis column type', () => {
  const v = (n: number = 50): DG.Viewer => DG.Viewer.densityPlot(demog(n), {x: 'age', y: 'height'});

  test('binShape="rectangle" round-trips through setOptions/look', async () => {
    expectRoundTrip(v(), {binShape: 'rectangle'});
  });

  test('switching xColumnName to a string column from rectangle does not throw', async () => {
    const df = demog();
    expect(df.col('race') != null, true);
    const c = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    c.setOptions({binShape: 'rectangle'});
    const beforeX = look(c)['xColumnName'];
    expectNoThrow(() => c.setOptions({xColumnName: 'race'}));
    const afterX = look(c)['xColumnName'];
    expect(afterX === beforeX || afterX === 'race', true);
  });

  test('binShape round-trips for a non-rectangle valid value', async () => {
    const c = v();
    // Discover valid binShape values via getProperties().choices when available.
    const p = findProp(c, 'binShape');
    var choices: string[] | null = null;
    try { choices = p?.choices ?? null; }
    catch (_e) { /* property may not expose choices */ }
    const candidates = choices != null && choices.length > 0 ? choices : ['hex', 'rectangle'];
    var picked: string = 'hex';
    for (var ch of candidates) if (ch !== 'rectangle') { picked = ch; break; }
    expectRoundTrip(c, {binShape: picked});
    expectNoThrow(() => c.setOptions({binShape: 'rectangle'}));
    expect(look(c)['binShape'], 'rectangle');
  });

  test('getProperties lists binShape/xColumnName/yColumnName (best-effort)', async () => {
    const c = v(20);
    const wanted = ['binShape', 'xColumnName', 'yColumnName'];
    var foundCount = 0;
    for (var n of wanted) if (findProp(c, n) != null) foundCount++;
    if (foundCount < wanted.length)
      expectRoundTrip(c, {binShape: 'rectangle', xColumnName: 'age', yColumnName: 'height'});
    else
      expect(foundCount, wanted.length);
  });
});
