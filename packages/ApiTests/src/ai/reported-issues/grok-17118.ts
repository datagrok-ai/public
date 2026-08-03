import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, expectRoundTrip, findProp, look} from '../helpers';

// Regression coverage for GROK-17118: density plot binShape + axis column type compatibility.
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
    const p = findProp(c, 'binShape');
    let choices: string[] | null = null;
    try {
      choices = p?.choices ?? null;
    } catch (_e) {/* property may not expose choices */}
    const candidates = choices != null && choices.length > 0 ? choices : ['hex', 'rectangle'];
    let picked: string = 'hex';
    for (const ch of candidates) {
      if (ch !== 'rectangle') {
        picked = ch;
        break;
      }
    }
    expectRoundTrip(c, {binShape: picked});
    expectNoThrow(() => c.setOptions({binShape: 'rectangle'}));
    expect(look(c)['binShape'], 'rectangle');
  });

  test('getProperties lists binShape/xColumnName/yColumnName', async () => {
    const c = v(20);
    for (const n of ['binShape', 'xColumnName', 'yColumnName'])
      expect(findProp(c, n) != null, true);
  });
});
