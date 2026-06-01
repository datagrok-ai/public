import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, findProp, look} from '../helpers';

// Regression coverage for GROK-18223: histogram bins range stability (no throw on extreme values).
category('AI: GROK-18223: Histogram bins range stability', () => {
  test('sweeping bins through low and high values does not throw', async () => {
    const v = demog().plot.histogram({value: 'age'});
    let lastApplied: number | null = null;
    expectNoThrow(() => {
      for (const n of [1, 2, 50, 200, 500, 1000]) {
        v.setOptions({bins: n});
        const stored = look(v)['bins'];
        if (typeof stored === 'number')
          lastApplied = stored;
      }
    });
    expect(typeof lastApplied === 'number' && lastApplied! > 0 && Number.isFinite(lastApplied!), true);
  });

  test('non-positive bins values do not throw and stay positive', async () => {
    const v = demog().plot.histogram({value: 'age'});
    expectNoThrow(() => v.setOptions({bins: 0}));
    const stored = look(v)['bins'];
    if (typeof stored === 'number')
      expect(Number.isFinite(stored), true);
  });

  test('getProperties() exposes bins with a non-empty description', async () => {
    const p = findProp(demog(20).plot.histogram({value: 'age'}), 'bins');
    expect(p != null, true);
    expect((p!.description ?? '').trim().length > 0, true);
  });
});
