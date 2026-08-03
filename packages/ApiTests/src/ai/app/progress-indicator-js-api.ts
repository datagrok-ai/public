import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectFiresWithin, expectNoThrow} from '../helpers';

category('AI: App: Progress Indicator', () => {
  // Create a TaskBarProgressIndicator, run body, always close (close is no-throw).
  async function withPI(body: (pi: DG.TaskBarProgressIndicator) => Promise<void> | void): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('x');
    try {
      await body(pi);
    } finally {
      expectNoThrow(() => pi.close());
    }
  }

  test('create returns a TaskBarProgressIndicator', async () => {
    await withPI((pi) => {
      expect(pi instanceof DG.TaskBarProgressIndicator, true);
      expect(pi instanceof DG.ProgressIndicator, true);
    });
  });

  // A fresh indicator's percent default is implementation-defined (may be null until the first
  // update), so there is no deterministic default to assert; the update() test below is the real
  // percent/description contract.

  test('update(pct, desc) reflects in percent + description', async () => {
    await withPI((pi) => {
      pi.update(50, 'half');
      expect(pi.percent, 50);
      expect(pi.description, 'half');
    });
  });

  test('onProgressUpdated fires on update', async () => {
    await withPI((pi) => expectFiresWithin(pi.onProgressUpdated, () => pi.update(70, 'more')));
  });

  test('boundary update values are not clamped', async () => {
    await withPI((pi) => {
      pi.update(0, '');
      expect(pi.percent, 0);
      pi.update(100, '');
      expect(pi.percent, 100);
      pi.update(150, ''); // un-clamped
      expect(pi.percent, 150);
      pi.update(null as any, 'keep'); // update(null) preserves the prior percent
      expect(pi.percent, 150);
    });
  });

  test('close() no-throw and idempotent', async () => {
    const pi = DG.TaskBarProgressIndicator.create('x');
    expectNoThrow(() => pi.close());
    expectNoThrow(() => pi.close());
  });

  test('base ProgressIndicator.create + canceled getter', async () => {
    const pi = DG.ProgressIndicator.create();
    expect(pi instanceof DG.ProgressIndicator, true);
    expect(pi.canceled, false);
  });
}, {owner: 'agolovko@datagrok.ai'});
