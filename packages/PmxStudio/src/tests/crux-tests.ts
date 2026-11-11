import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {cruxPool} from '../crux';
import type {CruxModules, CruxResult} from '../crux';

/** Sets of workers that answer when told: `go()` resolves a set's pending calls; `ended` counts its terminations,
 *  `loads` the modules the pool had it load. From the `failFrom`-th set on, opening fails. `loaded` answers the
 *  pool's question of what a set has loaded (the pause waits on it). */
function fakeSets(failFrom = Infinity, loaded: () => Promise<string[]> = async () => ['stats']) {
  const made: {ended: number, loads: string[], go: () => void}[] = [];
  let opened = 0;
  const open = (): Promise<CruxModules> => {
    if (opened++ >= failFrom)
      return Promise.reject(new Error('no workers'));
    const waiting: (() => void)[] = [];
    const s = {ended: 0, loads: [] as string[], go: () => waiting.splice(0).forEach((f) => f())};
    made.push(s);
    return Promise.resolve({
      call: () => new Promise<CruxResult>((r) => waiting.push(() => r({} as CruxResult))),
      loaded: loaded as any, load: async (m: string) => void s.loads.push(m), terminate: async () => void s.ended++,
    });
  };
  return {made, open, opens: () => opened};
}

const IDLE = 20;

category('Crux pool', () => {
  test('sets open for overlapping calls, are not closed while a call runs, and shrink to one after a pause',
    async () => {
      const {made, open} = fakeSets();
      const pool = cruxPool(open, 4, IDLE);
      const one = pool.call('f', {});
      const two = pool.call('f', {});
      await delay(0);
      expect(pool.sets(), 2);
      made[0].go();
      await one;
      await delay(IDLE * 3); // a pause with a call pending
      expect(made.every((s) => s.ended === 0), true, 'a set was closed with a call running');
      made[1].go();
      await two;
      await delay(IDLE * 3);
      expect(made[0].ended + made[1].ended, 2, 'the old sets are not both closed');
      expect(pool.sets(), 1);
    });

  test('a call during the countdown keeps the sets', async () => {
    const {made, open} = fakeSets();
    const pool = cruxPool(open, 4, IDLE);
    const first = pool.call('f', {});
    await delay(0);
    made[0].go();
    await first;
    await delay(IDLE / 2);
    const next = pool.call('f', {});
    await delay(IDLE * 2);
    expect(made[0].ended, 0, 'closed under a call');
    made[0].go();
    await next;
  });

  test('a call that starts and ends while the pause reads the sets keeps them', async () => {
    let answer = () => {};
    const {made, open} = fakeSets(Infinity, () => new Promise((r) => answer = () => r(['stats'])));
    const pool = cruxPool(open, 4, IDLE);
    const first = pool.call('f', {});
    await delay(0);
    made[0].go();
    await first;
    await delay(IDLE * 2); // the pause has begun, and waits on loaded()
    const quick = pool.call('f', {});
    await delay(0);
    made[0].go();
    await quick;
    answer();
    await delay(0);
    expect(made[0].ended, 0, 'swapped although a call came in during the pause');
  });

  test('wake reopens the sets a pause closed, with their modules, and never opens what was not', async () => {
    const {made, open} = fakeSets();
    const pool = cruxPool(open, 4, IDLE);
    pool.wake(); // never opened: nothing to reopen
    expect(pool.sets(), 0);
    const calls = [pool.call('f', {}), pool.call('f', {}), pool.call('f', {})];
    await delay(0);
    made.forEach((s) => s.go());
    await Promise.all(calls);
    pool.wake(); // open: nothing opened
    expect(made.length, 3, 'a wake opened sets on an open pool');
    await delay(IDLE * 3);
    expect(pool.sets(), 1, 'the pause kept more than one set');
    pool.wake(); // at rest: back to three, each with the modules the sets had
    expect(pool.sets(), 3);
    expect(made.length, 4 + 2);
    await delay(0);
    expect(made.slice(3).map((s) => s.loads.join()).join('|'), 'stats|stats|stats', 'reopened without the modules');
    pool.wake(); // just woken: nothing
    expect(pool.sets(), 3);
    await delay(IDLE * 3);
    expect(pool.sets(), 1, 'no pause after a wake with no call');
    // at rest, a call, then a wake: the call took the warm set, the wake opens nothing
    const after = pool.call('f', {});
    await delay(0);
    made[made.length - 1].go();
    await after;
    const opened = made.length;
    pool.wake();
    expect(made.length, opened, 'a wake after a call reopened sets');
  });

  test('a wake with nothing to reopen does not cycle the set; open, it starts the countdown again', async () => {
    const {made, open} = fakeSets();
    const pool = cruxPool(open, 4, 100);
    const one = pool.call('f', {});
    await delay(0);
    made[0].go();
    await one;
    await delay(70);
    pool.wake(); // a gesture: someone is there
    await delay(60);
    expect(made.length, 1, 'rested although a gesture came 60 ms ago');
    await delay(150);
    expect(made.length, 2, 'never rested after the wake'); // (the pause swaps the set for a fresh one)
    pool.wake(); // at rest with the one set it had: nothing to reopen, no countdown
    await delay(250);
    expect(made.length, 2, 'a wake at rest with one set swapped it again');
  });

  test('a call that fails with no set open leaves no set opened after the pause', async () => {
    const {open, opens} = fakeSets(0); // no set opens
    const pool = cruxPool(open, 4, IDLE);
    await pool.call('f', {}).catch(() => null);
    await delay(IDLE * 3);
    expect(opens(), 1, 'a set was opened that no call asked for');
    expect(pool.sets(), 0);
  });

  test('a replacement set that fails to open leaves the pool, and the next call opens another', async () => {
    const {made, open} = fakeSets(1); // the first set opens; every later one fails
    const pool = cruxPool(open, 4, IDLE);
    const first = pool.call('f', {});
    await delay(0);
    made[0].go();
    await first;
    await delay(IDLE * 3); // swapped for a set that failed to open
    expect(made[0].ended, 1);
    expect(pool.sets(), 0, 'a dead set stayed in the pool');
    let failed = false;
    await pool.call('f', {}).catch(() => failed = true); // a new set is tried (and fails here), not the dead one
    expect(failed, true);
    expect(pool.sets(), 0);
  });
});
