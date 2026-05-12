// Sync NM ↔ Async NM parity gate.
//
// optimizer-nelder-mead-sync.ts is a mechanical mirror of
// optimizer-nelder-mead.ts. Any drift between the two algorithms — even an
// arithmetic re-ordering that's IEEE-equivalent on paper — produces a
// different bit pattern and trips this test. Closure objectives only, no
// FuncCall plumbing, so the test is fast and the failure points straight
// at the NM code.

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {optimizeNM, optimizeNMSync} from './imports';
import type {Extremum} from './imports';
import {defaultNmSettings} from './utils';

function normalizeExtremum(e: Extremum): any {
  return {
    point: Array.from(e.point),
    cost: e.cost,
    iterCount: e.iterCount,
    iterCosts: e.iterCosts.slice(0, e.iterCount),
  };
}

const objectives: ReadonlyArray<{
  name: string;
  f: (x: Float64Array) => number;
  seed: Float64Array;
}> = [
  {
    name: 'quadratic_2d',
    f: (x) => (x[0] - 3) ** 2 + (x[1] + 1) ** 2,
    seed: new Float64Array([0, 0]),
  },
  {
    name: 'rosenbrock_2d',
    f: (x) => (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2,
    seed: new Float64Array([-1.2, 1]),
  },
  {
    name: 'himmelblau_2d',
    f: (x) => (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2,
    seed: new Float64Array([0.5, 0.5]),
  },
];

category('ComputeUtils: Fitting / NM sync↔async parity', () => {
  for (const {name, f, seed} of objectives) {
    test(`nm_sync_vs_async_${name}`, async () => {
      const settings = defaultNmSettings({maxIter: 200, tolerance: 1e-9});
      const syncExt = optimizeNMSync(f, seed, settings);
      const asyncExt = await optimizeNM(async (x) => f(x), seed, settings);
      expectDeepEqual(normalizeExtremum(syncExt), normalizeExtremum(asyncExt));
    });
  }

  test('nm_sync_vs_async_threshold_short_circuit', async () => {
    // Threshold parameter must short-circuit identically on both paths.
    const f = (x: Float64Array) => (x[0] - 3) ** 2 + (x[1] + 1) ** 2;
    const seed = new Float64Array([2.9, -1.1]);
    const settings = defaultNmSettings({maxIter: 200, tolerance: 1e-9});
    const threshold = 1e-2;
    const syncExt = optimizeNMSync(f, seed, settings, threshold);
    const asyncExt = await optimizeNM(async (x) => f(x), seed, settings, threshold);
    expectDeepEqual(normalizeExtremum(syncExt), normalizeExtremum(asyncExt));
    expect(syncExt.cost <= threshold, true, 'threshold not honored');
  });
});
