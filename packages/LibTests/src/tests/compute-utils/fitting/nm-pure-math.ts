// Pure-math fixtures — lock down Nelder-Mead and the bounds/sampler
// invariants without any FuncCall plumbing. The fastest tests in the suite
// and the strongest correctness signal for a worker rewrite.

import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
// `executors.ts` selects the executor; tests with closure-only objectives
// fall back to main on the worker iteration.
import {ALL_EXECUTORS, runFitting} from './executors';
import {defaultNmSettings, reproSettings, noEarlyStopping, earlyStop,
  rangeBound, formulaBound, isMonotoneNonIncreasing} from './utils';
import {makeBoundsChecker} from './imports';
import type {Extremum, ValueBoundsData} from './imports';

const wrap = (
  f: (x: Float64Array) => number,
): ((x: Float64Array) => Promise<number | undefined>) => async (x) => f(x);

// expectDeepEqual doesn't natively recognize TypedArrays — coerce to plain
// arrays + truncate variable-length tail so structural diffing covers each
// extremum end-to-end with one call.
function normalizeExtremum(e: Extremum): any {
  return {
    point: Array.from(e.point),
    cost: e.cost,
    iterCount: e.iterCount,
    iterCosts: e.iterCosts.slice(0, e.iterCount),
  };
}

for (const executor of ALL_EXECUTORS) {
  category(`ComputeUtils: Fitting / NM pure math (${executor})`, async () => {
    test('nm_quadratic_2d', async () => {
      // f(x,y) = (x-3)^2 + (y+1)^2, minimum at (3,-1), cost 0.
      const objective = wrap((x) => (x[0] - 3) ** 2 + (x[1] + 1) ** 2);
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {x: rangeBound(-10, 10, 'x'), y: rangeBound(-10, 10, 'y')},
        samplesCount: 4,
        settings: defaultNmSettings({maxIter: 200}),
        reproSettings: reproSettings(7),
        earlyStoppingSettings: noEarlyStopping(),
      });
      expect(result.extremums.length > 0, true, 'no extremums returned');
      const best = result.extremums.reduce((a, b) => (a.cost < b.cost ? a : b));
      expectFloat(best.cost, 0, 1e-5);
      expectFloat(best.point[0], 3, 1e-2);
      expectFloat(best.point[1], -1, 1e-2);
    });

    test('nm_rosenbrock_2d', async () => {
      // (1-x)^2 + 100(y - x^2)^2, minimum at (1,1).
      const objective = wrap((x) => (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2);
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {x: rangeBound(-2, 2, 'x'), y: rangeBound(-2, 3, 'y')},
        samplesCount: 16,
        settings: defaultNmSettings({maxIter: 500, tolerance: 1e-9}),
        reproSettings: reproSettings(11),
        earlyStoppingSettings: noEarlyStopping(),
      });
      const best = result.extremums.reduce((a, b) => (a.cost < b.cost ? a : b));
      // Rosenbrock is hard for NM from a single start; with 16 starts at least
      // one should land near the global minimum.
      expect(best.cost < 1e-2, true, `best cost ${best.cost} > 1e-2`);
      expectFloat(best.point[0], 1, 0.1);
      expectFloat(best.point[1], 1, 0.1);
    });

    test('nm_himmelblau_finds_multiple_minima', async () => {
      // Himmelblau has 4 minima, all with cost 0. Coverage threshold (>=3 of 4)
      // is calibrated for seed=3 with 32 samples — re-tune both if the seed or
      // sample count changes.
      const objective = wrap((x) =>
        (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2);
      const minima = [
        [3, 2], [-2.805118, 3.131312], [-3.779310, -3.283186], [3.584428, -1.848126],
      ];
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {x: rangeBound(-5, 5, 'x'), y: rangeBound(-5, 5, 'y')},
        samplesCount: 32,
        settings: defaultNmSettings({maxIter: 300}),
        reproSettings: reproSettings(3),
        earlyStoppingSettings: noEarlyStopping(),
      });
      const found = new Set<number>();
      for (const e of result.extremums) {
        if (e.cost > 1e-3) continue;
        for (let i = 0; i < minima.length; ++i) {
          const m = minima[i];
          if (Math.hypot(e.point[0] - m[0], e.point[1] - m[1]) < 0.1)
            found.add(i);
        }
      }
      expect(found.size >= 3, true, `found ${found.size}/4 Himmelblau minima`);
    });

    test('nm_with_bounds_clip', async () => {
      // f(x) = (x-10)^2 with x in [0, 5]. NM itself is unconstrained — bounds
      // are enforced by the cost function returning undefined out-of-box (see
      // cost-functions.ts:33,62). Mirror that here: wrap the raw objective
      // with makeBoundsChecker. With the box hard-clipped, the optimum sits
      // at the boundary x=5.
      const bounds: Record<string, ValueBoundsData> = {x: rangeBound(0, 5, 'x')};
      const checker = makeBoundsChecker(bounds, ['x']);
      const objective = async (x: Float64Array): Promise<number | undefined> => {
        if (!checker(x)) return;
        return (x[0] - 10) ** 2;
      };
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: bounds,
        samplesCount: 8,
        settings: defaultNmSettings({maxIter: 200}),
        reproSettings: reproSettings(5),
        earlyStoppingSettings: noEarlyStopping(),
      });
      for (const e of result.extremums) {
        expect(e.iterCount > 0, true, 'iterCount must be positive for monotonicity check');
        expect(isMonotoneNonIncreasing(e.iterCosts, e.iterCount), true, 'iterCosts not monotone');
      }
      const best = result.extremums.reduce((a, b) => (a.cost < b.cost ? a : b));
      // Stronger §1 invariant: every reported point satisfies the bound checker.
      for (const e of result.extremums)
        expect(checker(e.point), true, `reported point ${Array.from(e.point)} fails boundsChecker`);
      // The unconstrained minimum is at x=10, so the constrained best must
      // sit at the upper edge.
      expectFloat(best.point[0], 5, 1e-2);
    });

    test('nm_formula_bounds', async () => {
      // For a single varied input x with formula bottom `a-1` and top `a+1`
      // (a is a fixed const at 5), starting points must satisfy 4 <= x <= 6.
      // f(x) = (x-5)^2.
      const objective = wrap((x) => (x[0] - 5) ** 2);
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {
          a: {type: 'const', value: 5},
          x: formulaBound('a-1', 'a+1', 'x'),
        },
        samplesCount: 4,
        settings: defaultNmSettings({maxIter: 100}),
        reproSettings: reproSettings(13),
        earlyStoppingSettings: noEarlyStopping(),
      });
      const best = result.extremums.reduce((a, b) => (a.cost < b.cost ? a : b));
      expectFloat(best.cost, 0, 1e-6);
      expectFloat(best.point[0], 5, 1e-3);
    });

    test('nm_formula_bounds_clips_optimum', async () => {
      // Constrained-vs-unconstrained mismatch: f(x) = (x - 10)^2 has its
      // unconstrained minimum at x=10, but the formula bound 0 <= x <= a-1
      // (a fixed at 5) clips x to [0, 4]. The constrained minimum is x=4
      // with cost 36. Closure-based tests don't gate bounds inside the
      // executor (see nm_with_bounds_clip), so we wrap the objective with
      // the shared `makeBoundsChecker` — directly exercising the formula
      // path of the bounds checker.
      const bounds: Record<string, ValueBoundsData> = {
        a: {type: 'const', value: 5},
        x: formulaBound('0', 'a - 1', 'x'),
      };
      const checker = makeBoundsChecker(bounds, ['x']);
      const objective = async (x: Float64Array): Promise<number | undefined> => {
        if (!checker(x)) return;
        return (x[0] - 10) ** 2;
      };
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: bounds,
        samplesCount: 8,
        settings: defaultNmSettings({maxIter: 200}),
        reproSettings: reproSettings(17),
        earlyStoppingSettings: noEarlyStopping(),
      });
      for (const e of result.extremums)
        expect(checker(e.point), true, `reported point ${Array.from(e.point)} fails boundsChecker`);
      const best = result.extremums.reduce((a, b) => (a.cost < b.cost ? a : b));
      expectFloat(best.point[0], 4, 1e-3);
      expectFloat(best.cost, 36, 1e-2);
    });

    test('nm_seeded_reproducibility', async () => {
      // Two runs with the same seed must produce byte-identical Extremum.point
      // and iterCosts. This is the load-bearing invariant for dual-executor
      // parity: any divergence between main and worker shows up here.
      const objective = wrap((x) => (x[0] - 3) ** 2 + (x[1] + 1) ** 2);
      const args = {
        objectiveFunc: objective,
        inputsBounds: {x: rangeBound(-10, 10, 'x'), y: rangeBound(-10, 10, 'y')},
        samplesCount: 5,
        settings: defaultNmSettings({maxIter: 200}),
        reproSettings: reproSettings(99),
        earlyStoppingSettings: noEarlyStopping(),
      };
      const a = await runFitting(executor, args);
      const b = await runFitting(executor, args);

      // Layer 1 — structural diff. floatChecker uses strict `<`, so tolerance
      // 0 rejects byte-identical floats; Number.MIN_VALUE is the smallest
      // positive that still admits equality.
      expectDeepEqual(
        a.extremums.map(normalizeExtremum),
        b.extremums.map(normalizeExtremum),
        {floatTolerance: Number.MIN_VALUE},
      );

      // Layer 2 — strict byte-identity via Object.is, catching ULP-level drift
      // that a positive tolerance would miss (e.g. reordered FMAs in a worker
      // rewrite producing IEEE-equivalent but bit-different results).
      expect(a.extremums.length, b.extremums.length, 'extremum count differs');
      for (let i = 0; i < a.extremums.length; ++i) {
        const ea = a.extremums[i];
        const eb = b.extremums[i];
        expect(Object.is(ea.cost, eb.cost), true, `cost[${i}] not bit-identical`);
        for (let j = 0; j < ea.point.length; ++j) {
          expect(Object.is(ea.point[j], eb.point[j]), true,
            `point[${i}][${j}] not bit-identical`);
        }
        for (let j = 0; j < ea.iterCount; ++j) {
          expect(Object.is(ea.iterCosts[j], eb.iterCosts[j]), true,
            `iterCosts[${i}][${j}] not bit-identical`);
        }
      }
    });

    test('nm_early_stopping_threshold', async () => {
      // f(x) = (x-3)^2 — easy. Threshold 1e-3, stopAfter 3, samplesCount 5:
      // every start converges easily so the loop must actually break early
      // at the third valid point (samples 4 and 5 never run). Five > three
      // guarantees the stopAfter boundary is exercised; if production
      // silently stopped honouring stopAfter, this test would fail.
      const objective = wrap((x) => (x[0] - 3) ** 2);
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {x: rangeBound(-10, 10, 'x')},
        samplesCount: 5,
        settings: defaultNmSettings({maxIter: 200}),
        reproSettings: reproSettings(2),
        earlyStoppingSettings: earlyStop(1e-3, 3, false),
      });
      expect(result.extremums.length, 3, 'expected exactly 3 extremums (early-break)');
      for (const e of result.extremums)
        expect(e.cost <= 1e-3, true, `extremum cost ${e.cost} above threshold`);
    });

    test('nm_above_threshold_fill', async () => {
      // With threshold=0 no point passes; with useAboveThresholdPoints=true
      // production code (optimizer.ts:96-97) appends ALL above-threshold
      // extremums — it does NOT cap at stopAfter. So result count == samples.
      // (Plan §1 entry initially claimed `=== stopAfter` — that was wrong;
      // production behaviour is "dump the whole pool".)
      //
      // Objective is (x-3)^2 + 1 so cost is provably >= 1 > 0 everywhere —
      // no seed can ever land in validExtremums, even by floating-point
      // accident. This eliminates a theoretical flake on lucky starts.
      const objective = wrap((x) => (x[0] - 3) ** 2 + 1);
      const samplesCount = 8;
      const stopAfter = 5;
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {x: rangeBound(-10, 10, 'x')},
        samplesCount,
        settings: defaultNmSettings({maxIter: 5}),
        reproSettings: reproSettings(17),
        earlyStoppingSettings: earlyStop(0, stopAfter, true),
      });
      expect(result.extremums.length, samplesCount,
        `expected all ${samplesCount} extremums (uncapped backfill), got ${result.extremums.length}`);
    });

    test('nm_iter_costs_monotone', async () => {
      // For any successful extremum, iterCosts[0..iterCount] must never
      // increase across iterations. This is a structural NM invariant that
      // any worker rewrite must preserve.
      const objective = wrap((x) => (x[0] - 1) ** 2 + (x[1] + 2) ** 2 + (x[2] - 0.5) ** 4);
      const result = await runFitting(executor, {
        objectiveFunc: objective,
        inputsBounds: {
          x: rangeBound(-5, 5, 'x'), y: rangeBound(-5, 5, 'y'), z: rangeBound(-5, 5, 'z'),
        },
        samplesCount: 4,
        settings: defaultNmSettings({maxIter: 150}),
        reproSettings: reproSettings(23),
        earlyStoppingSettings: noEarlyStopping(),
      });
      for (const e of result.extremums) {
        expect(e.iterCount > 0, true, 'iterCount must be positive for monotonicity check');
        expect(isMonotoneNonIncreasing(e.iterCosts, e.iterCount), true,
          `iterCosts not monotone for extremum cost=${e.cost}`);
      }
    });
  });
}
