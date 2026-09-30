// Fitting with the sci-comp methods (L-BFGS-B, PSO, L-BFGS, Adam): the IOptimizer adapter,
// bounds passing (native / smooth penalty / undefined cutoff), the method registry, and
// the Diff Studio worker path. Ground truth: known minima of analytic objectives and models.

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {getIVP, getIvp2WebWorker, getPipelineCreator} from 'diff-grok';
import {runFitting} from './executors';
import {LOSS, METHOD, OPTIMIZERS, getFittedParamsFinalized, getOptimizerBounds, makeBoundsChecker,
  runOptimizerFinalized} from './imports';
import type {Extremum, OptimizerInputsConfig, OptimizerOutputsConfig, ValueBoundsData} from './imports';
import {makeExpDecayFunc} from './script-fixtures';
import {formulaBound, isMonotoneNonIncreasing, noEarlyStopping, rangeBound, reproSettings} from './utils';

const SCI_COMP_METHODS = [METHOD.LBFGSB, METHOD.PSO, METHOD.LBFGS, METHOD.ADAM];
const GRADIENT_METHODS = [METHOD.LBFGSB, METHOD.LBFGS, METHOD.ADAM];

function methodSettings(method: METHOD, overrides: Record<string, number> = {}): Map<string, number> {
  const m = new Map<string, number>();
  OPTIMIZERS.get(method)!.settingsOpts.forEach((opts, key) => m.set(key, overrides[key] ?? opts.default));
  return m;
}

const best = (extremums: Extremum[]) => extremums.reduce((a, b) => (a.cost < b.cost ? a : b));

/** The adapter maps sci-comp's cost history to the "Fitting profile" data */
function expectHistory(e: Extremum, method: METHOD): void {
  expect(e.iterCount > 0, true, `${method}: empty cost history`);
  expect(e.iterCosts.length, e.iterCount, `${method}: iterCosts length differs from iterCount`);
  expect(isMonotoneNonIncreasing(e.iterCosts, e.iterCount), true, `${method}: iterCosts not monotone`);
  expect(e.cost <= e.iterCosts[e.iterCount - 1], true, `${method}: cost above the last history value`);
}

// y = a * exp(-b * t) at t in [0, 5], scaled RMSE against a* = 2.5, b* = 0.7
const EXP_T = Array.from({length: 20}, (_, i) => i / 19 * 5);
const EXP_TARGET = EXP_T.map((t) => 2.5 * Math.exp(-0.7 * t));
const EXP_BOUNDS = {lower: Float64Array.of(0.5, 0.1), upper: Float64Array.of(5, 2)};
const expDecayLoss = (x: Float64Array) =>
  Math.sqrt(EXP_T.reduce((s, t, i) => s + ((x[0] * Math.exp(-x[1] * t) - EXP_TARGET[i]) / 2.5) ** 2, 0) / 20);

category('ComputeUtils: Fitting / sci-comp methods', () => {
  test('bounds_formula_dims_unbounded', async () => {
    const bounds = getOptimizerBounds({
      a: {type: 'const', value: 5},
      x: rangeBound(0, 5, 'x'),
      y: formulaBound('a - 1', 'a + 1', 'y'),
      z: rangeBound(-1, 1, 'z'),
    });
    const expected = [[0, 5], [-Infinity, Infinity], [-1, 1]];
    expect(bounds.lower.length, expected.length, 'bounds dimension');
    expected.forEach(([lo, hi], i) => {
      expect(bounds.lower[i], lo, `lower[${i}]`);
      expect(bounds.upper[i], hi, `upper[${i}]`);
    });
  });

  for (const method of SCI_COMP_METHODS) {
    test(`${method}: quadratic_2d`, async () => {
      // f(x,y) = (x-3)^2 + (y+1)^2, minimum at (3,-1), cost 0.
      const result = await runFitting('main', {
        objectiveFunc: async (x) => (x[0] - 3) ** 2 + (x[1] + 1) ** 2,
        inputsBounds: {x: rangeBound(-10, 10, 'x'), y: rangeBound(-10, 10, 'y')},
        samplesCount: 3,
        method,
        settings: methodSettings(method),
        reproSettings: reproSettings(7),
        earlyStoppingSettings: noEarlyStopping(),
      });
      expect(result.extremums.length, 3, 'extremums count');
      result.extremums.forEach((e) => expectHistory(e, method));
      const b = best(result.extremums);
      expect(b.cost < 1e-3, true, `best cost ${b.cost} >= 1e-3`);
      expectFloat(b.point[0], 3, 0.05);
      expectFloat(b.point[1], -1, 0.05);
    });

    test(`${method}: bounds_clip_undefined_cost`, async () => {
      // f(x) = (x-10)^2, x in [0, 5]: out-of-box cost is undefined (as in the model cost
      // functions) and the adapter penalizes it; the constrained optimum is x = 5.
      const bounds: Record<string, ValueBoundsData> = {x: rangeBound(0, 5, 'x')};
      const checker = makeBoundsChecker(bounds, ['x']);
      const result = await runFitting('main', {
        objectiveFunc: async (x) => checker(x) ? (x[0] - 10) ** 2 : undefined,
        inputsBounds: bounds,
        samplesCount: 4,
        method,
        settings: methodSettings(method),
        reproSettings: reproSettings(5),
        earlyStoppingSettings: noEarlyStopping(),
      });
      for (const e of result.extremums)
        expect(checker(e.point), true, `reported point ${Array.from(e.point)} fails boundsChecker`);
      expectFloat(best(result.extremums).point[0], 5, 1e-2);
    });

    test(`${method}: formula_bounds_clip_optimum`, async () => {
      // Formula bound 0 <= x <= a-1 (a = 5) is unbounded for the optimizer: the undefined cutoff
      // alone keeps x in [0, 4]. f(x) = (x-10)^2, the constrained optimum is x = 4.
      const bounds: Record<string, ValueBoundsData> = {
        a: {type: 'const', value: 5},
        x: formulaBound('0', 'a - 1', 'x'),
      };
      const checker = makeBoundsChecker(bounds, ['x']);
      const result = await runFitting('main', {
        objectiveFunc: async (x) => checker(x) ? (x[0] - 10) ** 2 : undefined,
        inputsBounds: bounds,
        samplesCount: 8,
        method,
        settings: methodSettings(method),
        reproSettings: reproSettings(17),
        earlyStoppingSettings: noEarlyStopping(),
      });
      for (const e of result.extremums)
        expect(checker(e.point), true, `reported point ${Array.from(e.point)} fails boundsChecker`);
      expectFloat(best(result.extremums).point[0], 4, 1e-2);
    });

    test(`${method}: bounds_are_passed`, async () => {
      // No undefined cutoff: only the bounds passed to the optimizer keep x <= 5.
      const {optimizer} = OPTIMIZERS.get(method)!;
      const objective = async (x: Float64Array) => (x[0] - 10) ** 2;
      const x0 = Float64Array.of(1);
      const bounds = {lower: Float64Array.of(0), upper: Float64Array.of(5)};

      const bounded = await optimizer(objective, x0, methodSettings(method), undefined, bounds);
      // Penalty methods settle at 5 + 5 / (1 + mu), L-BFGS-B exactly on the bound
      expectFloat(bounded.point[0], 5, 1e-2, 'bounded optimum');
      if (method === METHOD.LBFGSB)
        expect(bounded.point[0] <= 5, true, `L-BFGS-B left the box: ${bounded.point[0]}`);

      const free = await optimizer(objective, x0, methodSettings(method));
      expectFloat(free.point[0], 10, 0.1, 'unbounded optimum');
    });

    test(`${method}: threshold_stops_early`, async () => {
      const {optimizer} = OPTIMIZERS.get(method)!;
      const x0 = Float64Array.of(1, 1.5);
      const objective = async (x: Float64Array) => expDecayLoss(x);
      const full = await optimizer(objective, x0, methodSettings(method), undefined, EXP_BOUNDS);
      const stopped = await optimizer(objective, x0, methodSettings(method), 1e-2, EXP_BOUNDS);
      expect(full.cost < 1e-2, true, `full run cost ${full.cost} does not reach the threshold`);
      expect(stopped.cost <= 1e-2, true, `stopped at cost ${stopped.cost} above the threshold`);
      expect(stopped.iterCount < full.iterCount, true,
        `threshold did not stop early: ${stopped.iterCount} vs ${full.iterCount} iterations`);
    });

    test(`${method}: cancel_stops_after_first_iteration`, async () => {
      const {optimizer} = OPTIMIZERS.get(method)!;
      const res = await optimizer(async (x) => expDecayLoss(x), Float64Array.of(1, 1.5), methodSettings(method),
        undefined, EXP_BOUNDS, () => true);
      expect(res.iterCount <= 2, true, `canceled run made ${res.iterCount} iterations`);
    });

    test(`${method}: exp_decay_model`, async () => {
      // Real JS-script model through the public API
      const func = makeExpDecayFunc();
      const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
      await targetCall.call();
      const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

      const fin = await runOptimizerFinalized({
        lossType: LOSS.RMSE,
        func,
        inputBounds: {a: rangeBound(0.5, 5, 'a'), b: rangeBound(0.1, 2, 'b'), N: {type: 'const', value: 20}},
        outputTargets: [{
          propName: 'simulation',
          type: DG.TYPE.DATA_FRAME,
          target: targetDf,
          argName: 't',
          cols: [targetDf.col('y')!],
        }],
        samplesCount: 3,
        similarity: 10,
        method,
        reproSettings: {reproducible: true, seed: 42},
      });

      expect(fin.fails, null, 'fitting failures');
      fin.allExtremums.forEach((e) => expectHistory(e, method));
      expectFloat(fin.allExtremums[0].point[0], 2.5, 0.02, 'recovered a*');
      expectFloat(fin.allExtremums[0].point[1], 0.7, 0.01, 'recovered b*');
    }, {timeout: 120000});
  }

  for (const method of GRADIENT_METHODS) {
    test(`${method}: max_evaluations_cap`, async () => {
      const {optimizer} = OPTIMIZERS.get(method)!;
      const cap = 30;
      let evals = 0;
      const objective = async (x: Float64Array) => {
        ++evals;
        return expDecayLoss(x);
      };
      await optimizer(objective, Float64Array.of(1, 1.5),
        methodSettings(method, {maxIterations: 10000, maxFunctionEvaluations: cap}), undefined, EXP_BOUNDS);
      // The cap is checked between iterations: one iteration (a gradient plus a line search) may overshoot it
      expect(evals <= cap + 25, true, `${evals} evaluations with the cap ${cap}`);
    });
  }

  test('PSO: reproducible_from_initial_point', async () => {
    const {optimizer} = OPTIMIZERS.get(METHOD.PSO)!;
    const run = () => optimizer(async (x) => expDecayLoss(x), Float64Array.of(1, 1.5), methodSettings(METHOD.PSO),
      undefined, EXP_BOUNDS);
    const [a, b] = [await run(), await run()];
    expect(Object.is(a.cost, b.cost), true, 'cost not bit-identical');
    a.point.forEach((v, i) => expect(Object.is(v, b.point[i]), true, `point[${i}] not bit-identical`));
  });

  test('forced_worker_executor_runs_on_main', async () => {
    // Only Nelder-Mead has the worker arm: a forced 'worker' must not run it for another method
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;
    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 5, 'a'), b: rangeBound(0.1, 2, 'b'), N: {type: 'const', value: 20}};
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation', type: DG.TYPE.DATA_FRAME, target: targetDf, argName: 't', cols: [targetDf.col('y')!]}];
    const run = (executor: 'main' | 'worker') => runOptimizerFinalized({
      lossType: LOSS.RMSE, func, inputBounds, outputTargets, samplesCount: 2, similarity: 10,
      method: METHOD.LBFGSB, reproSettings: {reproducible: true, seed: 42}, executor,
    });
    const [main, worker] = [await run('main'), await run('worker')];
    expect(worker.allExtremums.length, main.allExtremums.length, 'extremums count');
    main.allExtremums.forEach((e, i) =>
      expect(Object.is(e.cost, worker.allExtremums[i].cost), true, `cost[${i}] differs from the main arm`));
  }, {timeout: 120000});

  test('unknown_method_rejected', async () => {
    let message = '';
    try {
      await runOptimizerFinalized({
        lossType: LOSS.RMSE,
        func: makeExpDecayFunc(),
        inputBounds: {a: rangeBound(0.5, 5, 'a'), b: rangeBound(0.1, 2, 'b'), N: {type: 'const', value: 20}},
        outputTargets: [],
        method: 'Gradient descent' as METHOD,
      });
    } catch (e) {
      message = e instanceof Error ? e.message : String(e);
    }
    expect(message.startsWith('Unknown fitting method'), true, `unexpected error: '${message}'`);
  });
});

// y' = -k * y, y(0) = y0: fitted in the Diff Studio worker, where the optimizers come from the registry
const DECAY_MODEL = `#name: Decay
#equations:
  dy/dt = -k * y
#inits:
  y = 2.5
#parameters:
  k = 0.7
#argument: t
  start = 0
  finish = 5
  step = 0.05`;

/** Analytic solution of DECAY_MODEL with the same inputs & output as its Diff Studio script, which calls
 *  the DiffStudio package. Used only to materialize the fitted calls. */
function makeDecaySolutionFunc(): DG.Func {
  return DG.Script.create([
    '//name: DecaySolution',
    '//language: javascript',
    '//input: double _t0 = 0',
    '//input: double _t1 = 5',
    '//input: double _h = 0.05',
    '//input: double y = 2.5',
    '//input: double k = 0.7',
    '//output: dataframe df',
    '',
    'const n = Math.round((_t1 - _t0) / _h) + 1;',
    'const tArr = new Float64Array(n);',
    'const yArr = new Float64Array(n);',
    'for (let i = 0; i < n; ++i) {',
    '  tArr[i] = _t0 + i * _h;',
    '  yArr[i] = y * Math.exp(-k * tArr[i]);',
    '}',
    'df = DG.DataFrame.fromColumns([',
    '  DG.Column.fromFloat64Array(\'t\', tArr),',
    '  DG.Column.fromFloat64Array(\'y\', yArr),',
    ']);',
  ].join('\n'));
}

category('ComputeUtils: Fitting / Diff Studio methods', () => {
  for (const method of OPTIMIZERS.keys()) {
    test(`${method}: decay_model_in_worker`, async () => {
      const ivp = getIVP(DECAY_MODEL);
      const t = Float64Array.from({length: 21}, (_, i) => i * 0.25);
      const target = DG.DataFrame.fromColumns([
        DG.Column.fromFloat64Array('t', t),
        DG.Column.fromFloat64Array('y', t.map((ti) => 2.5 * Math.exp(-0.7 * ti))),
      ]);
      const fixedInputs = {_t0: 0, _t1: 5, _h: 0.05};
      const inputBounds: OptimizerInputsConfig = {
        _t0: {type: 'const', value: fixedInputs._t0},
        _t1: {type: 'const', value: fixedInputs._t1},
        _h: {type: 'const', value: fixedInputs._h},
        y: rangeBound(0.5, 5, 'y'),
        k: rangeBound(0.1, 2, 'k'),
      };

      const fin = await getFittedParamsFinalized({
        loss: LOSS.RMSE,
        ivp,
        ivp2ww: getIvp2WebWorker(ivp),
        pipelineCreator: getPipelineCreator(ivp),
        method,
        settings: methodSettings(method),
        bounds: inputBounds,
        variedInputNames: ['y', 'k'],
        fixedInputs,
        argColName: 't',
        funcCols: [target.col('y')!],
        target,
        samplesCount: 4,
        reproSettings: reproSettings(5),
        earlyStoppingSettings: noEarlyStopping(),
        func: makeDecaySolutionFunc(),
        inputBounds,
        outputTargets: [{propName: 'df', type: DG.TYPE.DATA_FRAME, target, argName: 't', cols: [target.col('y')!]}],
        similarity: 10,
      });

      // A registry that fails to load in the worker rejects the whole run; a per-point failure lands in fails
      expect(fin.fails, null, 'in-worker fitting failures');
      expect(fin.allExtremums.length, 4, 'extremums count');
      expectFloat(fin.allExtremums[0].point[0], 2.5, 0.01, 'recovered y0*');
      expectFloat(fin.allExtremums[0].point[1], 0.7, 0.01, 'recovered k*');
    }, {timeout: 120000});
  }
});
