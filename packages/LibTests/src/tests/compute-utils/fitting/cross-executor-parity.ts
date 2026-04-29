// Main↔worker byte-identity gate.
//
// Same fixed seed, same NM settings, same fixtures: any divergence between
// the two arms shows up here. The canary for silent numerical divergence —
// a worker rewrite that produces IEEE-equivalent but bit-different results
// trips the test instead of shipping unnoticed.

import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {LOSS, runOptimizer, OptimizerInputsConfig, OptimizerOutputsConfig} from './imports';
import {makeExpDecayFunc, makeMultiOutputFunc, makeRefDfPassthroughFunc,
  makeDayjsFormatFunc, makeDateGetTimeFunc} from './script-fixtures';
import {rangeBound, formulaBound} from './utils';
import {assertResultParity} from './parity-assertions';

category('ComputeUtils: Fitting / Cross-executor parity', () => {
  test('parity_exp_decay', async () => {
    // Same seed-list (reproSettings.seed) feeds both arms; results must agree.
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 5, 'a'),
      b: rangeBound(0.1, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 42},
    };

    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});

    assertResultParity(mainRes, wkrRes, 'exp_decay');
  });

  test('parity_formula_bound_clips', async () => {
    // Formula bound that strictly clips the unconstrained optimum: target
    // params are (a=2.5, b=0.7), but `b` is bounded to [0.05*a, 0.1*a],
    // i.e. at most 0.5 across the whole feasible `a` range. The
    // constrained optimum lives somewhere on the upper b-edge; both arms
    // must converge to byte-identical points. A worker that ignores the
    // formula bound walks toward (2.5, 0.7) and trips assertResultParity.
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 5, 'a'),
      b: formulaBound('0.05 * a', '0.1 * a', 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 23},
    };

    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});

    assertResultParity(mainRes, wkrRes, 'formula_bound_clips');
  });

  test('parity_const_dataframe_input', async () => {
    // Const DataFrame in inputBounds: the script body folds refDf.rowCount
    // into the simulation, so a worker that loses the DataFrame produces
    // a different cost. Pre-fix this also throws DataCloneError before
    // the worker runs at all.
    const func = makeRefDfPassthroughFunc();
    const refDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'k', [1, 2, 3, 4, 5, 6, 7]),
    ]);
    const targetCall = func.prepare({refDf, a: 0.4, N: 10});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      refDf: {type: 'const', value: refDf},
      a: rangeBound(0.05, 1, 'a'),
      N: {type: 'const', value: 10},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 31},
    };
    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});
    assertResultParity(mainRes, wkrRes, 'const_dataframe_input');
  });

  test('parity_const_dayjs_format', async () => {
    // Const Dayjs in inputBounds: body calls t0.year() — a method that only
    // exists on Dayjs, not on the ISO string the worker would see pre-fix.
    const func = makeDayjsFormatFunc();
    const t0 = dayjs('2024-06-15T00:00:00Z');
    const targetCall = func.prepare({t0, a: 1.5});
    await targetCall.call();
    const target = targetCall.getParamValue('y') as number;

    const inputBounds: OptimizerInputsConfig = {
      t0: {type: 'const', value: t0},
      a: rangeBound(0, 5, 'a'),
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'y',
      type: DG.TYPE.FLOAT,
      target,
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 41},
    };
    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});
    assertResultParity(mainRes, wkrRes, 'const_dayjs_format');
  });

  test('parity_const_date_gettime', async () => {
    // Const raw JS Date in inputBounds: body calls t0.getTime(), which only
    // works on Date (not Dayjs, not ISO string). Exercises the 'date'
    // reification branch.
    const func = makeDateGetTimeFunc();
    const t0 = new Date('2024-06-15T00:00:00Z');
    const targetCall = func.prepare({t0, a: 1.5});
    await targetCall.call();
    const target = targetCall.getParamValue('y') as number;

    const inputBounds: OptimizerInputsConfig = {
      t0: {type: 'const', value: t0},
      a: rangeBound(0, 5, 'a'),
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'y',
      type: DG.TYPE.FLOAT,
      target,
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 53},
    };
    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});
    assertResultParity(mainRes, wkrRes, 'const_date_gettime');
  });

  test('parity_formula_bound_references_dataframe', async () => {
    // Formula bound that calls a method on a const DataFrame: the upper
    // bound for `b` is `refDf.rowCount * 0.1`. Both arms must evaluate
    // the formula identically.
    const func = makeExpDecayFunc();
    const refDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'k', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
    ]);
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      refDf: {type: 'const', value: refDf},
      a: rangeBound(0.5, 5, 'a'),
      b: formulaBound('0.1', 'refDf.rowCount * 0.1', 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 67},
    };
    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});
    assertResultParity(mainRes, wkrRes, 'formula_bound_dataframe');
  });

  test('parity_formula_bound_references_dayjs', async () => {
    // Formula bound that calls a method on a const Dayjs: upper bound for
    // `b` is `t0.year() / 10000` ≈ 0.2024. Both arms must agree.
    const func = makeExpDecayFunc();
    const t0 = dayjs('2024-06-15T00:00:00Z');
    const targetCall = func.prepare({a: 2.5, b: 0.15, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      t0: {type: 'const', value: t0},
      a: rangeBound(0.5, 5, 'a'),
      b: formulaBound('0.1', 't0.year() / 10000', 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 71},
    };
    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});
    assertResultParity(mainRes, wkrRes, 'formula_bound_dayjs');
  });

  test('parity_multi_output', async () => {
    const func = makeMultiOutputFunc();
    const targetCall = func.prepare({a: 1.5, b: 0.75, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0, 5, 'a'),
      b: rangeBound(0, 5, 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y1')!, targetDf.col('y2')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 4,
      reproSettings: {reproducible: true, seed: 7},
    };

    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});

    assertResultParity(mainRes, wkrRes, 'multi_output');
  });

  test('parity_early_stop_index_order', async () => {
    // Many candidate seeds, generous threshold so most/all converge to
    // cost ≤ threshold. With early stopping and stopAfter=3, the main arm
    // picks the *first 3 valid seeds in index order*; pre-fix the worker
    // arm picked 3 in completion order, so the chosen extremum SET (and
    // thus their costs) diverged. This is the canary for the index-ordered
    // selection rule.
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 5, 'a'),
      b: rangeBound(0.1, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 10,
      reproSettings: {reproducible: true, seed: 101},
      // Threshold deliberately loose so all 10 seeds count as "valid".
      earlyStoppingSettings: {useEarlyStopping: true, costFuncThreshold: 1.0,
        stopAfter: 3, useAboveThresholdPoints: false},
    };

    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});

    expect(mainRes.extremums.length, 3, 'early_stop_index_order: main should yield 3 extremums');
    assertResultParity(mainRes, wkrRes, 'early_stop_index_order');
  });

  test('parity_early_stop_with_drain', async () => {
    // stopAfter=2 fires while higher-index seeds may still be in flight.
    // The worker arm must record those in-flight replies in their index
    // slots (not drop them) and the finalize loop must select the first
    // 2 valid by index. samplesCount=8 + small stopAfter forces the
    // drain-and-pick path under load.
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 5, 'a'),
      b: rangeBound(0.1, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 8,
      reproSettings: {reproducible: true, seed: 207},
      earlyStoppingSettings: {useEarlyStopping: true, costFuncThreshold: 1.0,
        stopAfter: 2, useAboveThresholdPoints: false},
    };

    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});

    expect(mainRes.extremums.length, 2, 'early_stop_with_drain: main should yield 2 extremums');
    assertResultParity(mainRes, wkrRes, 'early_stop_with_drain');
  });

  test('parity_above_threshold_spillover', async () => {
    // Tight threshold + useAboveThresholdPoints=true exercises the
    // spillover branch in finalize: the worker arm must concatenate
    // valid + above-threshold the same way EarlyStopTracker.finalize does
    // for the main arm.
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 5, 'a'),
      b: rangeBound(0.1, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const outputTargets: OptimizerOutputsConfig = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];

    const baseArgs = {
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 6,
      reproSettings: {reproducible: true, seed: 313},
      earlyStoppingSettings: {useEarlyStopping: true, costFuncThreshold: 1e-12,
        stopAfter: 4, useAboveThresholdPoints: true},
    };

    const [mainRes] = await runOptimizer({...baseArgs, executor: 'main'});
    const [wkrRes] = await runOptimizer({...baseArgs, executor: 'worker'});

    assertResultParity(mainRes, wkrRes, 'above_threshold_spillover');
  });
});
