// §4.3 of WORKERS_PLAN.md — end-to-end recovery tests via runOptimizer.
// Uses synthetic JS-language scripts with known (a*, b*, ...) and asserts
// the optimizer recovers them within tolerance.

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {LOSS, runOptimizer, OptimizerInputsConfig, OptimizerOutputsConfig} from './imports';
import {makeExpDecayFunc, makeLinearMixtureFunc, makeMultiOutputFunc, makeThrowingFunc}
  from './script-fixtures';
import {rangeBound} from './utils';

category('ComputeUtils: Fitting / End-to-end (main)', () => {
  test('e2e_exponential_decay', async () => {
    const func = makeExpDecayFunc();

    // Generate noiseless target with (a*, b*) = (2.5, 0.7).
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

    const [result, calls] = await runOptimizer({
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 8,
      reproSettings: {reproducible: true, seed: 42},
    });

    expect(result.extremums.length > 0, true, 'no extremums returned');
    const best = result.extremums.reduce((a, b) => a.cost < b.cost ? a : b);
    expectFloat(best.point[0], 2.5, 0.01, 'recovered a*');
    expectFloat(best.point[1], 0.7, 0.01, 'recovered b*');

    expect(calls.length > 0, true, 'no materialized FuncCalls');
    const sim = calls[0].getParamValue('simulation') as DG.DataFrame;
    expect(sim instanceof DG.DataFrame, true, 'best call did not produce a DataFrame');
    expect(sim.col('t') !== null && sim.col('y') !== null, true,
      'best call simulation missing expected columns');
  });

  test('e2e_linear_mixture', async () => {
    // Three-param recovery: y = a*t + b*sin(c*t).
    const func = makeLinearMixtureFunc();

    const targetCall = func.prepare({a: 1.5, b: 2.0, c: 3.0, N: 30});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;
    const targetMean = targetCall.getParamValue('mean') as number;

    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(0.5, 3, 'a'),
      b: rangeBound(0.5, 4, 'b'),
      c: rangeBound(1, 5, 'c'),
      N: {type: 'const', value: 30},
    };
    const outputTargets: OptimizerOutputsConfig = [
      {
        propName: 'simulation',
        type: DG.TYPE.DATA_FRAME,
        target: targetDf,
        argName: 't',
        cols: [targetDf.col('y')!],
      },
      {propName: 'mean', type: DG.TYPE.FLOAT, target: targetMean},
    ];

    const [result] = await runOptimizer({
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 16,
      reproSettings: {reproducible: true, seed: 42},
    });

    expect(result.extremums.length > 0, true, 'no extremums returned');
    const best = result.extremums.reduce((a, b) => a.cost < b.cost ? a : b);
    // Loose tolerance — 3-param NM with random restarts is finicky.
    expectFloat(best.point[0], 1.5, 0.05, 'recovered a*');
    expectFloat(best.point[1], 2.0, 0.05, 'recovered b*');
    expectFloat(best.point[2], 3.0, 0.05, 'recovered c*');
  });

  test('e2e_multi_output', async () => {
    // simulation has [t, y1, y2]. Target uses both y1 and y2 columns.
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

    const [result] = await runOptimizer({
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 8,
      reproSettings: {reproducible: true, seed: 42},
    });

    const best = result.extremums.reduce((a, b) => a.cost < b.cost ? a : b);
    expectFloat(best.point[0], 1.5, 0.01, 'recovered a*');
    expectFloat(best.point[1], 0.75, 0.01, 'recovered b*');
  });

  test('e2e_failing_initial_points', async () => {
    // makeThrowingFunc throws when a < 0. With a search range that includes
    // negatives, ~half the seeds fail. Verify fails DF is populated AND
    // good seeds still converge.
    const func = makeThrowingFunc();
    const inputBounds: OptimizerInputsConfig = {
      a: rangeBound(-5, 5, 'a'),
      b: rangeBound(0, 5, 'b'),
    };
    const outputTargets: OptimizerOutputsConfig = [
      {propName: 'y', type: DG.TYPE.FLOAT, target: 5},
    ];

    const [result] = await runOptimizer({
      lossType: LOSS.MAD,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 16,
      reproSettings: {reproducible: true, seed: 42},
    });

    // Some seeds threw → fails DF non-null with at least one row.
    expect(result.fails !== null, true, 'expected non-null fails DF');
    expect(result.fails!.rowCount > 0, true, 'expected at least one failed seed');
    // Good seeds still converge — extremums non-empty.
    expect(result.extremums.length > 0, true, 'expected at least one good extremum');
    const best = result.extremums.reduce((a, b) => a.cost < b.cost ? a : b);
    // Cost can be 0 only if a + b == 5 exactly; relax to "low".
    expect(best.cost < 1, true, `expected low cost, got ${best.cost}`);
  });
});
