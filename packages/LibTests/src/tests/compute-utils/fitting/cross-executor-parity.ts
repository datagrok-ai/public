// Main↔worker byte-identity gate.
//
// Same fixed seed, same NM settings, same fixtures: any divergence between
// the two arms shows up here. The canary for silent numerical divergence —
// a worker rewrite that produces IEEE-equivalent but bit-different results
// trips the test instead of shipping unnoticed.

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {LOSS, runOptimizer, OptimizerInputsConfig, OptimizerOutputsConfig} from './imports';
import {makeExpDecayFunc, makeMultiOutputFunc} from './script-fixtures';
import {rangeBound} from './utils';
import type {Extremum, OptimizationResult} from './imports';

function assertExtremumByteIdentical(a: Extremum, b: Extremum, label: string): void {
  expect(Object.is(a.cost, b.cost), true, `${label}: cost not bit-identical (${a.cost} vs ${b.cost})`);
  expect(a.iterCount, b.iterCount, `${label}: iterCount differs`);
  expect(a.point.length, b.point.length, `${label}: point.length differs`);
  for (let j = 0; j < a.point.length; ++j) {
    expect(Object.is(a.point[j], b.point[j]), true,
      `${label}: point[${j}] not bit-identical (${a.point[j]} vs ${b.point[j]})`);
  }
  const len = Math.min(a.iterCosts.length, b.iterCosts.length, a.iterCount);
  for (let j = 0; j < len; ++j) {
    expect(Object.is(a.iterCosts[j], b.iterCosts[j]), true,
      `${label}: iterCosts[${j}] not bit-identical`);
  }
}

function assertResultParity(main: OptimizationResult, wkr: OptimizationResult, label: string): void {
  expect(main.extremums.length, wkr.extremums.length, `${label}: extremum count differs`);
  // Sort by point lexicographically so we compare matching seed outcomes
  // even if the worker reordered them (work-stealing changes order).
  const sortKey = (e: Extremum) => Array.from(e.point).join(',');
  const m = [...main.extremums].sort((a, b) => sortKey(a).localeCompare(sortKey(b)));
  const w = [...wkr.extremums].sort((a, b) => sortKey(a).localeCompare(sortKey(b)));
  for (let i = 0; i < m.length; ++i)
    assertExtremumByteIdentical(m[i], w[i], `${label}[${i}]`);
}

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
});
