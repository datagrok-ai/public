// Tests for runOptimizerFinalized — the consolidated post-processing API.
//
// Asserts:
//  - selectedExtremums is sorted ascending by cost.
//  - selectedExtremums.length <= allExtremums.length and == calls.length.
//  - legacy runOptimizer returns equivalent (extremums, calls).

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {LOSS, runOptimizer, runOptimizerFinalized,
  OptimizerInputsConfig, OptimizerOutputsConfig} from './imports';
import {makeExpDecayFunc} from './script-fixtures';
import {rangeBound} from './utils';

category('ComputeUtils: Fitting / Finalized API', () => {
  test('finalized result is sorted and consistent', async () => {
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

    const fin = await runOptimizerFinalized({
      lossType: LOSS.RMSE,
      func,
      inputBounds,
      outputTargets,
      samplesCount: 5,
      similarity: 10,
      reproSettings: {reproducible: true, seed: 42},
    });

    expect(fin.allExtremums.length > 0, true, 'no extremums returned');
    for (let i = 1; i < fin.allExtremums.length; ++i) {
      expect(fin.allExtremums[i - 1].cost <= fin.allExtremums[i].cost, true,
        `allExtremums not sorted at index ${i}`);
    }
    for (let i = 1; i < fin.selectedExtremums.length; ++i) {
      expect(fin.selectedExtremums[i - 1].cost <= fin.selectedExtremums[i].cost, true,
        `selectedExtremums not sorted at index ${i}`);
    }
    expect(fin.selectedExtremums.length <= fin.allExtremums.length, true,
      'selectedExtremums longer than allExtremums');
    expect(fin.selectedExtremums.length === fin.calls.length, true,
      'calls length does not match selectedExtremums');

    const best = fin.allExtremums[0];
    expectFloat(best.point[0], 2.5, 0.01, 'recovered a*');
    expectFloat(best.point[1], 0.7, 0.01, 'recovered b*');
  });

  test('legacy runOptimizer matches runOptimizerFinalized', async () => {
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
      similarity: 10,
      reproSettings: {reproducible: true, seed: 42},
    };

    const fin = await runOptimizerFinalized(baseArgs);
    const [legacyResult, legacyCalls] = await runOptimizer(baseArgs);

    expect(legacyResult.extremums.length === fin.allExtremums.length, true,
      'legacy extremums length differs from allExtremums');
    for (let i = 0; i < legacyResult.extremums.length; ++i) {
      expectFloat(legacyResult.extremums[i].cost, fin.allExtremums[i].cost, 1e-9,
        `legacy extremum cost differs at ${i}`);
    }
    expect(legacyCalls.length === fin.calls.length, true,
      'legacy calls length differs');
  });
});
