import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {runOptimizerFinalized} from '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';
import type {OptimizerParams} from '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';
import {LOSS} from '@datagrok-libraries/compute-utils/function-views/src/fitting/constants';

function makeLineFunc(): DG.Func {
  return DG.Script.create([
    '//name: FitLine',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//output: dataframe line',
    '',
    'const tArr = new Float32Array([0, 1, 2, 3, 4]);',
    'line = DG.DataFrame.fromColumns([',
    '  DG.Column.fromFloat32Array(\'t\', tArr),',
    '  DG.Column.fromFloat32Array(\'y\', tArr.map((t) => a * t + b)),',
    ']);',
  ].join('\n'));
}

category('Fitting: RunOptimizer', () => {
  test('Returns the finalized calls through the function boundary', async () => {
    const func = makeLineFunc();
    const targetCall = func.prepare({a: 2, b: -1});
    await targetCall.call();
    const target = targetCall.getParamValue('line') as DG.DataFrame;

    const params: OptimizerParams = {
      lossType: LOSS.RMSE,
      func,
      inputBounds: {
        a: {type: 'changing', bottom: {name: 'a', type: 'value', value: -5}, top: {name: 'a', type: 'value', value: 5}},
        b: {type: 'changing', bottom: {name: 'b', type: 'value', value: -5}, top: {name: 'b', type: 'value', value: 5}},
      },
      outputTargets: [{propName: 'line', type: DG.TYPE.DATA_FRAME, target, argName: 't', cols: [target.col('y')!]}],
      samplesCount: 5,
      reproSettings: {reproducible: true, seed: 42},
    };

    const call = DG.Func.byName('Compute2:RunOptimizer').prepare({params});
    await call.call();
    const calls: DG.FuncCall[] = call.getOutputParamValue();

    expect(calls.length > 0, true, 'no calls returned');
    expectFloat(calls[0].inputs['a'], 2, 0.01, 'recovered a*');
    expectFloat(calls[0].inputs['b'], -1, 0.01, 'recovered b*');

    const direct = await runOptimizerFinalized(params);
    expect(calls.length, direct.calls.length, 'call count differs from runOptimizerFinalized');
    calls.forEach((c, i) => {
      expect(c.inputs['a'], direct.calls[i].inputs['a'], `call ${i} a`);
      expect(c.inputs['b'], direct.calls[i].inputs['b'], `call ${i} b`);
    });
  });
});
