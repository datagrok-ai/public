// Cost-function unit fixtures — probe makeConstFunction directly (no NM);
// assert cost values for known inputs, error paths, and bounds-checker
// semantics.

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {LOSS, makeConstFunction, OptimizerInputsConfig, OutputTargetItem} from './imports';
import {makeExpDecayFunc, makeScalarPairFunc, makeThrowingFunc} from './script-fixtures';
import {rangeBound} from './utils';

category('ComputeUtils: Fitting / Cost functions (main)', () => {
  test('cost_mad_scalar_targets', async () => {
    // f(a, b) = (a + b, a * b). Targets y1 = 7, y2 = 12.
    // At (3, 4): y1=7, y2=12 → MAD = max(0, 0) = 0.
    // At (2, 4): y1=6, y2=8  → MAD = max(|7-6|, |12-8|) = 4.
    const func = makeScalarPairFunc();
    const bounds: OptimizerInputsConfig = {
      a: rangeBound(-10, 10, 'a'),
      b: rangeBound(-10, 10, 'b'),
    };
    const targets: OutputTargetItem[] = [
      {propName: 'y1', type: DG.TYPE.FLOAT, target: 7},
      {propName: 'y2', type: DG.TYPE.FLOAT, target: 12},
    ];
    const cost = makeConstFunction(LOSS.MAD, func, bounds, targets);

    const c1 = await cost(new Float64Array([3, 4]));
    expectFloat(c1!, 0, 1e-6);

    const c2 = await cost(new Float64Array([2, 4]));
    expectFloat(c2!, 4, 1e-6);
  });

  test('cost_rmse_scalar_targets', async () => {
    // f(a, b) = (a + b, a * b). Targets y1 = 7, y2 = 12.
    // RMSE uses relative scaling: ((cur - target)/target)^2 per scalar.
    // At (2, 4): err1 = (7-6)/7, err2 = (12-8)/12.
    //           sumSq = (1/7)^2 + (4/12)^2 = 0.0204 + 0.1111 = 0.1315
    //           rmse  = sqrt(0.1315 / 2) = 0.2566.
    const func = makeScalarPairFunc();
    const bounds: OptimizerInputsConfig = {
      a: rangeBound(-10, 10, 'a'),
      b: rangeBound(-10, 10, 'b'),
    };
    const targets: OutputTargetItem[] = [
      {propName: 'y1', type: DG.TYPE.FLOAT, target: 7},
      {propName: 'y2', type: DG.TYPE.FLOAT, target: 12},
    ];
    const cost = makeConstFunction(LOSS.RMSE, func, bounds, targets);

    const c1 = await cost(new Float64Array([3, 4]));
    expectFloat(c1!, 0, 1e-6);

    const expected = Math.sqrt(((1 / 7) ** 2 + (4 / 12) ** 2) / 2);
    const c2 = await cost(new Float64Array([2, 4]));
    expectFloat(c2!, expected, 1e-4);
  });

  test('cost_mad_dataframe_target', async () => {
    // Exp-decay sim with (a*, b*) = (2.5, 0.7). Target DF is the same exp-decay
    // sampled at fewer arg points; getIndices nearest-match aligns sim to target.
    // At (a*, b*) the cost should be 0 (sim matches target exactly at the
    // aligned indices).
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const bounds: OptimizerInputsConfig = {
      a: rangeBound(0, 5, 'a'),
      b: rangeBound(0, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const targets: OutputTargetItem[] = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];
    const cost = makeConstFunction(LOSS.MAD, func, bounds, targets);

    const c0 = await cost(new Float64Array([2.5, 0.7]));
    expectFloat(c0!, 0, 1e-3);

    const c1 = await cost(new Float64Array([3.0, 0.7]));
    expect(c1! > 0, true, `expected positive cost at perturbed point, got ${c1}`);
  });

  test('cost_rmse_dataframe_target', async () => {
    // Same setup, RMSE — at the exact (a*, b*) cost ≈ 0.
    const func = makeExpDecayFunc();
    const targetCall = func.prepare({a: 1.0, b: 0.5, N: 20});
    await targetCall.call();
    const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;

    const bounds: OptimizerInputsConfig = {
      a: rangeBound(0, 5, 'a'),
      b: rangeBound(0, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const targets: OutputTargetItem[] = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: targetDf,
      argName: 't',
      cols: [targetDf.col('y')!],
    }];
    const cost = makeConstFunction(LOSS.RMSE, func, bounds, targets);

    const c0 = await cost(new Float64Array([1.0, 0.5]));
    expectFloat(c0!, 0, 1e-3);

    const c1 = await cost(new Float64Array([1.5, 0.5]));
    expect(c1! > 0, true, 'expected positive cost at perturbed point');
  });

  test('cost_inconsistent_tables', async () => {
    // Target DataFrame is missing the arg column 't'. getErrors throws
    // InconsistentTables, which propagates out of the cost function.
    const func = makeExpDecayFunc();
    const badTarget = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('y', new Float32Array([1, 2, 3])),
    ]);

    const bounds: OptimizerInputsConfig = {
      a: rangeBound(0, 5, 'a'),
      b: rangeBound(0, 2, 'b'),
      N: {type: 'const', value: 20},
    };
    const targets: OutputTargetItem[] = [{
      propName: 'simulation',
      type: DG.TYPE.DATA_FRAME,
      target: badTarget,
      argName: 't',
      cols: [badTarget.col('y')!],
    }];
    const cost = makeConstFunction(LOSS.MAD, func, bounds, targets);

    let threw = false;
    try {
      await cost(new Float64Array([1, 0.5]));
    } catch (e: any) {
      threw = true;
      expect(/argument column|inconsistent/i.test(String(e.message ?? e)), true,
        `expected InconsistentTables, got: ${e.message ?? e}`);
    }
    expect(threw, true, 'cost did not throw on missing arg column');
  });

  test('cost_throws_in_funccall', async () => {
    // Body throws when a < 0. Cost function rejects (NM-side will catch).
    const func = makeThrowingFunc();
    const bounds: OptimizerInputsConfig = {
      a: rangeBound(-10, 10, 'a'),
      b: rangeBound(-10, 10, 'b'),
    };
    const targets: OutputTargetItem[] = [
      {propName: 'y', type: DG.TYPE.FLOAT, target: 5},
    ];
    const cost = makeConstFunction(LOSS.MAD, func, bounds, targets);

    // a >= 0: works
    const cOk = await cost(new Float64Array([2, 3]));
    expectFloat(cOk!, 0, 1e-6);

    // a < 0: rejects
    let threw = false;
    try {
      await cost(new Float64Array([-1, 3]));
    } catch (e: any) {
      threw = true;
      expect(/negative a/i.test(String(e.message ?? e)), true,
        `expected body throw, got: ${e.message ?? e}`);
    }
    expect(threw, true, 'cost did not propagate funcCall throw');
  });

  test('cost_returns_undefined_on_oob', async () => {
    // Bounds [0, 10]; call with a = -1. boundsChecker rejects → cost
    // resolves to undefined (not throws), and NM treats as costOutside.
    const func = makeScalarPairFunc();
    const bounds: OptimizerInputsConfig = {
      a: rangeBound(0, 10, 'a'),
      b: rangeBound(0, 10, 'b'),
    };
    const targets: OutputTargetItem[] = [
      {propName: 'y1', type: DG.TYPE.FLOAT, target: 7},
    ];
    const cost = makeConstFunction(LOSS.MAD, func, bounds, targets);

    // expectDeepEqual handles undefined/null directly via its nullPredicate;
    // the regular `expect(actual, undefined)` falls through to default `true`.
    const cOob = await cost(new Float64Array([-1, 5]));
    expectDeepEqual(cOob, undefined);

    const cIn = await cost(new Float64Array([2, 5]));
    expect(cIn !== undefined, true, 'expected defined cost for in-bounds input');
  });
});
