// canHandle() — direct unit tests for the worker-arm gating predicate.
//
// canHandle decides whether `executor: 'auto'` routes a fit through the
// WorkerExecutor or falls back to MainExecutor. The contract is an explicit
// header annotation `//meta.workerSafe: true` on the JS-language DG.Script;
// no annotation means main-arm. These tests pin that contract.

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {canHandle, LOSS} from './imports';
import type {ExecutorArgs, OptimizerOutputsConfig} from './imports';
import {makeExpDecayFunc} from './script-fixtures';
import {rangeBound, defaultNmSettings, noEarlyStopping, reproSettings} from './utils';

function buildScript(lines: string[]): DG.Func {
  return DG.Script.create(lines.join('\n'));
}

function makeArgs(overrides: Partial<ExecutorArgs> = {}): ExecutorArgs {
  const targetDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', [0, 0.5, 1]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [1, 0.5, 0.25]),
  ]);
  const outputTargets: OptimizerOutputsConfig = [{
    propName: 'simulation',
    type: DG.TYPE.DATA_FRAME,
    target: targetDf,
    argName: 't',
    cols: [targetDf.col('y')!],
  }];
  return {
    objectiveFunc: async () => 0,
    inputsBounds: {a: rangeBound(0, 5, 'a'), b: rangeBound(0, 5, 'b')},
    samplesCount: 1,
    settings: defaultNmSettings(),
    reproSettings: reproSettings(),
    earlyStoppingSettings: noEarlyStopping(),
    func: makeExpDecayFunc(),
    outputTargets,
    lossType: LOSS.RMSE,
    ...overrides,
  };
}

category('ComputeUtils: Fitting / canHandle', () => {
  test('accepts JS script with meta.workerSafe: true', async () => {
    expect(canHandle(makeArgs()), true);
  });

  test('rejects JS script missing meta.workerSafe', async () => {
    const func = buildScript([
      '//name: NoAnnotation',
      '//language: javascript',
      '//input: double a',
      '//output: double y',
      '',
      'y = a;',
    ]);
    expect(canHandle(makeArgs({func})), false);
  });

  test('rejects JS script with meta.workerSafe: false', async () => {
    const func = buildScript([
      '//name: ExplicitFalse',
      '//language: javascript',
      '//meta.workerSafe: false',
      '//input: double a',
      '//output: double y',
      '',
      'y = a;',
    ]);
    expect(canHandle(makeArgs({func})), false);
  });

  test('accepts annotated script even when body textually mentions grok./ui.', async () => {
    // Pre-fix this body would be rejected by the regex even though the script
    // is opt-in. The annotation is now the sole gate — substring matches in
    // string literals or comments don't influence the decision.
    const func = buildScript([
      '//name: TextualMention',
      '//language: javascript',
      '//meta.workerSafe: true',
      '//input: double a',
      '//output: double y',
      '',
      '// note: do not call grok.shell or ui.div here',
      'y = a;',
    ]);
    expect(canHandle(makeArgs({func})), true);
  });

  test('rejects non-javascript language', async () => {
    const func = buildScript([
      '//name: PythonScript',
      '//language: python',
      '//meta.workerSafe: true',
      '//input: double a',
      '//output: double y',
      '',
      'y = a',
    ]);
    expect(canHandle(makeArgs({func})), false);
  });

  test('rejects when func is missing', async () => {
    expect(canHandle(makeArgs({func: undefined})), false);
  });

  test('rejects when outputTargets missing', async () => {
    expect(canHandle(makeArgs({outputTargets: undefined})), false);
  });

  test('rejects when lossType missing', async () => {
    expect(canHandle(makeArgs({lossType: undefined})), false);
  });
});
