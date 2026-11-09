import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import * as DG from 'datagrok-api/dg';
import {PipelineInstanceConfig, PipelineStateStatic} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {MockRun} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {
  HeadlessResult, RunState, StepRunState, WorkflowRunState, runWorkflowConfig, saveWorkflowConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {loadStateSnapshot} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import {loadInstanceState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, deleteRuns, savedCallIds} from '../../../test-utils';

const config: PipelineConfiguration = {
  id: 'pipeline',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links: [{id: 'link1', from: 'in:step1/res', to: 'out:step2/a'}],
};

const values = (step1: Record<string, any>, step2: Record<string, any>): PipelineInstanceConfig =>
  ({steps: [{id: 'step1', initialValues: step1}, {id: 'step2', initialValues: step2}]});

const steps = (state: RunState) => (state as WorkflowRunState).steps as StepRunState[];

const calculator = (calls: [string, Record<string, any>][]): MockRun => (configId, state) => {
  calls.push([configId, {a: state.a, b: state.b}]);
  if (state.a < 0)
    throw new Error('negative');
  return {res: configId === 'step2' ? state.a * state.b : state.a + state.b};
};

category('ComputeUtils: Driver headless', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  async function run(conf: PipelineConfiguration, instanceConfig: PipelineInstanceConfig, mockRun?: MockRun) {
    const pconf = await getProcessedConfig(conf);
    let result: HeadlessResult | undefined;
    testScheduler.run(() => {
      runWorkflowConfig({config: pconf, instanceConfig, mockMode: true, mockRun}).subscribe((res) => result = res);
    });
    return result!;
  }

  test('Runs every step with linked values', async () => {
    const calls: [string, Record<string, any>][] = [];
    const result = await run(config, values({a: 1, b: 2}, {b: 10}), calculator(calls));
    expectDeepEqual(result.isRun, true);
    expectDeepEqual(calls, [['step1', {a: 1, b: 2}], ['step2', {a: 3, b: 10}]]);
    expectDeepEqual(steps(result.state).map((step) => [step.id, step.status, step.initialValues]),
      [['step1', 'succeeded', {a: 1, b: 2, res: 3}], ['step2', 'succeeded', {a: 3, b: 10, res: 30}]]);
    expectDeepEqual(result.summary.map((item) => [item.kind, item.name, item.status]),
      [['workflow', 'pipeline', undefined], ['step', 'step1', 'succeeded'], ['step', 'step2', 'succeeded']]);
    expectDeepEqual(result.summary[0].rollup,
      {steps: 2, notLoaded: 0, failed: 0, outdated: 0, withErrors: 0, withWarnings: 0, inconsistent: 0});
  });

  test('Nothing runs when a ready step has validation errors', async () => {
    const calls: [string, Record<string, any>][] = [];
    const result = await run(config, values({a: 1}, {b: 10}), calculator(calls));
    expectDeepEqual(result.isRun, false);
    expectDeepEqual(calls, []);
    const [step1, step2] = steps(result.state);
    expectDeepEqual(step1.status, 'next error');
    expectDeepEqual(!!step1.validations.b?.errors?.length, true, {prefix: 'error on b'});
    expectDeepEqual(step2.status, 'pending');
    expectDeepEqual(result.summary[1].errors.length > 0, true, {prefix: 'summary errors'});
    expectDeepEqual(result.summary[0].rollup?.withErrors, 2, {prefix: 'step2 waits for a'});
  });

  test('A step waiting for an earlier one is checked during the run', async () => {
    const calls: [string, Record<string, any>][] = [];
    const result = await run(config, values({a: 1, b: 2}, {}), calculator(calls));
    expectDeepEqual(result.isRun, true);
    expectDeepEqual(calls.map(([id]) => id), ['step1']);
    const [step1, step2] = steps(result.state);
    expectDeepEqual(step1.status, 'succeeded');
    expectDeepEqual(step2.status, 'next error');
    expectDeepEqual(!!step2.validations.b?.errors?.length, true, {prefix: 'error on b'});
  });

  test('A failed step leaves the steps after it pending', async () => {
    const calls: [string, Record<string, any>][] = [];
    const result = await run(config, values({a: -1, b: 2}, {b: 10}), calculator(calls));
    expectDeepEqual(result.isRun, true);
    const [step1, step2] = steps(result.state);
    expectDeepEqual([step1.status, step1.runError], ['failed', 'Error: negative']);
    expectDeepEqual(step2.status, 'pending');
    expectDeepEqual(result.summary[1].runError, 'Error: negative');
    expectDeepEqual(result.summary[0].rollup?.failed, 1);
  });

  test('Pipeline validations sit on their workflow node', async () => {
    const conf: PipelineConfiguration = {
      ...config,
      links: [...config.links!, {
        id: 'check',
        type: 'pipelineValidator',
        from: 'res:step2/res',
        to: 'check',
        handler({controller}) {
          const res = controller.getFirst<number>('res');
          controller.setValidation('check', (res ?? 0) > 20 ? {errors: [{description: 'res too high'}]} : undefined);
        },
      }],
    };
    const result = await run(conf, values({a: 1, b: 2}, {b: 10}), calculator([]));
    expectDeepEqual(result.isRun, true);
    expectDeepEqual((result.state as WorkflowRunState).validations?.errors, [{description: 'res too high'}]);
    expectDeepEqual(result.summary[0].errors, ['res too high']);
  });

  test('Dynamic workflow steps come from the instance config', async () => {
    const conf: PipelineConfiguration = {
      id: 'pipeline',
      type: 'dynamic',
      stepTypes: [{id: 'step1', nqName: 'LibTests:TestAdd2'}, {id: 'step2', nqName: 'LibTests:TestMul2'}],
    };
    const calls: [string, Record<string, any>][] = [];
    const result = await run(conf, {steps: [
      {id: 'step2', initialValues: {a: 2, b: 3}},
      {id: 'step1', initialValues: {a: 4, b: 5}},
      {id: 'step2', initialValues: {a: 6, b: 7}},
    ]}, calculator(calls));
    expectDeepEqual(calls.map(([id]) => id), ['step2', 'step1', 'step2']);
    expectDeepEqual(steps(result.state).map((step) => step.status), ['succeeded', 'succeeded', 'succeeded']);
    expectDeepEqual(result.summary[0].rollup?.steps, 3);
  });

  test('A dynamic workflow item waiting for an earlier one does not stop the run', async () => {
    const conf: PipelineConfiguration = {
      id: 'pipeline',
      type: 'dynamic',
      stepTypes: [{id: 'step1', nqName: 'LibTests:TestAdd2'}, {id: 'step2', nqName: 'LibTests:TestMul2'}],
      links: [{id: 'link1', from: 'in:step1/res', to: 'out:step2/a', defaultRestrictions: {out: 'restricted'}}],
    };
    const calls: [string, Record<string, any>][] = [];
    const result = await run(conf, values({a: 1, b: 2}, {b: 10}), calculator(calls));
    expectDeepEqual(result.isRun, true);
    expectDeepEqual(calls, [['step1', {a: 1, b: 2}], ['step2', {a: 3, b: 10}]]);
  });

  test('A failed last step has the validations of the settled workflow', async () => {
    const conf: PipelineConfiguration = {
      ...config,
      links: [...config.links!, {
        id: 'checkA',
        type: 'validator',
        from: 'in:step2/a',
        to: 'out:step2/a',
        handler({controller}) {
          const a = controller.getFirst<number>('in');
          controller.setValidation('out', (a ?? 0) > 2 ? {warnings: [{description: 'a high'}]} : undefined);
        },
      }, {
        id: 'checkRes',
        type: 'pipelineValidator',
        from: 'a:step2/a',
        to: 'check',
        handler({controller}) {
          controller.setValidation('check', {notifications: [{description: `a is ${controller.getFirst('a')}`}]});
        },
      }],
    };
    const result = await run(conf, values({a: 1, b: 2}, {b: 10}), (configId, state) => {
      if (configId === 'step2')
        throw new Error('fail');
      return {res: state.a + state.b};
    });
    const [, step2] = steps(result.state);
    expectDeepEqual([step2.status, step2.runError], ['failed', 'Error: fail']);
    expectDeepEqual(step2.validations.a?.warnings, [{description: 'a high'}]);
    expectDeepEqual((result.state as WorkflowRunState).validations?.notifications, [{description: 'a is 3'}]);
  });

  test('A run state starts the workflow again with the steps run', async () => {
    const first = await run(config, values({a: 1, b: 2}, {b: 10}), calculator([]));
    const _instanceConfig: PipelineInstanceConfig = first.state;
    const calls: [string, Record<string, any>][] = [];
    const again = await run(config, first.state, calculator(calls));
    expectDeepEqual(calls, []);
    expectDeepEqual(steps(again.state).map((step) => [step.status, step.initialValues.res]),
      [['succeeded', 3], ['succeeded', 30]]);
  });

  test('A step editing a table in place leaves the caller table as it was', async () => {
    const conf: PipelineConfiguration =
      {id: 'pipeline', type: 'static', steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}]};
    const df = DG.DataFrame.fromCsv('x\n1');
    const result = await run(conf, {steps: [{id: 'step1', initialValues: {a: df, b: 2}}]}, (_id, state) => {
      state.a.set('x', 0, 5);
      return {res: 1};
    });
    expectDeepEqual(result.isRun, true);
    expectDeepEqual(df.get('x', 0), 1);
    expectDeepEqual(steps(result.state)[0].initialValues.a.get('x', 0), 5);
  });
});

category('ComputeUtils: Driver headless saving', async () => {
  test('A saved run state loads with its outputs', async () => {
    const pconf = await getProcessedConfig({...config, nqName: 'LibTests:MockWrapper1'} as PipelineConfiguration);
    const result = await runWorkflowConfig({config: pconf, instanceConfig: values({a: 1, b: 2}, {b: 10})}).toPromise();
    let metaCallId: string | undefined;
    let callIds: string[] = [];
    try {
      const metaCall = await saveWorkflowConfig({config: pconf, state: result!.state, metadata: {title: 'headless'}})
        .toPromise();
      metaCallId = metaCall.id;
      const [state] = await loadInstanceState(metaCall.id);
      callIds = savedCallIds(state);
      const loaded = await loadStateSnapshot(metaCall);
      const [step1, step2] = (loaded as PipelineStateStatic<any, {}>).steps;
      expectDeepEqual([step1.funcCall.outputs.res, step2.funcCall.outputs.res], [3, 30]);
      expectDeepEqual([step1.isOutputOutdated, step2.isOutputOutdated], [false, false]);
      expectDeepEqual(metaCall.friendlyName, 'headless');
    } finally {
      await deleteRuns([metaCallId, ...callIds]);
    }
  });
});
