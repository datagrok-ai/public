import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {
  HeadlessResult, StepRunState, WorkflowRunState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {loadStateSnapshot} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import {deleteSavedWorkflow} from './utils';

// through the function boundary, as compute-api calls them
async function callFunc(nqName: string, params: Record<string, any>) {
  const call = DG.Func.byName(nqName).prepare(params);
  await call.call();
  return call.getOutputParamValue();
}

category('Workflows: headless run', () => {
  test('A workflow run without UI returns its values and saves to history', async () => {
    const workflow = {nqName: 'Compute2:MockPipeline2', version: '1.0'};
    const result: HeadlessResult = await callFunc('Compute2:RunWorkflow', {...workflow,
      instanceConfig: {steps: [
        {id: 'stepAdd', initialValues: {a: 1, b: 2}},
        {id: 'stepMul', initialValues: {a: 10, b: 4}},
      ]},
    });
    expectDeepEqual(result.isRun, true);
    const steps = (result.state as WorkflowRunState).steps as StepRunState[];
    expectDeepEqual(steps.map((step) => [step.id, step.status, step.initialValues.result]),
      [['stepAdd', 'succeeded', 3], ['stepMul', 'succeeded', 40]]);
    expectDeepEqual(result.summary[0].rollup?.steps, 2);

    const run: DG.FuncCall = await callFunc('Compute2:SaveWorkflowRun',
      {...workflow, state: result.state, metadata: {title: 'headless'}});
    try {
      expectDeepEqual(run.friendlyName, 'headless');
      const saved = (await loadStateSnapshot(run)) as any;
      expectDeepEqual(saved.steps.map((step: any) => [step.funcCall.outputs.result, step.isOutputOutdated]),
        [[3, false], [40, false]]);
    } finally {
      await deleteSavedWorkflow(run.id);
    }
  });
});
