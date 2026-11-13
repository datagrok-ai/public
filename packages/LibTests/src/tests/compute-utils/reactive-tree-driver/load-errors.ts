import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration, historyUtils} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {loadInstanceState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {PipelineStateStatic, StepFunCallSerializedState} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {deleteRuns, savedCallIds} from '../../../test-utils';

async function errorOf(fn: () => Promise<unknown>) {
  try {
    await fn();
  } catch (e: any) {
    return String(e?.message ?? e);
  }
  return undefined;
}

category('ComputeUtils: Driver load errors', async () => {
  test('A missing step function names the step', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline', type: 'static', steps: [{id: 'gone', nqName: 'LibTests:NoSuchFunction'}],
    };
    const error = await errorOf(() => getProcessedConfig(config));
    expectDeepEqual(error, `Step 'gone': function 'LibTests:NoSuchFunction' not found, it was removed or renamed`);
  });

  test('A missing provider is named', async () => {
    const error = await errorOf(() => callHandler('LibTests:NoSuchProvider', {}).toPromise());
    expectDeepEqual(error, `Function 'LibTests:NoSuchProvider' not found`);
  });

  test('A saved step run that does not match its function names the step and the missing params', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const metaCall = await tree.save().toPromise();
    let otherRun: DG.FuncCall | undefined;
    try {
      otherRun = await historyUtils.saveRun(await DG.Func.byName('LibTests:TestIONamesB').prepare({x: 1, y: 2}).call());
      const [state] = await loadInstanceState(metaCall!.id);
      const steps = (state as PipelineStateStatic<StepFunCallSerializedState, {}>).steps;
      (steps[0] as StepFunCallSerializedState).funcCallId = otherRun!.id;
      const loaded = StateTree.fromInstanceState({state, config: pconf, isReadonly: false});
      const error = await errorOf(() => loaded.init().toPromise());
      expectDeepEqual(error, `Step 'step1': saved run ${otherRun!.id} does not match 'LibTests:TestAdd2', ` +
        `missing inputs a, b. The function was changed or removed after the run was saved`);
    } finally {
      await deleteRuns([metaCall?.id, otherRun?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });
});
