import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig, PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler} from '../../../test-utils';

const config: PipelineConfiguration = {
  id: 'pipeline',
  type: 'static',
  steps: [{id: 'step', nqName: 'LibTests:TestToleranceInputs'}],
};

category('ComputeUtils: Driver consistency tolerance', async () => {
  let testScheduler: TestScheduler;
  let pconf: PipelineConfigurationProcessed;

  before(async () => {
    testScheduler = createTestScheduler();
    pconf = await getProcessedConfig(config);
  });

  // a link assigns `assigned`, then the user edits the input to `edited`
  function isInconsistent(input: string, assigned: any, edited: any) {
    let inconsistent: boolean | undefined;
    testScheduler.run(({cold}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const step = () => tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      cold('-a').subscribe(() => step().getStateStore().setState(input, assigned));
      cold('--a').subscribe(() => step().getStateStore().editState(input, edited));
      cold('----a').subscribe(() => inconsistent = step().consistencyInfo$.value[input]?.inconsistent);
    });
    return inconsistent;
  }

  const df = (value: number) => DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', new Float32Array([value]))]);

  test('Without an annotation small changes stay within the default tolerance', async () => {
    expectDeepEqual(isInconsistent('plain', 1e-6, 5e-5), false);
    expectDeepEqual(isInconsistent('plain', 1e-6, 1e-3), true);
  });

  test('An absolute tolerance annotation replaces the default', async () => {
    expectDeepEqual(isInconsistent('abs', 1e-6, 5e-5), true);
    expectDeepEqual(isInconsistent('abs', 1e-6, 1e-6 + 1e-10), false);
  });

  test('A relative tolerance annotation scales with the value', async () => {
    expectDeepEqual(isInconsistent('plain', 1e6, 1e6 + 0.5), true, {prefix: 'Default flags float noise'});
    expectDeepEqual(isInconsistent('rel', 1e6, 1e6 + 0.5), false);
    expectDeepEqual(isInconsistent('rel', 1e6, 1e6 + 5), true);
  });

  test('A tolerance annotation applies to dataframe cells', async () => {
    expectDeepEqual(isInconsistent('df', df(1e-6), df(5e-5)), true);
    expectDeepEqual(isInconsistent('df', df(1e-6), df(1e-6)), false);
  });
});
