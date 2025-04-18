import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/utils/src/test';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';

category('ComputeUtils: Driver structure check hook running', async () => {
  let testScheduler: TestScheduler;
  before(async () => {
    testScheduler = new TestScheduler((actual, expected) => {
      // console.log(actual, expected);
      expectDeepEqual(actual, expected);
    });
  });

  const config1: PipelineConfiguration = {
    id: 'pipeline1',
    type: 'sequential',
    stepTypes: [
      {
        id: 'step1',
        nqName: 'LibTests:TestAdd2',
      },
      {
        id: 'step2',
        nqName: 'LibTests:TestMul2',
      },
    ],
    structureCheck(state) {
      if (state.type === 'funccall')
        return undefined;
      if (state.steps.length !== 2)
        return {errors: [{description: 'Expecting exactly 2 steps'}]};
      if (state.steps[0].type === state.steps[1].type)
        return {errors: [{description: 'Expecting different steps'}]};
    }
  };

  const config2: PipelineConfiguration = {
    id: 'pipeline2',
    type: 'sequential',
    stepTypes: [{...config1, initialSteps: [{id: 'step1'}, {id: 'step2'}]}],
    initialSteps: [{id: 'pipeline1'}],
  };

  test('Run StructureCheck', async () => {
    const pconf = await getProcessedConfig(config1);
    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const state = tree.toState({skipFuncCalls: true});
      expectDeepEqual((state as any).structureCheckResults, {
        "errors": [
          {
            "description": "Expecting exactly 2 steps"
          }
        ]
      })
    });
  });

  test('Run StructureCheck in nested items', async () => {
    const pconf = await getProcessedConfig(config2);
    testScheduler.run((helpers) => {
      const {expectObservable} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const state = tree.toState({skipFuncCalls: true});
      expectDeepEqual((state as any).steps[0].structureCheckResults, {
        "errors": [
          {
            "description": "Expecting different steps"
          }
        ]
      })
    });
  });
});
