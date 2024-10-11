import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {snapshotCompare} from '../../../test-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {PipelineInstanceConfig, PipelineStateStatic, StepFunCallState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {LoadedPipeline} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';

category('ComputeUtils: Driver state tree init', async () => {
  test('Process simple initial config', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Process simple initial config');
  });

  test('Process initial config with dynamic pipelines', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'pipelineSeq',
          type: 'sequential',
          stepTypes: [
            {
              id: 'stepMul',
              nqName: 'LibTests:TestMul2',
              friendlyName: 'mul',
            },
            {
              id: 'stepAdd',
              nqName: 'LibTests:TestAdd2',
              friendlyName: 'add',
            },
          ],
          initialSteps: [
            {
              id: 'stepMul',
            },
            {
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'stepAdd',
            },
          ],
        },
        {
          id: 'pipelinePar',
          type: 'parallel',
          stepTypes: [
            {
              id: 'stepMul',
              nqName: 'LibTests:TestMul2',
              friendlyName: 'mul',

            },
            {
              id: 'stepAdd',
              nqName: 'LibTests:TestAdd2',
              friendlyName: 'add',
            },
          ],
          initialSteps: [
            {
              id: 'stepAdd',
            },
            {
              id: 'stepAdd',
            },
            {
              id: 'stepMul',
            },
            {
              id: 'stepMul',
            },
          ],
        },
      ],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Process initial config with dynamic pipelines');
  });

  test('Process initial config with ref', async () => {
    const config: LoadedPipeline = {
      id: 'pipelineSeq',
      nqName: 'mockNqName',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
          friendlyName: 'mul',
        },
        {
          type: 'ref',
          provider: async (_p: any) => config,
        },
      ],
      initialSteps: [{
        id: 'stepMul',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Process initial config with ref');
  });
});

category('ComputeUtils: Driver init calls', async () => {
  test('Init function calls', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const state = tree.toState();
    (state as PipelineStateStatic<any>).steps.map(
      (x) => {
        if (!((x as StepFunCallState).funcCall instanceof DG.FuncCall))
          throw new Error(`funcCall is not an instance of DG.FuncCall`);

      },
    );
  });

  test('Init function calls options', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }],
    };
    const instanceConfig: PipelineInstanceConfig = {
      id: 'pipeline1',
      steps: [
        {
          id: 'step1',
        },
        {
          id: 'step2',
          values: {
            'a': 1,
            'b': 2,
          },
          inputRestrictions: {
            'a': 'disabled',
          },
        },
      ],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromInstanceConfig({instanceConfig, config: pconf});
    await tree.init().toPromise();
    const state = tree.toState();
    const fc = ((state as PipelineStateStatic<any>).steps[1] as StepFunCallState).funcCall!;
    expectDeepEqual(fc.inputs.a, 1);
    expectDeepEqual(fc.inputs.b, 2);
  });

});
