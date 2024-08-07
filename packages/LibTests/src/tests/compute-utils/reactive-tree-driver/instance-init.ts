import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {snapshotCompare} from '../../../test-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';

category('ComputeUtils: Driver state tree init', async () => {
  test('Process simple initial config', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2'
        },
        {
          id: 'step2',
          nqName: 'LibTests:TestMul2',
        }
      ],
      links: [{
        id: 'link1',
        from: 'step1/res',
        to: 'step2/a',
      }]
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig(pconf);
    const state = tree.toSerializedState(true);
    await snapshotCompare(state, 'simple initial config');
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
            }
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
          ]
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
            }
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
          ]
        }
      ],
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig(pconf);
    const state = tree.toSerializedState(true);
    await snapshotCompare(state, 'initial config with dynamic pipelines');
  });

  test('Process initial config with ref', async () => {
    const config: PipelineConfiguration = {
      id: 'pipelineSeq',
      globalId: 'MyPipelineGobalName',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
          friendlyName: 'mul',
        },
        {
          id: 'stepRef',
          type: 'ref',
          provider: async (_p: any) => ({globalId: 'MyPipelineGobalName', ...config}),
        }
      ],
      initialSteps: [{
        id: 'stepMul',
      }]
    };
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig(pconf);
    const state = tree.toSerializedState(true);
    await snapshotCompare(state, 'initial config with ref');
  });
});
