import * as DG from 'datagrok-api/dg';
import {category, test, before, delay} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import { snapshotCompare } from '../../../test-utils';

category('ComputeUtils: Driver config processing', async () => {
  before(async () => {});

  test('Process simple static config', async () => {
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
    await snapshotCompare(pconf, 'simple static config');
  });

  test('Process static config with dynamic pipelines', async () => {
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
          ]
        }
      ],
    }
    const pconf = await getProcessedConfig(config);
    await snapshotCompare(pconf, 'static config with dynamic pipelines');
  });

  test('Process config with globalId ref', async () => {
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
      ]
    };
    const pconf = await getProcessedConfig(config);
    await snapshotCompare(pconf, 'config with globalId ref');
  });
});
