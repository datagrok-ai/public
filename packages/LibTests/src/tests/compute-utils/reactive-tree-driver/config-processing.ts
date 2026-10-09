import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig, PipelineConfigurationStaticProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {snapshotCompare} from '../../../test-utils';
import {LoadedPipeline} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {PipelineStepConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-utils';

category('ComputeUtils: Driver config processing', async () => {
  before(async () => {});

  test('Process simple static config', async () => {
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
    await snapshotCompare(pconf, 'Process simple static config');
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
        },
      ],
    };
    const pconf = await getProcessedConfig(config);
    await snapshotCompare(pconf, 'Process static config with dynamic pipelines');
  });

  test('Normalize states shorthand', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'LibTests:TestAdd2',
          states: ['stepMeta1', {id: 'stepMeta2'}],
        },
      ],
      states: ['meta1', {id: 'meta2'}],
    };
    const pconf = await getProcessedConfig(config) as PipelineConfigurationStaticProcessed;
    expectDeepEqual(pconf.states, [{id: 'meta1'}, {id: 'meta2'}]);
    const step = pconf.steps[0] as PipelineStepConfigurationProcessed;
    expectDeepEqual(step.states, [{id: 'stepMeta1'}, {id: 'stepMeta2'}]);
  });

  test('Process config with globalId ref', async () => {
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
    };
    const pconf = await getProcessedConfig(config);
    await snapshotCompare(pconf, 'Process config with globalId ref');
  });

  test('A processing run prepares each function once to read its annotations', async () => {
    const config: PipelineConfiguration = {
      id: 'p',
      type: 'static',
      steps: [{id: 's1', nqName: 'LibTests:TestAdd2'}, {id: 's2', nqName: 'LibTests:TestAdd2'}],
      links: [{
        id: 'l', from: 'in_(template):s1/inputs(LibTests:TestAdd2)', to: 'out_(template):s2/inputs(LibTests:TestAdd2)',
      }],
    };
    const original = DG.Func.prototype.prepare;
    let prepared = 0;
    DG.Func.prototype.prepare = function(this: DG.Func, ...args: any[]) {
      if (this.nqName.toLowerCase() === 'libtests:testadd2')
        prepared++;
      return (original as any).apply(this, args);
    };
    try {
      const pconf = await getProcessedConfig(config);
      expectDeepEqual(prepared, 1, {prefix: 'One run'});
      const steps = (pconf as any).steps;
      expectDeepEqual(steps.map((step: any) => step.io.length), [3, 3], {prefix: 'Both steps get the io'});
      expectDeepEqual(steps[0].io !== steps[1].io, true, {prefix: 'Own io per step'});
      await getProcessedConfig(config);
      expectDeepEqual(prepared, 2, {prefix: 'Next run reads again'});
    } finally {
      DG.Func.prototype.prepare = original;
    }
  });
});
