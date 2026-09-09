import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {buildRefMap, getConfigByInstancePath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

const config: PipelineConfiguration = {
  id: 'root',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {
      id: 'nested',
      type: 'static',
      steps: [
        {id: 'stepA', nqName: 'LibTests:TestMul2'},
      ],
    },
    {
      id: 'dyn',
      type: 'dynamic',
      stepTypes: [
        {id: 'typeA', nqName: 'LibTests:TestAdd2'},
      ],
      initialSteps: [{id: 'typeA'}],
    },
  ],
};

async function resolve(path: string[]) {
  const pconf = await getProcessedConfig(config);
  const refMap = buildRefMap(pconf);
  return getConfigByInstancePath(path, pconf, refMap);
}

async function resolveError(path: string[]): Promise<string | undefined> {
  try {
    await resolve(path);
    return undefined;
  } catch (e: any) {
    return e.message;
  }
}

category('ComputeUtils: Driver config path resolution', async () => {
  test('Empty path resolves to the root config', async () => {
    expectDeepEqual((await resolve([])).id, 'root');
  });

  test('Top-level step path resolves', async () => {
    expectDeepEqual((await resolve(['step1'])).id, 'step1');
  });

  test('Nested pipeline and step paths resolve', async () => {
    expectDeepEqual((await resolve(['nested'])).id, 'nested');
    expectDeepEqual((await resolve(['nested', 'stepA'])).id, 'stepA');
  });

  test('Dynamic pipeline step type path resolves', async () => {
    expectDeepEqual((await resolve(['dyn', 'typeA'])).id, 'typeA');
  });

  test('Unknown segment throws', async () => {
    const msg = await resolveError(['nope']);
    expectDeepEqual(msg != null && msg.includes('nope'), true, {prefix: 'Error names the segment'});
  });

  test('One trailing segment past a leaf step throws', async () => {
    const msg = await resolveError(['step1', 'garbage']);
    expectDeepEqual(msg != null, true, {prefix: 'Path past a leaf rejected'});
    expectDeepEqual(msg != null && msg.includes('garbage'), true, {prefix: 'Error names the segment'});
  });

  test('Multiple trailing segments past a leaf throw with all segments listed', async () => {
    const msg = await resolveError(['step1', 'garbage', 'more']);
    expectDeepEqual(msg != null, true, {prefix: 'Path past a leaf rejected'});
    expectDeepEqual(msg != null && msg.includes('garbage') && msg.includes('more'), true, {prefix: 'Error names all segments'});
  });

  test('Trailing segment past a nested leaf throws', async () => {
    const msg = await resolveError(['nested', 'stepA', 'extra']);
    expectDeepEqual(msg != null, true, {prefix: 'Nested path past a leaf rejected'});
  });
});
