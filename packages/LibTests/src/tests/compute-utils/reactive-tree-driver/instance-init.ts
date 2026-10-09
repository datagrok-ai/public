import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {snapshotCompare} from '../../../test-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {isFuncCallSerializedState, normalizePipelineInstanceConfig, PipelineInstanceConfig, PipelineInstanceConfigInput, PipelineSerializedState, PipelineStateStatic, StepFunCallState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {LoadedPipeline} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';

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

  test('Self-referencing initialSteps throw a cycle error', async () => {
    const config: LoadedPipeline = {
      id: 'pipelineSeq',
      nqName: 'mockNqName',
      type: 'sequential',
      stepTypes: [
        {
          id: 'stepMul',
          nqName: 'LibTests:TestMul2',
        },
        {
          type: 'ref',
          provider: async (_p: any) => config,
        },
      ],
      initialSteps: [{id: 'stepMul'}, {id: 'pipelineSeq'}],
    };
    const pconf = await getProcessedConfig(config);
    let msg: string | undefined;
    try {
      StateTree.fromPipelineConfig({config: pconf});
    } catch (e) {
      msg = String(e);
    }
    expectDeepEqual(msg?.includes('Initial config cycle') ?? false, true);
  });

  test('Process initial config with additional data', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper5', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Process initial config with additional data');
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
    (state as PipelineStateStatic<any, {}>).steps.map(
      (x) => {
        if (!((x as StepFunCallState).funcCall instanceof DG.FuncCall))
          throw new Error(`funcCall is not an instance of DG.FuncCall`);
      },
    );
  });

  test('Step config initial values with an output start the step as run', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2', initialValues: {a: 1}},
        {id: 'step2', nqName: 'LibTests:TestAdd2', initialValues: {a: 1, res: 5}},
      ],
    };
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config)});
    await tree.init().toPromise();
    const step1 = tree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
    const step2 = tree.nodeTree.getItem([{idx: 1}]) as FuncCallNode;
    expectDeepEqual(step1.getStateStore().getState('a'), 1, {prefix: 'step1 a'});
    expectDeepEqual(step1.funcCallState$.value?.isOutputOutdated, true, {prefix: 'step1 not run'});
    expectDeepEqual(step2.getStateStore().getState('a'), 1, {prefix: 'step2 a'});
    expectDeepEqual(step2.getStateStore().getState('res'), 5, {prefix: 'step2 res'});
    expectDeepEqual(step2.funcCallState$.value?.isOutputOutdated, false, {prefix: 'step2 run'});
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
          initialValues: {
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
    const fc = ((state as PipelineStateStatic<any, {}>).steps[1] as StepFunCallState).funcCall!;
    expectDeepEqual(fc.inputs.a, 1);
    expectDeepEqual(fc.inputs.b, 2);
  });
});

category('ComputeUtils: Driver instance config defaults', async () => {
  const config: PipelineConfiguration = {
    id: 'root',
    type: 'static',
    steps: [
      {id: 'step1', nqName: 'LibTests:TestAdd2'},
      {
        id: 'nestedStatic',
        type: 'static',
        steps: [
          {id: 'inner1', nqName: 'LibTests:TestAdd2'},
          {
            id: 'innerDyn',
            type: 'parallel',
            stepTypes: [{id: 'add', nqName: 'LibTests:TestAdd2'}, {id: 'mul', nqName: 'LibTests:TestMul2'}],
            initialSteps: ['mul'],
          },
        ],
      },
      {
        id: 'nestedDyn',
        type: 'parallel',
        stepTypes: [{id: 'add', nqName: 'LibTests:TestAdd2'}, {id: 'mul', nqName: 'LibTests:TestMul2'}],
        initialSteps: ['add', 'mul'],
      },
    ],
  };

  const shape = (s: PipelineSerializedState): any =>
    isFuncCallSerializedState(s) ? s.configId : [s.configId, s.steps.map(shape)];

  async function buildShape(instanceConfig: PipelineInstanceConfigInput) {
    const tree = StateTree.fromInstanceConfig({
      instanceConfig: normalizePipelineInstanceConfig(instanceConfig),
      config: await getProcessedConfig(config),
    });
    return shape(tree.toSerializedState({disableNodesUUID: true}));
  }

  test('An empty instance config starts with the provider defaults', async () => {
    const config: PipelineConfiguration = {
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
    };
    const driver = new Driver(true);
    await driver.sendCommand({
      event: 'initPipeline',
      provider: '',
      config: await getProcessedConfig(config),
      instanceConfig: normalizePipelineInstanceConfig({}),
    });
    const steps = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}>).steps.map((s) => s.configId);
    driver.close();
    expectDeepEqual(steps, ['step1', 'step2']);
  });

  test('An empty instance config builds the whole default tree', async () => {
    expectDeepEqual(await buildShape({}), ['root', [
      'step1',
      ['nestedStatic', ['inner1', ['innerDyn', ['mul']]]],
      ['nestedDyn', ['add', 'mul']],
    ]]);
  });

  test('A nested static workflow without steps gets its config steps', async () => {
    expectDeepEqual(await buildShape({id: 'root', steps: [{id: 'nestedStatic'}]}), ['root', [
      ['nestedStatic', ['inner1', ['innerDyn', ['mul']]]],
    ]]);
  });

  test('A nested dynamic workflow without steps gets its initialSteps', async () => {
    expectDeepEqual(await buildShape({id: 'root', steps: [{id: 'nestedDyn'}]}), ['root', [
      ['nestedDyn', ['add', 'mul']],
    ]]);
  });

  test('Empty steps stay empty for static and dynamic workflows', async () => {
    const instanceConfig = {id: 'root', steps: [{id: 'nestedStatic', steps: []}, {id: 'nestedDyn', steps: []}]};
    expectDeepEqual(await buildShape(instanceConfig), ['root', [
      ['nestedStatic', []],
      ['nestedDyn', []],
    ]]);
  });

  test('Listed steps are built as given for static and dynamic workflows', async () => {
    const instanceConfig = {id: 'root', steps: [
      {id: 'nestedStatic', steps: ['inner1']},
      {id: 'nestedDyn', steps: ['mul', 'mul']},
    ]};
    expectDeepEqual(await buildShape(instanceConfig), ['root', [
      ['nestedStatic', ['inner1']],
      ['nestedDyn', ['mul', 'mul']],
    ]]);
  });

  test('A workflow left out of its parent steps is not built', async () => {
    expectDeepEqual(await buildShape({id: 'root', steps: ['step1']}), ['root', ['step1']]);
  });

  test('A script step entry sets its initial values', async () => {
    const tree = StateTree.fromInstanceConfig({
      instanceConfig: normalizePipelineInstanceConfig({id: 'root', steps: [{id: 'step1', initialValues: {a: 3}}]}),
      config: await getProcessedConfig(config),
    });
    await tree.init().toPromise();
    const fc = ((tree.toState() as PipelineStateStatic<any, {}>).steps[0] as StepFunCallState).funcCall!;
    expectDeepEqual(fc.inputs.a, 3);
  });

  async function runSetPipelineState(state: PipelineInstanceConfigInput) {
    const actionConfig: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [{
        id: 'nestedDyn',
        type: 'parallel',
        stepTypes: [
          {id: 'add', nqName: 'LibTests:TestAdd2'},
          {id: 'block', type: 'static', steps: [
            {id: 'inner1', nqName: 'LibTests:TestAdd2'},
            {id: 'inner2', nqName: 'LibTests:TestMul2'},
          ]},
        ],
        initialSteps: ['add'],
      }],
      actions: [{
        id: 'replace',
        from: [],
        position: 'none',
        to: 'out1:nestedDyn',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', state);
        },
      }],
    };
    const tree = StateTree.fromInstanceConfig({
      instanceConfig: {id: 'root', steps: [{id: 'nestedDyn', steps: []}]},
      config: await getProcessedConfig(actionConfig),
      mockMode: true,
    });
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    return shape(tree.toSerializedState({disableNodesUUID: true}));
  }

  test('setPipelineState without steps gets the initialSteps', async () => {
    expectDeepEqual(await runSetPipelineState({id: 'nestedDyn'}), ['root', [['nestedDyn', ['add']]]]);
  });

  async function runReplaceInPlace(state: PipelineInstanceConfigInput) {
    const replaceConfig: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {
          id: 'nestedDyn',
          type: 'parallel',
          stepTypes: [{id: 'add', nqName: 'LibTests:TestAdd2'}],
          initialSteps: ['add'],
        },
        {id: 'step2', nqName: 'LibTests:TestMul2'},
      ],
      actions: [{
        id: 'replace',
        from: [],
        position: 'none',
        to: 'out1:nestedDyn',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', state);
        },
      }],
    };
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(replaceConfig), mockMode: true});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    return tree;
  }

  test('setPipelineState replaces a workflow in place between its siblings', async () => {
    const tree = await runReplaceInPlace({id: 'nestedDyn', steps: ['add', 'add', 'add']});
    expectDeepEqual(shape(tree.toSerializedState({disableNodesUUID: true})), ['root', [
      'step1',
      ['nestedDyn', ['add', 'add', 'add']],
      'step2',
    ]]);
  });

  test('setPipelineState fills the initial values of the new steps', async () => {
    const tree = await runReplaceInPlace({id: 'nestedDyn', steps: [{id: 'add', initialValues: {a: 5}}]});
    const step = tree.nodeTree.getItem([{idx: 1}, {idx: 0}]) as FuncCallNode;
    expectDeepEqual(step.getStateStore().getState('a'), 5);
  });

  test('setPipelineState works when a step type references an outer workflow', async () => {
    const outerRefConfig: LoadedPipeline = {
      id: 'root',
      nqName: 'mockOuterRefRoot',
      type: 'static',
      steps: [{
        id: 'items',
        type: 'parallel',
        stepTypes: [
          {id: 'add', nqName: 'LibTests:TestAdd2'},
          {type: 'ref', provider: async (_p: any) => outerRefConfig},
        ],
        initialSteps: ['add'],
      }],
      actions: [{
        id: 'replace',
        from: [],
        position: 'none',
        to: 'out1:items',
        type: 'pipeline',
        handler({controller}) {
          controller.setPipelineState('out1', {id: 'items', steps: ['add', 'add']});
        },
      }],
    };
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(outerRefConfig), mockMode: true});
    await tree.init().toPromise();
    const action = [...tree.linksState.actions.values()][0];
    await tree.runAction(action.uuid).toPromise();
    expectDeepEqual(shape(tree.toSerializedState({disableNodesUUID: true})), ['root', [['items', ['add', 'add']]]]);
  });

  test('setPipelineState with a nested static workflow without steps gets its config steps', async () => {
    expectDeepEqual(await runSetPipelineState({id: 'nestedDyn', steps: [{id: 'block'}]}), ['root', [
      ['nestedDyn', [['block', ['inner1', 'inner2']]]],
    ]]);
  });
});
