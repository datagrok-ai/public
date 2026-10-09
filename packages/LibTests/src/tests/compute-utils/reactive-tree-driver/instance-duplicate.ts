import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {FuncCallNode, PipelineNodeBase} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {PipelineSerializedState, PipelineStateStatic, StepFunCallState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {LoadedPipeline} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {NodePath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/BaseTree';

type Shape = string | [string, Shape[]];

function shape(state: PipelineSerializedState): Shape {
  if (!('steps' in state))
    return state.configId;
  return [state.configId, state.steps.map(shape)];
}

function stepAt(tree: StateTree, path: number[]) {
  return tree.nodeTree.getItem(path.map((idx) => ({idx})) as NodePath) as FuncCallNode;
}

function pipelineAt(tree: StateTree, path: number[]) {
  return tree.nodeTree.getItem(path.map((idx) => ({idx})) as NodePath) as PipelineNodeBase;
}

async function makeTree(config: PipelineConfiguration, mockMode = false) {
  const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config), mockMode});
  await tree.init().toPromise();
  return tree;
}

category('ComputeUtils: Driver state tree duplicate', async () => {
  test('Duplicate copies unlinked inputs and inserts after the source', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [{
        id: 'seq',
        type: 'sequential',
        stepTypes: [
          {id: 'add', nqName: 'LibTests:TestAdd2'},
          {id: 'mul', nqName: 'LibTests:TestMul2'},
        ],
        initialSteps: ['add', 'mul'],
      }],
    });
    const source = stepAt(tree, [0, 0]);
    source.getStateStore().editState('a', 5);
    source.getStateStore().editState('b', 7);
    await tree.duplicateSubtree(source.uuid).toPromise();
    expectDeepEqual(shape(tree.toSerializedState()), ['root', [['seq', ['add', 'add', 'mul']]]]);
    const copy = stepAt(tree, [0, 1]);
    expectDeepEqual(copy.uuid !== source.uuid, true, {prefix: 'New node'});
    expectDeepEqual(copy.instancesWrapper.id !== source.instancesWrapper.id, true, {prefix: 'New call'});
    expectDeepEqual(copy.getStateStore().getState('a'), 5, {prefix: 'a'});
    expectDeepEqual(copy.getStateStore().getState('b'), 7, {prefix: 'b'});
    expectDeepEqual(copy.funcCallState$.value?.isOutputOutdated, true, {prefix: 'Outdated'});
  });

  test('Duplicate takes linked inputs from links', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {
          id: 'analyses',
          type: 'dynamic',
          stepTypes: [{id: 'regression', nqName: 'LibTests:TestMul2'}],
          initialSteps: ['regression'],
        },
      ],
      links: [{
        id: 'data-link',
        from: 'in1:step1/b',
        to: 'out1:analyses/all(regression)/a',
      }],
    });
    stepAt(tree, [0]).getStateStore().editState('b', 3);
    await tree.linksState.waitForLinks().toPromise();
    const source = stepAt(tree, [1, 0]);
    source.getStateStore().editState('a', 99);
    source.getStateStore().editState('b', 7);
    await tree.duplicateSubtree(source.uuid).toPromise();
    const copy = stepAt(tree, [1, 1]);
    expectDeepEqual(copy.getStateStore().getState('a'), 3, {prefix: 'Linked a'});
    expectDeepEqual(copy.getStateStore().getState('b'), 7, {prefix: 'Unlinked b'});
    expectDeepEqual(copy.consistencyInfo$.value.a?.inconsistent, false, {prefix: 'Consistent a'});
  });

  test('Duplicate of a run step keeps its run and marks link writes inconsistent', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {
          id: 'analyses',
          type: 'dynamic',
          stepTypes: [{id: 'regression', nqName: 'LibTests:TestMul2'}],
          initialSteps: ['regression'],
        },
      ],
      links: [{
        id: 'data-link',
        from: 'in1:step1/b',
        to: 'out1:analyses/all(regression)/a',
      }],
    }, true);
    stepAt(tree, [0]).getStateStore().editState('b', 3);
    await tree.linksState.waitForLinks().toPromise();
    const source = stepAt(tree, [1, 0]);
    source.getStateStore().editState('a', 99);
    await tree.runStep(source.uuid, {res: 21}).toPromise();
    await tree.duplicateSubtree(source.uuid).toPromise();
    const copy = stepAt(tree, [1, 1]);
    expectDeepEqual(copy.funcCallState$.value?.isOutputOutdated, false, {prefix: 'Run copy'});
    expectDeepEqual(copy.getStateStore().getState('res'), 21, {prefix: 'Output'});
    expectDeepEqual(copy.getStateStore().getState('a'), 99, {prefix: 'Overridden a'});
    expectDeepEqual(copy.consistencyInfo$.value.a?.inconsistent, true, {prefix: 'Inconsistent a'});
    expectDeepEqual(source.getStateStore().getState('a'), 99, {prefix: 'Source a'});
  });

  test('Output values in an instance config mark the step as run', async () => {
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2'},
        {id: 'step2', nqName: 'LibTests:TestAdd2'},
      ],
    };
    const tree = StateTree.fromInstanceConfig({
      config: await getProcessedConfig(config),
      instanceConfig: {id: 'root', steps: [
        {id: 'step1', initialValues: {a: 1, b: 2, res: 3}},
        {id: 'step2', initialValues: {a: 1, b: 2}},
      ]},
    });
    await tree.init().toPromise();
    const ran = stepAt(tree, [0]);
    expectDeepEqual(ran.funcCallState$.value?.isOutputOutdated, false, {prefix: 'Run step'});
    expectDeepEqual(ran.getStateStore().getState('res'), 3, {prefix: 'Output'});
    expectDeepEqual(ran.getStateStore().getState('a'), 1, {prefix: 'Input'});
    expectDeepEqual(stepAt(tree, [1]).funcCallState$.value?.isOutputOutdated, true, {prefix: 'Not run step'});
  });

  test('Duplicate copies dataframe inputs', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [{
        id: 'seq',
        type: 'sequential',
        stepTypes: [{id: 'df', nqName: 'LibTests:TestDF1'}],
        initialSteps: ['df'],
      }],
    });
    const source = stepAt(tree, [0, 0]);
    const df = DG.DataFrame.fromCsv('x\n1\n2');
    source.getStateStore().editState('df', df);
    await tree.duplicateSubtree(source.uuid).toPromise();
    const copied = stepAt(tree, [0, 1]).getStateStore().getState<DG.DataFrame>('df')!;
    expectDeepEqual(copied !== df, true, {prefix: 'Separate dataframe'});
    expectDeepEqual(copied.rowCount, 2, {prefix: 'Rows'});
  });

  test('Duplicate copies a workflow item with nested items', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [{
        id: 'items',
        type: 'parallel',
        stepTypes: [{
          id: 'block',
          type: 'static',
          steps: [
            {id: 'inner1', nqName: 'LibTests:TestAdd2'},
            {
              id: 'innerDyn',
              type: 'dynamic',
              stepTypes: [{id: 'mul', nqName: 'LibTests:TestMul2'}],
              initialSteps: ['mul'],
            },
          ],
          states: ['meta1'],
        }],
        initialSteps: ['block'],
      }],
    });
    const block = pipelineAt(tree, [0, 0]);
    const innerDyn = pipelineAt(tree, [0, 0, 1]);
    await tree.addSubTree(innerDyn.uuid, 'mul', 1).toPromise();
    block.getStateStore().setState('meta1', 'm');
    stepAt(tree, [0, 0, 0]).getStateStore().editState('a', 1);
    stepAt(tree, [0, 0, 1, 1]).getStateStore().editState('b', 2);
    await tree.duplicateSubtree(block.uuid).toPromise();
    expectDeepEqual(shape(tree.toSerializedState()), ['root', [['items', [
      ['block', ['inner1', ['innerDyn', ['mul', 'mul']]]],
      ['block', ['inner1', ['innerDyn', ['mul', 'mul']]]],
    ]]]]);
    expectDeepEqual(pipelineAt(tree, [0, 1]).getStateStore().getState('meta1'), 'm', {prefix: 'State'});
    expectDeepEqual(stepAt(tree, [0, 1, 0]).getStateStore().getState('a'), 1, {prefix: 'inner1 a'});
    expectDeepEqual(stepAt(tree, [0, 1, 1, 1]).getStateStore().getState('b'), 2, {prefix: 'mul b'});
  });

  test('Duplicate does not run onInit', async () => {
    let initRuns = 0;
    const tree = await makeTree({
      id: 'root',
      type: 'parallel',
      stepTypes: [{
        id: 'block',
        type: 'static',
        steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
        states: ['meta1'],
        onInit: {
          id: 'init',
          from: 'in1:step1/b',
          to: 'out1:meta1',
          handler({controller}) {
            initRuns++;
            controller.setAll('out1', initRuns);
          },
        },
      }],
    });
    await tree.addSubTree(tree.nodeTree.root.getItem().uuid, 'block', 0).toPromise();
    expectDeepEqual(initRuns, 1, {prefix: 'Init on add'});
    const block = pipelineAt(tree, [0]);
    block.getStateStore().setState('meta1', 42);
    await tree.duplicateSubtree(block.uuid).toPromise();
    expectDeepEqual(initRuns, 1, {prefix: 'Init on duplicate'});
    expectDeepEqual(pipelineAt(tree, [1]).getStateStore().getState('meta1'), 42, {prefix: 'Copied state'});
  });

  test('Duplicate copies a self-referencing workflow item', async () => {
    const config: LoadedPipeline = {
      id: 'pipelineSeq',
      nqName: 'mockNqName',
      type: 'sequential',
      stepTypes: [
        {id: 'stepMul', nqName: 'LibTests:TestMul2'},
        {type: 'ref', provider: async (_p: any) => config},
      ],
    };
    const tree = StateTree.fromInstanceConfig({
      config: await getProcessedConfig(config),
      instanceConfig: {id: 'pipelineSeq', steps: [
        {id: 'pipelineSeq', steps: [{id: 'stepMul', initialValues: {a: 4}}]},
      ]},
    });
    await tree.init().toPromise();
    await tree.duplicateSubtree(pipelineAt(tree, [0]).uuid).toPromise();
    expectDeepEqual(shape(tree.toSerializedState()), ['pipelineSeq', [
      ['pipelineSeq', ['stepMul']],
      ['pipelineSeq', ['stepMul']],
    ]]);
    expectDeepEqual(stepAt(tree, [1, 0]).getStateStore().getState('a'), 4, {prefix: 'a'});
  });

  test('Duplicate copies an item that references an outer workflow', async () => {
    const config: LoadedPipeline = {
      id: 'root',
      nqName: 'mockDuplicateOuterRef',
      type: 'static',
      steps: [{
        id: 'items',
        type: 'parallel',
        stepTypes: [{
          id: 'block',
          type: 'static',
          steps: [
            {id: 'inner1', nqName: 'LibTests:TestAdd2'},
            {
              id: 'nested',
              type: 'parallel',
              stepTypes: [{type: 'ref', provider: async (_p: any) => config}],
            },
          ],
        }],
        initialSteps: ['block'],
      }],
    };
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config)});
    await tree.init().toPromise();
    await tree.duplicateSubtree(pipelineAt(tree, [0, 0]).uuid).toPromise();
    expectDeepEqual(shape(tree.toSerializedState()), ['root', [['items', [
      ['block', ['inner1', ['nested', []]]],
      ['block', ['inner1', ['nested', []]]],
    ]]]]);
  });

  test('Instance config sets workflow states and skips onInit', async () => {
    let initRuns = 0;
    const config: PipelineConfiguration = {
      id: 'root',
      type: 'static',
      steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
      states: ['meta1'],
      onInit: {
        id: 'init',
        from: 'in1:step1/b',
        to: 'out1:meta1',
        handler({controller}) {
          initRuns++;
          controller.setAll('out1', 'init');
        },
      },
    };
    const tree = StateTree.fromInstanceConfig({
      config: await getProcessedConfig(config),
      instanceConfig: {id: 'root', skipOnInit: true, initialValues: {meta1: 'given'}},
    });
    await tree.init().toPromise();
    expectDeepEqual(initRuns, 0, {prefix: 'Init runs'});
    expectDeepEqual(pipelineAt(tree, []).getStateStore().getState('meta1'), 'given', {prefix: 'State'});
  });

  test('Duplicate of a static workflow step fails', async () => {
    const driver = new Driver(true);
    await driver.sendCommand({
      event: 'initPipeline',
      provider: '',
      config: await getProcessedConfig({
        id: 'root',
        type: 'static',
        steps: [{id: 'step1', nqName: 'LibTests:TestAdd2'}],
      }),
    });
    const step = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}>).steps[0];
    const res = await driver.sendCommand({event: 'duplicateDynamicItem', uuid: step.uuid});
    expectDeepEqual(res, null, {prefix: 'Resolved value'});
    expectDeepEqual(driver.logger.errors[0]?.context, 'command:duplicateDynamicItem', {prefix: 'Error context'});
    driver.close();
  });

  test('Duplicate command adds a dynamic item', async () => {
    const driver = new Driver(true);
    await driver.sendCommand({
      event: 'initPipeline',
      provider: '',
      config: await getProcessedConfig({
        id: 'root',
        type: 'sequential',
        stepTypes: [{id: 'add', nqName: 'LibTests:TestAdd2'}],
        initialSteps: ['add'],
      }),
    });
    const step = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}>).steps[0];
    await driver.sendCommand({event: 'duplicateDynamicItem', uuid: step.uuid});
    const steps = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}>).steps;
    expectDeepEqual(steps.map((s) => s.configId), ['add', 'add']);
    driver.close();
  });
});
