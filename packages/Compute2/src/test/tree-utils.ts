import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {
  PipelineState, StepDynamicDescription, ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {
  ConsistencyInfo, FuncCallStateInfo,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {RestrictionType} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {
  couldBeSaved, findTreeNodeByPath, findTreeNodeParrent, getRelevantGlobalActions, hasAddControls,
  hasSubtreeAnyInconsistencies, hasSubtreeFixableInconsistencies, isDeletable, isDuplicable, isEachDraggable,
  statesToStatus,
} from '../utils';
import {AugmentedStat} from '../components/TreeWizard/types';
import {mockDynamicPipeline, mockFuncCall, mockStaticPipeline} from './state-fixtures';

const call = (state: Partial<FuncCallStateInfo> = {}): FuncCallStateInfo => ({
  isRunning: false, isRunnable: true, isOutputOutdated: true, runError: undefined, pendingDependencies: [], ...state,
});
const errors = {a: {errors: ['bad']}};
const warnings = {b: {warnings: ['careful']}};
const inconsistent = (restriction: RestrictionType): Record<string, ConsistencyInfo> =>
  ({c: {restriction, inconsistent: true, assignedValue: 1}});
const action = (id: string, position: ViewAction['position']): ViewAction => ({id, uuid: id, position, visible: true});

category('TreeWizard: step status', () => {
  test('Running, failed and pending come first', async () => {
    expect(statesToStatus(call({isRunning: true, runError: 'x'})), 'running');
    expect(statesToStatus(call({runError: 'x', pendingDependencies: ['s']})), 'failed');
    expect(statesToStatus(call({pendingDependencies: ['s']}), errors), 'pending');
    expect(statesToStatus(call({pendingDependencies: ['s'], isOutputOutdated: false})), 'pending executed');
  });

  test('A run step reports inconsistencies, then warnings, then changes', async () => {
    const done = call({isOutputOutdated: false});
    expect(statesToStatus(done, warnings, inconsistent('restricted')), 'succeeded inconsistent');
    expect(statesToStatus(done, {}, inconsistent('disabled')), 'succeeded inconsistent');
    expect(statesToStatus(done, errors, inconsistent('info')), 'succeeded warn');
    expect(statesToStatus(done, warnings), 'succeeded warn');
    expect(statesToStatus(done, {}, inconsistent('info')), 'succeeded info');
    expect(statesToStatus(done), 'succeeded');
  });

  test('A step to run reports errors before warnings', async () => {
    expect(statesToStatus(call(), {...errors, ...warnings}), 'next error');
    expect(statesToStatus(call(), warnings), 'next warn');
    expect(statesToStatus(call(), {}, inconsistent('restricted')), 'next warn');
    expect(statesToStatus(call(), {}, inconsistent('info')), 'next');
    expect(statesToStatus(call()), 'next');
  });
});

category('TreeWizard: tree helpers', () => {
  const buildTree = () => mockStaticPipeline('root', [
    mockFuncCall('s1'),
    mockDynamicPipeline('dyn', [mockFuncCall('s2')], {
      stepTypes: [{configId: 's2'}],
      actions: [action('dynGlobal', 'globalmenu'), action('dynButton', 'buttons')],
    }),
  ], {nqName: 'Pkg:Root', actions: [action('rootGlobal', 'globalmenu')]});

  test('Parents and paths are found', async () => {
    const tree = buildTree();
    expect(findTreeNodeParrent('s2', tree)?.uuid, 'dyn');
    expect(findTreeNodeParrent('s1', tree)?.uuid, 'root');
    expect(findTreeNodeParrent('root', tree) === undefined, true);
    expect(findTreeNodeByPath([0, 1, 0], tree)?.state.uuid, 's2');
    expect(findTreeNodeByPath([0], tree)?.state.uuid, 'root');
    expect(findTreeNodeByPath([0, 5], tree) === undefined, true);
  });

  test('Global actions are collected from the root down to the step', async () => {
    const tree = buildTree();
    expectDeepEqual(getRelevantGlobalActions(tree, 's2').map((a) => a.id), ['rootGlobal', 'dynGlobal']);
    expectDeepEqual(getRelevantGlobalActions(tree, 's1').map((a) => a.id), ['rootGlobal']);
  });

  test('Add and save controls follow the workflow state', async () => {
    const stepTypes: StepDynamicDescription[] = [{configId: 's2'}];
    expect(hasAddControls(mockDynamicPipeline('d', [], {stepTypes})), true);
    expect(hasAddControls(mockDynamicPipeline('d', [], {stepTypes, isReadonly: true})), false);
    expect(hasAddControls(mockDynamicPipeline('d', [], {stepTypes: [{configId: 's2', disableUIAdding: true}]})), false);
    expect(hasAddControls(mockStaticPipeline('s', [])), false);
    expect(couldBeSaved(mockStaticPipeline('p', [], {nqName: 'Pkg:P'})), true);
    expect(couldBeSaved(mockStaticPipeline('p', [], {nqName: 'Pkg:P', disableHistory: true})), false);
    expect(couldBeSaved(mockStaticPipeline('p', [])), false);
    expect(couldBeSaved(mockFuncCall('s')), false);
  });
});

category('TreeWizard: subtree inconsistencies', () => {
  const tree = mockStaticPipeline('root', [mockFuncCall('info'), mockFuncCall('locked', {isReadonly: true}),
    mockFuncCall('waiting'), mockFuncCall('ok')]);
  const consistency = {
    info: inconsistent('info'), locked: inconsistent('restricted'), waiting: inconsistent('restricted'),
  };
  const calls = {waiting: call({pendingDependencies: ['ok']})};

  test('Fixable inconsistencies skip info, read-only and waiting steps', async () => {
    expect(hasSubtreeFixableInconsistencies(tree, calls, consistency) === undefined, true);
    expect(hasSubtreeFixableInconsistencies(tree, {}, consistency)?.state.uuid, 'waiting');
  });

  test('Any inconsistency counts info but still skips read-only and waiting steps', async () => {
    expect(hasSubtreeAnyInconsistencies(tree, calls, consistency)?.state.uuid, 'info');
    expect(hasSubtreeAnyInconsistencies(tree, calls, {locked: inconsistent('info')}) === undefined, true);
  });
});

category('TreeWizard: step permissions', () => {
  const stat = (data: PipelineState, parent?: PipelineState) => ({
    data, children: [], parent: parent ? {data: parent, parent: null, children: []} : null,
  }) as unknown as AugmentedStat;
  const permissions = (s: AugmentedStat) => [!!isEachDraggable(s), isDeletable(s), isDuplicable(s)];
  const inDynamic = (stepType: StepDynamicDescription, opts: {isReadonly?: boolean} = {}) =>
    stat(mockFuncCall('s'), mockDynamicPipeline('dyn', [], {stepTypes: [stepType], ...opts}));

  test('Only items of an editable dynamic workflow can be moved, removed or duplicated', async () => {
    expectDeepEqual(permissions(stat(mockFuncCall('s'))), [false, false, false], {prefix: 'Root'});
    expectDeepEqual(permissions(stat(mockFuncCall('s'), mockStaticPipeline('p', []))), [false, false, false],
      {prefix: 'Static parent'});
    expectDeepEqual(permissions(inDynamic({configId: 's'}, {isReadonly: true})), [false, false, false],
      {prefix: 'Read-only parent'});
    expectDeepEqual(permissions(inDynamic({configId: 's'})), [true, true, true], {prefix: 'Editable'});
  });

  test('Each disableUI flag turns off its own action for that step type', async () => {
    expectDeepEqual(permissions(inDynamic({configId: 's', disableUIDragging: true})), [false, true, true]);
    expectDeepEqual(permissions(inDynamic({configId: 's', disableUIRemoving: true})), [true, false, true]);
    expectDeepEqual(permissions(inDynamic({configId: 's', disableUIAdding: true})), [true, true, false]);
    expectDeepEqual(permissions(inDynamic({configId: 'other', disableUIRemoving: true})), [true, true, true],
      {prefix: 'Other type'});
  });
});
