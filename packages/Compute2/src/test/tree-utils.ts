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
  stepMoveIndex, dropIndex, focusCandidates,
} from '../components/TreeWizard/tree-utils';
import {statesToStatus} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/run-summary';
import {AugmentedStat} from '../components/TreeWizard/types';
import {mockDynamicPipeline, mockFuncCall, mockStaticPipeline} from './state-fixtures';

const call = (state: Partial<FuncCallStateInfo> = {}): FuncCallStateInfo => ({
  isRunning: false, isRunnable: true, isOutputOutdated: true, runError: undefined, pendingDependencies: [], ...state,
});
const errors = {a: {errors: ['bad']}};
const warnings = {b: {warnings: ['careful']}};
const inconsistent = (restriction: RestrictionType): Record<string, ConsistencyInfo> =>
  ({c: {restriction, inconsistent: true, assignedValue: 1}});
const action = (id: string, position: ViewAction['position']): ViewAction => ({id, uuid: id, position});

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

  test('A read-only item of an editable workflow can be moved and removed but not duplicated', async () => {
    const parent = mockDynamicPipeline('dyn', [], {stepTypes: [{configId: 's'}]});
    expectDeepEqual(permissions(stat(mockFuncCall('s', {isReadonly: true}), parent)), [true, true, false]);
  });

  test('Each disableUI flag turns off its own action for that step type', async () => {
    expectDeepEqual(permissions(inDynamic({configId: 's', disableUIDragging: true})), [false, true, true]);
    expectDeepEqual(permissions(inDynamic({configId: 's', disableUIRemoving: true})), [true, false, true]);
    expectDeepEqual(permissions(inDynamic({configId: 's', disableUIAdding: true})), [true, true, false]);
    expectDeepEqual(permissions(inDynamic({configId: 'other', disableUIRemoving: true})), [true, true, true],
      {prefix: 'Other type'});
  });
});

category('TreeWizard: keyboard moves and focus', () => {
  // a parent stat with its children linked both ways, like he-tree builds them
  const tree = (parent: PipelineState, children: PipelineState[]) => {
    const root = {data: parent, parent: null, children: []} as unknown as AugmentedStat;
    root.children = children.map((data) => ({data, parent: root, children: []}) as unknown as AugmentedStat);
    return root;
  };
  const PINNED: StepDynamicDescription = {configId: 'pinned', disableUIDragging: true};
  const list = (stepTypes: StepDynamicDescription[] = [{configId: 's'}, PINNED]) =>
    tree(mockDynamicPipeline('dyn', [], {stepTypes}), [mockFuncCall('s'), mockFuncCall('s'), mockFuncCall('pinned')]);

  test('A step moves one place within its list and never out of it', async () => {
    const [first, second] = list().children;
    expectDeepEqual([stepMoveIndex(first, -1), stepMoveIndex(first, 1)], [undefined, 1], {prefix: 'First'});
    expectDeepEqual([stepMoveIndex(second, -1), stepMoveIndex(second, 1)], [0, 2], {prefix: 'Past a pinned step'});
  });

  test('A step that cannot be dragged does not move by keys either', async () => {
    const pinned = list().children[2];
    expectDeepEqual(stepMoveIndex(pinned, -1), undefined, {prefix: 'Pinned'});
    const fixed = tree(mockStaticPipeline('p', []), [mockFuncCall('a'), mockFuncCall('b')]).children[0];
    expectDeepEqual(stepMoveIndex(fixed, 1), undefined, {prefix: 'Static parent'});
    const readonly = tree(mockDynamicPipeline('dyn', [], {stepTypes: [{configId: 'a'}], isReadonly: true}),
      [mockFuncCall('a'), mockFuncCall('a')]).children[0];
    expectDeepEqual(stepMoveIndex(readonly, 1), undefined, {prefix: 'Read-only parent'});
  });

  test('A drop counts its position without the dragged step', async () => {
    expectDeepEqual([dropIndex(0, 2), dropIndex(2, 0), dropIndex(1, 1)], [1, 0, 1]);
  });

  test('The focus falls back to the next sibling, the previous one, then the parent', async () => {
    const root = tree(mockDynamicPipeline('dyn', []), [mockFuncCall('a'), mockFuncCall('b'), mockFuncCall('c')]);
    expectDeepEqual(focusCandidates(root.children[1]), ['b', 'c', 'a', 'dyn'], {prefix: 'Middle'});
    expectDeepEqual(focusCandidates(root.children[2]), ['c', 'b', 'dyn'], {prefix: 'Last'});
    expectDeepEqual(focusCandidates(root), ['dyn'], {prefix: 'Root'});
  });
});
