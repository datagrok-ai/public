import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration, historyUtils} from '@datagrok-libraries/compute-utils';
import {
  getProcessedConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {
  CONFIG_PATH, OUTPUT_OUTDATED_PATH, loadFuncCall, readOutputOutdated,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {
  StateSnapshot, StepSnapshot, loadStateSnapshot, saveStateSnapshot, snapshotCall, snapshotState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import type {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {deleteRuns, savedCallIds} from '../../../test-utils';

// step1 ran, step2 got its a from the link and never ran
async function liveTree() {
  const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
  const pconf = await getProcessedConfig(config);
  const tree = StateTree.fromPipelineConfig({config: pconf});
  await tree.init().toPromise();
  const step1 = tree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
  step1.getStateStore().editState('a', 1);
  step1.getStateStore().editState('b', 2);
  await tree.runStep(step1.uuid).toPromise();
  await tree.linksState.waitForLinks().toPromise();
  return {tree, pconf};
}

// step1 failed with a run error, step2 has the restricted a its link set, overridden
async function errorTree() {
  const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper4', {version: '1.0'}).toPromise();
  const pconf = await getProcessedConfig(config);
  const tree = StateTree.fromPipelineConfig({config: pconf});
  await tree.init().toPromise();
  const step1 = tree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
  step1.getStateStore().editState('a', -1);
  step1.getStateStore().editState('b', -2);
  await step1.getStateStore().run().toPromise();
  await tree.linksState.waitForLinks().toPromise();
  (tree.nodeTree.getItem([{idx: 1}]) as FuncCallNode).getStateStore().editState('a', 5);
  return {tree, pconf};
}

const values = <T, >(states: Record<string, {value: T}>) =>
  Object.fromEntries(Object.entries(states).map(([uuid, state$]) => [uuid, state$.value]));

const liveStates = (tree: StateTree) =>
  ({calls: values(tree.getFuncCallStates()), consistency: values(tree.getConsistency())});

const snapshotOf = (tree: StateTree) => {
  const {calls, consistency} = liveStates(tree);
  return snapshotState(tree.toState(), calls, consistency);
};

// saving gives each copied call its run id, so this also finds the runs of a save that failed partway
const callIds = (snapshot: StateSnapshot): string[] => 'steps' in snapshot ?
  snapshot.steps.flatMap((step) => callIds(step as StateSnapshot)) :
  [(snapshot as StepSnapshot).funcCall?.id].filter((id): id is string => !!id);

const stepStates = (tree: StateTree) => (tree.toState() as {steps: PipelineState[]}).steps;

const steps = (snapshot: StateSnapshot) => (snapshot as {steps: StepSnapshot[]}).steps;

const withoutCalls = (snapshot: StateSnapshot) => ({
  ...snapshot,
  steps: steps(snapshot).map(({funcCall: _c, isOutputOutdated: _o, runError: _e, restrictions: _r, ...step}) => step),
});

const stepValues = (tree: StateTree) => [0, 1].map((idx) => {
  const step = tree.nodeTree.getItem([{idx}]) as FuncCallNode;
  const store = step.getStateStore();
  return {a: store.getState('a'), b: store.getState('b'), res: store.getState('res'),
    restriction: step.instancesWrapper.inputRestrictions$.value.a,
    outdated: step.funcCallState$.value?.isOutputOutdated, runError: step.funcCallState$.value?.runError};
});

category('ComputeUtils: Driver state snapshots', async () => {
  test('A snapshot of live state has the serialized structure and copies of the calls', async () => {
    const {tree} = await liveTree();
    const snapshot = snapshotOf(tree);
    expectDeepEqual(withoutCalls(snapshot), tree.toSerializedState(), {prefix: 'Structure'});
    const [step1, step2] = steps(snapshot);
    const liveCall = (tree.nodeTree.getItem([{idx: 0}]) as FuncCallNode).instancesWrapper.getInstance()!.getFuncCall();
    expectDeepEqual(step1.funcCall !== liveCall, true, {prefix: 'Copied call'});
    liveCall.inputs['a'] = 100;
    expectDeepEqual(step1.funcCall!.inputs['a'], 1, {prefix: 'Kept value'});
    expectDeepEqual([step1.isOutputOutdated, step2.isOutputOutdated], [false, true], {prefix: 'Statuses'});
    expectDeepEqual(step2.restrictions, {a: {type: 'restricted', assignedValue: 3}}, {prefix: 'Restrictions'});
    tree.close();
  });

  test('Saving a snapshot gives the run that saving the tree gives', async () => {
    const {tree, pconf} = await liveTree();
    const snapshot = snapshotOf(tree);
    let metaCall: DG.FuncCall | undefined;
    let treeCall: DG.FuncCall | undefined;
    try {
      metaCall = await saveStateSnapshot(snapshot, {title: 'Snapshot'});
      treeCall = await tree.save(undefined, {title: 'Snapshot'}).toPromise();
      const [fromSnapshot, fromTree] = await Promise.all([metaCall, treeCall!].map(async (call) => {
        const loaded = await StateTree.load({dbId: call.id, config: pconf}).toPromise();
        return loaded.init().toPromise();
      }));
      const options = {disableNodesUUID: true, disableCallsUUID: true};
      expectDeepEqual(fromSnapshot!.toSerializedState(options), fromTree!.toSerializedState(options),
        {prefix: 'Structure'});
      expectDeepEqual(stepValues(fromSnapshot!), stepValues(fromTree!), {prefix: 'Values'});
      expectDeepEqual(metaCall.options['title'], treeCall!.options['title'], {prefix: 'Title'});
    } finally {
      await deleteRuns([metaCall?.id, ...callIds(snapshot), treeCall?.id, ...savedCallIds(tree.toSerializedState())]);
      tree.close();
    }
  });

  test('A saved run loads with its step runs, a step that fails to load keeps the error', async () => {
    const {tree} = await liveTree();
    const metaCall = (await tree.save().toPromise())!;
    try {
      const loaded = (await loadStateSnapshot(await historyUtils.loadRun(metaCall.id)))!;
      const [step1, step2] = steps(loaded);
      expectDeepEqual([step1.funcCall?.inputs['a'], step1.isOutputOutdated, step2.isOutputOutdated], [1, false, true]);
      const broken = DG.Script.create(
        '//name: snapshotProbe\n//language: javascript\n//output: int x\nx = 1;').prepare();
      broken.options[CONFIG_PATH] = serialize({...withoutCalls(loaded), steps: [
        steps(withoutCalls(loaded) as StateSnapshot)[0], {...steps(withoutCalls(loaded) as StateSnapshot)[1],
          funcCallId: 'snapshot-run-that-does-not-exist'},
      ]}, {useJsonDF: false});
      const [ok, failed] = steps((await loadStateSnapshot(broken))!);
      expectDeepEqual([!!ok.funcCall, !!ok.loadError, !!failed.funcCall, !!failed.loadError, failed.isOutputOutdated],
        [true, false, false, true, true]);
      expectDeepEqual(await loadStateSnapshot(ok.funcCall!), undefined, {prefix: 'Plain run'});
    } finally {
      await deleteRuns([metaCall.id, ...savedCallIds(tree.toSerializedState())]);
      tree.close();
    }
  });

  test('A step snapshot saves its call as a new run with the metadata', async () => {
    const call = DG.Func.byName('LibTests:TestAdd2').prepare({a: 1, b: 2});
    const snapshot = snapshotCall(call, {isOutputOutdated: true});
    call.inputs['a'] = 100;
    try {
      const saved = await saveStateSnapshot(snapshot, {title: 'Step'});
      expectDeepEqual(saved.id !== call.id, true, {prefix: 'New run'});
      const loaded = await historyUtils.loadRun(saved.id);
      expectDeepEqual([loaded.options['title'], loaded.inputs['a'], readOutputOutdated(loaded)], ['Step', 1, true]);
      expectDeepEqual(loaded.options[OUTPUT_OUTDATED_PATH] != null, true, {prefix: 'Outdated flag saved'});
    } finally {
      await deleteRuns(callIds(snapshot));
    }
  });

  test('A workflow snapshot saves and restores run errors and overridden restrictions', async () => {
    const {tree, pconf} = await errorTree();
    const snapshot = snapshotOf(tree);
    let metaCall: DG.FuncCall | undefined;
    try {
      metaCall = await saveStateSnapshot(snapshot);
      const loaded = await StateTree.load({dbId: metaCall.id, config: pconf}).toPromise();
      await loaded.init().toPromise();
      expectDeepEqual(stepValues(loaded), stepValues(tree), {prefix: 'Restored'});
      expectDeepEqual(stepValues(loaded).map(({runError, restriction}) => ({runError, restriction})), [
        {runError: 'Error: Test error', restriction: undefined},
        {runError: undefined, restriction: {type: 'restricted', assignedValue: 2}},
      ], {prefix: 'Run error and restriction'});
      const step2 = loaded.nodeTree.getItem([{idx: 1}]) as FuncCallNode;
      expectDeepEqual(step2.consistencyInfo$.value.a?.inconsistent, true, {prefix: 'Overridden a'});
    } finally {
      await deleteRuns([metaCall?.id, ...callIds(snapshot)]);
      tree.close();
    }
  });

  test('A workflow step snapshot saves as a run that restores its run error, status and restrictions', async () => {
    const {tree} = await errorTree();
    const {calls, consistency} = liveStates(tree);
    const [failed, overridden] = stepStates(tree).map((state) => snapshotState(state, calls, consistency));
    try {
      for (const step of [failed, overridden])
        await saveStateSnapshot(step, {title: 'Step'});
      const [failedRun, overriddenRun] =
        await Promise.all([failed, overridden].map((step) => loadFuncCall(callIds(step)[0], false)));
      expectDeepEqual([failedRun.runError, failedRun.isOutputOutdated, failedRun.adapter.getState('a')],
        ['Error: Test error', true, -1], {prefix: 'Failed step'});
      expectDeepEqual([overriddenRun.restrictions, overriddenRun.isOutputOutdated, overriddenRun.adapter.getState('a')],
        [{a: {type: 'restricted', assignedValue: 2}}, true, 5], {prefix: 'Overridden step'});
      expectDeepEqual(overriddenRun.adapter.getFuncCall().options['title'], 'Step', {prefix: 'Title'});
    } finally {
      await deleteRuns([...callIds(failed), ...callIds(overridden)]);
      tree.close();
    }
  });

  test('A nested item snapshot saves as a run of its own and loads back into its slot', async () => {
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockWrapper3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const tree = StateTree.fromPipelineConfig({config: outerPconf});
    await tree.init().toPromise();
    const add = tree.nodeTree.getItem([{idx: 0}, {idx: 0}]) as FuncCallNode;
    add.getStateStore().editState('a', 1);
    add.getStateStore().editState('b', 2);
    await tree.runStep(add.uuid).toPromise();
    const {calls, consistency} = liveStates(tree);
    const snapshot = snapshotState(stepStates(tree)[0], calls, consistency);
    const restored = StateTree.fromPipelineConfig({config: outerPconf});
    let metaCall: DG.FuncCall | undefined;
    try {
      metaCall = await saveStateSnapshot(snapshot);
      expectDeepEqual(metaCall.func.nqName, 'LibTests:MockWrapper2', {prefix: 'Own run'});
      await restored.init().toPromise();
      const root = restored.nodeTree.getItem([]);
      await restored.loadSubTree(root.uuid, metaCall.id, 'pipelinePar', 1, false).toPromise();
      const options = {disableNodesUUID: true, disableCallsUUID: true};
      expectDeepEqual(StateTree.toSerializedStateRec(restored.nodeTree.getNode([{idx: 1}]), options),
        {...StateTree.toSerializedStateRec(tree.nodeTree.getNode([{idx: 0}]), options), sourceRunId: metaCall.id},
        {prefix: 'Structure'});
      const loadedAdd = restored.nodeTree.getItem([{idx: 1}, {idx: 0}]) as FuncCallNode;
      expectDeepEqual([loadedAdd.getStateStore().getState('a'), loadedAdd.getStateStore().getState('res'),
        loadedAdd.funcCallState$.value?.isOutputOutdated], [1, 3, false], {prefix: 'Values'});
    } finally {
      await deleteRuns([metaCall?.id, ...callIds(snapshot)]);
      tree.close();
      restored.close();
    }
  });

  test('A dataframe restriction is copied into the snapshot and saved with it', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapperDF', {version: '1.0'}).toPromise();
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config)});
    await tree.init().toPromise();
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('x', new Int32Array([1, 2, 3]))]);
    const step = tree.nodeTree.getItem([{idx: 1}]) as FuncCallNode;
    step.getStateStore().setState('df', df, 'restricted');
    const {calls, consistency} = liveStates(tree);
    const snapshot = snapshotState(stepStates(tree)[1], calls, consistency) as StepSnapshot;
    try {
      const copied = snapshot.restrictions!.df!.assignedValue as DG.DataFrame;
      const live = step.instancesWrapper.inputRestrictions$.value.df!.assignedValue as DG.DataFrame;
      expectDeepEqual([copied instanceof DG.DataFrame, copied !== live, copied.rowCount], [true, true, 3],
        {prefix: 'Copied'});
      live.set('x', 0, 100);
      expectDeepEqual(copied.get('x', 0), 1, {prefix: 'Independent'});
      await saveStateSnapshot(snapshot);
      const {restrictions} = await loadFuncCall(callIds(snapshot)[0], false);
      const saved = restrictions.df!.assignedValue as DG.DataFrame;
      expectDeepEqual([restrictions.df!.type, saved.rowCount, saved.get('x', 0)], ['restricted', 3, 1],
        {prefix: 'Saved'});
    } finally {
      await deleteRuns(callIds(snapshot));
      tree.close();
    }
  });
});
