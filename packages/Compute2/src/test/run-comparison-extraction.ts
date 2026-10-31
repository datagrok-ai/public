import * as DG from 'datagrok-api/dg';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import type {PipelineConfiguration} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {
  getProcessedConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import type {
  PipelineSerializedState, PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  CONFIG_PATH, OUTPUT_OUTDATED_PATH, readOutputOutdated,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {snapshotState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {
  parseComparisonDefaults, entryFromDataFrame, collectWorkflowSteps, entryFromFuncCall, entryFromSnapshot,
  entryFromLiveCall, loadRunEntries,
} from '../components/RunComparison/entry-extraction';
import {deleteSavedRun, deleteSavedWorkflow} from './utils';
import {STEP_HISTORY_OPTION} from '../components/RunStorage/run-storage';

category('RunComparison: annotation defaults', () => {
  test('parses the comparison JSON annotation', async () => {
    const defaults = parseComparisonDefaults({
      comparison: '{"index": "time", "split": "species", "mode": "timeseries", "units": "s"}',
    });
    expect(defaults.index, 'time');
    expect(defaults.split, 'species');
    expect(defaults.mode, 'timeseries');
    expect(defaults.units, 's');
  });

  test('ignores unknown mode and units values', async () => {
    const defaults = parseComparisonDefaults({comparison: '{"mode": "spiral", "units": "weeks"}'});
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });

  test('malformed JSON produces no defaults', async () => {
    const defaults = parseComparisonDefaults({comparison: '{index: time'});
    expect(defaults.index == null, true);
    expect(defaults.split == null, true);
  });

  test('null annotation produces no defaults', async () => {
    // JSON.parse('null') succeeds with null — must not crash on property access
    const defaults = parseComparisonDefaults({comparison: 'null'});
    expect(defaults.index == null, true);
    expect(defaults.mode == null, true);
  });

  test('missing annotations produce no defaults', async () => {
    const defaults = parseComparisonDefaults(undefined);
    expect(defaults.index == null, true);
    expect(defaults.split == null, true);
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });
});

category('RunComparison: raw table entries', () => {
  const makeDf = () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'time', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20, 30]),
    ]);
    df.name = 'data';
    return df;
  };

  test('distinct tables with equal names and row counts get distinct ids', async () => {
    expect(entryFromDataFrame(makeDf()).id === entryFromDataFrame(makeDf()).id, false);
  });

  test('the same table maps to the same id', async () => {
    const df = makeDf();
    expect(entryFromDataFrame(df).id, entryFromDataFrame(df).id);
  });
});

category('RunComparison: workflow steps', () => {
  const step = (configId: string, friendlyName: string, funcCallId: string) =>
    ({type: 'funccall', configId, friendlyName, funcCallId});
  // only the fields the step walk reads
  const state = {
    type: 'static', configId: 'wf', friendlyName: 'Workflow', steps: [
      step('config', 'Configuration', 'c'),
      {type: 'dynamic', configId: 'days', friendlyName: 'Daily cultivation', steps: [
        step('day', 'Day', 'd1'), step('day', 'Day', 'd2'), step('other', 'Other', 'o'),
      ]},
    ],
  } as unknown as PipelineSerializedState;

  test('repeated steps get positional paths and labels', async () => {
    const steps = collectWorkflowSteps(state);
    expectArray(steps.map((s) => s.path), ['wf/config', 'wf/days/day[1]', 'wf/days/day[2]', 'wf/days/other']);
    expectArray(steps.map((s) => s.friendlyPath),
      ['Configuration', 'Daily cultivation · Day 1', 'Daily cultivation · Day 2', 'Daily cultivation · Other']);
    expectArray(steps.map((s) => s.step.funcCallId), ['c', 'd1', 'd2', 'o']);
  });
});

category('RunComparison: workflow step failures', () => {
  test('steps that cannot be loaded are reported, not thrown', async () => {
    const call = DG.Script.create(
      '//name: rcStepFailureProbe\n//language: javascript\n//output: int x\nx = 1;').prepare();
    const missing = (configId: string, friendlyName: string) =>
      ({type: 'funccall', configId, friendlyName, funcCallId: `${configId}-run-that-does-not-exist`});
    call.options[CONFIG_PATH] = serialize({
      type: 'static', configId: 'wf', friendlyName: 'Workflow',
      steps: [missing('config', 'Configuration'), missing('result', 'Result')],
    });
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'workflow');
    expectArray(entry.failedSteps ?? [], ['Configuration', 'Result']);
    expect(entry.nodes.tables.length + entry.nodes.scalars.length, 0);
  });

  test('saved steps that never ran are skipped, not reported', async () => {
    const call = DG.Script.create(
      '//name: rcPendingStepProbe\n//language: javascript\n//output: int x\nx = 1;').prepare();
    call.options[CONFIG_PATH] = serialize({
      type: 'static', configId: 'wf', steps: [{type: 'funccall', configId: 'pending', friendlyName: 'Pending'}],
    });
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'workflow');
    expect(entry.failedSteps === undefined, true);
    expect(entry.nodes.scalars.length + entry.nodes.tables.length, 0);
  });
});

category('RunComparison: workflow step paths', () => {
  test('nested pipelines join ids and names, unnamed steps fall back to their id', async () => {
    const state = {
      type: 'static', configId: 'wf', friendlyName: 'Workflow', steps: [
        {type: 'static', configId: 'stage', friendlyName: 'Stage', steps: [
          {type: 'static', configId: 'inner', steps: [
            {type: 'funccall', configId: 'calc', friendlyName: 'Calc', funcCallId: 'x'},
          ]},
          {type: 'funccall', configId: 'pending', friendlyName: 'Pending'},
        ]},
      ],
    } as unknown as PipelineSerializedState;
    const steps = collectWorkflowSteps(state);
    expectArray(steps.map((s) => s.path), ['wf/stage/inner/calc', 'wf/stage/pending']);
    expectArray(steps.map((s) => s.friendlyPath), ['Stage · inner · Calc', 'Stage · Pending']);
  });

  test('a root that is a single step uses its own id and name', async () => {
    const solo = {type: 'funccall', configId: 'solo', funcCallId: 's'} as unknown as PipelineSerializedState;
    const [step] = collectWorkflowSteps(solo);
    expect(step.path, 'solo');
    expect(step.friendlyPath, 'solo');
    expect(step.step.funcCallId, 's');
  });
});

category('RunComparison: function run extraction', () => {
  const probeCall = () => DG.Script.create([
    '//name: rcExtractionProbe',
    '//language: javascript',
    '//input: double a {caption: Data}',
    '//input: double x',
    '//input: string label',
    '//output: double b {caption: Data}',
    '//output: dataframe result',
    '//output: dataframe empty',
  ].join('\n')).prepare({a: 1, x: 2, label: 'l'});

  test('scalars and tables come from inputs and outputs, shared captions get the property name', async () => {
    const call = probeCall();
    call.inputs['x'] = NaN;
    call.outputs['b'] = 3;
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'time', [0, 1])]);
    call.outputs['result'] = df;
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'function');
    const scalars = new Map(entry.nodes.scalars.map((s) => [s.path, s]));
    expectArray([...scalars.keys()], ['a', 'x', 'b']);
    expect(scalars.get('a')!.name, 'Data (a)');
    expect(scalars.get('b')!.name, 'Data (b)');
    expect(scalars.get('a')!.value, 1);
    expect(scalars.get('x')!.value === null, true);
    expect(scalars.get('b')!.value, 3);
    expectArray(entry.nodes.tables.map((t) => t.path), ['result']);
    expect(entry.dataFrames.get('result')!.dart === df.dart, true);
  });

  test('a malformed workflow config degrades to an entry without steps', async () => {
    const call = probeCall();
    call.options[CONFIG_PATH] = '{not json';
    const entry = await entryFromFuncCall(call);
    expect(entry.sourceKind, 'workflow');
    expect(entry.nodes.scalars.length + entry.nodes.tables.length, 0);
    expect(entry.failedSteps === undefined, true);
  });
});

const snapshotProbe = () => {
  const call = DG.Script.create([
    '//name: rcSnapshotProbe',
    '//language: javascript',
    '//input: double a',
    '//output: double b',
    '//output: dataframe result',
  ].join('\n')).prepare({a: 1});
  call.outputs['b'] = 2;
  call.outputs['result'] = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'time', [0, 1])]);
  return call;
};

const callState = (isOutputOutdated: boolean) =>
  ({isOutputOutdated, isRunning: false, isRunnable: true, runError: undefined, pendingDependencies: []});

const nodePaths = (entry: {nodes: {scalars: {path: string}[], tables: {path: string}[]}}) =>
  [...entry.nodes.scalars, ...entry.nodes.tables].map((node) => node.path);

const outdatedPaths = (entry: {nodes: {scalars: {path: string, outdated?: boolean}[],
  tables: {path: string, outdated?: boolean}[]}}) =>
  [...entry.nodes.scalars, ...entry.nodes.tables].filter((node) => node.outdated).map((node) => node.path);

category('RunComparison: snapshots', () => {
  // the same tree twice, as saved (funcCallId) and as live (funcCall)
  const tree = (leaf: (configId: string, name: string) => object) => ({
    type: 'static', configId: 'wf', friendlyName: 'Workflow', steps: [
      leaf('config', 'Configuration'),
      {type: 'dynamic', configId: 'days', friendlyName: 'Days', steps: [leaf('day', 'Day'), leaf('day', 'Day')]},
    ],
  });

  test('live and saved states get the same step paths', async () => {
    const saved = tree((configId, friendlyName) =>
      ({type: 'funccall', configId, friendlyName, funcCallId: configId})) as unknown as PipelineSerializedState;
    const live = tree((configId, friendlyName) =>
      ({type: 'funccall', configId, friendlyName, uuid: configId})) as unknown as PipelineState;
    const savedSteps = collectWorkflowSteps(saved);
    const liveSteps = collectWorkflowSteps(live);
    expectArray(liveSteps.map((s) => s.path), savedSteps.map((s) => s.path));
    expectArray(liveSteps.map((s) => s.friendlyPath), savedSteps.map((s) => s.friendlyPath));
  });

  test('an edited step snapshot has the paths of its saved run', async () => {
    const call = snapshotProbe();
    const saved = await entryFromFuncCall(call);
    const step = {type: 'funccall', uuid: 'u', configId: 'step', funcCall: call} as unknown as PipelineState;
    const snapshot = entryFromSnapshot(snapshotState(step, {u: callState(false)}, {}), 'Current');
    expectArray(nodePaths(snapshot), nodePaths(saved));
    expect(snapshot.sourceKind, 'function');
    expect(snapshot.snapshot !== undefined, true);
    expect(saved.snapshot === undefined, true);
  });

  test('a snapshot keeps its tables when the source changes', async () => {
    const call = snapshotProbe();
    const snapshot = entryFromLiveCall(call, {isOutputOutdated: false}, 'Current');
    const source = call.outputs['result'] as DG.DataFrame;
    source.col('time')!.set(0, 100);
    const copy = snapshot.dataFrames.get('result')!;
    expect(copy.dart === source.dart, false);
    expect(copy.col('time')!.get(0), 0);
  });

  test('outputs of outdated steps are flagged, inputs are not', async () => {
    const state = {
      type: 'static', configId: 'wf', steps: [
        {type: 'funccall', configId: 'fresh', friendlyName: 'Fresh', uuid: 'f', funcCall: snapshotProbe()},
        {type: 'funccall', configId: 'stale', friendlyName: 'Stale', uuid: 's', funcCall: snapshotProbe()},
        {type: 'funccall', configId: 'loading', friendlyName: 'Loading', uuid: 'l'},
      ],
    } as unknown as PipelineState;
    const entry = entryFromSnapshot(snapshotState(state, {f: callState(false), s: callState(true)}, {}), 'Current');
    expect(entry.sourceKind, 'workflow');
    expectArray(outdatedPaths(entry), ['wf/stale/b', 'wf/stale/result']);
    expectArray(entry.outdatedSteps ?? [], ['Stale']);
    expect(nodePaths(entry).some((path) => path.startsWith('wf/loading')), false);
  });

  test('a saved run with outdated outputs flags them', async () => {
    const call = snapshotProbe();
    expect((await entryFromFuncCall(call)).outdatedSteps === undefined, true);
    call.options[OUTPUT_OUTDATED_PATH] = serialize(true);
    const entry = await entryFromFuncCall(call);
    expectArray(outdatedPaths(entry), ['b', 'result']);
    expectArray(entry.outdatedSteps ?? [], [entry.modelName]);
  });
});

category('RunComparison: driver states', () => {
  test('a live tree and its serialized state get the same step paths', async () => {
    const config: PipelineConfiguration = {
      id: 'wf', type: 'static', steps: [
        {id: 'config', nqName: 'Compute2:SimpleInputs2', friendlyName: 'Configuration'},
        {
          id: 'days', type: 'dynamic', friendlyName: 'Days',
          stepTypes: [{id: 'day', nqName: 'Compute2:SimpleInputs2', friendlyName: 'Day'}],
          initialSteps: ['day', 'day'],
        },
      ],
    };
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config), mockMode: true});
    await tree.init().toPromise();
    const live = tree.toState() as PipelineState & {steps: PipelineState[]};
    const saved = tree.toSerializedState() as PipelineSerializedState & {steps: PipelineSerializedState[]};
    const paths = (steps: {path: string, friendlyPath: string}[]) =>
      steps.map(({path, friendlyPath}) => `${path} | ${friendlyPath}`);
    expectArray(paths(collectWorkflowSteps(live)), paths(collectWorkflowSteps(saved)));
    expectArray(collectWorkflowSteps(live).map((s) => s.path), ['wf/config', 'wf/days/day[1]', 'wf/days/day[2]']);
    // a nested item is the root of its own saved run
    expectArray(paths(collectWorkflowSteps(live.steps[1])), paths(collectWorkflowSteps(saved.steps[1])));
  });
});

category('RunComparison: saved original', () => {
  test('a saved run loads with its outdated outputs and the paths of its snapshot', async () => {
    const call = DG.Func.byName('Compute2:SimpleInputs2').prepare({a: 1, b: 2, c: 'x'});
    call.outputs['aout'] = 1;
    call.outputs['bout'] = 2;
    call.outputs['cout'] = 'x';
    call.options[OUTPUT_OUTDATED_PATH] = serialize(true);
    const saved = await historyUtils.saveRun(call, {newId: true});
    try {
      const missing = 'rc-run-that-does-not-exist';
      const {entries, failed} = await loadRunEntries([saved.id, missing]);
      expectArray(failed, [missing]);
      expect(entries.length, 1);
      const snapshot = entryFromLiveCall(call, {isOutputOutdated: true}, 'Current');
      expectArray(nodePaths(entries[0]), nodePaths(snapshot));
      expectArray(outdatedPaths(entries[0]), ['aout', 'bout']);
      expectArray(outdatedPaths(snapshot), ['aout', 'bout']);
    } finally {
      await deleteSavedRun(saved.id);
    }
  });
});

// MockPipeline1 with a run first step and a never run second one, snapshotted like TreeWizard does
const liveSnapshot = async () => {
  const config = await callHandler<PipelineConfiguration>('Compute2:MockPipeline1', {}).toPromise();
  const tree = StateTree.fromInstanceConfig({
    config: await getProcessedConfig(config),
    instanceConfig: {steps: [{id: 'step1', initialValues: {a: 1, b: 2, res: 3}}, {id: 'step2', initialValues: {b: 4}}]},
  });
  await tree.init().toPromise();
  const values = <T, >(states: Record<string, {value: T}>) =>
    Object.fromEntries(Object.entries(states).map(([uuid, state$]) => [uuid, state$.value]));
  const snapshotOf = (state: PipelineState) =>
    snapshotState(state, values(tree.getFuncCallStates()), values(tree.getConsistency()));
  return {tree, snapshot: snapshotOf(tree.toState()), snapshotOf};
};

category('RunComparison: saving snapshots', () => {
  test('a workflow snapshot keeps its values when the source changes', async () => {
    const {tree, snapshot} = await liveSnapshot();
    const entry = entryFromSnapshot(snapshot, 'Current');
    const live = (tree.toState() as PipelineState & {steps: {funcCall?: DG.FuncCall}[]}).steps[0].funcCall!;
    live.inputs['a'] = 100;
    tree.close();
    expect(entry.nodes.scalars.find((node) => node.path === 'pipeline1/step1/a')!.value, 1);
    expectArray(entry.outdatedSteps ?? [], ['step2']);
  });

  test('a workflow snapshot saves and loads with its values and paths', async () => {
    const {tree, snapshot} = await liveSnapshot();
    tree.close();
    const entry = entryFromSnapshot(snapshot, 'Current');
    const runId = (await entry.snapshot!.save!({title: 'rc saved snapshot'})).id;
    try {
      const {entries, failed} = await loadRunEntries([runId]);
      expectArray(failed, []);
      expect(entries[0].name, 'rc saved snapshot');
      expectArray(nodePaths(entries[0]), nodePaths(entry));
      expectArray(entries[0].nodes.scalars.map((node) => node.value), entry.nodes.scalars.map((node) => node.value));
      expectArray(entries[0].outdatedSteps ?? [], ['step2']);
    } finally {
      await deleteSavedWorkflow(runId);
    }
  });

  test('a workflow step snapshot saves to the step history with the paths of its saved run', async () => {
    const {tree, snapshotOf} = await liveSnapshot();
    const step = (tree.toState() as PipelineState & {steps: PipelineState[]}).steps[0];
    const entry = entryFromSnapshot(snapshotOf(step), 'Current', {stepHistory: true});
    tree.close();
    const runId = (await entry.snapshot!.save!({title: 'rc saved step'})).id;
    try {
      const saved = await historyUtils.loadRun(runId);
      expect(saved.options[STEP_HISTORY_OPTION], 'true');
      expect(getRunName(saved), 'rc saved step');
      const reloaded = await entryFromFuncCall(saved);
      expect(reloaded.sourceKind, 'function');
      expectArray(nodePaths(reloaded), nodePaths(entry));
      expectArray(reloaded.nodes.scalars.map((node) => node.value), entry.nodes.scalars.map((node) => node.value));
    } finally {
      await deleteSavedRun(runId);
    }
  });

  test('a function snapshot RFV would not save has no save', async () => {
    const save = (state: {isOutputOutdated: boolean, runError?: string}) =>
      entryFromLiveCall(snapshotProbe(), state, 'Current').snapshot!.save;
    expect(save({isOutputOutdated: true}) === undefined, true);
    expect(save({isOutputOutdated: false, runError: 'Failed'}) === undefined, true);
    expect(save({isOutputOutdated: false}) !== undefined, true);
  });

  test('a function snapshot saves as a run of its own', async () => {
    const call = snapshotProbe();
    const snapshot = entryFromLiveCall(call, {isOutputOutdated: false}, 'Current');
    call.inputs['a'] = 100;
    const runId = (await snapshot.snapshot!.save!({title: 'rc saved call'})).id;
    try {
      expect(runId === call.id, false);
      const saved = await historyUtils.loadRun(runId);
      expect(getRunName(saved), 'rc saved call');
      expect(saved.inputs['a'], 1);
      expect(readOutputOutdated(saved), false);
    } finally {
      await deleteSavedRun(runId);
    }
  });
});
