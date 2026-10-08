import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig, PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {isFuncCallSerializedState, PipelineSerializedState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

// These tests run in real time: in-place dataframe edits reach the driver through platform
// onDataChanged events, which the RxJS test scheduler cannot drive.

const makeDf = (v: number) => DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', new Float32Array([v]))]);
const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

function makeConfig(restriction?: string): PipelineConfiguration {
  return {
    id: 'root',
    type: 'static',
    states: ['src', 'copy'],
    steps: [{id: 't', nqName: 'LibTests:TestToleranceInputs'}],
    links: [
      {id: 'feed', from: 'in:src', to: 'out:t/df', ...(restriction ? {defaultRestrictions: restriction as any} : {})},
      {id: 'to-state', from: 'in:src', to: 'out:copy'},
    ],
  };
}

async function countClones(fn: () => Promise<void>) {
  const original = DG.DataFrame.prototype.clone;
  let clones = 0;
  DG.DataFrame.prototype.clone = function(this: DG.DataFrame, ...args: any[]) {
    clones++;
    return (original as any).apply(this, args);
  };
  try {
    await fn();
  } finally {
    DG.DataFrame.prototype.clone = original;
  }
  return clones;
}

category('ComputeUtils: Driver dataframe copies: platform', async () => {
  test('A dataframe clone keeps the id', async () => {
    const df = makeDf(1);
    expectDeepEqual(df.clone().id, df.id);
  });

  test('A clone of a re-ided dataframe gets the original id back', async () => {
    const df = makeDf(1);
    const copy = df.clone();
    copy.id = crypto.randomUUID();
    expectDeepEqual(copy.id !== df.id, true, {prefix: 'Copy re-ided'});
    expectDeepEqual(copy.clone().id, df.id);
  });

  test('Uploading a dataframe whose id is taken stores it under a new table', async () => {
    const df = makeDf(1);
    const other = df.clone();
    other.col('x')!.set(0, 2);
    const ids: string[] = [];
    try {
      ids.push(await grok.dapi.tables.uploadDataFrame(df));
      ids.push(await grok.dapi.tables.uploadDataFrame(other));
      expectDeepEqual(ids[0] !== ids[1], true, {prefix: 'Second upload got its own table'});
      expectDeepEqual((await grok.dapi.tables.getTable(ids[0])).col('x')!.get(0), 1, {prefix: 'First table kept'});
    } finally {
      for (const id of ids) {
        const info = await grok.dapi.tables.find(id);
        if (info)
          await grok.dapi.tables.delete(info);
      }
    }
  });
});

category('ComputeUtils: Driver dataframe copies', async () => {
  let pconf: PipelineConfigurationProcessed;
  let pconfNone: PipelineConfigurationProcessed;

  before(async () => {
    pconf = await getProcessedConfig(makeConfig());
    pconfNone = await getProcessedConfig(makeConfig('none'));
  });

  async function build(config: PipelineConfigurationProcessed, readonly = false) {
    let tree: StateTree;
    if (readonly) {
      const state = StateTree.fromPipelineConfig({config, mockMode: true}).toSerializedState();
      const steps = (state as any).steps.map((s: PipelineSerializedState) => isFuncCallSerializedState(s) ? {...s, isReadonly: true} : s);
      tree = StateTree.fromInstanceState({state: {...state, steps} as PipelineSerializedState, config, isReadonly: false, mockMode: true});
    } else
      tree = StateTree.fromPipelineConfig({config, mockMode: true});
    await tree.init().toPromise();
    return tree;
  }

  const step = (tree: StateTree) => tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
  const input = (tree: StateTree) => step(tree).getStateStore().getState('df') as DG.DataFrame | undefined;
  const snapshot = (tree: StateTree) => step(tree).instancesWrapper.inputRestrictions$.value.df?.assignedValue as DG.DataFrame | undefined;
  const inconsistent = (tree: StateTree) => step(tree).consistencyInfo$.value.df?.inconsistent;

  async function settle(tree: StateTree) {
    await sleep(200);
    await tree.linksState.waitForLinks().toPromise();
  }

  async function write(tree: StateTree, src: DG.DataFrame) {
    tree.nodeTree.root.getItem().getStateStore().setState('src', src);
    await settle(tree);
  }

  async function runStep(tree: StateTree) {
    await tree.runStep(step(tree).uuid, {res: 1}).toPromise();
    await settle(tree);
  }

  async function editInPlace(tree: StateTree, df: DG.DataFrame) {
    df.col('x')!.set(0, 99);
    await settle(tree);
  }

  test('Not run step gets its own input and a separate snapshot', async () => {
    const tree = await build(pconf);
    const src = makeDf(1);
    const clones = await countClones(() => write(tree, src));
    expectDeepEqual(clones, 3, {prefix: 'Clones: input, snapshot, workflow state'});
    expectDeepEqual(input(tree) !== src, true, {prefix: 'Input is not the source'});
    expectDeepEqual(snapshot(tree) !== src && snapshot(tree) !== input(tree), true, {prefix: 'Snapshot is separate'});
    expectDeepEqual(input(tree)!.id !== src.id, true, {prefix: 'Input has its own id'});
    await editInPlace(tree, input(tree)!);
    expectDeepEqual(inconsistent(tree), true, {prefix: 'In-place edit is flagged'});
  });

  test('Already run step keeps its input and stores one snapshot', async () => {
    const tree = await build(pconf);
    await runStep(tree);
    const src = makeDf(1);
    const clones = await countClones(() => write(tree, src));
    expectDeepEqual(clones, 2, {prefix: 'Clones: snapshot, workflow state'});
    expectDeepEqual(input(tree) === undefined, true, {prefix: 'Input untouched'});
    expectDeepEqual(snapshot(tree) !== src, true, {prefix: 'Snapshot is not the source'});
    expectDeepEqual(inconsistent(tree), true, {prefix: 'Step is inconsistent'});
  });

  test('Reset to consistent gives the input its own copy and id', async () => {
    const tree = await build(pconf);
    await runStep(tree);
    const src = makeDf(1);
    await write(tree, src);
    const clones = await countClones(async () => {
      step(tree).instancesWrapper.setToConsistent('df');
      await settle(tree);
    });
    expectDeepEqual(clones, 1, {prefix: 'Clones'});
    expectDeepEqual(input(tree) !== src && input(tree) !== snapshot(tree), true, {prefix: 'Input is a separate copy'});
    expectDeepEqual(input(tree)!.id !== src.id, true, {prefix: 'Input has its own id'});
    await editInPlace(tree, input(tree)!);
    expectDeepEqual(inconsistent(tree), true, {prefix: 'In-place edit is flagged'});
  });

  test('Update gives the input its own copy and id', async () => {
    const tree = await build(pconf);
    await runStep(tree);
    const src = makeDf(1);
    await write(tree, src);
    const clones = await countClones(async () => {
      await step(tree).instancesWrapper.overrideToConsistent().toPromise();
      await settle(tree);
    });
    expectDeepEqual(clones, 1, {prefix: 'Clones'});
    expectDeepEqual(input(tree) !== src && input(tree) !== snapshot(tree), true, {prefix: 'Input is a separate copy'});
    expectDeepEqual(input(tree)!.id !== src.id, true, {prefix: 'Input has its own id'});
    await editInPlace(tree, input(tree)!);
    expectDeepEqual(inconsistent(tree), true, {prefix: 'In-place edit is flagged'});
  });

  test('Read-only step stores one snapshot', async () => {
    const tree = await build(pconf, true);
    const src = makeDf(1);
    const clones = await countClones(() => write(tree, src));
    expectDeepEqual(clones, 2, {prefix: 'Clones: snapshot, workflow state'});
    expectDeepEqual(input(tree) === undefined, true, {prefix: 'Input untouched'});
    expectDeepEqual(snapshot(tree) !== src, true, {prefix: 'Snapshot is not the source'});
    expectDeepEqual(inconsistent(tree), true, {prefix: 'Step is inconsistent'});
  });

  test('Restriction none copies only the input', async () => {
    const tree = await build(pconfNone);
    const src = makeDf(1);
    const clones = await countClones(() => write(tree, src));
    expectDeepEqual(clones, 2, {prefix: 'Clones: input, workflow state'});
    expectDeepEqual(snapshot(tree) === undefined, true, {prefix: 'No snapshot'});
    expectDeepEqual(input(tree) !== src, true, {prefix: 'Input is not the source'});
    expectDeepEqual(input(tree)!.id !== src.id, true, {prefix: 'Input has its own id'});
  });

  test('A workflow state gets its own copy', async () => {
    const tree = await build(pconf);
    const src = makeDf(1);
    await write(tree, src);
    const copy = tree.nodeTree.root.getItem().getStateStore().getState('copy') as DG.DataFrame;
    expectDeepEqual(copy !== src, true, {prefix: 'State is not the source'});
    expectDeepEqual(copy.id !== src.id, true, {prefix: 'State has its own id'});
  });
});
