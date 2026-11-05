import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {Observable, of, throwError} from 'rxjs';
import {catchError, filter, switchMap, take, timeout} from 'rxjs/operators';

// the first add step drives a message, meta, a restricted value and a title on the first mul step,
// a message on the workflow node and the visibility of the workflow's report action
const config: PipelineConfiguration = {
  id: 'scenarios',
  type: 'dynamic',
  stepTypes: [
    {id: 'add', nqName: 'LibTests:TestAdd2'},
    {id: 'mul', nqName: 'LibTests:TestMul2'},
  ],
  initialSteps: ['add', 'mul'],
  actions: [{id: 'report', position: 'buttons', from: [], to: 'out:first(mul)/a', handler() {}}],
  links: [{
    id: 'r',
    type: 'rule',
    from: 'x:first(add)/a',
    to: ['y:first(mul)/a', 'z:first(mul)/b', 'title:first(mul)/title'],
    effects: [
      'hide(z, when: gt(x, 5))',
      'error(y, "Too big", when: gt(x, 5))',
      'set(y, x, restriction: "restricted")',
      'describe(title, cat("From ", x))',
    ],
  }, {
    id: 'msg',
    type: 'rule',
    from: 'x:first(add)/a',
    to: 'scenarios(node)',
    effects: ['error(scenarios, "Too big", when: gt(x, 5))'],
  }, {
    id: 'act',
    type: 'rule',
    from: 'x:first(add)/a',
    to: [],
    actions: 'scenarios',
    effects: ['hideAction(scenarios, "report", when: gt(x, 5))'],
  }],
} as PipelineConfiguration;

// without links, so only the tree commands and input edits mark the workflow edited
const plainConfig = {...config, actions: [], links: []} as PipelineConfiguration;

// waits for a published value and fails with `what` when it does not come
const until = <T>(obs$: Observable<T>, pred: (value: T) => boolean, what: string) => obs$.pipe(
  filter(pred),
  take(1),
  timeout(5000),
  catchError(() => throwError(new Error(`timed out waiting for ${what}`))),
).toPromise();

const initPipeline = async (driver: Driver, pipeline = config) => await driver.sendCommand(
  {event: 'initPipeline', provider: '', config: await getProcessedConfig(pipeline)}) as StateTree;
const steps = (driver: Driver) => (driver.currentState$.value as any).steps.map((step: any) => step.uuid) as string[];
const editInput = (tree: StateTree, idx: number, io: string, value: any) =>
  (tree.nodeTree.getNode([{idx}]).getItem() as FuncCallNode).getStateStore().editState(io, value);
const settled = (driver: Driver) => until(driver.treeMutationsLocked$, (locked) => !locked, 'the tree to unlock');

category('ComputeUtils: Driver published states', async () => {
  test('The driver publishes the states of every node and republishes them after a mutation', async () => {
    const driver = new Driver(true);
    const tree = await initPipeline(driver);
    const root = (driver.currentState$.value as any).uuid as string;
    const [, mul] = steps(driver);
    editInput(tree, 0, 'a', 10);
    await until(driver.currentValidations$.value[mul], (v) => !!v.a?.errors?.length, 'the validation');
    await until(driver.currentMeta$.value[mul].pipe(switchMap((meta) => meta.b ?? of(undefined))),
      (meta) => meta?.hidden === true, 'the meta');
    await until(driver.currentConsistency$.value[mul], (c) => c.a?.restriction === 'restricted', 'the consistency');
    await until(driver.currentCallsState$.value[mul], (state) => state != null, 'the call state');
    await until(driver.nodesDescriptions$.value[mul], (d) => d?.title === 'From 10', 'the description');
    await until(driver.currentPipelineValidations$.value[root], (v) => !!v?.errors?.length, 'the node message');
    await until(driver.hiddenActions$, (hidden) => hidden.size === 1, 'the hidden action');
    const linkIds = [...new Set(driver.currentLinks$.value.filter((l) => !l.isAction).map((l) => l.id))].sort();
    expectDeepEqual(linkIds, ['act::nodemeta', 'msg::pipelineValidator', 'r::data', 'r::meta', 'r::nodemeta',
      'r::validator'], {prefix: 'links'});

    await settled(driver);
    await driver.sendCommand({event: 'addDynamicItem', parentUuid: root, itemId: 'mul', position: 2});
    const all = steps(driver);
    expectDeepEqual(all.length, 3, {prefix: 'steps after adding'});
    const published = {
      validations: driver.currentValidations$.value,
      meta: driver.currentMeta$.value,
      consistency: driver.currentConsistency$.value,
      calls: driver.currentCallsState$.value,
      descriptions: driver.nodesDescriptions$.value,
    };
    for (const [name, byUuid] of Object.entries(published))
      expectDeepEqual(all.every((uuid) => uuid in byUuid), true, {prefix: `${name} after adding`});
    expectDeepEqual(root in driver.currentPipelineValidations$.value, true, {prefix: 'node messages after adding'});
    driver.close();
  });

  test('Unsaved edits follow input edits and tree commands and wait for locks', async () => {
    const driver = new Driver(true);
    const restart = async () => {
      await settled(driver);
      const tree = await initPipeline(driver, plainConfig);
      expectDeepEqual(driver.hasNotSavedEdits$.value, false, {prefix: 'after init'});
      return tree;
    };
    const edited = (what: string) => until(driver.hasNotSavedEdits$, (value) => value, what);

    let tree = await restart();
    editInput(tree, 0, 'b', 1);
    await edited('an input edit');

    for (const command of [
      () => ({event: 'addDynamicItem', parentUuid: (driver.currentState$.value as any).uuid, itemId: 'mul',
        position: 2}),
      () => ({event: 'moveDynamicItem', uuid: steps(driver)[1], position: 0}),
      () => ({event: 'duplicateDynamicItem', uuid: steps(driver)[0]}),
      () => ({event: 'removeDynamicItem', uuid: steps(driver)[1]}),
    ]) {
      tree = await restart();
      await settled(driver);
      const msg = command() as any;
      await driver.sendCommand(msg);
      await edited(msg.event);
    }

    tree = await restart();
    await settled(driver);
    tree.globalROLocked$.next(true);
    editInput(tree, 0, 'b', 2);
    expectDeepEqual(driver.hasNotSavedEdits$.value, false, {prefix: 'while locked'});
    tree.globalROLocked$.next(false);
    await edited('the edit made while locked');
    driver.close();
  });

  test('Hidden actions follow the current tree', async () => {
    const driver = new Driver(true);
    let tree = await initPipeline(driver);
    const root = (driver.currentState$.value as any).uuid as string;
    const report = tree.linksState.getNodeActionsData(root)![0].uuid;
    expectDeepEqual(driver.hiddenActions$.value.size, 0, {prefix: 'before the edit'});
    editInput(tree, 0, 'a', 10);
    await until(driver.hiddenActions$, (hidden) => hidden.has(report), 'the report action to hide');
    await settled(driver);
    tree = await initPipeline(driver);
    expectDeepEqual(driver.hiddenActions$.value.size, 0, {prefix: 'after a new init'});
    editInput(tree, 0, 'a', 10);
    const newRoot = (driver.currentState$.value as any).uuid as string;
    const newReport = tree.linksState.getNodeActionsData(newRoot)![0].uuid;
    await until(driver.hiddenActions$, (hidden) => hidden.has(newReport) && !hidden.has(report),
      'the new tree\'s report action to hide');
    driver.close();
  });
});
