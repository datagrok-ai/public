import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {compileEffect} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/rule-formula';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrows, expectThrowsAsync} from '../../../test-utils';

// raises cooling/a above cooling/b, by `step` when the button passes it
const raise = (extra: Record<string, any> = {}) => ({
  id: 'raise',
  position: 'none' as const,
  from: 'b:cooling/b',
  to: 'a:cooling/a',
  handler({controller}: any) {
    controller.setAll('a', controller.getFirst('b') + (controller.getAdditionalParam('step') ?? 10));
  },
  ...extra,
});

const dayAction = {
  id: 'raise',
  position: 'none' as const,
  from: 'in:all(day)/a',
  to: 'out:all(day)/b',
  handler({controller}: any) {
    controller.setAll('out', 1);
  },
};

const workflow = (links: any[], rootActions: any[] = [raise()], dayActions: any[] = [], days = 0) => ({
  id: 'pipeline1',
  type: 'static',
  steps: [
    {id: 'cooling', nqName: 'LibTests:TestAdd2'},
    {
      id: 'days',
      type: 'parallel',
      stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}],
      initialSteps: Array.from({length: days}, () => ({id: 'day'})),
      actions: dayActions,
    },
  ],
  links,
  actions: rootActions,
}) as PipelineConfiguration;

const tooCold = (effects: any[], extra: Record<string, any> = {}) => ({
  id: 'tooCold',
  type: 'rule',
  debounce: 0,
  from: ['a:cooling/a', 'b:cooling/b'],
  to: 't:cooling/a',
  actions: 'fixes',
  effects,
  ...extra,
});

// extra keys fail too, so a plain message is told apart from one with actions
const expectExact = (actual: any, expected: any, options: {prefix?: string} = {}) =>
  expectDeepEqual(actual, expected, {...options, forbidAdditionalProps: true});

const OFFER = 'actions: [action(fixes, "raise", "Raise it")]';

category('ComputeUtils: Driver rule message actions', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  const cooling = (tree: StateTree) => tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
  const messages = (tree: StateTree) => cooling(tree).validationInfo$.value['a'];
  const rootAction = (tree: StateTree) =>
    tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid)![0].uuid;
  const coolingAction = (tree: StateTree) =>
    tree.linksState.getNodeActionsData(cooling(tree).uuid)![0].uuid;

  function run(pconf: any, body: (tree: StateTree, flush: () => void) => void) {
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      body(tree, flush);
    });
  }

  const setInputs = (tree: StateTree, flush: () => void, values: Record<string, any>) => {
    for (const [io, value] of Object.entries(values))
      cooling(tree).getStateStore().setState(io, value);
    flush();
  };

  test('A message offers the action and running it applies the action', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([`error(t, "Too cold", when: lt(a, b), ${OFFER})`])]));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1, b: 20});
      const uuid = rootAction(tree);
      expectExact(messages(tree)?.errors,
        [{description: 'Too cold', actions: [{actionName: 'Raise it', action: uuid}]}], {prefix: 'offered'});
      tree.runAction(uuid);
      flush();
      expectDeepEqual([cooling(tree).getStateStore().getState('a'), messages(tree)?.errors ?? []], [30, []],
        {prefix: 'after running'});
    });
  });

  test('Only the messages with actions offer them', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([
      `error(t, "Too cold", when: lt(a, b), ${OFFER})`,
      'error(t, "Not a model input")',
      'warning(t, "Check b", when: gt(b, 10))',
    ])]));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1, b: 20});
      const uuid = rootAction(tree);
      expectExact(messages(tree)?.errors, [
        {description: 'Too cold', actions: [{actionName: 'Raise it', action: uuid}]},
        {description: 'Not a model input'},
      ], {prefix: 'errors'});
      expectExact(messages(tree)?.warnings, [{description: 'Check b'}], {prefix: 'warnings'});
    });
  });

  test('Every item of an array message offers the actions', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([`warning(t, ["One", "Two"], ${OFFER})`])]));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1});
      const offered = [{actionName: 'Raise it', action: rootAction(tree)}];
      expectExact(messages(tree)?.warnings,
        [{description: 'One', actions: offered}, {description: 'Two', actions: offered}]);
    });
  });

  test('The button name and parameters are formulas and reach the action', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([
      'error(t, "Too cold", when: lt(a, b), ' +
        'actions: [action(fixes, "raise", cat("Raise above ", b), additionalParams: obj(step: 5))])',
    ])]));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1, b: 20});
      const [item] = (messages(tree)?.errors as any[])[0].actions;
      expectExact(item, {actionName: 'Raise above 20', action: rootAction(tree), additionalParams: {step: 5}},
        {prefix: 'b is 20'});
      setInputs(tree, flush, {b: 40});
      const [next] = (messages(tree)?.errors as any[])[0].actions;
      expectDeepEqual(next.actionName, 'Raise above 40', {prefix: 'b is 40'});
      tree.runAction(next.action, next.additionalParams);
      flush();
      expectDeepEqual(cooling(tree).getStateStore().getState('a'), 45, {prefix: 'raised by the step'});
    });
  });

  test('A message whose action is not found keeps the plain message', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([
      'error(t, "Too cold", actions: [action(fixes, "missing", "Raise it")])',
    ])]));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1});
      expectExact(messages(tree)?.errors, [{description: 'Too cold'}]);
    });
  });

  test('Actions are found on the node they are shown on, as for validators', async () => {
    const shown = [raise({visibleOn: 'cooling'})];
    const onCooling = await getProcessedConfig(workflow([tooCold([`error(t, "Too cold", ${OFFER})`],
      {actions: 'fixes:cooling'})], shown));
    run(onCooling, (tree, flush) => {
      setInputs(tree, flush, {a: 1});
      expectExact(messages(tree)?.errors,
        [{description: 'Too cold', actions: [{actionName: 'Raise it', action: coolingAction(tree)}]}],
        {prefix: 'node where it is shown'});
    });
    const onRoot = await getProcessedConfig(workflow([tooCold([`error(t, "Too cold", ${OFFER})`])], shown));
    run(onRoot, (tree, flush) => {
      setInputs(tree, flush, {a: 1});
      expectExact(messages(tree)?.errors, [{description: 'Too cold'}], {prefix: 'node where it is defined'});
    });
  });

  test('A message offers an action recreated by a tree mutation', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([`error(t, "Too cold", ${OFFER})`],
      {actions: 'fixes:days'})], [], [dayAction], 1));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1});
      const offered = () => (messages(tree)?.errors as any[])[0].actions[0].action as string;
      expectDeepEqual(tree.linksState.actions.has(offered()), true, {prefix: 'before the mutation'});
      tree.addSubTree(tree.nodeTree.getNode([{idx: 1}]).getItem().uuid, 'day', 1).subscribe();
      flush();
      expectDeepEqual(tree.linksState.actions.has(offered()), true, {prefix: 'after the mutation'});
    });
  });

  test('A hidden action is still offered', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([`error(t, "Too cold", ${OFFER})`])],
      [raise({showWhen: 'guard: missingStep'})]));
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1});
      const uuid = rootAction(tree);
      expectDeepEqual(tree.linksState.isActionVisible(uuid), false, {prefix: 'hidden'});
      expectExact(messages(tree)?.errors,
        [{description: 'Too cold', actions: [{actionName: 'Raise it', action: uuid}]}], {prefix: 'offered'});
    });
  });

  test('A separate message offers its actions', async () => {
    const pconf = await getProcessedConfig(workflow([tooCold([
      'error(t, "Too cold", when: lt(a, b), separate: true, actions: [action(fixes, "raise", cat("Raise above ", b))])',
      'hide(t, when: gt(a, 100))',
    ])]));
    const link = pconf.links!.find((l) => l.id === 'tooCold::validator:1')! as any;
    expectDeepEqual([link.from.map((io: any) => io.name), link.actions.map((io: any) => io.name)],
      [['a', 'b'], ['fixes']], {prefix: 'link'});
    run(pconf, (tree, flush) => {
      setInputs(tree, flush, {a: 1, b: 20});
      expectExact(messages(tree)?.errors,
        [{description: 'Too cold', actions: [{actionName: 'Raise above 20', action: rootAction(tree)}]}]);
    });
  });
});

category('ComputeUtils: Driver rule message actions expansion', async () => {
  const rule = (effects: any[], extra: Record<string, any> = {}) => getProcessedConfig(workflow([tooCold(effects, {
    actions: ['fixes', 'other:days'],
    ...extra,
  })], [raise()], [dayAction], 1));

  test('A validator link gets only the actions queries its messages use', async () => {
    const pconf = await rule([`error(t, "Too cold", ${OFFER})`, 'hideAction(other, "raise")']);
    const actions = Object.fromEntries(pconf.links!.map((link: any) =>
      [link.id, (link.actions ?? []).map((io: any) => io.name)]));
    expectExact(actions, {'tooCold::validator': ['fixes'], 'tooCold::nodemeta': ['other']});
  });

  test('Reject invalid message actions', async () => {
    await expectThrowsAsync(() => rule(['error(t, "Too cold", actions: [action(zz, "raise", "Raise it")])',
      'hideAction(other, "raise")', 'hideAction(fixes, "raise")']),
    /message action references unknown actions alias zz/);
    await expectThrowsAsync(() => rule([{effect: 'error', targets: 't', message: 'x',
      actions: [{actions: 'fixes', action: '', actionName: 'x'}]}, 'hideAction(other, "raise")']),
    /message action on fixes needs an action id/);
    await expectThrowsAsync(() => rule([{effect: 'hide', targets: 't',
      actions: [{actions: 'fixes', action: 'raise', actionName: 'x'}]}, 'hideAction(other, "raise")']),
    /on error, warning or notification only/);
    await expectThrowsAsync(() => rule([`error(p, "Too cold", ${OFFER})`, 'hide(t)', 'hideAction(other, "raise")'],
      {to: ['t:cooling/a', 'p(node):days']}),
    /\(node\) alias p messages take no actions/);
  });

  test('Formulas take actions on messages only', async () => {
    expectExact(compileEffect(
      'error(t, "m", actions: [action(f, "a", "A"), action(g, "b", cat("B", x), additionalParams: obj(k: 1))])'), {
      effect: 'error', targets: 't', message: 'm', actions: [
        {actions: 'f', action: 'a', actionName: 'A'},
        {actions: 'g', action: 'b', actionName: {cat: ['B', {var: 'x'}]}, additionalParams: {literal: {k: 1}}},
      ],
    });
    expectThrows(() => compileEffect('error(t, "m", actions: action(f, "a", "A"))'),
      /actions are a list of action\(alias, "id", name\) calls/);
    expectThrows(() => compileEffect('error(t, "m", actions: [cat(f, "a")])'),
      /actions are a list of action\(alias, "id", name\) calls/);
    expectThrows(() => compileEffect('error(t, "m", actions: [action(f.x, "a", "A")])'),
      /The actions alias is an alias/);
    expectThrows(() => compileEffect('error(t, "m", actions: [action(f, a, "A")])'), /The action id must be a string/);
    expectThrows(() => compileEffect('error(t, "m", actions: [action(f, "a")])'), /action needs actionName/);
    expectThrows(() => compileEffect('error(t, "m", actions: [action(f, "a", "A", params: 1)])'),
      /action has no option params \(expected additionalParams\)/);
    expectThrows(() => compileEffect('hide(t, actions: [action(f, "a", "A")])'), /hide has no option actions/);
    expectThrows(() => compileEffect('error(t, action(f, "a", "A"))'),
      /action is only allowed in the actions option of a message/);
  });
});
