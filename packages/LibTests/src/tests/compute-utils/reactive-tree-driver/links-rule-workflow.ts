import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {PipelineNodeBase} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrowsAsync, twoSteps} from '../../../test-utils';

type TreeMode = {batchLinks?: boolean, annotationChecks?: boolean};
const driverMode: TreeMode = {batchLinks: true, annotationChecks: true};

const withAction = (rule: any, action: Record<string, any> = {}): PipelineConfiguration => ({
  ...twoSteps([rule]),
  actions: [{id: 'addDay', position: 'none', from: [], to: 'out:step2/a', handler() {}, ...action}],
} as PipelineConfiguration);

const daysPipeline = (rule: any, initialDays = 0): PipelineConfiguration => ({
  id: 'pipeline1',
  type: 'static',
  steps: [{
    id: 'days',
    type: 'parallel',
    stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}],
    initialSteps: Array.from({length: initialDays}, () => ({id: 'day'})),
    links: [rule],
  }],
});

const hasErrors = (node: PipelineNodeBase) =>
  Object.values(node.pipelineValidations$.value).some((result) => !!result?.errors?.length);

category('ComputeUtils: Driver links rule: workflow level', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Workflow effects expand into node meta and pipeline validator links', async () => {
    const pconf = await getProcessedConfig(withAction({
      id: 'r', type: 'rule', from: 'm:step1/a', to: ['t:step2/a', 'title:step2/title', 'p(node):step2'],
      actions: 'act',
      effects: [
        {effect: 'hide', targets: 't'},
        {effect: 'describe', targets: 'title', value: {var: 'm'}},
        {effect: 'hideAction', targets: 'act', action: 'addDay'},
        {effect: 'error', targets: 'p', message: 'bad'},
      ],
    }));
    const links = pconf.links!;
    expectDeepEqual(links.map((l) => l.id), ['r::meta', 'r::nodemeta', 'r::pipelineValidator']);
    expectDeepEqual(links.map((l) => l.type), ['meta', 'nodemeta', 'pipelineValidator']);
    expectDeepEqual(links.map((l) => l.to.map((io) => io.name)), [['t'], ['title'], ['p']]);
    expectDeepEqual(links[1].actions!.map((io) => io.name), ['act']);
    expectDeepEqual(links[2].to[0].flags, ['node']);
    expectDeepEqual(links.map((l) => l.params!.effects.map((e: any) => e.effect)),
      [['hide'], ['describe', 'hideAction'], ['error']]);
  });

  test('Reject invalid workflow rules', async () => {
    const rule = (extra: any) => getProcessedConfig(withAction({id: 'r', type: 'rule', from: 'm:step1/a', ...extra}));
    await expectThrowsAsync(() => rule({to: 't:step2/a', effects: [{effect: 'describe', targets: 't', value: 1}]}),
      /must end with one of/);
    await expectThrowsAsync(() => rule({actions: 'act', effects: [{effect: 'hide', targets: 'act'}]}),
      /takes only hideAction or showAction/);
    await expectThrowsAsync(() => rule({to: 't:step2/a', effects: [
      {effect: 'hideAction', targets: 't', action: 'addDay'},
    ]}), /targets unknown actions alias t/);
    await expectThrowsAsync(() => rule({to: 'p(node):step2', effects: [{effect: 'set', targets: 'p', value: 1}]}),
      /takes only error, warning or notification/);
    await expectThrowsAsync(() => rule({to: 'title:step2/title', effects: [
      {effect: 'describe', targets: 'title', value: 1}, {effect: 'hide', targets: 'title'},
    ]}), /targeted as both/);
    await expectThrowsAsync(() => rule({to: ['t:step2/a', 'p(node):step2'], effects: [
      {effect: 'error', targets: ['t', 'p'], message: 'x'},
    ]}), /mixes targets of different kinds/);
    await expectThrowsAsync(() => rule({to: 'title:step2/title', effects: [
      {effect: 'describe', targets: 'title', value: {var: '$position'}},
    ]}), /\$position needs a base/);
    await expectThrowsAsync(() => rule({actions: 'act', to: 't:step2/a', effects: ['hide(t)']}),
      /actions alias act is not targeted by any effect/);
    await expectThrowsAsync(() => rule({actions: 'm', effects: ['hideAction(m, "addDay")']}),
      /actions alias m collides with an io alias/);
    await expectThrowsAsync(() => rule({
      actions: 'act', effects: ['hideAction(act, "addDay")', 'showAction(act, "addDay")'],
    }), /hide\/show applied twice to action act\/addDay/);
    await expectThrowsAsync(() => rule({actions: 'act_(template):step2', effects: ['hideAction(act_step2, "x")']}),
      /\(template\) flag is not allowed in rule queries/);
    await expectThrowsAsync(() => rule({actions: 'act', effects: [{effect: 'hideAction', targets: 'act'}]}),
      /hideAction needs an action id/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{id: 'l', from: 'in:step1/a', to: 'out(node):step2'}])),
      /output out uses \(node\); only rules and pipeline validators allow it/);
    const flagged = await getProcessedConfig(twoSteps([{
      id: 'pv', type: 'pipelineValidator', from: 'in:step1/a', to: 'out(node):step2', handler() {},
    }]));
    expectDeepEqual(flagged.links![0].to[0].flags, ['node'], {prefix: 'pipeline validator accepts (node)'});
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{id: 'l', from: 'in(node):step1', to: 'out:step2/a'}])),
      /\(node\) flag is not allowed in input queries/);
    await expectThrowsAsync(() => rule({from: 'm(node):step1', to: 't:step2/a', effects: ['hide(t)']}),
      /\(node\) flag is not allowed in input queries/);
    const hiddenByNode = withAction({id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a', effects: ['hide(t)']},
      {hideWhen: 'guard(node): step1'});
    await expectThrowsAsync(() => getProcessedConfig(hiddenByNode), /\(node\) flag is not allowed in hideWhen queries/);
  });

  test('Describe sets a description field and clears it when off', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 'title:step2/title', when: 'gt(m, 0)',
      effects: ['describe(title, cat("m=", m))'],
    }]));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const description = tree.nodeTree.getNode([{idx: 1}]).getItem().nodeDescription;
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      expectDeepEqual(description.getState('title'), 'm=1', {prefix: 'on'});
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', -1);
      flush();
      expectDeepEqual(description.getState('title'), undefined, {prefix: 'off'});
    });
  });

  test('Describe names base-expanded steps by $position', async () => {
    const pconf = await getProcessedConfig(daysPipeline({
      id: 'r', type: 'rule', base: 'base:expand(day)', from: 'v:same(@base, day)/a', to: 'title:same(@base, day)/title',
      effects: [{effect: 'describe', targets: 'title', value: {cat: ['Day ', {'+': [{var: '$position'}, 1]}]}}],
    }, 3));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const days = [0, 1, 2].map((idx) => tree.nodeTree.getNode([{idx: 0}, {idx}]).getItem());
      days.forEach((day) => day.getStateStore().setState('a', 1));
      flush();
      expectDeepEqual(days.map((day) => day.nodeDescription.getState('title')), ['Day 1', 'Day 2', 'Day 3']);
    });
  });

  const hidesAction = (mode: TreeMode) => async () => {
    const pconf = await getProcessedConfig(withAction({
      id: 'r', type: 'rule', from: 'm:step1/a', actions: 'act', when: 'gt(m, 0)',
      effects: ['hideAction(act, "addDay")'],
    }));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, ...mode});
      tree.init().subscribe();
      flush();
      const uuid = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid)![0].uuid;
      const input = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      input.setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), false, {prefix: 'on'});
      input.setState('a', -1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), true, {prefix: 'off'});
    });
  };

  test('Rule hides an action and shows it again when off', hidesAction({}));

  test('Rule hides an action and shows it again when off in Driver mode', hidesAction(driverMode));

  test('showAction shows an action only while the rule holds', async () => {
    const pconf = await getProcessedConfig(withAction({
      id: 'r', type: 'rule', from: 'm:step1/a', actions: 'act', when: 'gt(m, 0)',
      effects: ['showAction(act, "addDay")'],
    }));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const uuid = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid)![0].uuid;
      const input = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      input.setState('a', -1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), false, {prefix: 'off'});
      input.setState('a', 1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), true, {prefix: 'on'});
    });
  });

  test('A rule that is off shows an action hidden by hideWhen', async () => {
    const pconf = await getProcessedConfig(withAction({
      id: 'r', type: 'rule', from: 'm:step1/a', actions: 'act', when: 'gt(m, 0)',
      effects: [{effect: 'hideAction', targets: 'act', action: 'addDay'}],
    }, {hideWhen: 'guard: step1'}));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const uuid = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid)![0].uuid;
      expectDeepEqual(tree.linksState.isActionVisible(uuid), false, {prefix: 'hidden by hideWhen'});
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', -1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), true, {prefix: 'shown by the rule'});
    });
  });

  test('Two describe effects on one field do not erase each other', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 'title:step2/title',
      effects: ['describe(title, "Hot", when: gt(m, 30))', 'describe(title, "Cold", when: lte(m, 30))'],
    }]));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const description = tree.nodeTree.getNode([{idx: 1}]).getItem().nodeDescription;
      const input = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      input.setState('a', 40);
      flush();
      expectDeepEqual(description.getState('title'), 'Hot', {prefix: 'hot'});
      input.setState('a', 20);
      flush();
      expectDeepEqual(description.getState('title'), 'Cold', {prefix: 'cold'});
    });
  });

  const sourceAndDays = (rule: any, rootRule = true): PipelineConfiguration => ({
    id: 'pipeline1',
    type: 'static',
    steps: [
      {id: 'src', nqName: 'LibTests:TestAdd2'},
      {
        id: 'days', type: 'parallel', stepTypes: [{id: 'day', nqName: 'LibTests:TestAdd2'}],
        initialSteps: [{id: 'day'}], ...(rootRule ? {} : {links: [rule]}),
      },
    ],
    ...(rootRule ? {links: [rule]} : {}),
  });

  test('A separate data effect not reading the count leaves user values alone after tree changes', async () => {
    const pconf = await getProcessedConfig(sourceAndDays({
      id: 'r', type: 'rule', from: ['m:src/a', 'n:days/count()'], to: ['t:src/b', 'title:title'],
      effects: ['set(t, m, separate: true)', 'describe(title, cat("n=", n))'],
    }));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const src = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      src.setState('a', 1);
      flush();
      src.setState('b', 5);
      flush();
      tree.addSubTree(tree.nodeTree.getNode([{idx: 1}]).getItem().uuid, 'day', 1).subscribe();
      flush();
      expectDeepEqual(src.getState('b'), 5, {prefix: 'user value kept'});
    });
  });

  test('A data effect reading a count runs after tree changes, not at init', async () => {
    const pconf = await getProcessedConfig(sourceAndDays({
      id: 'r', type: 'rule', from: 'n:count()', to: 't:first(day)/b', effects: ['set(t, n)'],
    }, false));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const day = tree.nodeTree.getNode([{idx: 1}, {idx: 0}]).getItem().getStateStore();
      expectDeepEqual(day.getState('b'), undefined, {prefix: 'not at init'});
      tree.addSubTree(tree.nodeTree.getNode([{idx: 1}]).getItem().uuid, 'day', 1).subscribe();
      flush();
      expectDeepEqual(day.getState('b'), 2, {prefix: 'after adding a day'});
    });
  });

  test('hideWhen that starts to hide an action wins over a rule that showed it before', async () => {
    const pconf = await getProcessedConfig({
      ...sourceAndDays({id: 'r', type: 'rule', from: 'm:src/a', actions: 'act', when: 'gt(m, 0)',
        effects: ['hideAction(act, "addDay")']}),
      actions: [{
        id: 'addDay', position: 'none', from: [], to: 'out:src/b', hideWhen: 'guard: days/first(day)', handler() {},
      }],
    } as PipelineConfiguration);
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const uuid = tree.linksState.getNodeActionsData(tree.nodeTree.root.getItem().uuid)![0].uuid;
      tree.removeSubtree(tree.nodeTree.getNode([{idx: 1}, {idx: 0}]).getItem().uuid).subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', -1);
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), true, {prefix: 'shown by the off rule'});
      tree.addSubTree(tree.nodeTree.getNode([{idx: 1}]).getItem().uuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(tree.linksState.isActionVisible(uuid), false, {prefix: 'hidden by hideWhen again'});
    });
  });

  const validatesSteps = (mode: TreeMode) => async () => {
    const pconf = await getProcessedConfig(daysPipeline({
      id: 'r', type: 'rule', from: 'n:count()', to: 'self(node)', when: 'eq(n, 0)',
      effects: ['error(self, "Add at least one day")'],
    }));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, ...mode});
      tree.init().subscribe();
      flush();
      const days = tree.nodeTree.getNode([{idx: 0}]).getItem() as PipelineNodeBase;
      expectDeepEqual(hasErrors(days), true, {prefix: 'no days'});
      tree.addSubTree(days.uuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(hasErrors(days), false, {prefix: 'one day'});
      tree.removeSubtree(tree.nodeTree.getNode([{idx: 0}, {idx: 0}]).getItem().uuid).subscribe();
      flush();
      expectDeepEqual(hasErrors(days), true, {prefix: 'removed again'});
    });
  };

  // tree mode runs meta links on input changes only, the Driver also at init
  const describesSteps = (mode: TreeMode, initTitle?: string) => async () => {
    const pconf = await getProcessedConfig(daysPipeline({
      id: 'r', type: 'rule', from: 'n:count()', to: 'title:title', effects: ['describe(title, cat("Days: ", n))'],
    }, 1));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, ...mode});
      tree.init().subscribe();
      flush();
      const days = tree.nodeTree.getNode([{idx: 0}]).getItem();
      expectDeepEqual(days.nodeDescription.getState('title'), initTitle, {prefix: 'init'});
      tree.addSubTree(days.uuid, 'day', 1).subscribe();
      flush();
      expectDeepEqual(days.nodeDescription.getState('title'), 'Days: 2', {prefix: 'after adding a day'});
    });
  };

  test('A describe reading a count follows every tree change', describesSteps({}));

  test('A describe reading a count follows every tree change in Driver mode', describesSteps(driverMode, 'Days: 1'));

  test('Describe on tags writes this rule\'s tags', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 'tags:step2/tags', effects: ['describe(tags, ["fit", "pk"])'],
    }]));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore().setState('a', 1);
      flush();
      const description = tree.nodeTree.getNode([{idx: 1}]).getItem().nodeDescription;
      const tags = description.getState<Record<string, string[]>>('tags');
      expectDeepEqual(Object.values(tags ?? {}), [['fit', 'pk']]);
    });
  });

  test('$position follows inserted and removed base steps', async () => {
    const pconf = await getProcessedConfig(daysPipeline({
      id: 'r', type: 'rule', base: 'base:expand(day)', from: 'v:same(@base, day)/a', to: 'title:same(@base, day)/title',
      effects: [{effect: 'describe', targets: 'title', value: {cat: ['Day ', {'+': [{var: '$position'}, 1]}]}}],
    }, 2));
    testScheduler.run(({flush}) => {
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      flush();
      const daysNode = tree.nodeTree.getNode([{idx: 0}]);
      const titles = () => daysNode.getChildren().map(({item}) => item.getItem().nodeDescription.getState('title'));
      daysNode.getChildren().forEach(({item}) => item.getItem().getStateStore().setState('a', 1));
      flush();
      expectDeepEqual(titles(), ['Day 1', 'Day 2'], {prefix: 'initial'});
      tree.addSubTree(daysNode.getItem().uuid, 'day', 0).subscribe();
      flush();
      expectDeepEqual(titles(), ['Day 1', 'Day 2', 'Day 3'], {prefix: 'after inserting at the start'});
      tree.removeSubtree(daysNode.getChildren()[0].item.getItem().uuid).subscribe();
      flush();
      expectDeepEqual(titles(), ['Day 1', 'Day 2'], {prefix: 'after removing the first'});
    });
  });

  test('A count validates the workflow after every tree change', validatesSteps({}));

  test('A count validates the workflow after every tree change in Driver mode', validatesSteps(driverMode));
});
