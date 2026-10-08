import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {
  getProcessedConfig, PipelineConfigurationStaticProcessed,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {
  formatLinkIO, formatLinkSegment, parseLinkIO,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/LinkSpec';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {
  inspectConfig, inspectLinks, inspectNode, toInspectorJSON,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/inspection';
import {
  isLinkLogItem, LOG_EVENT_TYPES, LogItem,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {NodePath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/BaseTree';
import {
  explainLinkMatch, matchNodeLink,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/link-matching';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {BehaviorSubject} from 'rxjs';

const steps = [
  {id: 'step1', nqName: 'LibTests:TestAdd2'},
  {id: 'step2', nqName: 'LibTests:TestMul2'},
];

async function makeTree(config: PipelineConfiguration) {
  const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config), mockMode: true});
  await tree.init().toPromise();
  return tree;
}

async function notMatched(links: any[], extra: Partial<PipelineConfiguration> = {}) {
  const tree = await makeTree({id: 'root', type: 'static', steps, links, ...extra} as PipelineConfiguration);
  return inspectLinks(tree).notMatched;
}

category('ComputeUtils: Driver inspection', async () => {
  test('A step path typo reaches no nodes', async () => {
    const [entry, ...rest] = await notMatched([{id: 'l1', from: 'in:stpe1/res', to: 'out:step2/a'}]);
    expectDeepEqual(rest.length, 0, {prefix: 'Count'});
    expectDeepEqual([entry.id, entry.isAction, entry.node.configId], ['l1', false, 'root']);
    expectDeepEqual(entry.explanation.failedAlias, 'in');
    expectDeepEqual(entry.explanation.aliases[0], {name: 'in', kind: 'input', optional: false, nodes: 0, ios: 0});
  });

  test('An io name typo reaches the step but not the io', async () => {
    const [entry] = await notMatched([{id: 'l1', from: 'in:step1/rez', to: 'out:step2/a'}]);
    expectDeepEqual(entry.explanation.failedAlias, 'in');
    expectDeepEqual(entry.explanation.aliases[0], {name: 'in', kind: 'input', optional: false, nodes: 1, ios: 0});
    expectDeepEqual(entry.explanation.aliases[1], {name: 'out', kind: 'output', optional: false, nodes: 1, ios: 1});
  });

  test('An alias referring to an empty alias depends on it', async () => {
    const [entry] = await notMatched([{id: 'l1', from: 'in:stpe1/res', to: 'out:after(@in, step2)/a'}]);
    expectDeepEqual(entry.explanation.failedAlias, 'in');
    expectDeepEqual(entry.explanation.aliases[1],
      {name: 'out', kind: 'output', optional: false, nodes: 0, dependsOn: 'in'});
  });

  test('A matching not clause blocks the link', async () => {
    const [entry] = await notMatched([{id: 'l1', not: 'blk:step1', from: 'in:step1/res', to: 'out:step2/a'}]);
    expectDeepEqual(entry.explanation.blockedBy, 'blk');
    expectDeepEqual(entry.explanation.failedAlias, undefined);
    expectDeepEqual(entry.explanation.aliases[0], {name: 'blk', kind: 'not', optional: false, nodes: 1});
  });

  test('An unmatched action is listed as an action', async () => {
    const [entry] = await notMatched([], {actions: [{
      id: 'act1', from: 'in:step1/res', to: 'out:stpe2/a', position: 'none', handler() {},
    }]} as any);
    expectDeepEqual([entry.id, entry.isAction, entry.explanation.failedAlias], ['act1', true, 'out']);
  });

  test('An empty dynamic workflow leaves the base empty', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [{
        id: 'seq',
        type: 'sequential',
        stepTypes: steps,
        initialSteps: [],
        links: [{id: 'l1', base: 'base:expand(step1)', from: 'in:same(@base)/res', to: 'out:same(@base)/a'}],
      }],
    });
    const [entry] = inspectLinks(tree).notMatched;
    expectDeepEqual([entry.id, entry.node.configId, entry.node.path], ['l1', 'seq', 'seq[0]']);
    expectDeepEqual(entry.explanation,
      {aliases: [{name: 'base', kind: 'base', optional: false, nodes: 0}], failedAlias: 'base'});
  });

  test('A base link matched at some nodes is not listed as not matched', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'sequential',
      stepTypes: steps,
      initialSteps: ['step1', 'step2', 'step1'],
      links: [{id: 'l1', base: 'base:expand(step1)', from: 'in:same(@base)/res', to: 'out:after+(@base, step2)/a'}],
    } as PipelineConfiguration);
    const {matched, notMatched} = inspectLinks(tree);
    expectDeepEqual(notMatched, []);
    expectDeepEqual(matched.map((l) => [l.id, l.basePath, l.resolvedOutputs]),
      [['l1', 'step1[0]', {out: ['step2[1]::a']}]]);
  });

  test('An optional alias without matches is shown as empty', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps, links: [
      {id: 'l1', from: ['in:step1/res', 'extra(optional):step3/a'], to: 'out:step2/a'},
    ]});
    const {matched, notMatched} = inspectLinks(tree);
    expectDeepEqual(notMatched, []);
    expectDeepEqual(matched[0].resolvedInputs, {in: ['step1[0]::res']});
    expectDeepEqual(matched[0].emptyAliases, {extra: 'optional'});
  });

  test('Targets dropped by linked are shown', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps, links: [
      {id: 'l1', from: 'in:step1/res', to: 'out:step2/a'},
      {
        id: 'l2', type: 'meta', from: 'in:step1/res',
        to: '_(template):step2/inputs(LibTests:TestMul2, $linked)', handler() {},
      },
    ]});
    const l2 = inspectLinks(tree).matched.find((l) => l.id === 'l2')!;
    expectDeepEqual(l2.resolvedOutputs, {b: ['step2[1]::b']});
    expectDeepEqual(l2.emptyAliases, {a: 'linked'});
    expectDeepEqual(l2.linkedTargets, {a: ['step2[1]::a']});
  });

  test('Links with the same id in two steps keep their own spec', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps: [
        {id: 'p1', type: 'static', steps: [steps[0]], links: [{id: 'l', from: 'in:step1/res', to: 'out:step1/a'}]},
        {id: 'p2', type: 'static', steps: [steps[1]], links: [
          {id: 'l', type: 'validator', from: 'in:step2/res', to: 'out:step2/a', handler() {}},
        ]},
      ],
    } as PipelineConfiguration);
    const entries = inspectLinks(tree).matched.map((l) => [l.id, l.node.configId, l.type, !!l.hasHandler]);
    expectDeepEqual(entries, [['l', 'p1', 'data', false], ['l', 'p2', 'validator', true]]);
  });

  test('Explaining agrees with matching', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps, links: [
      {id: 'ok', from: 'in:step1/res', to: 'out:step2/a'},
      {id: 'pathTypo', from: 'in:stpe1/res', to: 'out:step2/a'},
      {id: 'ioTypo', from: 'in:step1/res', to: 'out:step2/z'},
      {id: 'optional', from: ['in:step1/res', 'extra(optional):step3/a'], to: 'out:step2/a'},
      {id: 'dependsOn', from: 'in:stpe1/res', to: 'out:after(@in, step2)/a'},
      {id: 'blocked', not: 'blk:step1', from: 'in:step1/res', to: 'out:step2/a'},
      {id: 'notClear', not: 'blk:step3', from: 'in:step1/res', to: 'out:step2/a'},
      {id: 'call', from: 'fc(call,optional):step1', to: 'out:step2/a'},
      {id: 'pipeline', type: 'pipeline', from: [], to: 'out:stpe2', handler() {}},
    ]} as PipelineConfiguration);
    const root = tree.nodeTree.root;
    const specs = (root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    const results = specs.map((spec) => {
      const {failedAlias, blockedBy} = explainLinkMatch(root, spec);
      return [spec.id, !!matchNodeLink(root, spec), failedAlias ?? blockedBy ?? null];
    });
    expectDeepEqual(results, [
      ['ok', true, null], ['pathTypo', false, 'in'], ['ioTypo', false, 'out'], ['optional', true, null],
      ['dependsOn', false, 'in'], ['blocked', false, 'blk'], ['notClear', true, null], ['call', true, null],
      ['pipeline', false, 'out'],
    ]);
    const aliases = (id: string) => explainLinkMatch(root, specs.find((s) => s.id === id)!).aliases;
    expectDeepEqual(aliases('notClear'), [
      {name: 'blk', kind: 'not', optional: false, nodes: 0},
      {name: 'in', kind: 'input', optional: false, nodes: 1, ios: 1},
      {name: 'out', kind: 'output', optional: false, nodes: 1, ios: 1},
    ], {prefix: 'notClear'});
    expectDeepEqual(aliases('call'), [
      {name: 'fc', kind: 'input', optional: true, nodes: 1},
      {name: 'out', kind: 'output', optional: false, nodes: 1, ios: 1},
    ], {prefix: 'call'});
  });

  test('Explaining a base link agrees with matching', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'sequential',
      stepTypes: steps,
      initialSteps: ['step1', 'step2', 'step1'],
      links: [
        {id: 'some', base: 'base:expand(step1)', from: 'in:same(@base)/res', to: 'out:after+(@base, step2)/a'},
        {id: 'none', base: 'base:expand(step1)', from: 'in:same(@base)/res', to: 'out:after+(@base, step3)/a'},
      ],
    } as PipelineConfiguration);
    const root = tree.nodeTree.root;
    const [some, none] = (root.getItem().config as PipelineConfigurationStaticProcessed).links!;
    expectDeepEqual([matchNodeLink(root, some)?.length, explainLinkMatch(root, some).failedAlias], [1, undefined]);
    expectDeepEqual(matchNodeLink(root, none), undefined);
    expectDeepEqual(explainLinkMatch(root, none), {
      aliases: [
        {name: 'base', kind: 'base', optional: false, nodes: 2},
        {name: 'in', kind: 'input', optional: false, nodes: 1, ios: 1},
        {name: 'out', kind: 'output', optional: false, nodes: 0, ios: 0},
      ],
      failedAlias: 'out',
    });
  });

  test('A repeated alias throws in matching but not in explaining', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps});
    const config = await getProcessedConfig({id: 'root', type: 'static', steps, links: [
      {id: 'dup', from: ['x:step1/res', 'x:step2/res'], to: 'out:step2/a'},
      {id: 'dupAfterTypo', from: ['a:stpe1/res', 'x:step1/res', 'x:step2/res'], to: 'out:step2/a'},
    ]}) as PipelineConfigurationStaticProcessed;
    const [dup, dupAfterTypo] = config.links!;
    const root = tree.nodeTree.root;
    let error = '';
    try {
      matchNodeLink(root, dup);
    } catch (e) {
      error = (e as Error).message;
    }
    expectDeepEqual(error.startsWith('Duplicate io name x'), true, {prefix: 'Matching'});
    expectDeepEqual(matchNodeLink(root, dupAfterTypo), undefined);
    expectDeepEqual(explainLinkMatch(root, dupAfterTypo).failedAlias, 'a');
  });

  test('A step node shows its io and states', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps});
    const step = tree.nodeTree.getItem([{idx: 0}] as NodePath);
    const node = inspectNode(tree, step.uuid)!;
    expectDeepEqual([node.configId, node.path, node.type, node.isReadonly], ['step1', 'step1[0]', 'funccall', false]);
    expectDeepEqual(node.inputs?.map((io) => io.name), ['a', 'b']);
    expectDeepEqual(node.outputs?.map((io) => io.name), ['res']);
    expectDeepEqual(node.callState?.isOutputOutdated, true, {prefix: 'Outdated'});
    expectDeepEqual([node.validations, node.consistency, node.meta], [{}, {}, {}]);
  });

  test('A workflow node shows its states', async () => {
    const tree = await makeTree(
      {id: 'root', type: 'static', steps, states: [{id: 'counter'}]} as PipelineConfiguration);
    tree.nodeTree.root.getItem().getStateStore().editState('counter', 3);
    const node = inspectNode(tree, tree.nodeTree.root.getItem().uuid)!;
    expectDeepEqual([node.configId, node.path, node.type, node.states, node.pipelineValidations],
      ['root', '', 'static', {counter: 3}, {}]);
    expectDeepEqual(inspectNode(tree, 'missing'), undefined);
  });

  test('Link queries are printed with all their parts', async () => {
    const [io] = parseLinkIO('in(optional):before(@base, step1|step2, stop)/#after(@base, t1&t2)/res', 'input');
    expectDeepEqual(formatLinkIO(io),
      {name: 'in', flags: ['optional'], path: ['before(@base, step1|step2, stop)', '#after(@base,t1&t2)', 'res']});
    const [plain] = parseLinkIO('in:all(step1)/step2/res', 'input');
    expectDeepEqual(plain.segments.map(formatLinkSegment), ['all(step1)', 'step2', 'res']);
    expectDeepEqual(formatLinkIO({...plain, unlinked: true}).unlinked, true);
  });

  test('Inspector JSON makes values readable', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', [1, 2])]);
    const data = toInspectorJSON(
      {df, map: new Map([['k', 1]]), set: new Set([1]), subj: new BehaviorSubject(5), fn: () => 1});
    expectDeepEqual(data, {df: '#DataFrame(2 rows, 1 cols)', map: {k: 1}, set: [1], subj: 5, fn: '#Handler'});
    expectDeepEqual(toInspectorJSON(undefined), {});
  });

  test('Config io is split into inputs and outputs', async () => {
    const config = await getProcessedConfig(
      {id: 'root', type: 'static', steps, links: [{id: 'l1', from: 'in:step1/res', to: 'out:step2/a'}]});
    const data = inspectConfig(config);
    expectDeepEqual([data.steps[0].inputs, data.steps[0].outputs, data.steps[0].io],
      [[{name: 'a', type: 'double'}, {name: 'b', type: 'double'}], [{name: 'res', type: 'double'}], undefined]);
    expectDeepEqual(data.links[0].from, [{name: 'in', path: ['step1', 'res']}]);
  });

  test('Link log items are recognized', async () => {
    const base = {uuid: '', timestamp: new Date()};
    const link = {...base, type: 'actionAdded', linkUUID: '', prefix: [], id: 'a'} as LogItem;
    expectDeepEqual([isLinkLogItem(link), isLinkLogItem({...base, type: 'treeUpdateStarted'})], [true, false]);
    expectDeepEqual(LOG_EVENT_TYPES.filter((t) => isLinkLogItem({...base, type: t.type} as LogItem)).length, 6);
  });

  test('The driver reads links and nodes of the current tree', async () => {
    const driver = new Driver(true);
    const config = await getProcessedConfig({id: 'root', type: 'static', steps, links: [
      {id: 'l1', from: 'in:step1/res', to: 'out:step2/a'},
      {id: 'l2', from: 'in:stpe1/res', to: 'out:step2/a'},
    ]});
    expectDeepEqual(driver.inspectLinks(), {matched: [], notMatched: []}, {prefix: 'Before init'});
    const tree = await driver.sendCommand({event: 'initPipeline', provider: '', config}) as StateTree;
    const {matched, notMatched} = driver.inspectLinks();
    expectDeepEqual([matched.map((l) => l.id), notMatched.map((l) => l.id)], [['l1'], ['l2']]);
    expectDeepEqual(driver.inspectNode(tree.nodeTree.root.getItem().uuid)?.configId, 'root');
    driver.close();
  });
});
