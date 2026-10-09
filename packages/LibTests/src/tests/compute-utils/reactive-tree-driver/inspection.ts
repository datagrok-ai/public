import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {
  getProcessedConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {
  inspectConfig, inspectLinks, toInspectorJSON,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/inspection';
import {
  isLinkLogItem, LOG_EVENT_TYPES, LogItem,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
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

category('ComputeUtils: Driver inspection', async () => {
  test('A base link lists the base it matched', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'sequential',
      stepTypes: steps,
      initialSteps: ['step1', 'step2', 'step1'],
      links: [{id: 'l1', base: 'base:expand(step1)', from: 'in:same(@base)/res', to: 'out:after+(@base, step2)/a'}],
    } as PipelineConfiguration);
    expectDeepEqual(inspectLinks(tree).map((l) => [l.id, l.basePath, l.resolvedOutputs]),
      [['l1', 'step1[0]', {out: ['step2[1]::a']}]]);
  });

  test('An optional alias without matches is shown as empty', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps, links: [
      {id: 'l1', from: ['in:step1/res', 'extra(optional):step3/a'], to: 'out:step2/a'},
    ]});
    const [link] = inspectLinks(tree);
    expectDeepEqual(link.resolvedInputs, {in: ['step1[0]::res']});
    expectDeepEqual(link.emptyAliases, {extra: 'optional'});
  });

  test('Targets dropped by linked are shown', async () => {
    const tree = await makeTree({id: 'root', type: 'static', steps, links: [
      {id: 'l1', from: 'in:step1/res', to: 'out:step2/a'},
      {
        id: 'l2', type: 'meta', from: 'in:step1/res',
        to: '_(template):step2/inputs(LibTests:TestMul2, $linked)', handler() {},
      },
    ]});
    const l2 = inspectLinks(tree).find((l) => l.id === 'l2')!;
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
    const entries = inspectLinks(tree).map((l) => [l.id, l.node.configId, l.type, !!l.hasHandler]);
    expectDeepEqual(entries, [['l', 'p1', 'data', false], ['l', 'p2', 'validator', true]]);
  });

  test('Actions show whether they are visible', async () => {
    const tree = await makeTree({
      id: 'root',
      type: 'static',
      steps,
      links: [{id: 'l1', from: 'in:step1/res', to: 'out:step2/a'}],
      actions: [
        {id: 'shown', from: 'in:step1/res', to: 'out:step2/a', position: 'buttons', handler() {}},
        {id: 'hidden', from: 'in:step1/res', to: 'out:step2/a', position: 'buttons', hideWhen: 'h:step1', handler() {}},
      ],
    } as PipelineConfiguration);
    expectDeepEqual(inspectLinks(tree).map((l) => [l.id, l.isAction, l.visible]),
      [['l1', false, undefined], ['shown', true, true], ['hidden', true, false]]);
  });

  test('Links show the link, rule or action as written', async () => {
    const link = {id: 'l1', from: 'in:step1/res', to: 'out:step2/a'};
    const rule = {
      id: 'r', type: 'rule', from: 'm:step1/res', to: ['t:step2/a', 'u:step2/b'],
      effects: [{effect: 'hide', targets: 't'}, {effect: 'set', targets: 'u', value: {var: 'm'}}],
    };
    const action = {id: 'act', from: 'in:step1/res', to: 'out:step2/a', position: 'none', handler() {}};
    const written = JSON.parse(JSON.stringify([link, rule, {...action, handler: '#Handler'}]));
    const tree = await makeTree(
      {id: 'root', type: 'static', steps, links: [link, rule], actions: [action]} as PipelineConfiguration);
    expectDeepEqual(inspectLinks(tree).map((l) => [l.id, l.original]), [
      ['l1', written[0]], ['r::meta', written[1]], ['r::data', written[1]], ['act', written[2]],
    ]);
  });

  test('Inspector JSON makes values readable', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', [1, 2])]);
    const data = toInspectorJSON(
      {df, map: new Map([['k', 1]]), set: new Set([1]), subj: new BehaviorSubject(5), fn: () => 1});
    expectDeepEqual(data, {df: '#DataFrame(2 rows, 1 cols)', map: {k: 1}, set: [1], subj: 5, fn: '#Handler'});
    expectDeepEqual(toInspectorJSON(undefined), {});
  });

  test('Config is shown as written with resolved step io', async () => {
    const config = await getProcessedConfig({
      id: 'root',
      type: 'static',
      steps: [
        steps[0],
        {id: 'seq', type: 'sequential', stepTypes: [{...steps[1], disableUIControlls: true}], initialSteps: []},
      ],
      links: [{id: 'l1', from: 'in:step1/res', to: 'out:seq/step2/a', handler() {}}],
    } as PipelineConfiguration);
    const inputs = [{name: 'a', type: 'double'}, {name: 'b', type: 'double'}];
    const outputs = [{name: 'res', type: 'double'}];
    expectDeepEqual(inspectConfig(config), {
      id: 'root',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestAdd2', inputs, outputs},
        {id: 'seq', type: 'sequential', initialSteps: [], stepTypes: [
          {id: 'step2', nqName: 'LibTests:TestMul2', disableUIControlls: true, inputs, outputs},
        ]},
      ],
      links: [{id: 'l1', from: 'in:step1/res', to: 'out:seq/step2/a', handler: '#Handler'}],
    });
  });

  test('Link log items are recognized', async () => {
    const base = {uuid: '', timestamp: new Date()};
    const link = {...base, type: 'actionAdded', linkUUID: '', prefix: [], id: 'a'} as LogItem;
    expectDeepEqual([isLinkLogItem(link), isLinkLogItem({...base, type: 'treeUpdateStarted'})], [true, false]);
    expectDeepEqual(LOG_EVENT_TYPES.filter((t) => isLinkLogItem({...base, type: t.type} as LogItem)).length, 6);
  });

  test('The driver reads links and the config of the current tree', async () => {
    const driver = new Driver(true);
    const config = await getProcessedConfig({id: 'root', type: 'static', steps, links: [
      {id: 'l1', from: 'in:step1/res', to: 'out:step2/a'},
      {id: 'l2', from: 'in:stpe1/res', to: 'out:step2/a'},
    ]});
    expectDeepEqual([driver.inspectLinks(), driver.inspectConfig()], [[], undefined], {prefix: 'Before init'});
    await driver.sendCommand({event: 'initPipeline', provider: '', config});
    expectDeepEqual(driver.inspectLinks().map((l) => l.id), ['l1']);
    expectDeepEqual(driver.inspectConfig().links.map((l: any) => l.from), ['in:step1/res', 'in:stpe1/res']);
    driver.close();
  });
});
