import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {resolveSources} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallInstancesBridge} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrowsAsync} from '../../../test-utils';

const twoSteps = (links: PipelineLinkConfigurationInput<string | string[]>[]): PipelineConfiguration => ({
  id: 'pipeline1',
  type: 'static',
  steps: [
    {id: 'step1', nqName: 'LibTests:TestAdd2'},
    {id: 'step2', nqName: 'LibTests:TestMul2'},
  ],
  links,
});

category('ComputeUtils: Driver rule js sources', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Rules validate js sources', async () => {
    const badRule = (rule: any) => expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      effects: [{effect: 'error', targets: 't', message: 'm'}], ...rule,
    }])));
    await badRule({sources: {v: {js: {args: ['nope'], fn: () => 1}}}});
    await badRule({sources: {v: {js: {args: 'x', fn: () => 1}}}});
    await badRule({sources: {v: {js: {args: ['x'], fn: 1}}}});
    await badRule({sources: {v: {js: {args: ['x'], fn: () => []}}},
      effects: [{effect: 'verdicts', targets: 't', source: 'v'}]});
    const pconf: any = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {v: {js: {args: ['x'], fn: (x: number) => x * 2}}},
      effects: [{effect: 'error', targets: 't', message: {var: 'v'}}],
    }]));
    expectDeepEqual(pconf.links.map((link: any) => link.id), ['r::validator']);
    expectDeepEqual(pconf.links[0].params.sources.v.js.args, ['x']);
  });

  test('A js source feeds items and meta on every run', async () => {
    let calls = 0;
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: ['m:step1/a', 'other:step1/b'],
      to: 't:step2/a',
      sources: {list: {js: {args: ['m'], fn: (m: number) => {
        calls++;
        return m > 0 ? ['x', 'y'] : ['z'];
      }}}},
      effects: [
        {effect: 'items', targets: 't', items: {var: 'list'}},
        {effect: 'meta', targets: 't', meta: {count: {len: {var: 'list'}}}},
      ],
    }]));
    const metas: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const store = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => store.setState('a', 1));
      cold('--a').subscribe(() => metas.push(outBridge.meta.a.value));
      cold('---a').subscribe(() => store.setState('b', 5));
      cold('----a').subscribe(() => metas.push(outBridge.meta.a.value));
      cold('-----a').subscribe(() => store.setState('a', -1));
      cold('------a').subscribe(() => metas.push(outBridge.meta.a.value));
    });
    expectDeepEqual(metas, [
      {items: ['x', 'y'], count: 2},
      {items: ['x', 'y'], count: 2},
      {items: ['z'], count: 1},
    ]);
    expect(calls, 3);
  });

  test('A js source may be async', async () => {
    let calls = 0;
    const sum = {js: {args: ['x', 'y'], fn: async (x: number, y: number) => {
      calls++;
      return x + y;
    }}};
    const controller = (values: Record<string, any>) => ({getFirst: (name: string) => values[name]}) as any;
    expectDeepEqual(await resolveSources(controller({x: 1, y: 2}), {sum}), {sum: 3});
    expect(resolveSources(controller({x: 1, y: 2}), {sum}) instanceof Promise, true);
    expectDeepEqual(await resolveSources(controller({x: 2, y: 2}), {sum}), {sum: 4});
    expect(calls, 3);
    const list = {js: {args: ['x'], fn: (x: number) => [x]}};
    expectDeepEqual(resolveSources(controller({x: 7}), {list}), {list: [7]});
    expectDeepEqual(resolveSources(controller({x: undefined}), {list}), {list: [undefined]});
  });
});
