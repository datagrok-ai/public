import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
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
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
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
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {x: {js: {args: ['x'], fn: () => 1}}}, effects: [{effect: 'error', targets: 't', message: 'm'}],
    }])), /collides with an input alias/);
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

  test('Rules validate func sources', async () => {
    const badRule = (rule: any) => expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      effects: [{effect: 'error', targets: 't', message: 'm'}], ...rule,
    }])));
    await badRule({sources: {v: {func: {args: {a: {var: 'x'}}}}}});
    await badRule({sources: {v: {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'nope'}}}}}});
    await badRule({sources: {v: {func: {name: 'LibTests:TestAdd2', args: [{var: 'x'}]}}}});
    await badRule({sources: {v: {query: {sql: 'select 1'}}}});
    await badRule({sources: {v: {query: {connection: 'System:Datagrok', sql: 'select 1', args: {a: {var: 'nope'}}}}}});
    const pconf: any = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {v: {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'x'}, b: 5}}}},
      effects: [{effect: 'error', targets: 't', message: {var: 'v'}}],
    }]));
    expectDeepEqual(pconf.links[0].params.sources.v.func.args, {a: {var: 'x'}, b: 5});
  });

  test('A func source calls a platform function', async () => {
    const controller = (values: Record<string, any>) => ({getFirst: (name: string) => values[name]}) as any;
    const sum = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'x'}, b: 5}}};
    const pending = resolveSources(controller({}), {sum}, {$all: {}, x: 1});
    expect(pending instanceof Promise, true);
    expectDeepEqual(await pending, {sum: 6});
    const presets = {func: {name: 'LibTests:TestPresets'}};
    const {presets: df} = await resolveSources(controller({}), {presets}) as any;
    expectDeepEqual(df.col('preset').toList(), ['fast', 'exact']);
  });

  test('A query source runs sql on a connection', async () => {
    const controller = () => ({getFirst: () => undefined}) as any;
    const me = await grok.dapi.users.current();
    const users = {query: {
      connection: 'System:Datagrok',
      sql: 'select login from users where login = @login',
      args: {login: {var: 'login'}},
    }};
    const {users: df} = await resolveSources(controller(), {users}, {$all: {}, login: me.login}) as any;
    expectDeepEqual(df.col('login').toList(), [me.login]);
    const declared = {query: {
      connection: 'System:Datagrok',
      sql: '--input: string login\nselect login from users where login = @login',
      args: {login: {var: 'login'}},
    }};
    const {declared: df2} = await resolveSources(controller(), {declared}, {$all: {}, login: me.login}) as any;
    expectDeepEqual(df2.rowCount, 1);
  });

  test('A file source loads a table', async () => {
    const controller = () => ({getFirst: () => undefined}) as any;
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {v: {file: ''}}, effects: [{effect: 'error', targets: 't', message: 'm'}],
    }])));
    const births = {file: 'System:DemoFiles/births.csv'};
    const {births: df} = await resolveSources(controller(), {births}) as any;
    expect(df.rowCount > 0, true);
  });

  test('A table source takes a dataframe or CSV text', async () => {
    const badRule = (table: any) => expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {v: {table}}, effects: [{effect: 'error', targets: 't', message: 'm'}],
    }])));
    await badRule('');
    await badRule(5);
    await badRule({csv: ''});
    await badRule({options: {delimiter: ';'}});
    const df = DG.DataFrame.fromCsv('preset\nfast\nexact');
    const cache = new Map<string, any>();
    const controller = () => ({getFirst: () => undefined, sourceCache: cache}) as any;
    const sources = {
      given: {table: df},
      text: {table: 'preset,v\nfast,1\nexact,2'},
      parsed: {table: {csv: 'preset;v\nfast;1', options: {delimiter: ';'}}},
    };
    const first: any = resolveSources(controller(), sources);
    expect(first instanceof Promise, false);
    expect(first.given === df, true);
    expectDeepEqual(first.text.col('v').toList(), [1, 2]);
    expectDeepEqual(first.parsed.columns.names(), ['preset', 'v']);
    const second: any = resolveSources(controller(), sources);
    expect(second.text === first.text, true);
    expect(second.parsed === first.parsed, true);
  });

  test('A table source feeds a rule in the tree', async () => {
    const df = DG.DataFrame.fromCsv('preset\nfast\nexact');
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      sources: {presets: {table: df}},
      effects: [{effect: 'items', targets: 't', items: {column: [{var: 'presets'}, 'preset']}}],
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
    });
    expectDeepEqual(metas, [{items: ['fast', 'exact']}]);
  });

  test('Sources that read no input resolve once per link', async () => {
    const original = grok.functions.call;
    let calls = 0;
    (grok.functions as any).call = (name: string, params: any) => {
      calls++;
      return original.call(grok.functions, name, params);
    };
    try {
      const cache = new Map<string, any>();
      const controller = () => ({getFirst: () => undefined, sourceCache: cache}) as any;
      const fixed = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: 2}}};
      const varying = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'x'}, b: 2}}};
      expectDeepEqual(await resolveSources(controller(), {fixed, varying}, {$all: {}, x: 1}), {fixed: 3, varying: 3});
      expectDeepEqual(await resolveSources(controller(), {fixed, varying}, {$all: {}, x: 5}), {fixed: 3, varying: 7});
      expect(calls, 3);
      expect(cache.has('fixed'), true);
      expect(cache.has('varying'), false);
    } finally {
      (grok.functions as any).call = original;
    }
  });

  test('A func source without args runs on every run', async () => {
    const original = grok.functions.call;
    let calls = 0;
    (grok.functions as any).call = (name: string, params: any) => {
      calls++;
      return original.call(grok.functions, name, params);
    };
    try {
      const cache = new Map<string, any>();
      const controller = () => ({getFirst: () => undefined, sourceCache: cache}) as any;
      const presets = {func: {name: 'LibTests:TestPresets'}};
      await resolveSources(controller(), {presets});
      await resolveSources(controller(), {presets});
      expect(calls, 2);
      expect(cache.has('presets'), false);
    } finally {
      (grok.functions as any).call = original;
    }
  });

  test('A failed source is not cached and is retried', async () => {
    const original = grok.functions.call;
    let fail = true;
    (grok.functions as any).call = (name: string, params: any) => {
      if (fail) {
        fail = false;
        return Promise.reject(new Error('Source down'));
      }
      return original.call(grok.functions, name, params);
    };
    try {
      const cache = new Map<string, any>();
      const controller = () => ({getFirst: () => undefined, sourceCache: cache}) as any;
      const fixed = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: 2}}};
      await expectThrowsAsync(async () => {
        await resolveSources(controller(), {fixed});
      }, /Source down/);
      expect(cache.has('fixed'), false);
      expectDeepEqual(await resolveSources(controller(), {fixed}), {fixed: 3});
      expect(cache.has('fixed'), true);
    } finally {
      (grok.functions as any).call = original;
    }
  });

  test('A failing source applies no effects and the next run recovers', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      sources: {list: {js: {args: ['m'], fn: (m: number) => {
        if (m < 0)
          throw new Error('Negative input');
        return [String(m)];
      }}}},
      effects: [{effect: 'items', targets: 't', items: {var: 'list'}}],
    }]));
    const logger = new DriverLogger();
    const metas: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, logger});
      tree.init().subscribe();
      const store = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => store.setState('a', 1));
      cold('--a').subscribe(() => metas.push(outBridge.meta.a.value));
      cold('---a').subscribe(() => store.setState('a', -1));
      cold('----a').subscribe(() => metas.push(outBridge.meta.a.value));
      cold('-----a').subscribe(() => store.setState('a', 2));
      cold('------a').subscribe(() => metas.push(outBridge.meta.a.value));
    });
    expectDeepEqual(metas, [{items: ['1']}, {items: ['1']}, {items: ['2']}]);
    expectDeepEqual(logger.errors.map((error) => [error.context, error.message]), [['link:r::meta', 'Negative input']]);
  });

  test('Sources run while the rule is off', async () => {
    let calls = 0;
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      when: {'>': [{var: 'm'}, 0]},
      sources: {list: {js: {args: ['m'], fn: () => {
        calls++;
        return ['x'];
      }}}},
      effects: [{effect: 'items', targets: 't', items: {var: 'list'}}],
    }]));
    const metas: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const store = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => store.setState('a', -1));
      cold('--a').subscribe(() => store.setState('a', -2));
      cold('---a').subscribe(() => metas.push(outBridge.meta.a.value));
    });
    expectDeepEqual(metas, [{}]);
    expect(calls, 2);
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
