import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {resolveSources} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-sources';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallInstancesBridge} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {DriverLogger} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrowsAsync, noCallController, twoSteps} from '../../../test-utils';

const valuesController = (values: Record<string, any>) => ({getFirst: (name: string) => values[name]}) as any;

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
    const sum = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'x'}, b: 5}}};
    const pending = resolveSources(valuesController({}), {sum}, {$all: {}, x: 1});
    expect(pending instanceof Promise, true);
    expectDeepEqual(await pending, {sum: 6});
    const presets = {func: {name: 'LibTests:TestPresets'}};
    const {presets: df} = await resolveSources(valuesController({}), {presets}) as any;
    expectDeepEqual(df.col('preset').toList(), ['fast', 'exact']);
  });

  test('A query source runs sql on a connection', async () => {
    const me = await grok.dapi.users.current();
    const users = {query: {
      connection: 'System:Datagrok',
      sql: 'select login from users where login = @login',
      args: {login: {var: 'login'}},
    }};
    const {users: df} = await resolveSources(noCallController(), {users}, {$all: {}, login: me.login}) as any;
    expectDeepEqual(df.col('login').toList(), [me.login]);
    const declared = {query: {
      connection: 'System:Datagrok',
      sql: '--input: string login\nselect login from users where login = @login',
      args: {login: {var: 'login'}},
    }};
    const {declared: df2} = await resolveSources(noCallController(), {declared}, {$all: {}, login: me.login}) as any;
    expectDeepEqual(df2.rowCount, 1);
  });

  test('A file source loads a table', async () => {
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'bad', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources: {v: {file: ''}}, effects: [{effect: 'error', targets: 't', message: 'm'}],
    }])));
    const compounds = {file: 'System:AppData/LibTests/compounds.csv'};
    const {compounds: df} = await resolveSources(noCallController(), {compounds}) as any;
    expectDeepEqual(df.columns.names(), ['compound', 'clearance', 'volume']);
    expectDeepEqual(df.col('compound').toList(), ['CMP-0001', 'CMP-0002']);
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
    const sources = {
      given: {table: df},
      text: {table: 'preset,v\nfast,1\nexact,2'},
      parsed: {table: {csv: 'preset;v\nfast;1', options: {delimiter: ';'}}},
    };
    const first: any = resolveSources(noCallController(cache), sources);
    expect(first instanceof Promise, false);
    expect(first.given === df, true);
    expectDeepEqual(first.text.col('v').toList(), [1, 2]);
    expectDeepEqual(first.parsed.columns.names(), ['preset', 'v']);
    const second: any = resolveSources(noCallController(cache), sources);
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
      const fixed = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: 2}}};
      const varying = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'x'}, b: 2}}};
      const controller = noCallController(cache);
      expectDeepEqual(await resolveSources(controller, {fixed, varying}, {$all: {}, x: 1}), {fixed: 3, varying: 3});
      expectDeepEqual(await resolveSources(controller, {fixed, varying}, {$all: {}, x: 5}), {fixed: 3, varying: 7});
      expect(calls, 3);
      expect(cache.has('fixed'), true);
      expect(cache.has('varying'), false);
    } finally {
      (grok.functions as any).call = original;
    }
  });

  test('A func source without args resolves once per link', async () => {
    const original = grok.functions.call;
    let calls = 0;
    (grok.functions as any).call = (name: string, params: any) => {
      calls++;
      return original.call(grok.functions, name, params);
    };
    try {
      const cache = new Map<string, any>();
      const presets = {func: {name: 'LibTests:TestPresets'}};
      await resolveSources(noCallController(cache), {presets});
      await resolveSources(noCallController(cache), {presets});
      expect(calls, 1);
      expect(cache.has('presets'), true);
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
      const fixed = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: 2}}};
      await expectThrowsAsync(async () => {
        await resolveSources(noCallController(cache), {fixed});
      }, /Source down/);
      expect(cache.has('fixed'), false);
      expectDeepEqual(await resolveSources(noCallController(cache), {fixed}), {fixed: 3});
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
    expectDeepEqual(await resolveSources(valuesController({x: 1, y: 2}), {sum}), {sum: 3});
    expect(resolveSources(valuesController({x: 1, y: 2}), {sum}) instanceof Promise, true);
    expectDeepEqual(await resolveSources(valuesController({x: 2, y: 2}), {sum}), {sum: 4});
    expect(calls, 3);
    const list = {js: {args: ['x'], fn: (x: number) => [x]}};
    expectDeepEqual(resolveSources(valuesController({x: 7}), {list}), {list: [7]});
    expectDeepEqual(resolveSources(valuesController({x: undefined}), {list}), {list: [undefined]});
  });
});

const countCalls = async (fn: (calls: () => number) => Promise<void>, fail = false) => {
  const original = grok.functions.call;
  let calls = 0;
  let failing = fail;
  (grok.functions as any).call = (name: string, params: any) => {
    calls++;
    if (failing) {
      failing = false;
      return Promise.reject(new Error('Source down'));
    }
    return original.call(grok.functions, name, params);
  };
  try {
    await fn(() => calls);
  } finally {
    (grok.functions as any).call = original;
  }
};

category('ComputeUtils: Driver rule source chains', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Sources read only inputs and the sources declared before them', async () => {
    const process = (sources: any) => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'x:step1/a', to: 't:step1/a',
      sources, effects: [{effect: 'error', targets: 't', message: 'm'}],
    }]));
    const later = /reads w, which is declared after it/;
    const fn = (w: number) => w;
    await expectThrowsAsync(() => process({v: {js: {args: ['w'], fn}}, w: {js: {args: ['x'], fn}}}), later);
    await expectThrowsAsync(() => process({v: {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'w'}}}},
      w: {js: {args: ['x'], fn}}}), later);
    await expectThrowsAsync(() => process({v: {query: {connection: 'System:Datagrok', sql: 'select 1',
      args: {a: {var: 'w'}}}}, w: {js: {args: ['x'], fn}}}), later);
    await expectThrowsAsync(() => process({v: 'func("LibTests:TestAdd2", a: w)', w: {js: {args: ['x'], fn}}}), later);
    await expectThrowsAsync(() => process({v: {js: {args: ['v'], fn}}}), /reads itself/);
    await expectThrowsAsync(() => process({v: {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'v'}}}}}),
      /reads itself/);
    await expectThrowsAsync(() => process({w: {js: {args: ['x'], fn}},
      v: {validators: {input: 'w', names: ['LibTests:CheckTolerance']}}}), /unknown input alias w/);
    await expectThrowsAsync(() => process({1: {js: {args: ['x'], fn}}}), /must not be a number/);
    const pconf: any = await process({
      w: {js: {args: ['x'], fn}},
      v: {js: {args: ['w', 'x'], fn}},
      s: 'func("LibTests:TestAdd2", a: v, b: x)',
      q: {query: {connection: 'System:Datagrok', sql: 'select 1', args: {a: {var: 's'}}}},
    });
    expectDeepEqual(Object.keys(pconf.links[0].params.sources), ['w', 'v', 's', 'q']);
  });

  test('A source reads the sources declared before it', async () => {
    const last = {js: {args: ['x'], fn: (x: number) => x * 10}};
    const half = {js: {args: ['last'], fn: (t: number) => t / 2}};
    const sum = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'half'}, b: {var: 'x'}}}};
    const sync = resolveSources(valuesController({x: 2}), {last, half}, {$all: {}, x: 2});
    expect(sync instanceof Promise, false);
    expectDeepEqual(sync, {last: 20, half: 10});
    expectDeepEqual(await resolveSources(valuesController({x: 2}), {last, half, sum}, {$all: {}, x: 2}),
      {last: 20, half: 10, sum: 12});
  });

  test('Independent async sources start together', async () => {
    const original = grok.functions.call;
    const started: number[] = [];
    const releases: (() => void)[] = [];
    (grok.functions as any).call = (_name: string, params: any) => {
      started.push(params.a);
      return new Promise((resolve) => releases.push(() => resolve(params.a + params.b)));
    };
    try {
      const one = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: {var: 'x'}}}};
      const two = {func: {name: 'LibTests:TestAdd2', args: {a: 2, b: {var: 'x'}}}};
      const pending = resolveSources(noCallController(), {one, two}, {$all: {}, x: 1});
      expectDeepEqual(started, [1, 2]);
      releases.forEach((release) => release());
      expectDeepEqual(await pending, {one: 2, two: 3});
    } finally {
      (grok.functions as any).call = original;
    }
  });

  test('A source waits only for the async sources it reads', async () => {
    const order: string[] = [];
    let release!: () => void;
    const gate = new Promise<void>((resolve) => release = resolve);
    const slow = {js: {args: ['x'], fn: async (x: number) => {
      await gate;
      order.push('slow');
      return x;
    }}};
    const after = {js: {args: ['slow'], fn: (s: number) => {
      order.push('after');
      return s * 2;
    }}};
    const fast = {js: {args: ['x'], fn: (x: number) => {
      order.push('fast');
      return x + 1;
    }}};
    const pending = resolveSources(valuesController({x: 1}), {slow, after, fast});
    expectDeepEqual(order, ['fast']);
    release();
    expectDeepEqual(await pending, {slow: 1, after: 2, fast: 2});
    expectDeepEqual(order, ['fast', 'slow', 'after']);
  });

  test('A chain of sources that read no input is loaded once per link', async () => {
    await countCalls(async (calls) => {
      const cache = new Map<string, any>();
      const base = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: 2}}};
      const next = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'base'}, b: 1}}};
      const varying = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'next'}, b: {var: 'x'}}}};
      const controller = noCallController(cache);
      expectDeepEqual(await resolveSources(controller, {base, next, varying}, {$all: {}, x: 1}),
        {base: 3, next: 4, varying: 5});
      expectDeepEqual(await resolveSources(controller, {base, next, varying}, {$all: {}, x: 5}),
        {base: 3, next: 4, varying: 9});
      expect(calls(), 4);
      expectDeepEqual([...cache.keys()], ['base', 'next']);
    });
  });

  test('A failed source fails the sources reading it and both are retried', async () => {
    await countCalls(async () => {
      const cache = new Map<string, any>();
      const base = {func: {name: 'LibTests:TestAdd2', args: {a: 1, b: 2}}};
      const next = {func: {name: 'LibTests:TestAdd2', args: {a: {var: 'base'}, b: 1}}};
      await expectThrowsAsync(async () => {
        await resolveSources(noCallController(cache), {base, next});
      }, /Source down/);
      expectDeepEqual([...cache.keys()], []);
      expectDeepEqual(await resolveSources(noCallController(cache), {base, next}), {base: 3, next: 4});
      expectDeepEqual([...cache.keys()], ['base', 'next']);
    }, true);
  });

  test('A separate effect brings the sources its sources read', async () => {
    const pconf: any = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: ['m:step1/a', 'o:step1/b'], to: 't:step2/a',
      sources: {
        last: {js: {args: ['m'], fn: (m: number) => m}},
        unused: {js: {args: ['o'], fn: (o: number) => o}},
        half: 'func("LibTests:TestAdd2", a: last, b: 1)',
      },
      effects: [
        {effect: 'items', targets: 't', items: {var: 'half'}, separate: true},
        {effect: 'meta', targets: 't', meta: {o: {var: 'o'}}},
      ],
    }]));
    const byId = Object.fromEntries(pconf.links.map((link: any) => [link.id, link]));
    expectDeepEqual(Object.keys(byId['r::meta:1'].params.sources), ['last', 'half']);
    expectDeepEqual(byId['r::meta:1'].from.map((io: any) => io.name), ['m']);
  });

  test('A source chain feeds a rule in the tree', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      sources: {
        last: {js: {args: ['m'], fn: (m: number) => m * 10}},
        half: {js: {args: ['last'], fn: (t: number) => t / 2}},
        list: {js: {args: ['half', 'm'], fn: (h: number, m: number) => [String(h), String(m)]}},
      },
      effects: [{effect: 'items', targets: 't', items: {var: 'list'}}],
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
      cold('---a').subscribe(() => store.setState('a', 4));
      cold('----a').subscribe(() => metas.push(outBridge.meta.a.value));
    });
    expectDeepEqual(metas, [{items: ['5', '1']}, {items: ['20', '4']}]);
  });
});

const cities = {c: {choices: {input: 'city', call: '$call'}}};
const callController = (call: DG.FuncCall) => ({
  hasCall: () => true,
  getFirst: (name: string) => name === '$call' ? call : call.inputs[name],
  getMatchedPositions: () => [{path: [], position: 0, ioName: 'city'}],
}) as any;
const items = async (call: DG.FuncCall): Promise<string[]> =>
  (await resolveSources(callController(call), cities)).c.items;

// evaluated choices need platform 1.28.0 or later
category('ComputeUtils: Driver rule choices sources', async () => {
  const prepare = (region: string) => DG.Func.byName('LibTests:TestCountingChoices').prepare({region, n: 1});

  test('Choices are evaluated again only when a dependency changes', async () => {
    const call = prepare('FR');
    const [first] = await items(call);
    call.inputs['n'] = 2;
    expectDeepEqual(await items(call), [first]);
    call.inputs['region'] = 'DE';
    const [changed] = await items(call);
    expect(changed.startsWith('DE-'), true);
    call.inputs['region'] = 'FR';
    const [back] = await items(call);
    expect(back.startsWith('FR-') && back !== first, true);
  });

  test('Steps of one function keep their own choices', async () => {
    const frCall = prepare('FR');
    const usCall = prepare('US');
    const [fr] = await items(frCall);
    const [us] = await items(usCall);
    expect(fr.startsWith('FR-') && us.startsWith('US-'), true);
    expectDeepEqual(await items(frCall), [fr]);
  });
});

type Deferred<T> = {promise: Promise<T>, resolve: (value: T) => void, reject: (error: any) => void};

function deferred<T>(): Deferred<T> {
  let resolve!: (value: T) => void;
  let reject!: (error: any) => void;
  const promise = new Promise<T>((res, rej) => {
    resolve = res;
    reject = rej;
  });
  return {promise, resolve, reject};
}

// a FuncCall stand-in whose evaluations wait until the test answers them
function stubChoicesCall(inputs: Record<string, any>) {
  const evaluations: {region: string, answer: Deferred<any>}[] = [];
  let started = deferred<void>();
  const call = {
    inputs,
    func: {inputs: Object.keys(inputs).map((name) => ({name}))},
    evalParamChoices: () => {
      const answer = deferred<any>();
      evaluations.push({region: inputs.region, answer});
      started.resolve();
      started = deferred<void>();
      return answer.promise;
    },
  } as unknown as DG.FuncCall;
  return {
    call,
    evaluations,
    nextEvaluation: () => started.promise,
    answer: (idx: number) => evaluations[idx].answer.resolve(
      {items: [`${evaluations[idx].region}-${idx}`], values: {}, lookup: null, dependsOn: ['region']}),
    fail: (idx: number) => evaluations[idx].answer.reject(new Error('Choices failed')),
  };
}

category('ComputeUtils: Driver rule choices in flight', async () => {
  test('A dependency change during an evaluation evaluates again', async () => {
    const {call, evaluations, nextEvaluation, answer} = stubChoicesCall({region: 'EU', n: 1});
    const inFlight = items(call);
    call.inputs['region'] = 'US';
    const latest = items(call);
    const second = nextEvaluation();
    answer(0);
    await second;
    answer(1);
    expectDeepEqual(await latest, ['US-1']);
    expectDeepEqual(await inFlight, ['EU-0']);
    expectDeepEqual(await items(call), ['US-1']);
    expect(evaluations.length, 2);
  });

  test('An unrelated change during an evaluation shares it', async () => {
    const {call, evaluations, answer} = stubChoicesCall({region: 'EU', n: 1});
    const inFlight = items(call);
    call.inputs['n'] = 2;
    const shared = items(call);
    answer(0);
    expectDeepEqual(await shared, ['EU-0']);
    expectDeepEqual(await inFlight, ['EU-0']);
    expect(evaluations.length, 1);
  });

  test('A failed evaluation fails its waiters and the next run evaluates again', async () => {
    const {call, evaluations, answer, fail} = stubChoicesCall({region: 'EU', n: 1});
    const runs = Promise.allSettled([items(call), items(call)]);
    fail(0);
    const errors = (await runs).map((run) => run.status === 'rejected' ? run.reason.message : run.status);
    expectDeepEqual(errors, ['Choices failed', 'Choices failed']);
    const retry = items(call);
    answer(1);
    expectDeepEqual(await retry, ['EU-1']);
    expect(evaluations.length, 2);
  });
});
