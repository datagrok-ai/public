import {category, test, before} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FuncCallInstancesBridge} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {PipelineLinkConfigurationInput} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import * as DG from 'datagrok-api/dg';
import {evaluate} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-expressions';
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

category('ComputeUtils: Driver links rule', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Expand rule into family links', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: ['m:step1/a', 'n:step1/b'],
      to: ['t1:step2/a', 't2:step2/b'],
      when: {'>': [{var: 'm'}, 0]},
      debounce: 0,
      effects: [
        {effect: 'hide', targets: 't1'},
        {effect: 'items', targets: ['t2'], items: ['x', 'y']},
        {effect: 'error', targets: 't1', message: {cat: ['bad ', {var: 'all.n'}]}},
        {effect: 'set', targets: 't2', value: {var: 'm'}},
      ],
    }]));
    const links = pconf.links!;
    expectDeepEqual(links.map((l) => l.id), ['r::meta', 'r::validator', 'r::data']);
    expectDeepEqual(links.map((l) => l.type), ['meta', 'validator', 'data']);
    expectDeepEqual(links.map((l) => l.from.map((io) => io.name)), [['m', 'n'], ['m', 'n'], ['m', 'n']]);
    expectDeepEqual(links.map((l) => l.to.map((io) => io.name)), [['t1', 't2'], ['t1'], ['t2']]);
    expectDeepEqual(links.map((l) => l.params!.effects.map((e: any) => e.effect)),
      [['hide', 'items'], ['error'], ['set']]);
    const when = {'>': [{var: 'm'}, 0]};
    expectDeepEqual(links.map((l) => l.params!.when), [when, when, when]);
    expectDeepEqual((links[1] as any).debounce, 0);
    expectDeepEqual(links.every((l) => typeof l.handler === 'function'), true);
  });

  test('Expand only families present', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/a',
      effects: [{effect: 'show', targets: 't'}],
    }]));
    expectDeepEqual(pconf.links!.map((l) => l.id), ['r::meta']);
  });

  test('Reject invalid rules', async () => {
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [{effect: 'hide', targets: 'zzz'}],
    }])), /unknown output alias zzz/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: ['t:step2/a', 'u:step2/b'],
      effects: [{effect: 'hide', targets: 't'}],
    }])), /output alias u is not targeted/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      when: {'>': [{var: 'q'}, 0]},
      effects: [{effect: 'hide', targets: 't'}],
    }])), /unknown input alias q/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [{effect: 'items', targets: 't', items: {var: 'nope'}}],
    }])), /unknown input alias nope/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [],
    }])), /effects list is empty/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [{effect: 'items', targets: 't', items: [{var: 'nested'}, 'y']}],
    }])), /unknown input alias nested/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      when: {'!!': [{var: ['dflt', 1]}]},
      effects: [{effect: 'hide', targets: 't'}],
    }])), /unknown input alias dflt/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      when: {missing: ['gone']},
      effects: [{effect: 'hide', targets: 't'}],
    }])), /unknown input alias gone/);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      effects: [{effect: 'hide', targets: 't'}, {effect: 'show', targets: ['t']}],
    }])), /hide\/show applied twice/);
    await getProcessedConfig(twoSteps([{
      id: 'ok', type: 'rule', from: 'm:step1/a', to: 't:step2/a',
      when: {missing_some: [1, ['m', 'all.m']]},
      effects: [{effect: 'meta', targets: 't', meta: {cfg: {literal: {var: 'not an alias'}}}}],
    }]));
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([{
      id: 'r', type: 'rule', from: 'm(template):step1/a|b', to: 't:step2/a',
      effects: [{effect: 'hide', targets: 't'}],
    }])), /\(template\) flag is not allowed/);
  });

  test('Hide and show toggle with the condition', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: ['t:step2/a', 'u:step2/b'],
      when: {'>': [{var: 'm'}, 0]},
      effects: [
        {effect: 'hide', targets: 't'},
        {effect: 'show', targets: 'u'},
      ],
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a-b').subscribe((v) => {
        inNode.getItem().getStateStore().setState('a', v === 'a' ? 1 : -1);
      });
      expectObservable(outBridge.meta.a).toBe('ab-c', {a: undefined, b: {hidden: true}, c: {hidden: false}});
      expectObservable(outBridge.meta.b).toBe('ab-c', {a: undefined, b: {hidden: false}, c: {hidden: true}});
    });
  });

  test('Items and meta effects evaluate expressions', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/a',
      effects: [
        {effect: 'items', targets: 't', items: {if: [{'>': [{var: 'm'}, 0]}, ['x', 'y'], ['z']]}},
        {effect: 'meta', targets: 't', meta: {units: 'K', twice: {'*': [{var: 'm'}, 2]}}},
      ],
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a-b').subscribe((v) => {
        inNode.getItem().getStateStore().setState('a', v === 'a' ? 1 : -1);
      });
      expectObservable(outBridge.meta.a).toBe('ab-c', {
        a: undefined,
        b: {items: ['x', 'y'], units: 'K', twice: 2},
        c: {items: ['z'], units: 'K', twice: -2},
      });
    });
  });

  test('Items and meta effects are dropped when the condition is off', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/a',
      when: {'>': [{var: 'm'}, 0]},
      effects: [
        {effect: 'items', targets: 't', items: ['x', 'y']},
        {effect: 'meta', targets: 't', meta: {units: 'K'}},
      ],
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a-b').subscribe((v) => {
        inNode.getItem().getStateStore().setState('a', v === 'a' ? 1 : -1);
      });
      expectObservable(outBridge.meta.a).toBe('ab-c', {
        a: undefined,
        b: {items: ['x', 'y'], units: 'K'},
        c: {},
      });
    });
  });

  test('Custom expression operations', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'x', [1, 2, 3]),
      DG.Column.fromList('int', 'n', [1, 2, 3]),
      DG.Column.fromList('string', 's', ['a', 'b', 'c']),
    ]);
    df.col('s')!.semType = 'Text';
    const ctx = {all: {}, df, s: 'abcd', list: [1, 2]};
    expectDeepEqual(evaluate({columns: [{var: 'df'}]}, ctx), ['x', 'n', 's']);
    expectDeepEqual(evaluate({columns: [{var: 'df'}, 'numerical']}, ctx), ['x', 'n']);
    expectDeepEqual(evaluate({columns: [{var: 'df'}, 'categorical']}, ctx), ['s']);
    expectDeepEqual(evaluate({columns: [{var: 'df'}, 'int']}, ctx), ['n']);
    expectDeepEqual(evaluate({columns: [{var: 'df'}, 'Text']}, ctx), ['s']);
    expectDeepEqual(evaluate({columns: [{var: 'missing'}]}, ctx), []);
    expectDeepEqual(evaluate({len: [{var: 'df'}]}, ctx), 3);
    expectDeepEqual(evaluate({len: [{var: 's'}]}, ctx), 4);
    expectDeepEqual(evaluate({len: [{var: 'list'}]}, ctx), 2);
    expectDeepEqual(evaluate({len: [{var: 'missing'}]}, ctx), 0);
    const spec = [['x', 'double'], ['n', 'string'], ['s', 'Text'], ['q'], 'x'];
    expectDeepEqual(evaluate({columnsMissing: [{var: 'df'}, spec]}, ctx), ['n (string)', 'q']);
    expectDeepEqual(evaluate({columnsMissing: [{var: 'missing'}, spec]}, ctx),
      ['x (double)', 'n (string)', 's (Text)', 'q', 'x']);
    expectDeepEqual(evaluate({columnsMissing: [{var: 'df'}, []]}, ctx), []);
  });

  test('Rules validate dataframe columns and feed dropdown items', async () => {
    const spec = [['x', 'double'], ['s', 'string']];
    const pconf = await getProcessedConfig({
      id: 'pipeline1',
      type: 'static',
      steps: [
        {id: 'step1', nqName: 'LibTests:TestDF1'},
        {id: 'step2', nqName: 'LibTests:TestAdd2'},
      ],
      links: [{
        id: 'schema',
        type: 'rule',
        from: 'df:step1/df',
        to: 't:step1/df',
        debounce: 0,
        when: {and: [{var: 'df'}, {'!!': {columnsMissing: [{var: 'df'}, spec]}}]},
        effects: [{
          effect: 'error', targets: 't', message: {cat: ['Missing: ', {columnsMissing: [{var: 'df'}, spec]}]},
        }],
      }, {
        id: 'choices',
        type: 'rule',
        from: 'table:step1/df',
        to: 'c:step2/a',
        effects: [{effect: 'items', targets: 'c', items: {columns: [{var: 'table'}, 'numerical']}}],
      }],
    });
    const partial = DG.DataFrame.fromColumns([DG.Column.fromList('double', 'x', [1])]);
    const full = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'x', [1]),
      DG.Column.fromList('int', 'n', [1]),
      DG.Column.fromList('string', 's', ['a']),
    ]);
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      const snap = () => snapshots.push([node.validationInfo$.value, outBridge.meta.a.value]);
      cold('-a').subscribe(() => node.getStateStore().setState('df', partial));
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => node.getStateStore().setState('df', full));
      cold('----a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [
      [{df: {errors: [{description: 'Missing: s (string)'}], warnings: [], notifications: []}}, {items: ['x']}],
      [{}, {items: ['x', 'n']}],
    ]);
  });

  test('Literal escape and null condition', async () => {
    const ctx = {all: {}, m: 1};
    expectDeepEqual(evaluate({literal: {foo: 1}}, ctx), {foo: 1});
    expectDeepEqual(evaluate({if: [{'>': [{var: 'm'}, 0]}, {literal: {var: 'kept'}}, 'no']}, ctx), {var: 'kept'});
    expectDeepEqual(evaluate([{literal: {a: 1}}, {literal: {b: 2}}], ctx), [{a: 1}, {b: 2}]);
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/a',
      when: null,
      effects: [{effect: 'hide', targets: 't'}, {effect: 'meta', targets: 't', meta: {cfg: {literal: {foo: 1}}}}],
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => inStore.setState('a', -5));
      expectObservable(outBridge.meta.a).toBe('ab', {a: undefined, b: {hidden: true, cfg: {foo: 1}}});
    });
  });

  test('Validation message follows the condition', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: ['a1:step1/a', 'b1:step1/b'],
      to: 't:step1/a',
      when: {'<': [{var: 'a1'}, {var: 'b1'}]},
      debounce: 0,
      effects: [{effect: 'error', targets: 't', message: {cat: ['a must be at least ', {var: 'b1'}]}}],
    }]));
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const store = node.getStateStore();
      cold('-a').subscribe(() => {
        store.setState('b', 2);
        store.setState('a', 1);
      });
      cold('--a').subscribe(() => snapshots.push(node.validationInfo$.value));
      cold('---a').subscribe(() => store.setState('a', 5));
      cold('----a').subscribe(() => snapshots.push(node.validationInfo$.value));
    });
    expectDeepEqual(snapshots, [
      {a: {errors: [{description: 'a must be at least 2'}], warnings: [], notifications: []}},
      {},
    ]);
  });

  test('Set writes a default value without consistency tracking', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/b',
      when: {'>': [{var: 'm'}, 0]},
      effects: [{effect: 'set', targets: 't', value: {'*': [{var: 'm'}, 2]}}],
    }]));
    const restrictions: any[] = [];
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => inNode.getItem().getStateStore().setState('a', 2));
      cold('--a').subscribe(() => restrictions.push(outBridge.inputRestrictions$.value['b']));
      cold('---a').subscribe(() => inNode.getItem().getStateStore().setState('a', -1));
      cold('----a').subscribe(() => restrictions.push(outBridge.inputRestrictions$.value['b']));
      expectObservable(outBridge.getStateChanges('b')).toBe('ab', {a: undefined, b: 4});
    });
    expectDeepEqual(restrictions, [undefined, undefined]);
  });

  test('Set with a restriction clears it when off', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/b',
      when: {'>': [{var: 'm'}, 0]},
      effects: [{effect: 'set', targets: 't', value: {'*': [{var: 'm'}, 2]}, restriction: 'restricted'}],
    }]));
    const restrictions: any[] = [];
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inNode = tree.nodeTree.getNode([{idx: 0}]);
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => inNode.getItem().getStateStore().setState('a', 2));
      cold('--a').subscribe(() => restrictions.push(outBridge.inputRestrictions$.value['b']));
      cold('---a').subscribe(() => inNode.getItem().getStateStore().setState('a', -1));
      cold('----a').subscribe(() => restrictions.push(outBridge.inputRestrictions$.value['b']));
      expectObservable(outBridge.getStateChanges('b')).toBe('ab', {a: undefined, b: 4});
    });
    expectDeepEqual(restrictions, [{assignedValue: 4, type: 'restricted'}, undefined]);
  });

  test('Data effects do not touch a step whose outputs are current', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: ['t:step2/a', 'u:step2/b'],
      effects: [
        {effect: 'set', targets: 't', value: 42},
        {effect: 'clear', targets: 'u'},
      ],
    }, {
      id: 'tracked',
      type: 'rule',
      from: 'n:step1/b',
      to: 'v:step2/a',
      effects: [{effect: 'set', targets: 'v', value: 7, restriction: 'restricted'}],
    }]));
    const restrictions: any[] = [];
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outNode = tree.nodeTree.getNode([{idx: 1}]).getItem() as FuncCallNode;
      const outBridge = outNode.getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => {
        outBridge.setState('a', 1, 'none');
        outBridge.setState('b', 2, 'none');
        outNode.setOutdatedStatus(false);
      });
      cold('--a').subscribe(() => inStore.setState('a', 5));
      cold('---a').subscribe(() => inStore.setState('b', 5));
      cold('----a').subscribe(() => restrictions.push(outBridge.inputRestrictions$.value['a']));
      expectObservable(outBridge.getStateChanges('a')).toBe('ab', {a: undefined, b: 1});
      expectObservable(outBridge.getStateChanges('b')).toBe('ab', {a: undefined, b: 2});
    });
    expectDeepEqual(restrictions, [{assignedValue: 7, type: 'restricted'}]);
  });

  test('Unmatched optional targets are skipped', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: ['t:step2/a', 'u(optional):step2/nonexistent'],
      effects: [{effect: 'hide', targets: ['t', 'u']}],
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => inStore.setState('a', 1));
      expectObservable(outBridge.meta.a).toBe('ab', {a: undefined, b: {hidden: true}});
    });
  });

  test('Set with runOnInit provides a default at init', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: [],
      to: 't:step2/b',
      runOnInit: true,
      effects: [{effect: 'set', targets: 't', value: 42}],
    }]));
    const values: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('--a').subscribe(() => values.push([outBridge.getState('b'), outBridge.inputRestrictions$.value['b']]));
    });
    expectDeepEqual(values, [[42, undefined]]);
  });

  test('Default validator treats only null and empty string as missing', async () => {
    const pconf = await getProcessedConfig(twoSteps([]));
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const store = node.getStateStore();
      const snap = () => snapshots.push(Object.keys(node.validationInfo$.value).sort());
      cold('-a').subscribe(() => {
        store.setState('a', 0);
        store.setState('b', false);
      });
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => {
        store.setState('a', '');
        store.setState('b', null);
      });
      cold('----a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [[], ['a', 'b']]);
  });

  test('Clear resets the target on every input change', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/b',
      effects: [{effect: 'clear', targets: 't'}],
    }]));
    const restrictions: any[] = [];
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outBridge = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore() as FuncCallInstancesBridge;
      cold('-a').subscribe(() => outBridge.setState('b', 5));
      cold('--a').subscribe(() => inStore.setState('a', 1));
      cold('---a').subscribe(() => outBridge.setState('b', 7));
      cold('----a').subscribe(() => inStore.setState('a', 2));
      cold('-----a').subscribe(() => restrictions.push(outBridge.inputRestrictions$.value['b']));
      expectObservable(outBridge.getStateChanges('b')).toBe('abcdc', {a: undefined, b: 5, c: null, d: 7});
    });
    expectDeepEqual(restrictions, [undefined]);
  });

  test('Clear resets the target only while the condition holds', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'm:step1/a',
      to: 't:step2/b',
      when: {'>': [{var: 'm'}, 0]},
      effects: [{effect: 'clear', targets: 't'}],
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outStore = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore();
      cold('-a').subscribe(() => outStore.setState('b', 5));
      cold('--a').subscribe(() => inStore.setState('a', -1));
      cold('---a').subscribe(() => inStore.setState('a', 1));
      cold('----a').subscribe(() => outStore.setState('b', 7));
      cold('-----a').subscribe(() => inStore.setState('a', 2));
      expectObservable(outStore.getStateChanges('b')).toBe('ab-cdc', {a: undefined, b: 5, c: null, d: 7});
    });
  });

  test('Hidden input suppresses default and custom validators', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'r',
      type: 'rule',
      from: 'bb:step1/b',
      to: 'ta:step1/a',
      when: {'>': [{var: 'bb'}, 0]},
      effects: [{effect: 'hide', targets: 'ta'}],
    }, {
      id: 'custom',
      type: 'validator',
      from: 'va:step1/a',
      to: 'vt:step1/a',
      debounce: 0,
      handler({controller}) {
        controller.setValidation('vt', {warnings: [{description: 'custom'}]});
      },
    }]));
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, defaultValidators: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const bridge = node.getStateStore() as FuncCallInstancesBridge;
      const snap = () => snapshots.push({
        validations: node.validationInfo$.value,
        runnable: bridge.isRunable$.value,
        meta: bridge.meta.a.value,
      });
      cold('-a').subscribe(snap);
      cold('--a').subscribe(() => node.getStateStore().setState('b', 1));
      cold('---a').subscribe(snap);
      cold('----a').subscribe(() => node.getStateStore().setState('b', -1));
      cold('-----a').subscribe(snap);
    });
    const missing = {errors: [{description: 'Missing value'}], warnings: [], notifications: []};
    const missingAndCustom = {
      errors: [{description: 'Missing value'}], warnings: [{description: 'custom'}], notifications: [],
    };
    expectDeepEqual(snapshots, [
      {validations: {a: missingAndCustom, b: missing}, runnable: false, meta: {hidden: false}},
      {validations: {}, runnable: true, meta: {hidden: true}},
      {validations: {a: missingAndCustom}, runnable: false, meta: {hidden: false}},
    ]);
  });

  test('Handlers read link params', async () => {
    const pconf = await getProcessedConfig(twoSteps([{
      id: 'l',
      from: 'in:step1/a',
      to: 'out:step2/a',
      params: {k: 3},
      handler({controller}) {
        controller.setAll('out', controller.getFirst<number>('in') + controller.getParam('k'));
      },
    }]));
    testScheduler.run((helpers) => {
      const {expectObservable, cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const inStore = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore();
      const outStore = tree.nodeTree.getNode([{idx: 1}]).getItem().getStateStore();
      cold('-a-b').subscribe((v) => inStore.setState('a', v === 'a' ? 1 : 2));
      expectObservable(outStore.getStateChanges('a')).toBe('ab-c', {a: undefined, b: 4, c: 5});
    });
  });
});
