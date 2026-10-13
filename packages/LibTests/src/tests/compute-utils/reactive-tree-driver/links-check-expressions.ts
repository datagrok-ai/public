import {category, test, before, expect} from '@datagrok-libraries/test/src/test';
import {getProcessedConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {expandChecks} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/checks';
import {evaluate} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/rule-expressions';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {FuncCallNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FuncCallInstancesBridge} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/FuncCallInstancesBridge';
import {TestScheduler} from 'rxjs/testing';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createTestScheduler, expectThrowsAsync, twoSteps, errors} from '../../../test-utils';

const EXPRESSIONS = 'LibTests:TestExpressionInputs';


// GrokScript expression checks need platform 1.28.0 or later
category('ComputeUtils: Driver links check expressions', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Annotations parse expressions', async () => {
    const pconf: any = await getProcessedConfig({id: 'p', type: 'static', steps: [{id: 's', nqName: EXPRESSIONS}]});
    const io = Object.fromEntries(pconf.steps[0].io.map((item: any) => [item.id, item]));
    expectDeepEqual(io.hv.checks, {visible: 'k > 1'});
    expectDeepEqual(io.foo.checks, {validator: 'bar > 3'});
    expectDeepEqual(io.code.checks, {validator: 'startsWith(value, "12")'});
    expect(io.k.checks === undefined, true);
  });

  test('Expand expression checks', async () => {
    const [visible] = expandChecks({visible: 'k > 1'});
    expect(visible.family, 'meta');
    expect(visible.needsInputs, true);
    expectDeepEqual(visible.params,
      {when: {'==': [{script: 'k > 1'}, false]}, effects: [{effect: 'hide', targets: ['$target']}]});
    const [validator] = expandChecks({validator: 'bar > 3'});
    expect(validator.family, 'validator');
    expect(validator.needsInputs, true);
    expectDeepEqual(validator.params, {
      when: {and: [{'!': {missing: ['value']}}, {'!!': {scriptVerdict: 'bar > 3'}}]},
      effects: [{effect: 'error', targets: ['$target'], message: {scriptVerdict: 'bar > 3'}}],
    });
    const [regex] = expandChecks({validator: '/^a/'});
    expect(regex.needsInputs, false);
    const [custom] = expandChecks({validator: 'bar > 3'}, {message: 'Bar too small'});
    expect((custom.params.effects[0] as any).message, 'Bar too small');
  });

  test('Script operations', async () => {
    const ctx = (data: Record<string, any>) => ({$all: {}, ...data});
    expect(evaluate({script: 'k > 1'}, ctx({k: 2})), true);
    expect(evaluate({script: 'k > 1'}, ctx({k: 1})), false);
    expect(evaluate({script: 'k + value'}, ctx({k: 1, value: 2})), 3);
    expect(evaluate({script: 'nope > 1'}, ctx({k: 1})) === undefined, true);
    expect(evaluate({scriptVerdict: 'bar > 3'}, ctx({bar: 2})), 'bar > 3');
    expect(evaluate({scriptVerdict: 'bar > 3'}, ctx({bar: 5})) === null, true);
    expect(evaluate({scriptVerdict: 'startsWith(value, "12")'}, ctx({value: '12ab'})) === null, true);
    expect(evaluate({scriptVerdict: 'startsWith(value, "12")'}, ctx({value: 'ab'})), 'startsWith(value, "12")');
    expect(evaluate({scriptVerdict: 'nope > 1'}, ctx({bar: 2})), 'Error during validation: "nope > 1"');
  });

  test('Enabled hides an input, combined with visible', async () => {
    const pconf =
      await getProcessedConfig({id: 'p', type: 'static', steps: [{id: 's', nqName: 'LibTests:TestEnabledInputs'}]});
    const hidden: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const bridge = tree.nodeTree.getNode([{idx: 0}]).getItem().getStateStore() as FuncCallInstancesBridge;
      const snap = () => hidden.push([bridge.meta.en.value?.hidden, bridge.meta.both.value?.hidden]);
      [2, 1, 0].forEach((k, idx) => {
        cold(`${'-'.repeat(2 * idx + 1)}a`).subscribe(() => bridge.setState('k', k, 'none'));
        cold(`${'-'.repeat(2 * idx + 2)}a`).subscribe(snap);
      });
    });
    expectDeepEqual(hidden, [[false, false], [true, true], [true, true]]);
  });

  test('Expression checks run as default links', async () => {
    const pconf = await getProcessedConfig({id: 'p', type: 'static', steps: [{id: 's', nqName: EXPRESSIONS}]});
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true, annotationChecks: true});
      StateTree.loadOrCreateCalls(tree, true).subscribe();
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const bridge = node.getStateStore() as FuncCallInstancesBridge;
      const store = node.getStateStore();
      const snap = () => snapshots.push({validation: node.validationInfo$.value, hv: bridge.meta.hv.value});
      cold('-a').subscribe(() => {
        store.setState('k', 2);
        store.setState('hv', 1);
        store.setState('foo', 5);
        store.setState('bar', 2);
        store.setState('code', '12ab');
      });
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => {
        store.setState('k', 1);
        store.setState('bar', 5);
        store.setState('code', 'ab');
      });
      cold('----a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [
      {validation: {foo: errors('bar > 3')}, hv: {hidden: false}},
      {validation: {code: errors('startsWith(value, "12")')}, hv: {hidden: true}},
    ]);
  });

  test('Check links with expressions read vars', async () => {
    const pconf: any = await getProcessedConfig(twoSteps([
      {id: 'gate', type: 'check', io: 'step1/a', check: {validator: 'other > 3'}, vars: {other: 'step1/b'}},
      {
        id: 'vis', type: 'check', io: 'step1/b', check: {visible: 'other > 3'}, vars: {other: 'step1/a'},
        message: 'ignored', severity: 'warning', debounce: 100,
      },
    ]));
    expectDeepEqual(pconf.links.map((link: any) => [link.id, link.type, link.from.map((io: any) => io.name)]), [
      ['gate::validator', 'validator', ['value', 'other']],
      ['vis::visible', 'meta', ['value', 'other']],
    ]);
    expect(pconf.links[1].debounce === 100, false);
    expectDeepEqual(pconf.links[1].params.effects, [{effect: 'hide', targets: ['$target']}]);
    const snapshots: any[] = [];
    testScheduler.run((helpers) => {
      const {cold} = helpers;
      const tree = StateTree.fromPipelineConfig({config: pconf, mockMode: true});
      tree.init().subscribe();
      const node = tree.nodeTree.getNode([{idx: 0}]).getItem() as FuncCallNode;
      const bridge = node.getStateStore() as FuncCallInstancesBridge;
      const store = node.getStateStore();
      const snap = () => snapshots.push({validation: node.validationInfo$.value, b: bridge.meta.b.value});
      cold('-a').subscribe(() => {
        store.setState('a', 1);
        store.setState('b', 2);
      });
      cold('--a').subscribe(snap);
      cold('---a').subscribe(() => store.setState('b', 5));
      cold('----a').subscribe(snap);
      cold('-----a').subscribe(() => store.setState('a', 4));
      cold('------a').subscribe(snap);
    });
    expectDeepEqual(snapshots, [
      {validation: {a: errors('other > 3')}, b: {hidden: true}},
      {validation: {}, b: {hidden: true}},
      {validation: {}, b: {hidden: false}},
    ]);
    await expectThrowsAsync(() => getProcessedConfig(twoSteps([
      {id: 'bad', type: 'check', io: 'step1/a', check: {visible: 'value > 3'}, vars: {value: 'step1/b'}},
    ])));
  });
});
