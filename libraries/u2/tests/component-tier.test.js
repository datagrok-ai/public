/* The second binding tier of Component (VP-1..VP-4): a declared property without a same-named
   signal becomes a cached, echo-suppressed bind step under `propertyTier`; `link` is the one
   bridge; `Component.is`/`Control.is` are structural; `propertyFields` copies a getter-backed
   property. The designer sites the tier reaches — the panel's live values and the async designer
   actions — are exercised here too, over `tests/dg-stub.mjs`. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal, isWritableSignal} from '../src/core/signals.js';
import {Component, Control} from '../src/core/component.js';
import {propertyFields} from '../src/core/property-like.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {StateSource} from '../src/sources/state.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecNodeRef} from '../src/dg/designer/node-ref.js';
import {propsFor} from '../src/dg/designer/prop-model.js';
import {Property, Stream} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {DesignerActions} = await import('../src/dg/designer/actions.js');

function tier(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** A property-tier control over a plain object: `level` writable (the object clamps it to 10 and,
 * like the platform, announces every write — P7), `unit` read-only, changes through a stream. */
class Gauge extends Control {
  constructor(target, changes = null) {
    super();
    this.target = target;
    this.changes = changes;
    this.reads = 0;
  }

  get propertyTier() {
    return true;
  }

  getProperties() {
    this.reads++;
    return [
      {name: 'level', type: 'int', get: () => this.target.level,
        set: (_t, v) => {
          this.target.level = Math.min(v, 10);
          this.changes?.fire('level');
        }},
      {name: 'unit', type: 'string', get: () => this.target.unit},
    ];
  }

  propertyChanges() {
    return this.changes;
  }
}

const GAUGE = {
  tag: 'u2-gauge',
  description: 'Property-tier control for the tier tests',
  props: [
    {name: 'level', type: 'int', bindable: true, twoWay: true},
    {name: 'unit', type: 'string', bindable: true},
  ],
  example: {tag: 'u2-gauge'},
};

tier('Component.is / Control.is are structural', () => {
  const scope = new Scope();
  const el = document.createElement('div');
  const plain = {scope, bindStep() {}, bindProps() {}, root: el};
  assert.equal(Control.is(plain), true);
  assert.equal(Component.is(plain), true);
  assert.equal(Control.is({...plain, root: undefined}), false);
  assert.equal(Component.is({bindStep() {}, bindProps() {}}), false, 'a bare bind source has no scope');
  const source = new StateSource();
  assert.equal(Component.is(source), true);
  assert.equal(Control.is(source), false);
  assert.equal(Component.is(el), false);
  assert.equal(Component.is(null), false);
  source.dispose();
  scope.dispose();
});

tier('sameValue: arrays element-wise, plain objects key-wise, everything else by identity', () => {
  assert.equal(Component.sameValue([1, 2], [1, 2]), true);
  assert.equal(Component.sameValue({a: [1]}, {a: [1]}), true);
  assert.equal(Component.sameValue([1], [1, 2]), false);
  assert.equal(Component.sameValue({}, []), false);
  assert.equal(Component.sameValue({a: 1}, {b: 1}), false);
  assert.equal(Component.sameValue(NaN, NaN), true);
  const date = new Date(0);
  assert.equal(Component.sameValue(date, new Date(0)), false, 'non-plain objects by identity');
  assert.equal(Component.sameValue(date, date), true);
});

tier('propertyTarget: a meta-stamped component is its own target, and its property reaches its signal', () => {
  const reg = new Registry();
  reg.register({tag: 'u2-ti', description: '', props: [{name: 'value', type: 'string'}],
    create: () => new TextInput({value: 'x'}), example: {tag: 'u2-ti'}});
  const input = reg.get('u2-ti').create({});
  assert.equal(input.propertyTier, false, 'the tier is off for a plain component');
  assert.equal(input.propertyTarget, input);
  const [value] = input.getProperties();
  assert.equal(value.get(input.propertyTarget), input.value.value);
  value.set(input.propertyTarget, 'y');
  assert.equal(input.value.value, 'y');
  assert.equal(input.getProperties().find((p) => p.name === 'nope'), undefined);
  assert.deepEqual(input.bindProps().map((p) => p.name), ['value'], 'bindProps unchanged with the tier off');
  input.dispose();
});

tier('property tier: cached steps, read-only as computed, notifications, no allocation in bindProps', async () => {
  const target = {level: 1, unit: 'mg'};
  const changes = new Stream();
  const gauge = new Gauge(target, changes);

  const props = gauge.bindProps();
  gauge.bindProps();
  assert.deepEqual(props.map((p) => [p.name, p.writable]), [['level', true], ['unit', false]]);
  assert.equal(gauge.reads, 1, 'getProperties() is indexed once');
  assert.equal(gauge._u2.propertySignals.size, 0, 'bindProps allocates no signal');

  const level = gauge.bindStep('level');
  assert.equal(gauge.bindStep('level'), level, 'cached');
  assert.equal(isWritableSignal(level), true);
  assert.equal(level.value, 1);
  level.value = 5;
  assert.equal(target.level, 5, 'a write goes through the property on propertyTarget');

  const unit = gauge.bindStep('unit');
  assert.equal(isWritableSignal(unit), false, 'no set: a computed');
  assert.equal(unit.value, 'mg');
  assert.equal(gauge.bindStep('nope'), null);
  assert.equal(gauge.reads, 1, 'steps read through the index');

  let runs = 0;
  gauge.effect(() => {
    level.value;
    runs++;
  });
  assert.equal(runs, 1);
  target.level = 7;
  changes.fire('level');
  assert.equal(level.value, 5, 'a notification never refreshes inside the event');
  await flush();
  assert.equal(level.value, 7, 'a named notification refreshes the step in a microtask');
  assert.equal(runs, 2);
  changes.fire('level');
  await flush();
  assert.equal(runs, 2, 'a same-value notification does not re-fire dependents');
  target.unit = 'kg';
  target.level = 8;
  changes.fire(null);
  await flush();
  assert.equal(unit.value, 'kg', 'null refreshes every cached step');
  assert.equal(level.value, 8);
  assert.equal(runs, 3);
  assert.equal(changes.count, 1, 'one change subscription');

  level.value = 50;
  assert.equal(target.level, 10, 'the object normalized the write');
  assert.equal(level.value, 50, 'the step converges on the platform\'s notification, never inside the effect');
  await flush();
  assert.equal(level.value, 10);
  assert.equal(runs, 5);

  gauge.dispose();
  assert.equal(changes.count, 0, 'the subscription dies with the scope');
});

tier('link: the one bridge, both ways, echo-suppressed by value', () => {
  const target = {level: 1, unit: 'mg'};
  const gauge = new Gauge(target);
  const source = signal([3]);
  const arrays = new Gauge({list: [1]});
  arrays.getProperties = () => [{name: 'list', type: 'object', get: () => arrays.target.list,
    set: (_t, v) => arrays.target.list = [...v]}];
  arrays.link('list', source, true);
  assert.deepEqual(arrays.target.list, [3]);
  const written = arrays.target.list;
  source.value = [3];
  assert.equal(arrays.target.list, written, 'an equal array is no write');
  arrays.bindStep('list').value = [4, 5];
  assert.deepEqual(source.value, [4, 5]);

  const lvl = signal(2);
  gauge.link('level', lvl, true);
  assert.equal(target.level, 2);
  lvl.value = 6;
  assert.equal(target.level, 6);
  gauge.bindStep('level').value = 9;
  assert.equal(lvl.value, 9);
  const unit = signal('g');
  const warnings = [];
  const warn = console.warn;
  console.warn = (m) => warnings.push(m);
  try {
    gauge.link('unit', unit, true);
    gauge.link('unit', unit, true);
  } finally {
    console.warn = warn;
  }
  assert.equal(gauge.bindStep('unit').value, 'mg', 'a read-only step is not written');
  assert.deepEqual(warnings, ['u2: component.unit is read-only — edits will not flow back'], 'said once');
  gauge.dispose();
  arrays.dispose();
});

tier('propertyFields copies the named fields of a getter-backed property', () => {
  const p = new Property('allowEdit', 'bool', {description: 'd', category: 'Data', nullable: true, choices: ['a']});
  assert.deepEqual(Object.keys({...p}), ['dart'], 'a spread sees nothing');
  assert.deepEqual(propertyFields(p), {name: 'allowEdit', propertyType: 'bool', type: 'bool',
    description: 'd', category: 'Data', nullable: true, choices: ['a']});
});

tier('propsFor: a property-tier node reads its live values through propertyTarget', () => {
  const reg = new Registry();
  const target = {level: 3, unit: 'mg'};
  reg.register({...GAUGE, create: () => new Gauge(target)});
  const node = {tag: 'u2-gauge', name: 'g', props: {level: 3}};
  const instance = renderSpec({$schema: 'dg-ui/1', root: node}, new SpecContext(), reg);
  const ref = new SpecNodeRef(instance, node);
  const [section] = propsFor(ref);
  assert.deepEqual(section.values, {level: 3, unit: 'mg'});
  target.level = 4;
  assert.equal(section.read().level, 4, 'live, not the document');
  instance.dispose();
});

tier('designer actions: an async produce applies on a patch, nothing on null, and sees the build', async () => {
  const reg = new Registry();
  const seen = [];
  reg.register({...GAUGE, create: () => new Gauge({level: 1, unit: 'mg'}), designerActions: [
    {name: 'Set level', produce: async (node, built) => {
      seen.push(built);
      return {op: 'set-prop', name: 'level', value: 5};
    }},
    {name: 'Cancelled', produce: () => Promise.resolve(null)},
    {name: 'Styled', icon: 'paint-brush', produce: () => null},
  ]});
  const node = {tag: 'u2-gauge', name: 'g'};
  const instance = renderSpec({$schema: 'dg-ui/1', root: node}, new SpecContext(), reg);
  const applied = [];
  const actions = new DesignerActions({
    editor: signal({canApply: () => null, uniqueNames: () => (name) => name}),
    instance: () => instance,
    selected: () => node,
    multi: () => [node],
    select() {},
    pending() {},
    apply: (patch) => applied.push(patch),
    refocus() {},
  });
  actions.run('Cancelled');
  await flush();
  assert.deepEqual(applied, []);
  actions.run('Set level');
  await flush();
  assert.deepEqual(applied, [{op: 'set-prop', node, name: 'level', value: 5}]);
  assert.equal(seen[0], instance.nodes().get(node), 'produce receives what the node built');
  const icons = Object.fromEntries(actions.actions().map((a) => [a.name, a.icon]));
  assert.equal(icons.Styled, 'paint-brush', 'a registry verb carries its own icon');
  assert.equal(icons['Set level'], undefined, 'and none is made up for one that does not');
  instance.dispose();
});

tier('Component.init on a foreign object: own state under _u2, nothing else touched, ambient adoption', () => {
  const foreign = Object.create(Component.prototype);
  foreign._properties = 'platform';
  foreign._functions = 'platform';
  const owner = new Scope();
  Scope.runWith(owner, () => Component.init(foreign));
  assert.equal(foreign._properties, 'platform');
  assert.equal(foreign._functions, 'platform');
  assert.equal(foreign.scope instanceof Scope, true);
  assert.equal(Component.is(foreign), true);
  assert.deepEqual(foreign.getFunctions(), [], 'functions live under _u2');
  assert.equal(foreign.aiDescription, null);
  owner.dispose();
  assert.equal(foreign.scope.isDisposed, true, 'adopted by the ambient scope');
});

tier('property tier: a list property answering a new array per read never loops', async () => {
  const target = {list: ['a']};
  const changes = new Stream();
  const gauge = new Gauge(target, changes);
  let reads = 0;
  gauge.getProperties = () => [{name: 'list', type: 'object',
    get: () => {
      reads++;
      return [...target.list];
    },
    set: (_t, v) => {
      target.list = [...v];
      changes.fire('list');
    }}];
  const list = gauge.bindStep('list');
  let runs = 0;
  gauge.effect(() => {
    list.value;
    runs++;
  });
  list.value = ['a', 'b'];
  await flush();
  assert.deepEqual(target.list, ['a', 'b']);
  assert.equal(runs, 2, 'one write, one dependent run — no echo');
  changes.fire('list');
  changes.fire(null);
  await flush();
  assert.equal(runs, 2, 'equal-content arrays refresh nothing');
  assert.ok(reads < 10, `bounded reads (${reads})`);
  gauge.dispose();
});

tier('property tier: a property answering a fresh wrapper per read writes once per user write, never at bind, never loops', async () => {
  class Level {
    constructor(n) {
      this.n = n;
    }
  }
  const target = {level: 1};
  const changes = new Stream();
  const gauge = new Gauge(target, changes);
  let writes = 0;
  gauge.getProperties = () => [{name: 'level', type: 'object',
    get: () => new Level(target.level),
    set: (_t, v) => {
      writes++;
      target.level = v.n;
      changes.fire('level');
    }}];
  const level = gauge.bindStep('level');
  await flush();
  assert.equal(writes, 0, 'the tier compares against what it last observed, not a fresh wrapper');
  level.value = new Level(5);
  for (let i = 0; i < 5; i++)
    await flush();
  assert.equal(target.level, 5);
  assert.equal(writes, 1, 'one write per user write — the notification it triggers writes nothing back');
  assert.equal(level.value.n, 5);
  changes.fire('level');
  await flush();
  assert.equal(writes, 1, 'a refresh answering a new wrapper is no write either');
  gauge.dispose();
});

/** A gauge whose property reads are counted per name, both steps cached, the counters zeroed. */
function countedGauge(target, changes) {
  const gauge = new Gauge(target, changes);
  const gets = {level: 0, unit: 0};
  gauge.getProperties = () => [
    {name: 'level', type: 'int', get: () => {
      gets.level++;
      return target.level;
    }, set: (_t, v) => target.level = v},
    {name: 'unit', type: 'string', get: () => {
      gets.unit++;
      return target.unit;
    }},
  ];
  const level = gauge.bindStep('level');
  const unit = gauge.bindStep('unit');
  gets.level = 0;
  gets.unit = 0;
  return {gauge, gets, level, unit};
}

tier('property tier: two notifications of one name in one tick are one read, in a microtask', async () => {
  const target = {level: 1, unit: 'mg'};
  const changes = new Stream();
  const {gauge, gets, level} = countedGauge(target, changes);
  target.level = 2;
  changes.fire('level');
  changes.fire('level');
  assert.deepEqual(gets, {level: 0, unit: 0}, 'nothing is read inside the event');
  assert.equal(level.value, 1);
  await flush();
  assert.deepEqual(gets, {level: 1, unit: 0});
  assert.equal(level.value, 2);
  changes.fire('nope');
  await flush();
  assert.deepEqual(gets, {level: 1, unit: 0}, 'a name without a cached step reads nothing');
  gauge.dispose();
});

tier('property tier: a null notification refreshes every cached step once, whatever else the tick marked', async () => {
  const target = {level: 1, unit: 'mg'};
  const changes = new Stream();
  const {gauge, gets, level, unit} = countedGauge(target, changes);
  target.level = 2;
  target.unit = 'kg';
  changes.fire('level');
  changes.fire(null);
  changes.fire('unit');
  changes.fire(null);
  await flush();
  assert.deepEqual(gets, {level: 1, unit: 1});
  assert.equal(level.value, 2);
  assert.equal(unit.value, 'kg');
  gauge.dispose();
});

tier('property tier: a notification in flight when the component is disposed reads nothing', async () => {
  const target = {level: 1, unit: 'mg'};
  const changes = new Stream();
  const gauge = new Gauge(target, changes);
  let gets = 0;
  gauge.getProperties = () => [{name: 'level', type: 'int', get: () => {
    gets++;
    return target.level;
  }}];
  gauge.bindStep('level');
  gets = 0;
  changes.fire('level');
  gauge.dispose();
  assert.equal(changes.count, 0, 'the subscription died with the scope');
  changes.fire('level');
  await flush();
  assert.equal(gets, 0);
});
