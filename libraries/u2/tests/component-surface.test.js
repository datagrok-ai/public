/* The Component introspection surface (DD7): what a component answers without the platform —
   registry-generated properties over live signals, declared functions, component and DOM events,
   the AI briefing and the status snapshot. The dg half (a real DG.Widget delegating to all of it)
   is tests/widget-host.test.js. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Component, Control} from '../src/core/component.js';
import {button} from '../src/core/elements.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';
import {renderSpec} from '../src/spec/spec.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function widget(name, body) {
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

/** A component of the shape a registration builds: signal-backed props, a callable, an event. */
class Gauge extends Control {
  constructor() {
    super();
    this.value = signal(10);
    this.caption = 'Pressure';
    this.registerFunction({
      name: 'reset',
      description: 'Sets the value back to zero.',
      inputs: [{name: 'to', type: 'int'}],
      apply: (params) => {
        this.value.value = params?.to ?? 0;
        return this.value.peek();
      },
    });
  }

  bump() {
    this.value.value = this.value.peek() + 1;
    this.fireEvent('bumped', this.value.peek());
  }
}

const GAUGE_META = {
  tag: 'u2-gauge',
  description: 'A gauge',
  props: [
    {name: 'value', type: 'int', description: 'The reading.'},
    {name: 'caption', type: 'string'},
  ],
  events: ['bumped'],
  example: {tag: 'u2-gauge'},
  create: () => new Gauge(),
};

function gauge() {
  const g = new Gauge();
  g.componentMeta = GAUGE_META;
  return g;
}

widget('properties are generated from the registry metadata, over live signals', () => {
  const g = gauge();
  const props = g.getProperties();
  assert.deepEqual(props.map((p) => p.name), ['value', 'caption']);
  assert.equal(props[0].description, 'The reading.');

  assert.equal(props[0].get(null), 10);
  g.bump();
  assert.equal(props[0].get(null), 11, 'reads the signal, not a snapshot');

  props[0].set(null, 42);
  assert.equal(g.value.peek(), 42, 'writes through to the signal');

  assert.equal(props[1].get(null), 'Pressure', 'a plain member reads as itself');
  assert.equal(props[1].set, undefined, 'and stays read-only');
  g.dispose();
});

widget('a component without metadata has no properties', () => {
  const g = new Gauge();
  g.name = 'gauge1';
  assert.deepEqual(g.getProperties(), []);
  g.dispose();
});

widget('the registry stamps the metadata of the tag it built', () => {
  const reg = new Registry();
  reg.register(GAUGE_META);
  const meta = reg.get('u2-gauge');
  const loose = meta.create({});
  assert.equal(loose.componentMeta, meta, 'every registry-driven construction is stamped');
  loose.dispose();

  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-gauge', name: 'g'}}, undefined, reg);
  const built = instance.node('g');
  assert.equal(built.componentMeta, meta);
  assert.deepEqual(built.getProperties().filter((p) => p.category !== 'Appearance').map((p) => p.name),
    ['value', 'caption']);
  instance.dispose();
});

widget('properties of a spec-rendered input read its live state, options included', () => {
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-text-input', name: 'q', props: {
    label: 'Search', placeholder: 'type here', value: 'aspirin', search: true,
  }}}, undefined, reg);
  const input = instance.node('q');
  const props = new Map(input.getProperties().map((p) => [p.name, p]));

  assert.equal(props.get('label').get(null), 'Search');
  assert.equal(props.get('placeholder').get(null), 'type here');
  assert.equal(props.get('search').get(null), true);
  assert.equal(props.get('value').get(null), 'aspirin');
  assert.equal(props.get('inline').set, undefined, 'a construction-only option renders read-only');
  assert.notEqual(props.get('label').set, undefined, 'an accessor-backed prop is writable');

  props.get('label').set(null, 'Filter');
  assert.equal(input.root.querySelector('.u2-input-label').textContent, 'Filter',
    'the accessor write reaches the live element');
  props.get('value').set(null, 'ibuprofen');
  assert.equal(input.value.peek(), 'ibuprofen', 'a signal-backed prop still writes through');
  instance.dispose();
});

widget('a component that keeps nothing still reports the props the spec set on it', () => {
  const reg = new Registry();
  reg.register({
    tag: 'u2-fake-button',
    create: () => new Control(button('', () => {})),
    description: 'A create that wraps a bare element and keeps no state of its own',
    props: [{name: 'text', type: 'string'}, {name: 'primary', type: 'bool'}],
    example: {tag: 'u2-fake-button'},
  });
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-fake-button', name: 'save',
    props: {text: 'Save', primary: true}}}, undefined, reg);
  const props = new Map(instance.node('save').getProperties().map((p) => [p.name, p]));

  assert.equal(props.get('text').get(null), 'Save');
  assert.equal(props.get('primary').get(null), true);
  assert.equal(props.get('text').set, undefined, 'a construction-time value renders read-only');
  instance.dispose();
});

widget('a prop the component stores under another name reports what the spec set', () => {
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-splitter', name: 'layout',
    props: {direction: 'vertical'}}}, undefined, reg);
  const direction = instance.node('layout').getProperties().find((p) => p.name === 'direction');

  assert.equal(direction.get(null), 'vertical');
  instance.dispose();
});

widget('a declared name prop is the form key, not the component identity', () => {
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-text-input', name: 'firstNameField',
    props: {label: 'First name', name: 'firstName'}}}, undefined, reg);
  const key = instance.node('firstNameField').getProperties().find((p) => p.name === 'name');

  assert.ok(key, 'the form key a registration declares IS a property');
  assert.equal(key.get(null), 'firstName');
  assert.notEqual(key.get(null), 'firstNameField', 'and never reads the spec identity');
  instance.dispose();
});

widget('declared functions are what the component registered', async () => {
  const g = gauge();
  const [reset] = g.getFunctions();
  assert.equal(reset.name, 'reset');
  assert.deepEqual(reset.inputs.map((p) => p.name), ['to']);
  assert.equal(await reset.apply({to: 7}), 7);
  assert.equal(g.value.peek(), 7);
  g.dispose();
});

widget('onEvent delivers component events, by id and unfiltered', () => {
  const g = gauge();
  const seen = [];
  const all = [];
  const sub = g.onEvent('bumped').subscribe((x) => seen.push(x));
  const every = g.onEvent().subscribe((x) => all.push(x));

  g.bump();
  assert.deepEqual(seen, [11]);
  assert.deepEqual(all, [11]);

  sub.unsubscribe();
  g.bump();
  assert.deepEqual(seen, [11], 'unsubscribed stops receiving');
  assert.deepEqual(all, [11, 12]);
  every.unsubscribe();
  g.dispose();
});

widget('onEvent on a control also delivers DOM events of that name from its root', () => {
  const g = gauge();
  const seen = [];
  const sub = g.onEvent('click').subscribe((e) => seen.push(e.type));

  fire(g.root, 'click');
  assert.deepEqual(seen, ['click']);

  sub.unsubscribe();
  fire(g.root, 'click');
  assert.deepEqual(seen, ['click'], 'the listener is removed with the subscription');
  g.dispose();
});

widget('a subscription left open is disposed with the component', () => {
  const g = gauge();
  const seen = [];
  g.onEvent('click').subscribe(() => seen.push(1));
  g.dispose();

  fire(g.root, 'click');
  assert.deepEqual(seen, []);
});

widget('a component-event subscription left open dies with the component too', () => {
  const g = gauge();
  const seen = [];
  g.onEvent('bumped').subscribe((x) => seen.push(x));
  g.onEvent().subscribe((x) => seen.push(x));
  g.dispose();

  g.bump();
  assert.deepEqual(seen, []);
});

widget('the status snapshot reports the declared events and description', () => {
  const g = gauge();
  const status = g.getWidgetStatus();
  assert.deepEqual(status.events, [{name: 'bumped', eventName: 'bumped', description: ''}]);
  assert.equal(status.description, 'A gauge');
  assert.deepEqual(status.parts, {});
  assert.equal(status.error, null);
  g.dispose();
});

widget('the AI briefing falls back to usage, then description, and is overridable', () => {
  const c = new Component();
  assert.equal(c.aiDescription, null);
  c.componentMeta = {props: [], description: 'A gauge'};
  assert.equal(c.aiDescription, 'A gauge');
  c.componentMeta = {props: [], description: 'A gauge', usage: 'Reach for it to show one number.'};
  assert.equal(c.aiDescription, 'Reach for it to show one number.');
  c.aiDescription = 'Set by hand';
  assert.equal(c.aiDescription, 'Set by hand');
  c.dispose();
});
