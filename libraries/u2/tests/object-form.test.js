/* Schema-driven form tests. datagrok-api cannot load in node, so the properties here are minimal
   fakes of the `PropertyLike` shape the generator consumes — a real DG.Property satisfies it
   structurally (checked at compile time, see src/dg/README.md). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextInput} from '../src/components/text-input.js';
import {propertyForm, objectForm} from '../src/dg/object-form.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function form(name, body) {
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

/** A `PropertyLike` over a plain field, as `DG.Property.js` builds one. */
function prop(name, propertyType, options = {}) {
  return {
    name, propertyType,
    get: (x) => x[name],
    set: (x, v) => x[name] = v,
    ...options,
  };
}

function mount(component) {
  document.body.append(component.root);
  return component;
}

function inputOf(f, name) {
  return f.inputs.find((i) => i.name === name);
}

function editorOf(f, name, selector = 'input') {
  return inputOf(f, name).root.querySelector(selector);
}

function kinds(f) {
  return f.inputs.map((i) => i.root.dataset.u2);
}

const PROPS = [
  prop('name', 'string', {caption: 'Compound', description: 'Registered name'}),
  prop('count', 'int', {min: 0, max: 10}),
  prop('weight', 'double'),
  prop('active', 'bool'),
  prop('series', 'string', {choices: ['a', 'b', 'c']}),
  prop('payload', 'map'),
];

function sample() {
  return {name: 'Aspirin', count: 3, weight: 1.5, active: true, series: 'b', payload: {x: 1}};
}

form('generates one input per property type, with a readonly fallback', () => {
  const f = mount(propertyForm(PROPS, sample()));
  assert.deepEqual(kinds(f), ['text-input', 'number-input', 'number-input', 'bool-input',
    'choice-input', 'text-input']);
  assert.equal(f.root.dataset.u2, 'object-form');

  const payload = inputOf(f, 'payload');
  assert.equal(payload.enabled, false, 'an uneditable type falls back to a disabled text input');
  assert.equal(payload.root.querySelector('input').disabled, true);
  assert.equal(payload.value.value, String({x: 1}));

  const series = editorOf(f, 'series', 'select');
  assert.deepEqual(series.children.map((o) => o.value), ['', 'a', 'b', 'c']);
  f.dispose();
});

form('choices win over the declared type, and a property with no setter is readonly', () => {
  const numeric = prop('grade', 'int', {choices: ['1', '2']});
  const computed = {name: 'total', propertyType: 'int', get: (x) => x.count * 2};
  const f = mount(propertyForm([numeric, computed], {grade: '2', count: 3}));
  assert.deepEqual(kinds(f), ['choice-input', 'text-input']);
  assert.equal(inputOf(f, 'grade').value.value, '2');
  assert.equal(inputOf(f, 'total').enabled, false);
  assert.equal(inputOf(f, 'total').value.value, '6');
  f.dispose();
});

form('maps caption, name and description, and reads the initial values', () => {
  const source = sample();
  const f = mount(propertyForm(PROPS, source));
  const name = inputOf(f, 'name');
  assert.equal(name.root.querySelector('.u2-input-label').textContent, 'Compound');
  assert.equal(name.root.title, 'Registered name');
  assert.equal(inputOf(f, 'count').root.querySelector('.u2-input-label').textContent, 'count');

  assert.deepEqual(f.getValues(), {name: 'Aspirin', count: 3, weight: 1.5, active: true,
    series: 'b', payload: String({x: 1})});
  assert.equal(editorOf(f, 'name').value, 'Aspirin');
  assert.equal(editorOf(f, 'active').checked, true);
  assert.equal(source.name, 'Aspirin', 'generation does not write back');
  f.dispose();
});

form('null and missing values fall back to the empty value of each editor', () => {
  const f = mount(propertyForm(PROPS, {name: null, series: null}));
  assert.deepEqual(f.getValues(), {name: '', count: null, weight: null, active: false,
    series: null, payload: ''});
  f.dispose();
});

form('an edit writes through prop.set and reports it once', () => {
  const source = sample();
  const changes = [];
  const f = mount(propertyForm(PROPS, source, {onChanged: (name, value) => changes.push([name, value])}));
  assert.deepEqual(changes, [], 'generation is not a change');

  const name = editorOf(f, 'name');
  name.value = 'Ibuprofen';
  fire(name, 'input');
  assert.equal(source.name, 'Ibuprofen');

  const active = editorOf(f, 'active');
  active.checked = false;
  fire(active, 'change');
  assert.equal(source.active, false);

  const series = editorOf(f, 'series', 'select');
  series.value = 'c';
  fire(series, 'change');
  assert.equal(source.series, 'c');

  assert.deepEqual(changes, [['name', 'Ibuprofen'], ['active', false], ['series', 'c']]);
  f.dispose();
});

form('refresh() re-reads the target without writing back or reporting', () => {
  const source = sample();
  const changes = [];
  const f = mount(propertyForm(PROPS, source, {onChanged: (name, value) => changes.push([name, value])}));

  source.name = 'Paracetamol';
  source.count = 7;
  source.active = false;
  assert.equal(inputOf(f, 'name').value.value, 'Aspirin', 'the form does not poll');

  f.refresh();
  assert.equal(inputOf(f, 'name').value.value, 'Paracetamol');
  assert.equal(editorOf(f, 'name').value, 'Paracetamol');
  assert.equal(inputOf(f, 'count').value.value, 7);
  assert.equal(inputOf(f, 'active').value.value, false);
  assert.deepEqual(changes, [], 'refresh is echo-suppressed');
  assert.equal(source.name, 'Paracetamol');
  f.dispose();
});

form('nullable: false becomes a required validator and drops the empty choice', () => {
  const props = [
    prop('name', 'string', {nullable: false}),
    prop('series', 'string', {choices: ['a', 'b'], nullable: false}),
  ];
  const f = mount(propertyForm(props, {name: '', series: null}));
  assert.equal(inputOf(f, 'name').validity.value, 'Value can\'t be empty');
  assert.equal(f.validity.value, 'Value can\'t be empty');
  assert.equal(f.validate(), false);
  assert.deepEqual(editorOf(f, 'series', 'select').children.map((o) => o.value), ['a', 'b']);

  const name = editorOf(f, 'name');
  name.value = 'Aspirin';
  fire(name, 'input');
  assert.equal(inputOf(f, 'name').validity.value, null);
  assert.equal(f.validity.value, 'Value can\'t be empty', 'the choice is still empty');

  const series = editorOf(f, 'series', 'select');
  series.value = 'b';
  fire(series, 'change');
  assert.equal(f.validity.value, null);
  assert.equal(f.validate(), true);
  f.dispose();
});

form('min and max clamp the generated number input', () => {
  const source = sample();
  const f = mount(propertyForm(PROPS, source));
  const count = editorOf(f, 'count');
  count.value = '42';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(source.count, 10);

  count.value = '-1';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(source.count, 0);

  const weight = editorOf(f, 'weight');
  weight.value = '-3.5';
  fire(weight, 'input');
  fire(weight, 'blur');
  assert.equal(source.weight, -3.5, 'a property with no bounds is not clamped');
  f.dispose();
});

form('include picks and orders the fields, exclude drops them', () => {
  const source = sample();
  const picked = mount(propertyForm(PROPS, source, {include: ['active', 'name', 'missing']}));
  assert.deepEqual(picked.inputs.map((i) => i.name), ['active', 'name']);
  picked.dispose();

  const rest = mount(propertyForm(PROPS, source, {exclude: ['payload', 'weight']}));
  assert.deepEqual(rest.inputs.map((i) => i.name), ['name', 'count', 'active', 'series']);
  rest.dispose();

  const both = mount(propertyForm(PROPS, source, {include: ['name', 'count'], exclude: ['count']}));
  assert.deepEqual(both.inputs.map((i) => i.name), ['name']);
  both.dispose();
});

form('an override replaces the input, other keys merge into the generated one', () => {
  const source = sample();
  const custom = new TextInput({label: 'Structure'});
  const f = mount(propertyForm(PROPS, source, {
    condensed: true,
    overrides: {
      name: {input: custom},
      count: {label: 'How many', inline: true},
    },
  }));
  assert.equal(f.root.classList.contains('u2-form-condensed'), true);
  assert.equal(inputOf(f, 'Structure'), custom, 'the caller-supplied input is laid out as is');
  assert.equal(custom.value.value, 'Aspirin', 'and seeded from the property');
  assert.equal(custom.root.parentNode, f.root.querySelector('.u2-form-rows'));

  const editor = custom.root.querySelector('input');
  editor.value = 'Ibuprofen';
  fire(editor, 'input');
  assert.equal(source.name, 'Ibuprofen', 'a custom input writes through the property too');

  const count = inputOf(f, 'count');
  assert.equal(count.root.querySelector('.u2-input-label'), null, 'inline drops the label');
  assert.equal(count.root.classList.contains('u2-input-inline'), true);
  f.dispose();
  assert.equal(custom.scope.isDisposed, false, 'the form never owns what it did not create');
  custom.dispose();
});

form('input() finds by property name, properties lists them in layout order', () => {
  const custom = new TextInput({label: 'Structure'});
  const f = mount(propertyForm(PROPS, sample(), {overrides: {name: {input: custom}}}));
  assert.deepEqual(f.properties.map((p) => p.name),
    ['name', 'count', 'weight', 'active', 'series', 'payload']);
  assert.equal(f.input('name'), custom, 'found by property name even when the input is captioned differently');
  assert.equal(f.input('count'), inputOf(f, 'count'));
  assert.equal(f.input('missing'), undefined);

  f.input('count').addValidator((v) => v === 3 ? 'not three' : null);
  assert.equal(f.validity.value, 'not three');
  f.dispose();
  custom.dispose();
});

form('objectForm enumerates the properties of a self-describing object', () => {
  const source = {
    name: 'Scatter plot', showAxes: true,
    getProperties: () => [prop('name', 'string'), prop('showAxes', 'bool')],
  };
  const f = mount(objectForm(source, {exclude: ['getProperties']}));
  assert.deepEqual(f.inputs.map((i) => i.name), ['name', 'showAxes']);
  assert.equal(f.target, source);

  const axes = editorOf(f, 'showAxes');
  axes.checked = false;
  fire(axes, 'change');
  assert.equal(source.showAxes, false);
  f.dispose();
});

form('disposing the form disposes every input it generated', () => {
  const f = mount(propertyForm(PROPS, sample()));
  const inputs = f.inputs.slice();
  assert.equal(inputs.length, 6);
  f.dispose();
  for (const input of inputs)
    assert.equal(input.scope.isDisposed, true, `${input.name} disposed`);
});
