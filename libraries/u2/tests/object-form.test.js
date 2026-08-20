/* Schema-driven form tests. datagrok-api cannot load in node, so the properties here are minimal
   fakes of the `PropertyLike` shape the generator consumes — a real DG.Property satisfies it
   structurally (checked at compile time, see src/dg/README.md). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {propertyForm, objectForm, inputForProperty, PlatformInputs} from '../src/dg/object-form.js';
import {QNum} from '../src/core/qnum.js';

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
  prop('payload', 'graphics'),
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

form('min and max flag the generated number input instead of clamping it', () => {
  const source = sample();
  const f = mount(propertyForm(PROPS, source));
  const count = editorOf(f, 'count');
  count.value = '42';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(source.count, 42, 'the out-of-range value reaches the target');
  assert.equal(f.input('count').validity.value, 'Value must be at most 10');
  assert.equal(f.validity.value, 'Value must be at most 10', 'and the form is invalid');

  count.value = '-1';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(source.count, -1);
  assert.equal(f.input('count').validity.value, 'Value must be at least 0');

  count.value = '5';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(f.validity.value, null);

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

function clickers(input) {
  return input.root.querySelectorAll('.u2-number-click');
}

function sliderOf(input) {
  return input.root.querySelector('.u2-number-slider');
}

form('number metadata rides the property into the generated editor', () => {
  const props = [
    prop('count', 'int', {min: 1, max: 10}),
    prop('free', 'int'),
    prop('dose', 'double', {min: 0, max: 1000, units: 'mg'}),
    prop('ratio', 'double'),
    prop('manual', 'int', {min: 1, max: 10, showPlusMinus: false, showSlider: true}),
  ];
  const f = mount(propertyForm(props, {count: 3, free: 1, dose: 250, ratio: 0.5, manual: 2}));
  assert.equal(clickers(inputOf(f, 'count')).length, 2, 'a bounded int gets the platform default clicker');
  assert.equal(clickers(inputOf(f, 'free')).length, 0, 'an unbounded one does not');
  assert.notEqual(sliderOf(inputOf(f, 'dose')), null, 'a bounded float gets the slider');
  assert.equal(inputOf(f, 'dose').root.querySelector('.u2-input-postfix').textContent, 'mg');
  assert.equal(sliderOf(inputOf(f, 'ratio')), null, 'no bounds, no slider');
  assert.equal(clickers(inputOf(f, 'manual')).length, 0, 'showPlusMinus: false wins over the default');
  assert.notEqual(sliderOf(inputOf(f, 'manual')), null, 'showSlider: true reaches an int too');
  f.dispose();
});

form('a format string goes through DG.format wherever the platform is loaded', () => {
  const props = [prop('dose', 'double', {format: '#0.00', step: 0.1})];
  try {
    globalThis.DG = {format: (value, format) => `${value} as ${format}`};
    const formatted = mount(propertyForm(props, {dose: 2.5}));
    assert.equal(editorOf(formatted, 'dose').value, '2.5 as #0.00');
    formatted.dispose();
  } finally {
    delete globalThis.DG;
  }
  const plain = mount(propertyForm(props, {dose: 2.5}));
  assert.equal(editorOf(plain, 'dose').value, '2.5', 'off-platform the input keeps its own precision');
  plain.dispose();
});

form('inputForProperty is that same mapping, for editors that learn the property late', () => {
  const unbound = mount(inputForProperty(null));
  assert.equal(unbound.root.dataset.u2, 'text-input', 'nothing bound, nothing to map');
  unbound.dispose();

  const series = mount(inputForProperty(prop('series', 'string', {choices: ['a', 'b'], nullable: false})));
  assert.equal(series.root.dataset.u2, 'choice-input', 'choices win over the type, as they do in Dart');
  assert.equal(series.nullable, false);
  series.dispose();

  const count = mount(inputForProperty(prop('count', 'int', {min: 1, max: 10})));
  assert.equal(count.root.dataset.u2, 'number-input');
  assert.equal(clickers(count).length, 2);
  count.value.value = 11;
  assert.equal(count.validity.value, 'Value must be at most 10');
  count.dispose();

  const labelled = mount(inputForProperty(prop('count', 'int'), {label: 'Count'}));
  assert.equal(labelled.root.querySelector('.u2-input-label').textContent, 'Count', 'options still win');
  labelled.dispose();
});

form('assumeWritable: a property with no setter still gets the editor its type asks for', () => {
  const param = {name: 'count', propertyType: 'int', min: 1, max: 10};
  const readonly = mount(inputForProperty(param));
  assert.equal(readonly.root.dataset.u2, 'text-input', 'nothing can be written back — a form shows it');
  assert.equal(readonly.enabled, false);
  readonly.dispose();

  const editor = mount(inputForProperty(param, {assumeWritable: true}));
  assert.equal(editor.root.dataset.u2, 'number-input', 'a value editor owns the value itself');
  assert.equal(editor.enabled, true);
  assert.equal(clickers(editor).length, 2, 'and gets the metadata with it');
  editor.dispose();
});

form('bigint and qnum get their own editors, not a number input', () => {
  const source = {id: 123456789012345678901234567890n, ic50: QNum.less(5.2)};
  const f = mount(propertyForm([prop('id', 'bigint'), prop('ic50', 'qnum')], source));
  assert.deepEqual(kinds(f), ['bigint-input', 'qnum-input']);
  assert.equal(editorOf(f, 'id').value, '123456789012345678901234567890');
  assert.equal(editorOf(f, 'ic50').value.startsWith('<5.2'), true);

  const id = inputOf(f, 'id');
  id.value.value = 42n;
  assert.equal(source.id, 42n, 'a bigint reaches the property as a bigint');

  const ic50 = inputOf(f, 'ic50');
  ic50.value.value = QNum.greater(7);
  assert.equal(QNum.getValue(source.ic50), 7, 'a qnum reaches it as the packed double');
  f.dispose();
});

form('a bigint property that hands over its digits as text is read all the same', () => {
  const f = mount(propertyForm([prop('id', 'bigint')], {id: '900719925474099101'}));
  assert.equal(inputOf(f, 'id').value.value, 900719925474099101n);
  f.dispose();

  const empty = mount(propertyForm([prop('id', 'bigint')], {}));
  assert.equal(inputOf(empty, 'id').value.value, null);
  empty.dispose();
});

form('a postfix override wins over the property units', () => {
  const dose = prop('dose', 'double', {units: 'mg'});
  const f = mount(propertyForm([dose], {dose: 1}, {overrides: {dose: {postfix: 'mg/kg'}}}));
  assert.equal(inputOf(f, 'dose').root.querySelector('.u2-input-postfix').textContent, 'mg/kg');
  f.dispose();

  const plain = mount(propertyForm([dose], {dose: 1}));
  assert.equal(inputOf(plain, 'dose').root.querySelector('.u2-input-postfix').textContent, 'mg');
  plain.dispose();
});

/** Every route the property metadata can ask for: `data-u2` names the editor that was built. */
function routes(cases) {
  for (const [property, kind, check] of cases) {
    const input = mount(inputForProperty(property));
    assert.equal(input.root.dataset.u2, kind, property.name);
    if (check)
      check(input);
    input.dispose();
  }
}

form('inputType picks the editor, ahead of every other hint', () => {
  routes([
    [prop('style', 'string', {inputType: 'Radio', choices: ['a', 'b']}), 'radio-input',
      (i) => assert.equal(i.root.querySelectorAll('input[type=radio]').length, 2)],
    [prop('shade', 'string', {inputType: 'Color'}), 'color-input'],
    [prop('face', 'string', {inputType: 'Font'}), 'font-input'],
    [prop('logo', 'string', {inputType: 'Image'}), 'image-input'],
    [prop('labels', 'list', {inputType: 'Tags', choices: ['a', 'b']}), 'tags-input'],
    [prop('note', 'string', {inputType: 'TextArea', editor: 'password'}), 'text-area'],
    [prop('grade', 'int', {inputType: 'Slider', min: 1, max: 5}), 'slider-input'],
    [prop('series', 'string', {inputType: 'Choice', choices: ['a'], editor: 'textarea'}), 'choice-input'],
  ]);
});

form('the editor hint picks it next, ahead of the choices and the type', () => {
  routes([
    [prop('note', 'string', {editor: 'textarea', choices: ['a', 'b']}), 'text-area'],
    // the platform's own `description` option carries the InputType spelling (property.ts:326)
    [prop('summary', 'string', {editor: 'TextArea'}), 'text-area'],
    [prop('secret', 'string', {editor: 'password'}), 'text-input',
      (i) => assert.equal(i.root.querySelector('input').type, 'password')],
    [prop('active', 'bool', {editor: 'switch'}), 'bool-input',
      (i) => assert.notEqual(i.root.querySelector('.u2-input-switch'), null)],
    [prop('dose', 'double', {editor: 'slider', min: 0, max: 10, step: 0.5}), 'slider-input',
      (i) => assert.deepEqual([i.root.querySelector('input').min, i.root.querySelector('input').max,
        i.root.querySelector('input').step], ['0', '10', '0.5'])],
  ]);
});

form('an editor hint the type does not accept is ignored, as the platform ignores it', () => {
  routes([
    // textarea and password are reached under `pt == Types.STRING` (input_base.dart:725-728)
    [prop('count', 'int', {editor: 'textarea'}), 'number-input'],
    [prop('count', 'int', {editor: 'password'}), 'number-input'],
    [prop('when', 'datetime', {editor: 'textarea'}), 'datetime-input'],
    // switch under `pt == Types.BOOL` (:702)
    [prop('count', 'int', {editor: 'switch'}), 'number-input'],
    [prop('label', 'string', {editor: 'switch'}), 'text-input'],
    // the slider (:705) is reached by everything the two bool branches above it did not take
    [prop('active', 'bool', {editor: 'slider'}), 'bool-input',
      (i) => assert.equal(i.root.querySelector('.u2-input-switch'), null)],
    [prop('grade', 'int', {editor: 'slider', min: 1, max: 5}), 'slider-input'],
  ]);
});

form('list, map and file route off the type, choices splitting the two list editors', () => {
  routes([
    [prop('labels', 'list'), 'list-input'],
    [prop('labels', 'list', {choices: ['a', 'b']}), 'multi-choice-input'],
    [prop('meta', 'map'), 'map-input'],
  ]);
});

form('a file property goes to the editor the platform layer registered', () => {
  const built = [];
  const unregister = PlatformInputs.register('file', (property, options) => {
    built.push(property.name);
    return new TextInput({...options, name: 'stub-file'});
  });
  try {
    const input = mount(inputForProperty(prop('report', 'file')));
    assert.deepEqual(built, ['report'], 'the property reaches the registered factory');
    assert.equal(input.name, 'stub-file');
    input.dispose();
  } finally {
    unregister();
  }

  const unwired = mount(inputForProperty(prop('report', 'file')));
  assert.equal(unwired.root.dataset.u2, 'text-input', 'off the platform, a read-only field');
  assert.equal(unwired.enabled, false);
  unwired.dispose();
});

form('lists and maps round-trip through the property without reporting a change', () => {
  const source = {labels: ['a', 'b'], meta: {x: '1'}, report: {path: 'demo.csv'}};
  const changes = [];
  const props = [prop('labels', 'list'), prop('meta', 'map'), prop('report', 'file')];
  const f = mount(propertyForm(props, source, {onChanged: (name, value) => changes.push([name, value])}));
  assert.deepEqual(kinds(f), ['list-input', 'map-input', 'text-input']);
  assert.deepEqual(f.getValues(), {labels: ['a', 'b'], meta: {x: '1'}, report: {path: 'demo.csv'}});
  assert.deepEqual(changes, [], 'generation is not a change, arrays and records included');

  inputOf(f, 'labels').value.value = ['c'];
  inputOf(f, 'meta').value.value = {y: '2'};
  assert.deepEqual(source.labels, ['c']);
  assert.deepEqual(source.meta, {y: '2'});
  assert.deepEqual(changes, [['labels', ['c']], ['meta', {y: '2'}]]);

  source.labels = ['d'];
  f.refresh();
  assert.deepEqual(inputOf(f, 'labels').value.value, ['d']);
  assert.deepEqual(changes.length, 2, 'refresh is echo-suppressed');
  f.dispose();
});

form('a list property that hands over its items as text is read all the same', () => {
  const f = mount(propertyForm([prop('labels', 'list')], {labels: 'a,"b,c"'}));
  assert.deepEqual(inputOf(f, 'labels').value.value, ['a', 'b,c']);
  f.dispose();

  const empty = mount(propertyForm([prop('labels', 'list'), prop('meta', 'map')], {}));
  assert.deepEqual(empty.getValues(), {labels: [], meta: {}});
  empty.dispose();
});
