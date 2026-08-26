/* The W1 FuncCallForm contract over the FuncCall doubles: routing, captions, categories,
   validation, literal defaults, and the two-way echo-suppressed binding. In an isolated copy
   where the module under test is not built, every test skips. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Func, FuncParam} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);

const mod = await import('../src/dg/funcs/func-form.js').then((m) => m, (e) => {
  if (e.code !== 'ERR_MODULE_NOT_FOUND')
    throw e;
  return null;
});
const skip = mod === null && 'src/dg/funcs/func-form.js not built in this copy; runs at the merge gate';
const funcForm = mod?.funcForm;

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function form(name, body) {
  test(name, {skip}, async () => {
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

function fp(name, propertyType, options = {}) {
  return new FuncParam(name, propertyType, options);
}

function callOf(...params) {
  return new Func('fceDemo', {inputs: params}).prepare({});
}

/** Counts every write and every report, so the two-way tests read plain ledgers. */
function wired(call, options = {}) {
  const sets = [];
  const setParamValue = call.setParamValue.bind(call);
  call.setParamValue = (name, value) => {
    sets.push([name, value]);
    setParamValue(name, value);
  };
  const changes = [];
  const f = mount(funcForm(call, {...options, onInputChanged: (name, value) => changes.push([name, value])}));
  return {f, sets, changes};
}

function mount(f) {
  document.body.append(f.root);
  return f;
}

function editorOf(f, name, selector = 'input') {
  return f.getInput(name).root.querySelector(selector);
}

function kinds(f) {
  return f.inputs.map((i) => i.root.dataset.u2);
}

function labelOf(f, name) {
  return f.getInput(name).root.querySelector('.u2-input-label')?.textContent;
}

function headers(f) {
  return [...f.root.querySelectorAll('.u2-form-category')].map((h) => h.textContent);
}

form('the W1 scalar set routes to the type editors, choices to a choice input', () => {
  const call = callOf(
    fp('name', 'string'),
    fp('count', 'int', {min: 0, max: 10}),
    fp('dose', 'double'),
    fp('active', 'bool'),
    fp('started', 'datetime'),
    fp('stage', 'string', {choices: ['ideation', 'running']}));
  const f = mount(funcForm(call));
  assert.deepEqual(kinds(f), ['text-input', 'number-input', 'number-input', 'bool-input',
    'datetime-input', 'choice-input']);
  assert.deepEqual(f.unsupported, []);
  assert.ok(f.inputs.every((i) => i.enabled !== false), 'a func param is editable despite carrying no set');
  f.dispose();
});

form('inputType picks the editor: Radio, TextArea, Color, Slider', () => {
  const call = callOf(
    fp('route', 'string', {inputType: 'Radio', choices: ['oral', 'iv']}),
    fp('notes', 'string', {inputType: 'TextArea'}),
    fp('shade', 'string', {inputType: 'Color'}),
    fp('progress', 'int', {inputType: 'Slider', min: 0, max: 100}));
  const f = mount(funcForm(call));
  assert.deepEqual(kinds(f), ['radio-input', 'text-area', 'color-input', 'slider-input']);
  assert.equal(f.getInput('route').root.querySelectorAll('input[type=radio]').length, 2);
  assert.deepEqual(f.unsupported, []);
  f.dispose();
});

form('a SERVICE_PARAM (editor: none) is skipped entirely — no field, not unsupported', () => {
  const call = callOf(fp('name', 'string'), fp('context', 'string', {options: {editor: 'none'}}));
  const f = mount(funcForm(call));
  assert.deepEqual(f.inputs.map((i) => i.name), ['name']);
  assert.equal(f.getInput('context'), undefined);
  assert.deepEqual(f.unsupported, []);
  f.dispose();
});

form('editorParam, layout, foreign editors and non-scalar types get no field and are named unsupported', () => {
  const call = callOf(
    fp('helper', 'string', {options: {editorParam: 'name'}}),
    fp('panel', 'string', {options: {editor: 'layout'}}),
    fp('molecule', 'string', {options: {editor: 'Chem:sketcher'}}),
    fp('table', 'dataframe'),
    fp('name', 'string'));
  const f = mount(funcForm(call));
  assert.deepEqual(f.inputs.map((i) => i.name), ['name']);
  assert.deepEqual([...f.unsupported], ['helper', 'panel', 'molecule', 'table']);
  f.dispose();
});

form('the four ib editor hints route to their editors — the documented W1 divergence', () => {
  const call = callOf(
    fp('notes', 'string', {options: {editor: 'textarea'}}),
    fp('memo', 'string', {options: {editor: 'TextArea'}}),
    fp('secret', 'string', {options: {editor: 'password'}}),
    fp('archived', 'bool', {options: {editor: 'switch'}}),
    fp('dose', 'double', {min: 0, max: 10, options: {editor: 'slider'}}));
  const f = mount(funcForm(call));
  assert.deepEqual(kinds(f), ['text-area', 'text-area', 'text-input', 'bool-input', 'slider-input']);
  assert.equal(editorOf(f, 'secret').type, 'password');
  assert.notEqual(f.getInput('archived').root.querySelector('.u2-input-switch'), null);
  assert.deepEqual(f.unsupported, []);
  f.dispose();
});

form('an editor hint the property type does not accept is ignored, and the type routes', () => {
  const call = callOf(
    fp('count', 'int', {options: {editor: 'textarea'}}),
    fp('label', 'string', {options: {editor: 'switch'}}));
  const f = mount(funcForm(call));
  assert.deepEqual(kinds(f), ['number-input', 'text-input']);
  assert.deepEqual(f.unsupported, []);
  f.dispose();
});

form('captions: camelCaseToWords over the friendlyName when set, else over the name', () => {
  const call = callOf(
    fp('rowCount', 'int'),
    fp('dose', 'double', {friendlyName: 'Dose (mg)'}),
    fp('stage', 'string', {friendlyName: 'stage'}),
    fp('dl', 'double', {friendlyName: 'doseLevel'}));
  const f = mount(funcForm(call));
  assert.equal(labelOf(f, 'rowCount'), 'Row Count');
  assert.equal(labelOf(f, 'dose'), 'Dose (mg)');
  assert.equal(labelOf(f, 'stage'), 'Stage');
  assert.equal(labelOf(f, 'dl'), 'Dose Level', 'a camelCase friendlyName is split too (ib:679)');
  assert.notEqual(f.getInput('rowCount'), undefined, 'inputs are keyed by param name');
  assert.equal(f.getInput('Row Count'), undefined, 'never by caption');
  f.dispose();
});

form('units come from the units field or the options map, as InputBase.getOption reads them', () => {
  const call = callOf(
    fp('dose', 'double', {options: {units: 'mg'}}),
    fp('volume', 'double', {units: 'ml'}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('dose').root.querySelector('.u2-input-postfix')?.textContent, 'mg');
  assert.equal(f.getInput('volume').root.querySelector('.u2-input-postfix')?.textContent, 'ml');
  f.dispose();
});

form('description becomes the tooltip on the input root', () => {
  const call = callOf(fp('dose', 'double', {description: 'Dose in milligrams'}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('dose').root.title, 'Dose in milligrams');
  f.dispose();
});

form('category headers render only when >1 category or the single one is not Misc', () => {
  const flat = mount(funcForm(callOf(fp('a', 'string'), fp('b', 'int'))));
  assert.deepEqual(headers(flat), [], 'one implicit Misc category — no headers');
  flat.dispose();

  const single = mount(funcForm(callOf(fp('a', 'string', {category: 'Advanced'}))));
  assert.deepEqual(headers(single), ['Advanced']);
  single.dispose();

  const grouped = mount(funcForm(callOf(
    fp('a', 'string'),
    fp('b', 'int', {category: 'Advanced'}),
    fp('c', 'bool'))));
  assert.deepEqual(headers(grouped), ['Misc', 'Advanced'], 'first-appearance order');
  assert.deepEqual(grouped.inputs.map((i) => i.name), ['a', 'c', 'b'],
    'grouped by category, declaration order within each');
  grouped.dispose();
});

form('a nullable choice gets the empty option; a non-nullable one does not, and requires a value', () => {
  const call = callOf(
    fp('stage', 'string', {choices: ['a', 'b'], nullable: true}),
    fp('route', 'string', {choices: ['oral', 'iv'], nullable: false}));
  const f = mount(funcForm(call));
  assert.deepEqual(editorOf(f, 'stage', 'select').children.map((o) => o.value), ['', 'a', 'b']);
  assert.deepEqual(editorOf(f, 'route', 'select').children.map((o) => o.value), ['oral', 'iv']);
  f.dispose();
});

form('min and max flag the input instead of clamping it — the value still reaches the call', () => {
  const call = callOf(fp('count', 'int', {min: 0, max: 10}));
  const f = mount(funcForm(call));
  const count = editorOf(f, 'count');
  count.value = '42';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(call.inputs.count, 42, 'not clamped');
  assert.equal(f.getInput('count').validity.value, 'Value must be at most 10');
  assert.equal(f.isValid, false);
  assert.equal(f.validateInputs(), false);

  count.value = '5';
  fire(count, 'input');
  fire(count, 'blur');
  assert.equal(f.isValid, true);
  assert.equal(f.validateInputs(), true);
  f.dispose();
});

form('nullable: false becomes a required validator', () => {
  const call = callOf(fp('name', 'string', {nullable: false}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('name').validity.value, 'Value can\'t be empty');
  assert.equal(f.isValid, false);

  const name = editorOf(f, 'name');
  name.value = 'Aspirin';
  fire(name, 'input');
  assert.equal(f.isValid, true);
  f.dispose();
});

form('a literal default is displayed without ever being written into the call', () => {
  const call = callOf(
    fp('count', 'int', {defaultValue: 3}),
    fp('stage', 'string', {isOptional: true, defaultValue: 'a'}));
  const {f, sets, changes} = wired(call);
  assert.equal(f.getInput('count').value.value, 3, 'seeded from defaultValue');
  assert.equal(f.getInput('stage').value.value, 'a', 'an optional param answers it through value already');
  assert.deepEqual([sets, changes], [[], []], 'display-only: no setParamValue, no onInputChanged');
  assert.equal(call.dart.params[0].dart.value, undefined);
  assert.equal(call.dart.params[1].dart.value, undefined);
  f.dispose();
});

form('an input edit is exactly one setParamValue and one onInputChanged, and the echo stops there', async () => {
  const call = callOf(fp('name', 'string'));
  const {f, sets, changes} = wired(call);
  assert.deepEqual([sets, changes], [[], []], 'generation is not a change');

  const name = editorOf(f, 'name');
  name.value = 'Aspirin';
  fire(name, 'input');
  await flush();
  assert.deepEqual(sets, [['name', 'Aspirin']]);
  assert.deepEqual(changes, [['name', 'Aspirin']]);
  assert.equal(call.inputs.name, 'Aspirin', 'the write went through the call');
  f.dispose();
});

form('a param fire refreshes the field with zero onInputChanged, and writes nothing back', async () => {
  const call = callOf(fp('name', 'string'));
  const {f, sets, changes} = wired(call);
  call.setParamValue('name', 'Ibuprofen');
  await flush();
  assert.equal(f.getInput('name').value.value, 'Ibuprofen', 'the field follows the call');
  assert.equal(editorOf(f, 'name').value, 'Ibuprofen');
  assert.deepEqual(changes, [], 'a param write is not an input change');
  assert.deepEqual(sets, [['name', 'Ibuprofen']], 'only the external write itself — no echo write-back');
  f.dispose();
});

form('a param fire with an unchanged value — the Dart List always-fire — writes nothing back', async () => {
  const call = callOf(fp('name', 'string'));
  const {f, sets, changes} = wired(call);
  call.setParamValue('name', 'Aspirin');
  await flush();
  sets.length = 0;
  const [p] = call.dart.params;
  p.dart.onChanged.fire(p);
  await flush();
  assert.deepEqual(sets, [], 'zero setParamValue calls — the echo terminates at the value comparison');
  assert.deepEqual(changes, []);
  assert.equal(f.getInput('name').value.value, 'Aspirin');
  f.dispose();
});

form('a datetime edit crosses the write boundary as dayjs when the platform global is there', async () => {
  globalThis.dayjs = (d) => ({$isDayjsObject: true, ms: d.getTime(), toDate: () => d});
  try {
    const call = callOf(fp('started', 'datetime'));
    const {f, sets} = wired(call);
    const date = new Date(2026, 7, 25, 10, 30);
    f.getInput('started').value.value = date;
    await flush();
    assert.equal(sets.length, 1);
    const [name, sent] = sets[0];
    assert.equal(name, 'started');
    assert.equal(sent.$isDayjsObject, true, 'the Date went out wrapped as dayjs');
    assert.equal(sent.ms, date.getTime());
    assert.equal(call.dart.params[0].dart.value, sent);
    assert.equal(f.getInput('started').value.value, date, 'the echo settles — coerce reads it back as the same Date');
    f.dispose();
  } finally {
    delete globalThis.dayjs;
  }
});

form('without the platform global a datetime edit passes the Date through', async () => {
  const call = callOf(fp('started', 'datetime'));
  const {f, sets} = wired(call);
  const date = new Date(2026, 0, 2);
  f.getInput('started').value.value = date;
  await flush();
  assert.deepEqual(sets, [['started', date]]);
  f.dispose();
});

form('replacing source rebinds: old subscriptions die, fields refresh and edit the new call', () => {
  const a = callOf(fp('name', 'string'));
  const b = callOf(fp('name', 'string'));
  a.setParamValue('name', 'one');
  b.setParamValue('name', 'two');
  const f = mount(funcForm(a));
  assert.equal(f.source, a);
  const [pa] = a.dart.params;
  const [pb] = b.dart.params;
  assert.equal(pa.onChanged.count, 1, 'one subscription per param — never accumulated');
  assert.equal(pb.onChanged.count, 0);

  f.source = b;
  assert.equal(f.source, b);
  assert.equal(pa.onChanged.count, 0, 'the old call\'s subscriptions are dead');
  assert.equal(pb.onChanged.count, 1, 'the new call is live');
  assert.equal(f.getInput('name').value.value, 'two', 'refreshed from the new call');

  a.setParamValue('name', 'stale');
  assert.equal(f.getInput('name').value.value, 'two', 'the old call no longer reaches the form');

  const name = editorOf(f, 'name');
  name.value = 'three';
  fire(name, 'input');
  assert.equal(b.inputs.name, 'three');
  assert.equal(a.inputs.name, 'stale', 'an edit lands in the new call only');
  f.dispose();
});

form('a field whose param the new call lacks goes dormant: disabled, unsubscribed, no writes', () => {
  const a = callOf(fp('name', 'string'), fp('count', 'int'));
  const b = callOf(fp('name', 'string'), fp('table', 'dataframe'));
  const f = mount(funcForm(a));
  const [, pCount] = a.dart.params;
  assert.equal(pCount.onChanged.count, 1);

  f.source = b;
  assert.equal(f.getInput('count').enabled, false, 'the stale field is disabled');
  assert.equal(pCount.onChanged.count, 0, 'its old-call subscription is dropped');
  assert.deepEqual([...f.unsupported], ['table'], 'unsupported recomputed over the new call');

  const count = editorOf(f, 'count');
  count.value = '7';
  fire(count, 'input');
  fire(count, 'blur');
  assert.deepEqual(b.inputs, {name: undefined, table: undefined}, 'the dead call took no write');

  f.source = a;
  assert.equal(f.getInput('count').enabled, true, 'a call that carries the param again revives the field');
  assert.equal(pCount.onChanged.count, 1);
  assert.deepEqual([...f.unsupported], []);
  const revived = editorOf(f, 'count');
  revived.value = '7';
  fire(revived, 'input');
  fire(revived, 'blur');
  assert.equal(a.inputs.count, 7);
  f.dispose();
});

form('twoWayBinding: false is one-way — edits reach the call, the call does not reach the form', () => {
  const call = callOf(fp('name', 'string'));
  const {f, changes} = wired(call, {twoWayBinding: false});
  call.setParamValue('name', 'Ibuprofen');
  assert.notEqual(f.getInput('name').value.value, 'Ibuprofen', 'no refresh from the call');

  const name = editorOf(f, 'name');
  name.value = 'Aspirin';
  fire(name, 'input');
  assert.equal(call.inputs.name, 'Aspirin', 'the write path stays');
  assert.deepEqual(changes, [['name', 'Aspirin']]);
  f.dispose();
});

form('disposal releases every param subscription', () => {
  const call = callOf(fp('name', 'string'), fp('count', 'int'));
  const f = mount(funcForm(call));
  assert.deepEqual(call.dart.params.map((p) => p.onChanged.count), [1, 1]);
  f.dispose();
  assert.deepEqual(call.dart.params.map((p) => p.onChanged.count), [0, 0]);
});
