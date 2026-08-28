/* The W1 FuncCallForm contract over the FuncCall doubles: routing, captions, categories,
   validation, literal defaults, and the two-way echo-suppressed binding. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Func, FuncParam} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);

// dynamic so the dg-stub loader hook above is registered first; a missing build is a hard failure
const {funcForm} = await import('../src/dg/funcs/func-form.js');
const {ParamState} = await import('../src/dg/funcs/param-sources.js');

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
    fp('frames', 'dataframe_list'),
    fp('name', 'string'));
  const f = mount(funcForm(call));
  assert.deepEqual(f.inputs.map((i) => i.name), ['name']);
  assert.deepEqual([...f.unsupported], ['helper', 'panel', 'molecule', 'frames']);
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

form('a literal default the call does not hold is never displayed — display == call, run included', () => {
  const call = callOf(
    fp('count', 'int', {defaultValue: 3}),
    fp('stage', 'string', {isOptional: true, defaultValue: 'a'}));
  const {f, sets, changes} = wired(call);
  assert.equal(f.getInput('count').value.value, null, 'a defaulted nullable param shows what the run will use — nothing');
  assert.equal(f.getInput('stage').value.value, 'a', 'an optional param answers it through value, so the run agrees');
  assert.deepEqual([sets, changes], [[], []], 'no setParamValue, no onInputChanged');
  assert.equal(call.dart.params[0].dart.value, undefined);
  assert.equal(call.dart.params[1].dart.value, undefined);
  f.dispose();
});

form('a source rebind seeds fields from the new call, never from defaultValue', () => {
  const {f} = wired(callOf(fp('count', 'int', {defaultValue: 3})));
  const filled = callOf(fp('count', 'int', {defaultValue: 3}));
  filled.setParamValue('count', 7);
  f.source = filled;
  assert.equal(f.getInput('count').value.value, 7, 'refreshed from the new call');
  f.source = callOf(fp('count', 'int', {defaultValue: 3}));
  assert.equal(f.getInput('count').value.value, null, 'display == call after a rebind too');
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

/* ------------------------------------------------------------------------------------------------
   W2: async choices / suggestions / computed defaults (plan-w2.md FP-W2-2…FP-W2-9). The doubles
   deliver every eval on a later tick; nothing below may rely on synchronous delivery. Timing:
   dep-triggered re-runs ride a 200 ms debounce, the initial run is immediate — states are polled
   with `until`, fixed sleeps appear only where a NEGATIVE outcome needs a bounded grace period. */

import {TextInput} from '../src/components/inputs/text-input.js';

/** The W1 `form` harness plus provider-registry hygiene for the W2 suites. */
function w2(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      Func.registry = [];
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** Doubles-only pins. */
function dbl(name, body) {
  test(name, async () => {
    try {
      await body();
    } finally {
      Func.registry = [];
      resetDom();
      await flush();
    }
  });
}

function provider(name, inputs, run) {
  const f = new Func(name, {inputs, run});
  Func.registry.push(f);
  return f;
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function until(cond, what, timeoutMs = 1500) {
  const start = Date.now();
  while (!cond()) {
    if (Date.now() - start > timeoutMs)
      assert.fail(`timed out waiting for ${what}`);
    await flush();
  }
}

function selectOf(f, name) {
  const select = editorOf(f, name, 'select');
  assert.ok(select, `${name} renders a choice editor`);
  return select;
}

/** The choice items on offer, the nullable empty option dropped. */
function items(f, name) {
  return selectOf(f, name).children.map((o) => o.value).filter((v) => v !== '');
}

function boxes(f, name) {
  const input = f.getInput(name);
  assert.ok(input, `${name} gets a field`);
  return input.root.querySelectorAll('.u2-multi-choice-item input').map((b) => b.value);
}

function stateOf(f, name) {
  return f.getInput(name).box.querySelector('.u2-param-source');
}

// --- the eval doubles themselves (green already: WO-3 does not touch them) ---

dbl('doubles: evalParamChoices resolves the command against the registry, never synchronously', async () => {
  provider('cities', [{name: 'country', propertyType: 'string'}],
    async (i) => i.country === 'DE' ? ['DE-1', 'DE-2'] : ['FR-1', 'FR-2']);
  const call = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities(country)'}}));
  call.setParamValue('country', 'DE');
  let landed = null;
  const pending = call.evalParamChoices('city').then((r) => landed = r);
  assert.equal(landed, null, 'delivery is never synchronous');
  await pending;
  assert.deepEqual(landed.items, ['DE-1', 'DE-2']);
  assert.deepEqual(landed.values, {'DE-1': 'DE-1', 'DE-2': 'DE-2'}, 'an array answer is identity values');
  assert.equal(landed.lookup, null);
  assert.deepEqual(landed.dependsOn, ['country'], 'the foreign bound param');
});

dbl('doubles: int items are stringified; a values/lookup answer passes through with no self dep', async () => {
  provider('ints', [], async () => [1, 2]);
  provider('cars', [{name: 'model', propertyType: 'string'}],
    async () => ({values: {Mazda: 'Mazda'}, lookup: {Mazda: {mpg: 21, cyl: 6}}}));
  const call = callOf(
    fp('n', 'string', {options: {choices: 'ints'}}),
    fp('model', 'string', {options: {choices: 'cars()'}}));
  const ints = await call.evalParamChoices('n');
  assert.deepEqual(ints.items, ['1', '2']);
  assert.deepEqual(ints.values, {'1': '1', '2': '2'});
  const cars = await call.evalParamChoices('model');
  assert.deepEqual(cars.items, ['Mazda']);
  assert.deepEqual(cars.lookup, {Mazda: {mpg: 21, cyl: 6}}, 'the propagate shape passes through');
  assert.deepEqual(cars.dependsOn, [], 'the param itself is never a dependency');
});

dbl('doubles: the single-input fallback binds the param itself — and the TYPED text for suggestions', async () => {
  const bound = [];
  provider('echo', [{name: 'q', propertyType: 'string'}], async (i) => {
    bound.push(i.q);
    return [String(i.q)];
  });
  const call = callOf(fp('city', 'string', {options: {choices: 'echo()', suggestions: 'echo()'}}));
  call.setParamValue('city', 'stored');
  const choices = await call.evalParamChoices('city');
  assert.deepEqual([bound, choices.dependsOn], [['stored'], []]);
  const suggestions = await call.evalParamSuggestions('city', 'typ');
  assert.deepEqual(bound, ['stored', 'typ'], 'suggestions bind the typed text, never the stored value');
  assert.deepEqual(suggestions, {items: ['typ'], tooltips: {}});
});

dbl('doubles: suggestions tooltips pass through; an unknown func rejects', async () => {
  provider('sug', [{name: 'q', propertyType: 'string'}],
    async (i) => ({items: [`${i.q}-1`], tooltips: {[`${i.q}-1`]: 'tip'}}));
  const call = callOf(
    fp('drug', 'string', {options: {suggestions: 'sug()'}}),
    fp('city', 'string', {options: {choices: 'nope()'}}));
  assert.deepEqual(await call.evalParamSuggestions('drug', 'a'), {items: ['a-1'], tooltips: {'a-1': 'tip'}});
  await assert.rejects(() => call.evalParamChoices('city'), /nope/);
  await assert.rejects(() => call.evalParamDefault('city'), /not found/i, 'no default command either');
});

dbl('doubles: evalParamDefault answers the provider value', async () => {
  provider('four', [], async () => 4);
  const call = callOf(fp('count', 'int', {options: {default: 'four()'}}));
  assert.equal(await call.evalParamDefault('count'), 4);
});

dbl('doubles: the evalDelayMs knob postpones delivery', async () => {
  provider('opts', [], async () => ['a']);
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()'}}));
  call.dart.evalDelayMs = 40;
  let landed = false;
  const pending = call.evalParamChoices('city').then(() => landed = true);
  await flush();
  assert.equal(landed, false, 'still in flight after a tick');
  await pending;
  assert.equal(landed, true);
});

dbl('doubles: a List value always fires onChanged, same reference included; scalars stay suppressed', async () => {
  const call = callOf(fp('tags', 'list'), fp('name', 'string'));
  const seen = [];
  call.inputParams['tags'].onChanged.subscribe(() => seen.push('tags'));
  call.inputParams['name'].onChanged.subscribe(() => seen.push('name'));
  const arr = ['a'];
  call.setParamValue('tags', arr);
  call.setParamValue('tags', arr);
  assert.deepEqual(seen, ['tags', 'tags'], 'the same-reference List write fires (func_call_param setter)');
  call.setParamValue('name', 'x');
  call.setParamValue('name', 'x');
  assert.deepEqual(seen, ['tags', 'tags', 'name'], 'the scalar suppression stays');
});

// --- routing (FP-W2-2) ---

w2('routing: the choices OPTION is the dynamic discriminator — annotation garbage in prop.choices never renders', async () => {
  provider('cities', [{name: 'city', propertyType: 'string'}], async () => ['Paris', 'Lyon']);
  const call = callOf(fp('city', 'string', {choices: ['ities('], options: {choices: 'cities()'}}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('city').root.dataset.u2, 'choice-input');
  await f.settled;
  assert.deepEqual(items(f, 'city'), ['Paris', 'Lyon'], 'the evaluator answers, the garbage never shows');
  f.dispose();
});

w2('routing: typed-setter choices without the option keep the W1 static path', async () => {
  const call = callOf(fp('stage', 'string', {choices: ['ideation', 'running']}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.deepEqual(items(f, 'stage'), ['ideation', 'running']);
  f.dispose();
});

w2('routing: a static list renders MultiChoiceInput, a choices option makes it dynamic, a bare list stays unsupported', async () => {
  provider('lists', [], async () => ['x', 'y']);
  const call = callOf(
    fp('tags', 'list', {choices: ['a', 'b']}),
    fp('cols', 'list', {options: {choices: 'lists()'}}),
    fp('other', 'list'));
  const f = mount(funcForm(call));
  assert.deepEqual([...f.unsupported], ['other']);
  assert.ok(f.getInput('tags'), 'a static list gets a field');
  assert.equal(f.getInput('tags').root.dataset.u2, 'multi-choice-input');
  assert.deepEqual(boxes(f, 'tags'), ['a', 'b']);
  assert.equal(f.getInput('cols').root.dataset.u2, 'multi-choice-input');
  await f.settled;
  assert.deepEqual(boxes(f, 'cols'), ['x', 'y']);
  f.dispose();
});

w2('routing: choices beats suggestions; suggestions alone builds a SuggestInput', async () => {
  provider('opts', [], async () => ['a']);
  provider('sug', [{name: 'q', propertyType: 'string'}], async () => ['s']);
  const call = callOf(
    fp('both', 'string', {options: {choices: 'opts()', suggestions: 'sug()'}}),
    fp('drug', 'string', {options: {suggestions: 'sug()'}}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('both').root.dataset.u2, 'choice-input');
  assert.equal(f.getInput('drug').root.dataset.u2, 'suggest-input');
  assert.equal(f.getInput('drug').root.classList.contains('u2-suggest-input'), true);
  await f.settled;
  f.dispose();
});

w2('a custom override editor owns the field: no dynamic wiring, no state element', async () => {
  let runs = 0;
  provider('opts', [{name: 'city', propertyType: 'string'}], async () => {
    runs++;
    return ['a'];
  });
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()'}}));
  const custom = new TextInput({name: 'city'});
  const f = mount(funcForm(call, {overrides: {city: {input: custom}}}));
  await f.settled;
  await sleep(50);
  assert.equal(runs, 0, 'the provider is never consulted');
  assert.equal(f.root.querySelector('.u2-param-source'), null);
  assert.equal(f.getInput('city'), custom);
  f.dispose();
  custom.dispose();
});

w2('a call without evalParamChoices renders a visible error state instead of throwing', async () => {
  const property = new FuncParam('city', 'string', {options: {choices: 'cities()'}});
  const param = {name: 'city', value: undefined, property,
    onChanged: {subscribe: () => ({unsubscribe() {}})}};
  const bare = {inputParams: {values: () => [param], city: param}, setParamValue() {}};
  const f = mount(funcForm(bare));
  await f.settled;
  const state = stateOf(f, 'city');
  assert.ok(state, 'a visible error state, never a silent skip');
  assert.ok(state.classList.contains('u2-param-source-error'));
  assert.ok(state.textContent.includes('Couldn\'t load choices:'), 'the operation prefixes the message');
  assert.equal(f.getInput('city').enabled, true, 'the field is not stuck loading');
  f.dispose();
});

// --- items application, value survival, prune-to-call (#9), dependsOn refresh ---

w2('a dep edit re-evals through the debounce; a vanished value is pruned INTO the call exactly once', async () => {
  let runs = 0;
  provider('cities', [{name: 'country', propertyType: 'string'}], async (i) => {
    runs++;
    return i.country === 'DE' ? ['DE-1', 'DE-2'] : ['FR-1', 'FR-2'];
  });
  const call = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities'}}));
  call.setParamValue('country', 'FR');
  call.setParamValue('city', 'FR-1');
  const {f, sets} = wired(call);
  await f.settled;
  assert.equal(runs, 1);
  assert.deepEqual(items(f, 'city'), ['FR-1', 'FR-2']);
  assert.equal(f.getInput('city').value.value, 'FR-1', 'a value still on offer survives setItems');
  assert.equal(call.inputs.city, 'FR-1');
  sets.length = 0;

  call.setParamValue('country', 'DE');
  await flush();
  assert.equal(runs, 1, 'the re-run waits out the 200 ms debounce');
  await until(() => items(f, 'city').includes('DE-1'), 'the dependent items to refresh');
  assert.equal(runs, 2, 'one dep edit is exactly one re-eval');
  await until(() => call.inputs.city === null, 'the prune to reach the call');
  assert.equal(f.getInput('city').value.value, null, 'the vanished value is pruned from the field');
  await sleep(250);
  assert.deepEqual(sets.filter(([n]) => n === 'city'), [['city', null]],
    'exactly one write — the refresh compares equal and the echo stops (divergence #9)');
  assert.equal(runs, 2, 'no echo storm');
  f.dispose();
});

w2('a multi-choice prune writes the filtered list into the call', async () => {
  provider('lists', [{name: 'country', propertyType: 'string'}],
    async (i) => i.country === 'DE' ? ['a'] : ['a', 'b', 'c']);
  const call = callOf(fp('country', 'string'), fp('tags', 'list', {options: {choices: 'lists()'}}));
  call.setParamValue('tags', ['a', 'c']);
  const {f, sets} = wired(call);
  await f.settled;
  assert.deepEqual(boxes(f, 'tags'), ['a', 'b', 'c']);
  assert.deepEqual(f.getInput('tags').value.value, ['a', 'c'], 'the checked set follows the call');
  sets.length = 0;

  call.setParamValue('country', 'DE');
  await until(() => boxes(f, 'tags').length === 1, 'the items to refresh');
  await until(() => Array.isArray(call.inputs.tags) && call.inputs.tags.length === 1, 'the prune to reach the call');
  assert.deepEqual(call.inputs.tags, ['a']);
  await sleep(250);
  assert.deepEqual(sets.filter(([n]) => n === 'tags'), [['tags', ['a']]],
    'one write despite the List always-fire — the refresh compares equal');
  f.dispose();
});

w2('a provider answering the same dep names causes zero re-subscriptions and zero extra evals', async () => {
  let runs = 0;
  provider('cities', [{name: 'country', propertyType: 'string'}], async (i) => {
    runs++;
    return [`${i.country ?? 'x'}-1`];
  });
  const call = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  const stream = call.inputParams['country'].onChanged;
  let resubs = 0;
  const subscribe = stream.subscribe.bind(stream);
  stream.subscribe = (next) => {
    resubs++;
    return subscribe(next);
  };
  const subCount = stream.count;

  call.setParamValue('country', 'DE');
  await until(() => runs === 2, 'the first re-eval');
  call.setParamValue('country', 'FR');
  await until(() => runs === 3, 'the second re-eval');
  await sleep(250);
  assert.equal(runs, 3, 'zero extra evals — reconciliation never bumps the version');
  assert.equal(resubs, 0, 'zero re-subscriptions for an unchanged dep set');
  assert.equal(stream.count, subCount);
  f.dispose();
});

w2('dep-set change reconciliation: the old dep stops triggering, the new one starts', async () => {
  let runsA = 0;
  let runsB = 0;
  provider('byA', [{name: 'a', propertyType: 'string'}], async () => {
    runsA++;
    return ['1'];
  });
  provider('byB', [{name: 'b', propertyType: 'string'}], async () => {
    runsB++;
    return ['2'];
  });
  const call = callOf(fp('a', 'string'), fp('b', 'string'),
    fp('city', 'string', {options: {choices: 'byA()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(runsA, 1);

  call.inputParams['city'].property.options['choices'] = 'byB()';
  call.setParamValue('a', 'x');
  await until(() => runsB === 1, 'the re-eval to read the new command');
  await sleep(250);
  assert.equal(runsA, 1, 'the old provider is done');

  call.setParamValue('a', 'y');
  await sleep(300);
  assert.equal(runsA + runsB, 2, 'the old dep no longer triggers');
  call.setParamValue('b', 'z');
  await until(() => runsB === 2, 'the new dep to trigger');
  f.dispose();
});

w2('no self-dependency: picking a value never re-runs the eval', async () => {
  let runs = 0;
  provider('opts', [{name: 'city', propertyType: 'string'}], async () => {
    runs++;
    return ['a', 'b'];
  });
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  const select = selectOf(f, 'city');
  select.value = 'a';
  fire(select, 'change');
  await until(() => call.inputs.city === 'a', 'the pick to write the call');
  await sleep(300);
  assert.equal(runs, 1, 'the write to the self-bound param never loops the source');
  f.dispose();
});

w2('a provider whose dependsOn names the param itself never self-subscribes', async () => {
  let runs = 0;
  provider('opts', [{name: 'city', propertyType: 'string'}], async () => {
    runs++;
    return ['a', 'b'];
  });
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()'}}));
  const evalChoices = call.evalParamChoices.bind(call);
  call.evalParamChoices = (name) => evalChoices(name).then((r) => ({...r, dependsOn: [name]}));
  const f = mount(funcForm(call));
  await f.settled;
  const select = selectOf(f, 'city');
  select.value = 'a';
  fire(select, 'change');
  await until(() => call.inputs.city === 'a', 'the pick to write the call');
  await sleep(300);
  assert.equal(runs, 1, 'the self-named dependency is filtered at reconcile');
  f.dispose();
});

w2('a dep param with no field still triggers the re-eval', async () => {
  let runs = 0;
  provider('byTable', [{name: 'table', propertyType: 'dataframe_list'}], async () => {
    runs++;
    return [`c${runs}`];
  });
  const call = callOf(fp('table', 'dataframe_list'), fp('city', 'string', {options: {choices: 'byTable()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.deepEqual([...f.unsupported], ['table'], 'the dep param renders no field');
  assert.equal(runs, 1);
  call.setParamValue('table', {name: 'demog'});
  await until(() => runs === 2, 'the fieldless dep to re-trigger');
  await until(() => items(f, 'city').includes('c2'), 'the refreshed items');
  f.dispose();
});

// --- abort / stale-drop ---

w2('stale-drop: a slow eval for dep A loses to a fast eval for dep B', async () => {
  const started = [];
  provider('cities', [{name: 'country', propertyType: 'string'}], async (i) => {
    started.push(i.country);
    await sleep(i.country === 'A' ? 150 : 1);
    return [`${i.country}-1`];
  });
  const call = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities()'}}));
  const f = mount(funcForm(call));
  await f.settled;

  call.setParamValue('country', 'A');
  await until(() => started.includes('A'), 'the slow eval to take off');
  call.setParamValue('country', 'B');
  await until(() => items(f, 'city').includes('B-1'), 'the fast eval to land');
  await sleep(200);
  assert.deepEqual(items(f, 'city'), ['B-1'], 'the aborted A landing never applies out of order');
  f.dispose();
});

w2('abort on rebind: a slow landing from the old call never applies to the new one', async () => {
  provider('cities', [{name: 'country', propertyType: 'string'}], async (i) => {
    await sleep(i.country === 'X' ? 120 : 1);
    return [`${i.country}-1`];
  });
  const make = (country) => {
    const c = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities()'}}));
    c.setParamValue('country', country);
    return c;
  };
  const callX = make('X');
  const callY = make('Y');
  const f = mount(funcForm(callX));
  f.source = callY;
  await f.settled;
  assert.deepEqual(items(f, 'city'), ['Y-1']);
  await sleep(180);
  assert.deepEqual(items(f, 'city'), ['Y-1'], 'the old generation\'s landing was dropped');
  assert.equal(callY.inputs.city, undefined, 'and nothing was written into the new call');
  f.dispose();
});

// --- propagate (FP-W2-4) ---

w2('propagate: a pick fills siblings through their input signals and the cascade terminates', async () => {
  let doseRuns = 0;
  provider('cars', [{name: 'model', propertyType: 'string'}], async () =>
    ({values: {Mazda: 'Mazda', Toyota: 'Toyota'},
      lookup: {Mazda: {mpg: 21, cyl: 6}, Toyota: {mpg: 28, cyl: 4}}}));
  provider('doses', [{name: 'mpg', propertyType: 'int'}], async (i) => {
    doseRuns++;
    return [`d${i.mpg}`];
  });
  const call = callOf(
    fp('model', 'string', {options: {choices: 'cars()', propagateChoice: 'all'}}),
    fp('mpg', 'int'),
    fp('cyl', 'int'),
    fp('dose', 'string', {options: {choices: 'doses()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.deepEqual(items(f, 'model'), ['Mazda', 'Toyota']);
  assert.equal(doseRuns, 1);

  const select = selectOf(f, 'model');
  select.value = 'Mazda';
  fire(select, 'change');
  await until(() => call.inputs.mpg === 21, 'the lookup row to reach the call through the sibling');
  assert.equal(call.inputs.model, 'Mazda');
  assert.equal(f.getInput('mpg').value.value, 21, 'the sibling input signal took the write');
  assert.equal(call.inputs.cyl, 6);
  await until(() => doseRuns === 2, 'the propagated mpg to re-eval its dependent source');
  await until(() => items(f, 'dose').includes('d21'), 'the dependent items with the propagated value');
  await sleep(250);
  assert.equal(doseRuns, 2, 'the cascade terminates');
  f.dispose();
});

w2('computed-default landings and external writes do not propagate', async () => {
  provider('cars', [{name: 'model', propertyType: 'string'}], async () =>
    ({values: {Mazda: 'Mazda', Toyota: 'Toyota'},
      lookup: {Mazda: {mpg: 21, cyl: 6}, Toyota: {mpg: 28, cyl: 4}}}));
  provider('pick', [], async () => 'Mazda');
  const call = callOf(
    fp('model', 'string', {options: {choices: 'cars()', propagateChoice: 'all', default: 'pick()'}}),
    fp('mpg', 'int'));
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(call.inputs.model, 'Mazda', 'the computed default landed');
  await sleep(250);
  assert.equal(call.inputs.mpg, undefined, 'a default landing never propagates');
  call.setParamValue('model', 'Toyota');
  await sleep(250);
  assert.equal(call.inputs.mpg, undefined, 'an external write echo never propagates');
  f.dispose();
});

w2('a system items-apply that changes the field value never propagates lookup values', async () => {
  provider('cars', [{name: 'model', propertyType: 'string'}, {name: 'country', propertyType: 'string'}],
    async (i) => i.country === 'B' ?
      {values: {Honda: 'Honda'}, lookup: {Honda: {mpg: 30}, 'null': {mpg: 99}}} :
      {values: {Mazda: 'Mazda'}, lookup: {Mazda: {mpg: 21}}});
  const call = callOf(
    fp('country', 'string'),
    fp('model', 'string', {options: {choices: 'cars()', propagateChoice: 'all'}}),
    fp('mpg', 'int'));
  const f = mount(funcForm(call));
  await f.settled;
  const select = selectOf(f, 'model');
  select.value = 'Mazda';
  fire(select, 'change');
  await until(() => call.inputs.mpg === 21, 'a user pick propagates');

  call.setParamValue('country', 'B');
  await until(() => call.inputs.model === null, 'the re-evaled items prune the vanished value into the call');
  await sleep(250);
  assert.equal(call.inputs.mpg, 21, 'the prune is a system apply — no lookup row is propagated');
  f.dispose();
});

w2('skipDefaultInit: items still apply, default writes and propagate wiring are suppressed', async () => {
  let defaultRuns = 0;
  provider('cars', [{name: 'model', propertyType: 'string'}], async () =>
    ({values: {Mazda: 'Mazda', Toyota: 'Toyota'},
      lookup: {Mazda: {mpg: 21, cyl: 6}, Toyota: {mpg: 28, cyl: 4}}}));
  provider('pick', [], async () => {
    defaultRuns++;
    return 'Mazda';
  });
  const call = callOf(
    fp('model', 'string', {options: {choices: 'cars()', propagateChoice: 'all', default: 'pick()'}}),
    fp('mpg', 'int'));
  const {f, sets} = wired(call, {skipDefaultInit: true});
  await f.settled;
  assert.deepEqual(items(f, 'model'), ['Mazda', 'Toyota'], 'items are not defaults (divergence #8)');
  assert.equal(defaultRuns, 0, 'the default eval is suppressed');
  assert.deepEqual(sets, [], 'nothing is written');
  const select = selectOf(f, 'model');
  select.value = 'Mazda';
  fire(select, 'change');
  await until(() => call.inputs.model === 'Mazda', 'the pick still writes the call');
  await sleep(250);
  assert.equal(call.inputs.mpg, undefined, 'propagate wiring is off under skipDefaultInit (fpe:585)');
  f.dispose();
});

// --- computed defaults (FP-W2-5, R6) ---

w2('a computed default is evaluated once and written through one setParamValue', async () => {
  let runs = 0;
  provider('four', [], async () => {
    runs++;
    return 4;
  });
  const call = callOf(fp('count', 'int', {options: {default: 'four()'}}), fp('name', 'string'));
  const {f, sets} = wired(call);
  await f.settled;
  assert.deepEqual(sets, [['count', 4]]);
  assert.equal(call.inputs.count, 4);
  assert.equal(f.getInput('count').value.value, 4, 'the write refreshed the field through the normal path');
  assert.equal(runs, 1);
  await sleep(250);
  assert.equal(runs, 1, 'defaults never re-run on their own');
  f.dispose();
});

w2('a param that already holds a value skips the default eval', async () => {
  let runs = 0;
  provider('four', [], async () => {
    runs++;
    return 4;
  });
  const call = callOf(fp('count', 'int', {options: {default: 'four()'}}));
  call.setParamValue('count', 9);
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(runs, 0);
  assert.equal(call.inputs.count, 9);
  f.dispose();
});

w2('a failing default shows the error state, writes nothing, and Retry re-runs the one eval', async () => {
  let fail = true;
  provider('four', [], async () => {
    if (fail)
      throw new Error('no default for you');
    return 4;
  });
  const call = callOf(fp('count', 'int', {options: {default: 'four()'}}), fp('name', 'string'));
  const {f, sets} = wired(call);
  await f.settled;
  const state = stateOf(f, 'count');
  assert.ok(state, 'the field carries a state element');
  assert.ok(state.classList.contains('u2-param-source-error'));
  assert.ok(state.textContent.includes('Couldn\'t compute the default: no default for you'),
    'the operation prefixes the raw message');
  assert.deepEqual(sets, [], 'a failed default writes nothing');
  assert.equal(f.getInput('count').enabled, true, 'the field stays usable');

  const name = editorOf(f, 'name');
  name.value = 'Aspirin';
  fire(name, 'input');
  await flush();
  assert.equal(call.inputs.name, 'Aspirin', 'the rest of the form is alive');

  fail = false;
  fire(state.querySelector('.u2-param-source-retry'), 'click');
  await until(() => call.inputs.count === 4, 'Retry to re-run the one default eval');
  f.dispose();
});

w2('a failed default stops claiming failure once the user fills the field by hand', async () => {
  provider('boom', [], async () => {
    throw new Error('down');
  });
  const call = callOf(fp('count', 'int', {options: {default: 'boom()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.ok(stateOf(f, 'count')?.classList.contains('u2-param-source-error'));

  const count = editorOf(f, 'count');
  count.value = '7';
  fire(count, 'input');
  await until(() => stateOf(f, 'count')?.classList.contains('u2-param-source-error') !== true,
    'the manual value to clear the error');
  assert.equal(call.inputs.count, 7);
  f.dispose();
});

w2('twoWayBinding: false — a user write still clears a failed default\'s error', async () => {
  provider('boom', [], async () => {
    throw new Error('down');
  });
  const call = callOf(fp('count', 'int', {options: {default: 'boom()'}}));
  const f = mount(funcForm(call, {twoWayBinding: false}));
  await f.settled;
  assert.ok(stateOf(f, 'count')?.classList.contains('u2-param-source-error'));

  f.getInput('count').value.value = 7;
  await until(() => stateOf(f, 'count')?.classList.contains('u2-param-source-error') !== true,
    'the user write to clear the error without a param subscription');
  assert.equal(call.inputs.count, 7);
  f.dispose();
});

w2('Retry never overwrites a value typed while the re-eval was in flight', async () => {
  let fail = true;
  provider('four', [], async () => {
    if (fail)
      throw new Error('down');
    return 4;
  });
  const call = callOf(fp('count', 'int', {options: {default: 'four()'}}));
  const {f, sets} = wired(call);
  await f.settled;
  const state = stateOf(f, 'count');
  assert.ok(state.classList.contains('u2-param-source-error'));

  fail = false;
  call.dart.evalDelayMs = 40;
  fire(state.querySelector('.u2-param-source-retry'), 'click');
  assert.equal(f.getInput('count').enabled, true, 'a retry is a refresh — the field stays editable');
  const count = editorOf(f, 'count');
  count.value = '7';
  fire(count, 'input');
  await until(() => stateOf(f, 'count')?.classList.contains('u2-param-source-error') !== true,
    'the error to clear');
  await sleep(80);
  assert.equal(call.inputs.count, 7, 'the landed retry never overwrote the typed value');
  assert.deepEqual(sets.filter(([, v]) => v === 4), [], 'no default write once the user has filled the field');
  f.dispose();
});

w2('rebind re-runs computed defaults only for the new call\'s still-null params', async () => {
  let runs = 0;
  provider('four', [], async () => {
    runs++;
    return 4;
  });
  const make = () => callOf(fp('count', 'int', {options: {default: 'four()'}}));
  const f = mount(funcForm(make()));
  await f.settled;
  assert.equal(runs, 1);

  const callB = make();
  f.source = callB;
  await f.settled;
  assert.equal(runs, 2, 'a still-null param on the new call re-runs');
  assert.equal(callB.inputs.count, 4);

  const callC = make();
  callC.setParamValue('count', 9);
  f.source = callC;
  await f.settled;
  assert.equal(runs, 2, 'a valued param does not');
  assert.equal(callC.inputs.count, 9);
  f.dispose();
});

w2('default-vs-items race: both landing orders settle on the ruled value', async () => {
  const build = async (itemsDelay, defaultDelay, def) => {
    Func.registry = [];
    provider('opts', [{name: 'city', propertyType: 'string'}], async () => {
      await sleep(itemsDelay);
      return ['a', 'b'];
    });
    provider('pick', [], async () => {
      await sleep(defaultDelay);
      return def;
    });
    const call = callOf(fp('city', 'string', {options: {choices: 'opts()', default: 'pick()'}}));
    const f = mount(funcForm(call));
    await f.settled;
    await flush();
    return {call, f};
  };

  const first = await build(1, 60, 'a');
  assert.equal(first.call.inputs.city, 'a', 'items first, default after — on offer, kept');
  assert.equal(first.f.getInput('city').value.value, 'a');
  first.f.dispose();

  const second = await build(60, 1, 'a');
  assert.equal(second.call.inputs.city, 'a', 'default first, items after — the value survives setItems');
  second.f.dispose();

  const third = await build(60, 1, 'zzz');
  assert.equal(third.call.inputs.city, null,
    'a default the items never offer is pruned to null (divergence #9)');
  third.f.dispose();
});

// --- settled (FP-W2-8) ---

w2('settled resolves once the initial loads land, immediately for a static form, and never rejects', async () => {
  const staticForm = mount(funcForm(callOf(fp('name', 'string'), fp('stage', 'string', {choices: ['a', 'b']}))));
  assert.ok(staticForm.settled instanceof Promise, 'settled is a class member (R5)');
  await staticForm.settled;
  staticForm.dispose();

  provider('slow', [{name: 'city', propertyType: 'string'}], async () => {
    await sleep(60);
    return ['a'];
  });
  provider('boom', [], async () => {
    throw new Error('down');
  });
  const call = callOf(
    fp('city', 'string', {options: {choices: 'slow()'}}),
    fp('count', 'int', {options: {default: 'boom()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.deepEqual(items(f, 'city'), ['a'], 'settled waited for the slow eval');
  assert.ok(stateOf(f, 'count')?.classList.contains('u2-param-source-error'),
    'and for the failed default — without rejecting');
  f.dispose();
});

w2('settled re-arms per rebind generation', async () => {
  provider('cities', [{name: 'country', propertyType: 'string'}], async (i) => {
    await sleep(20);
    return [`${i.country}-1`];
  });
  const make = (country) => {
    const c = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities()'}}));
    c.setParamValue('country', country);
    return c;
  };
  const f = mount(funcForm(make('A')));
  const first = f.settled;
  await first;
  assert.deepEqual(items(f, 'city'), ['A-1']);
  f.source = make('B');
  assert.notEqual(f.settled, first, 'a rebind re-arms settled');
  await f.settled;
  assert.deepEqual(items(f, 'city'), ['B-1']);
  f.dispose();
});

w2('settled settles when the form is disposed before the first landing', async () => {
  provider('slow', [{name: 'city', propertyType: 'string'}], async () => {
    await sleep(120);
    return ['a'];
  });
  const f = mount(funcForm(callOf(fp('city', 'string', {options: {choices: 'slow()'}}))));
  const settled = f.settled;
  f.dispose();
  assert.ok(await Promise.race([settled.then(() => true), sleep(400).then(() => false)]),
    'dispose resolves the pending settled instead of hanging it');
});

w2('a source rebind settles the old generation\'s settled', async () => {
  provider('slow', [{name: 'city', propertyType: 'string'}], async () => {
    await sleep(120);
    return ['a'];
  });
  provider('fast', [{name: 'city', propertyType: 'string'}], async () => ['b']);
  const f = mount(funcForm(callOf(fp('city', 'string', {options: {choices: 'slow()'}}))));
  const old = f.settled;
  f.source = callOf(fp('city', 'string', {options: {choices: 'fast()'}}));
  assert.ok(await Promise.race([old.then(() => true), sleep(400).then(() => false)]),
    'the replaced generation\'s settled still settles');
  await f.settled;
  assert.deepEqual(items(f, 'city'), ['b']);
  f.dispose();
});

// --- suggestions in the form (FP-W2-6) ---

w2('suggestions: the fetch receives the typed text, a pick writes through, rows carry tooltips', async () => {
  const texts = [];
  provider('sug', [{name: 'q', propertyType: 'string'}], async (i) => {
    texts.push(i.q);
    return {items: [`${i.q}-1`, `${i.q}-2`], tooltips: {[`${i.q}-1`]: 'first match'}};
  });
  const call = callOf(fp('drug', 'string', {options: {suggestions: 'sug()'}}));
  call.setParamValue('drug', 'stored');
  const f = mount(funcForm(call));
  assert.equal(f.getInput('drug').root.dataset.u2, 'suggest-input');
  await f.settled;
  assert.deepEqual(texts, [], 'suggestions are keystroke-driven — no initial query, no part in settled');

  const input = f.getInput('drug').root.querySelector('input');
  input.value = 'asp';
  fire(input, 'input');
  await until(() => texts.length > 0, 'the debounced fetch to fire');
  assert.deepEqual(texts, ['asp'], 'the TYPED text, never the stored value');
  assert.equal(call.inputs.drug, 'asp', 'free text writes the call as TextInput does');
  await until(() => document.body.querySelectorAll('.u2-typeahead-option').length === 2, 'the popup rows');
  const rows = document.body.querySelectorAll('.u2-typeahead-option');
  assert.ok(rows[0].getAttribute('title') === 'first match' ||
    rows[0].querySelector('[title="first match"]') !== null, 'the tooltip rides the row');
  fire(rows[0], 'pointerdown');
  await flush();
  assert.equal(f.getInput('drug').value.value, 'asp-1');
  assert.equal(call.inputs.drug, 'asp-1', 'the pick reaches the call');

  input.value = '';
  fire(input, 'input');
  await sleep(200);
  assert.deepEqual(texts, ['asp'], 'an empty field never queries the provider (minChars 1)');

  fire(input, 'keydown', {key: 'Escape'});
  fire(input, 'keydown', {key: 'Escape'});
  await flush();
  assert.equal(f.getInput('drug').value.value, 'asp-1',
    'Escape with the popup closed reverts the type-over to the committed pick');
  assert.equal(call.inputs.drug, 'asp-1', 'and the revert reaches the call');
  f.dispose();
});

// --- the state element and enabled reconciliation (FP-W2-7) ---

w2('the initial load disables the input and shows the grace-delayed spinner; ready re-enables, the orphan cycle ends enabled', async () => {
  const grace = ParamState.graceMs;
  ParamState.graceMs = 50;
  try {
    provider('slow', [{name: 'city', propertyType: 'string'}], async () => {
      await sleep(300);
      return ['a', 'b'];
    });
    const withCity = () => callOf(fp('city', 'string', {options: {choices: 'slow()'}}));
    const f = mount(funcForm(withCity()));
    const city = f.getInput('city');
    await flush();
    assert.equal(city.enabled, false, 'disabled while the initial load is in flight');
    const state = stateOf(f, 'city');
    assert.ok(state, 'the state element rides input.box');
    assert.equal(state.querySelector('.u2-loader-spinner'), null, 'no spinner before the grace elapses');
    await until(() => state.classList.contains('u2-param-source-loading'), 'the spinner after the grace');
    assert.ok(state.querySelector('.u2-loader-spinner'));
    await f.settled;
    assert.equal(city.enabled, true, 'ready re-enables');
    assert.equal(city.box.querySelector('.u2-param-source-loading'), null);
    assert.equal(city.box.querySelector('.u2-param-source-error'), null);

    f.source = callOf(fp('other', 'string'));
    assert.equal(city.enabled, false, 'orphaned');
    f.source = withCity();
    await f.settled;
    assert.equal(city.enabled, true, 'orphan → rebind → load ends enabled');
    f.dispose();
  } finally {
    ParamState.graceMs = grace;
  }
});

w2('the default error\'s Retry stays clickable while the other slot runs an initial load', async () => {
  let failChoices = true;
  provider('opts', [{name: 'city', propertyType: 'string'}], async () => {
    if (failChoices)
      throw new Error('c-down');
    await sleep(150);
    return ['a'];
  });
  provider('boom', [], async () => {
    throw new Error('d-down');
  });
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()', default: 'boom()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  const state = stateOf(f, 'city');
  assert.ok(state.classList.contains('u2-param-source-error'), 'both slots errored');

  failChoices = false;
  fire(state.querySelector('.u2-param-source-retry'), 'click');
  await until(() => f.getInput('city').enabled === false, 'the retried initial load to disable the input');
  assert.equal(state.querySelector('.u2-param-source-retry').disabled, false,
    'the disable sweep spares the default error\'s Retry');
  await until(() => f.getInput('city').enabled === true, 'the landing to re-enable');
  f.dispose();
});

w2('a fast initial eval never shows the spinner', async () => {
  provider('opts', [], async () => ['a', 'b']);
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()'}}));
  const f = mount(funcForm(call));
  assert.equal(stateOf(f, 'city')?.querySelector('.u2-loader-spinner') ?? null, null,
    'no spinner at the loading start');
  await f.settled;
  assert.equal(stateOf(f, 'city')?.classList.contains('u2-param-source-loading') ?? false, false);
  assert.deepEqual(items(f, 'city'), ['a', 'b']);
  f.dispose();
});

w2('a dependent refresh keeps the landed input enabled and swaps the items atomically', async () => {
  provider('cities', [{name: 'country', propertyType: 'string'}],
    async (i) => i.country === 'DE' ? ['DE-1', 'DE-2'] : ['FR-1', 'FR-2']);
  const call = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'cities()'}}));
  call.setParamValue('country', 'FR');
  const f = mount(funcForm(call));
  await f.settled;
  const city = f.getInput('city');
  assert.equal(city.enabled, true);

  const seen = [];
  call.setParamValue('country', 'DE');
  await until(() => {
    seen.push(city.enabled);
    return items(f, 'city').includes('DE-1');
  }, 'the refreshed items');
  assert.ok(seen.every((enabled) => enabled === true),
    'the input never went disabled through the debounce and eval window');
  assert.deepEqual(items(f, 'city'), ['DE-1', 'DE-2'], 'the items swapped in one go');
  f.dispose();
});

w2('an eval error keeps the stale items usable and the input enabled', async () => {
  let fail = false;
  provider('opts', [{name: 'country', propertyType: 'string'}], async () => {
    if (fail)
      throw new Error('backend down');
    return ['a', 'b'];
  });
  const call = callOf(fp('country', 'string'), fp('city', 'string', {options: {choices: 'opts()'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.deepEqual(items(f, 'city'), ['a', 'b']);

  fail = true;
  call.setParamValue('country', 'DE');
  await until(() => stateOf(f, 'city')?.classList.contains('u2-param-source-error') === true, 'the error state');
  assert.ok(stateOf(f, 'city').textContent.includes('Couldn\'t load choices: backend down'),
    'the operation prefixes the raw message');
  assert.deepEqual(items(f, 'city'), ['a', 'b'], 'the previous items stay, as the Dart select does');
  assert.equal(f.getInput('city').enabled, true, 'the input stays enabled on error');
  assert.ok(stateOf(f, 'city').querySelector('.u2-param-source-retry'), 'with a Retry');

  fail = false;
  fire(stateOf(f, 'city').querySelector('.u2-param-source-retry'), 'click');
  await until(() => stateOf(f, 'city')?.classList.contains('u2-param-source-error') !== true, 'the retry to clear it');
  f.dispose();
});

// --- descriptions, overrides bag, empty number ---

w2('descriptions: the selected item\'s description becomes the input tooltip', async () => {
  provider('fruits', [], async () => ['apple', 'pear']);
  const call = callOf(fp('fruit', 'string',
    {description: 'Pick one', options: {choices: 'fruits()', descriptions: {apple: 'An apple'}}}));
  const f = mount(funcForm(call));
  await f.settled;
  const select = selectOf(f, 'fruit');
  select.value = 'apple';
  fire(select, 'change');
  await flush();
  assert.equal(f.getInput('fruit').box.title, 'An apple', 'the tooltip rides the editor area, not the label');
  select.value = 'pear';
  fire(select, 'change');
  await flush();
  assert.equal(f.getInput('fruit').box.title, 'Pick one',
    'an undescribed item falls back to the property description');
  f.dispose();
});

w2('an option-bag override merges into a dynamic field without suppressing its wiring', async () => {
  provider('opts', [{name: 'city', propertyType: 'string'}], async () => ['a', 'b']);
  const call = callOf(fp('city', 'string', {options: {choices: 'opts()'}}));
  const f = mount(funcForm(call, {overrides: {city: {label: 'Town'}}}));
  await f.settled;
  assert.equal(labelOf(f, 'city'), 'Town');
  assert.deepEqual(items(f, 'city'), ['a', 'b'], 'the dynamic wiring stays');
  f.dispose();
});

w2('clearing a bounded number writes null into the call, never NaN or the empty string', async () => {
  const call = callOf(fp('count', 'int', {min: 0, max: 10}));
  const {f, sets} = wired(call);
  const count = editorOf(f, 'count');
  count.value = '5';
  fire(count, 'input');
  await flush();
  assert.equal(call.inputs.count, 5);
  sets.length = 0;

  count.value = '';
  fire(count, 'input');
  fire(count, 'blur');
  await flush();
  assert.deepEqual(sets, [['count', null]]);
  assert.equal(call.dart.params[0].dart.value, null, 'null, not NaN and not the empty string');
  assert.equal(f.isValid, true, 'an empty nullable number is not invalid');
  f.dispose();
});
