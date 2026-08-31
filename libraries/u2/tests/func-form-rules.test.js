/* The W4 FuncCallForm contract of plan-w4.md (FP-W4-2…FP-W4-9) over the FuncCall doubles and the
   scriptSync stub: visible:/enabled: expressions (keep-previous, composition, exemption #20),
   validator: regex/expression verdicts, named validators through evalParamValidators (stale-drop,
   inline failures #21, warnings on the notice channel), missingRequired, categoryGroups
   (#18/#19/#22) with header auto-hide — plus the verification-strategy traps: the echo audit,
   revalidate storms, the hide→show verdict round-trip and the W2 debounce interplay.
   The stub is installed per test and deleted in teardown; doubles stay getter-backed. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Func, FuncParam, dart, scriptSyncStub} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);

const {funcForm} = await import('../src/dg/funcs/func-form.js');

/** Every test runs with the stub installed and torn down, a clean document, an empty registry
 * and the live-scope count back where it was. */
function rules(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    globalThis.grok_ScriptSync = scriptSyncStub;
    try {
      await body();
    } finally {
      delete globalThis.grok_ScriptSync;
      dart.scriptSyncError = null;
      Func.registry = [];
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** No stub: the warn-once inert paths of a call/environment without the new members. */
function bare(name, body) {
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

function fp(name, propertyType, data = {}) {
  return new FuncParam(name, propertyType, data);
}

function callOf(...params) {
  return new Func('fceRules', {inputs: params}).prepare({});
}

function callWith(options, ...params) {
  return new Func('fceRules', {inputs: params, options}).prepare({});
}

function mount(f) {
  document.body.append(f.root);
  return f;
}

function editorOf(f, name, selector = 'input') {
  return f.getInput(name).root.querySelector(selector);
}

function typeInto(f, name, text) {
  const el = editorOf(f, name);
  el.value = text;
  fire(el, 'input');
}

function validityOf(f, name) {
  return f.getInput(name).validity.value;
}

function hiddenOf(f, name) {
  return f.getInput(name).root.hidden;
}

function headerEls(f) {
  return [...f.root.querySelectorAll('.u2-form-category')];
}

function headerOf(f, text) {
  const el = headerEls(f).find((h) => h.textContent === text);
  assert.ok(el, `a '${text}' header exists`);
  return el;
}

function stateOf(f, name) {
  return f.getInput(name).box.querySelector('.u2-param-source');
}

/** A registry-backed named validator: SYNC run over its single input. */
function validator(name, verdict, data = {}) {
  const f = new Func(name, {inputs: [{name: 'v', propertyType: 'string'}],
    run: (i) => verdict(i.v), ...data});
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

function countWarns(body) {
  const warn = console.warn;
  let count = 0;
  console.warn = () => count++;
  try {
    body();
  } finally {
    console.warn = warn;
  }
  return count;
}

/** The vehicle fixture of FP-W4-11, trimmed to the routes under test. `type` is preset in the
 * CALL, the way the live bench prepares it — expressions evaluate over what the call holds,
 * never over an uncommitted `defaultValue`. */
function vehicle() {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('cylinders', 'int', {options: {visible: 'type == "ICE"'}}),
    fp('tankVolume', 'double'),
    fp('tankExtension', 'bool', {options: {visible: 'type == "ICE"', enabled: 'tankVolume > 50'}}),
    fp('batteryCapacity', 'double', {options: {visible: 'type == "Electric"'}}));
  call.setParamValue('type', 'ICE');
  return call;
}

// --- visible: / enabled: expressions ---

rules('initial visibility comes from the call\'s values', () => {
  const f = mount(funcForm(vehicle()));
  assert.equal(hiddenOf(f, 'cylinders'), false, 'type preset to ICE');
  assert.equal(hiddenOf(f, 'batteryCapacity'), true);
  assert.equal(hiddenOf(f, 'type'), false, 'no expression, no hiding');
  f.dispose();
});

rules('a sibling edit flips visibility both ways, [hidden] on the roots', () => {
  const f = mount(funcForm(vehicle()));
  f.getInput('type').value.value = 'Electric';
  assert.equal(hiddenOf(f, 'cylinders'), true);
  assert.equal(hiddenOf(f, 'batteryCapacity'), false);
  f.getInput('type').value.value = 'ICE';
  assert.equal(hiddenOf(f, 'cylinders'), false);
  assert.equal(hiddenOf(f, 'batteryCapacity'), true);
  f.dispose();
});

rules('enabled: follows its expression both ways without touching visibility', () => {
  const f = mount(funcForm(vehicle()));
  assert.equal(f.getInput('tankExtension').enabled, false, '0 > 50 is false at open');
  assert.equal(hiddenOf(f, 'tankExtension'), false, 'visible under the ICE default');
  typeInto(f, 'tankVolume', '60');
  assert.equal(f.getInput('tankExtension').enabled, true);
  typeInto(f, 'tankVolume', '40');
  assert.equal(f.getInput('tankExtension').enabled, false);
  f.dispose();
});

rules('disabled-by-expression says why on the root title — checkbox and label hovers included — and clears it on re-enable (#23)', () => {
  const f = mount(funcForm(vehicle()));
  const input = f.getInput('tankExtension');
  assert.equal(input.root.title, 'Enabled when: tankVolume > 50', 'disabled at open (0 > 50)');
  assert.ok(input.root.querySelector('.u2-input-label'), 'the label hovers inside the titled root');
  assert.ok(input.root.querySelector('input[type="checkbox"]'), 'so does the disabled checkbox');
  assert.equal(input.box.title, '', 'never on the box — the label sits outside it');
  typeInto(f, 'tankVolume', '60');
  assert.equal(input.root.title, '', 'released on re-enable');
  typeInto(f, 'tankVolume', '40');
  assert.equal(input.root.title, 'Enabled when: tankVolume > 50', 'and back on re-disable');
  f.dispose();
});

rules('the expression title displaces a pre-existing root title and never clobbers a foreign one', () => {
  const call = callOf(
    fp('tankVolume', 'double'),
    fp('tankExtension', 'bool', {options: {enabled: 'tankVolume > 50'}}));
  call.setParamValue('tankVolume', 60);
  const f = mount(funcForm(call));
  const root = f.getInput('tankExtension').root;
  root.title = 'mine';
  typeInto(f, 'tankVolume', '40');
  assert.equal(root.title, 'Enabled when: tankVolume > 50');
  typeInto(f, 'tankVolume', '70');
  assert.equal(root.title, 'mine', 'the displaced title is restored');

  typeInto(f, 'tankVolume', '40');
  root.title = 'other';
  typeInto(f, 'tankVolume', '70');
  assert.equal(root.title, 'other', 'a title someone else wrote while disabled stays theirs');
  f.dispose();
});

rules('a script failure keeps the previous state and warns once (ib:328-331 parity)', () => {
  const f = mount(funcForm(vehicle()));
  assert.equal(hiddenOf(f, 'cylinders'), false);
  dart.scriptSyncError = true;
  const warns = countWarns(() => {
    f.getInput('type').value.value = 'Electric';
    typeInto(f, 'tankVolume', '10');
  });
  assert.equal(hiddenOf(f, 'cylinders'), false, 'previous visibility kept over two failures');
  assert.equal(warns, 4, 'once per field+key (4 expressions), never per keystroke');
  dart.scriptSyncError = null;
  typeInto(f, 'tankVolume', '20');
  assert.equal(hiddenOf(f, 'cylinders'), true, 'recovery follows the expression again');
  f.dispose();
});

rules('enabled composes with orphan: the expression term survives an orphan→rebind cycle', () => {
  const make = (type) => {
    const call = callOf(
      fp('type', 'string', {choices: ['ICE', 'Electric']}),
      fp('gears', 'int', {options: {enabled: 'type == "ICE"'}}));
    call.setParamValue('type', type);
    return call;
  };
  const f = mount(funcForm(make('ICE')));
  assert.equal(f.getInput('gears').enabled, true);

  f.source = callOf(fp('type', 'string', {choices: ['ICE', 'Electric']}));
  assert.equal(f.getInput('gears').enabled, false, 'orphaned disables regardless of the expression');

  f.source = make('Electric');
  assert.equal(f.getInput('gears').enabled, false, 'rebound, but the expression says disabled');
  f.getInput('type').value.value = 'ICE';
  assert.equal(f.getInput('gears').enabled, true, 'the cycle ends per expression');
  f.dispose();
});

rules('enabled composes with busy: a loading choices source disables regardless of the expression', async () => {
  Func.registry.push(new Func('cities', {inputs: [], run: async () => ['Boston', 'Dallas']}));
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('city', 'string', {options: {choices: 'cities()', enabled: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  call.dart.evalDelayMs = 40;
  const f = mount(funcForm(call));
  assert.equal(f.getInput('city').enabled, false, 'busy while the initial load is in flight');
  await f.settled;
  assert.equal(f.getInput('city').enabled, true, 'landed: the expression owns it again');
  f.getInput('type').value.value = 'Electric';
  assert.equal(f.getInput('city').enabled, false, 'and says disabled');
  f.dispose();
});

rules('a rebind re-applies the expressions over the new call\'s values immediately', () => {
  const make = (type) => {
    const call = callOf(
      fp('type', 'string', {choices: ['ICE', 'Electric']}),
      fp('cylinders', 'int', {options: {visible: 'type == "ICE"'}}));
    call.setParamValue('type', type);
    return call;
  };
  const f = mount(funcForm(make('ICE')));
  assert.equal(hiddenOf(f, 'cylinders'), false);
  f.source = make('Electric');
  assert.equal(hiddenOf(f, 'cylinders'), true, 'hidden at rebind, before any edit');
  f.source = make('ICE');
  assert.equal(hiddenOf(f, 'cylinders'), false);
  f.dispose();
});

rules('warn-once re-arms per bind: a rebind under a still-failing script warns again', () => {
  dart.scriptSyncError = true;
  const make = () => callOf(fp('cylinders', 'int', {options: {visible: 'type == "ICE"'}}),
    fp('type', 'string'));
  let f;
  const atBuild = countWarns(() => f = mount(funcForm(make())));
  assert.equal(atBuild, 1);
  const atRebind = countWarns(() => f.source = make());
  assert.equal(atRebind, 1, 'the form\'s warn map cleared on _armRules');
  f.dispose();
});

rules('two live forms over ONE call warn once EACH; one form\'s re-arm never resets the other', () => {
  dart.scriptSyncError = true;
  const call = callOf(fp('cylinders', 'int', {options: {visible: 'type == "ICE"'}}),
    fp('type', 'string'));
  let a, b;
  const atBuild = countWarns(() => {
    a = mount(funcForm(call));
    b = mount(funcForm(call));
  });
  assert.equal(atBuild, 2, 'per form, not one shared warn — the bench A/B pair scenario');
  const atRebind = countWarns(() => b.source = call);
  assert.equal(atRebind, 1, 'only the re-armed form warns again');
  const after = countWarns(() => a.getInput('type').value.value = 'x');
  assert.equal(after, 0, 'a stays warn-once for its key: b\'s clear did not reset it (b already warned)');
  a.dispose();
  b.dispose();
});

// --- exemption (#20): inactive fields and the value survival ---

rules('disabled-by-expression is inactive too: required suppressed while the field stays visible', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('gears', 'int', {nullable: false, options: {enabled: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  assert.equal(validityOf(f, 'gears'), 'Value can\'t be empty');
  assert.deepEqual(f.missingRequired.value, ['Gears']);

  f.getInput('type').value.value = 'Electric';
  assert.equal(hiddenOf(f, 'gears'), false, 'disabled, not hidden');
  assert.equal(f.getInput('gears').enabled, false);
  assert.equal(validityOf(f, 'gears'), null, 'a field the user cannot edit never blocks (#20)');
  assert.deepEqual(f.missingRequired.value, []);
  f.dispose();
});

rules('required is suppressed while hidden, and the value SURVIVES in the call', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('cylinders', 'int', {nullable: false, options: {visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  assert.equal(validityOf(f, 'cylinders'), 'Value can\'t be empty');
  typeInto(f, 'cylinders', '6');
  assert.equal(validityOf(f, 'cylinders'), null);
  assert.equal(call.inputs.cylinders, 6);

  f.getInput('type').value.value = 'Electric';
  assert.equal(hiddenOf(f, 'cylinders'), true);
  assert.equal(validityOf(f, 'cylinders'), null, 'hidden stays exempt');
  assert.equal(call.inputs.cylinders, 6, 'hiding never clears the call (Dart parity)');

  typeInto(f, 'cylinders', '');
  f.getInput('type').value.value = 'ICE';
  assert.equal(validityOf(f, 'cylinders'), 'Value can\'t be empty',
    'the verdict returns the moment the field shows again');
  f.dispose();
});

rules('an expression-invalid verdict clears when hidden and returns when shown, without an edit', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('minAge', 'int', {options: {visible: 'type == "ICE"', validator: 'minAge > 18'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  typeInto(f, 'minAge', '15');
  assert.equal(validityOf(f, 'minAge'), 'minAge > 18', 'false answers the expression text (parity)');

  f.getInput('type').value.value = 'Electric';
  assert.equal(validityOf(f, 'minAge'), null, 'validity null + gate open while hidden');
  f.getInput('type').value.value = 'ICE';
  assert.equal(validityOf(f, 'minAge'), 'minAge > 18', 'the SAME verdict returns on the flip');
  f.dispose();
});

// --- validator: regex and expressions ---

rules('a regex-literal validator maps to the Dart non-match message, empty-gated', () => {
  const f = mount(funcForm(callOf(fp('email', 'string', {options: {validator: '/^[^@]+@[^@]+$/'}}))));
  assert.equal(validityOf(f, 'email'), null, 'empty is the required validator\'s business alone');
  typeInto(f, 'email', 'nope');
  assert.equal(validityOf(f, 'email'), 'Value doesn\'t match /^[^@]+@[^@]+$/');
  assert.equal(f.getInput('email').root.querySelector('.u2-input-error').title,
    'Value doesn\'t match /^[^@]+@[^@]+$/', 'the clamped message\'s full text rides the title');
  typeInto(f, 'email', 'a@b.com');
  assert.equal(validityOf(f, 'email'), null);
  f.dispose();
});

rules('an expression validator: false → the expression text, string → the message, throw → the Dart wording', () => {
  const f = mount(funcForm(callOf(
    fp('minAge', 'int', {options: {validator: 'minAge > 18'}}),
    fp('note', 'string', {options: {validator: 'value == "ok" or "Only ok will do"'}}))));
  typeInto(f, 'minAge', '15');
  assert.equal(validityOf(f, 'minAge'), 'minAge > 18');
  typeInto(f, 'minAge', '21');
  assert.equal(validityOf(f, 'minAge'), null);

  typeInto(f, 'note', 'meh');
  assert.equal(validityOf(f, 'note'), 'Only ok will do');
  typeInto(f, 'note', 'ok');
  assert.equal(validityOf(f, 'note'), null);

  dart.scriptSyncError = true;
  typeInto(f, 'minAge', '15');
  assert.equal(validityOf(f, 'minAge'), 'Error during validation: "minAge > 18"');
  f.dispose();
});

rules('the verdict moves on a sibling change under an unchanged own value', () => {
  const f = mount(funcForm(callOf(
    fp('min', 'int'),
    fp('count', 'int', {options: {validator: 'count > min'}}))));
  typeInto(f, 'min', '3');
  typeInto(f, 'count', '5');
  assert.equal(validityOf(f, 'count'), null);
  typeInto(f, 'min', '10');
  assert.equal(validityOf(f, 'count'), 'count > min', 'the sibling edit revalidated count');
  typeInto(f, 'min', '2');
  assert.equal(validityOf(f, 'count'), null);
  f.dispose();
});

// --- missingRequired (FP-W4-7) ---

rules('missingRequired lists non-optional empty fields by label and reacts to edits', () => {
  const f = mount(funcForm(callOf(
    fp('firstName', 'string'),
    fp('age', 'int'),
    fp('nickname', 'string', {isOptional: true}))));
  assert.deepEqual(f.missingRequired.value, ['First Name', 'Age'],
    'labels, not param names; optional excluded');
  typeInto(f, 'firstName', 'Ada');
  assert.deepEqual(f.missingRequired.value, ['Age']);
  typeInto(f, 'age', '36');
  assert.deepEqual(f.missingRequired.value, []);
  f.dispose();
});

rules('missingRequired is visibility-aware: a hidden required field blocks nothing', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('cylinders', 'int', {options: {visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  assert.deepEqual(f.missingRequired.value, ['Cylinders'],
    'type holds a value in the call and does not read as missing');
  f.getInput('type').value.value = 'Electric';
  assert.deepEqual(f.missingRequired.value, [], 'hidden dropped from the gate');
  f.getInput('type').value.value = 'ICE';
  assert.deepEqual(f.missingRequired.value, ['Cylinders'], 'and back the moment it shows');
  f.dispose();
});

rules('missingRequired recomputes over a rebind (orphans excluded, values re-read)', () => {
  const f = mount(funcForm(callOf(fp('name', 'string'), fp('age', 'int'))));
  assert.deepEqual(f.missingRequired.value, ['Name', 'Age']);

  const next = callOf(fp('name', 'string'));
  next.setParamValue('name', 'Ada');
  f.source = next;
  assert.deepEqual(f.missingRequired.value, [], 'age orphaned out, name filled');

  const back = callOf(fp('name', 'string'), fp('age', 'int'));
  f.source = back;
  assert.deepEqual(f.missingRequired.value, ['Name', 'Age']);
  f.dispose();
});

rules('invalidFields lists invalid rendered fields by label; a hidden invalid field is exempt (#20)', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('minAge', 'int', {options: {visible: 'type == "ICE"', validator: 'minAge > 18'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  assert.deepEqual(f.invalidFields.value, []);
  typeInto(f, 'minAge', '15');
  assert.deepEqual(f.invalidFields.value, ['Min Age'], 'labels, not param names');

  f.getInput('type').value.value = 'Electric';
  assert.deepEqual(f.invalidFields.value, [], 'hidden dropped from the list');
  f.getInput('type').value.value = 'ICE';
  assert.deepEqual(f.invalidFields.value, ['Min Age'], 'and back the moment it shows');
  typeInto(f, 'minAge', '25');
  assert.deepEqual(f.invalidFields.value, []);
  f.dispose();
});

rules('a parse-invalid field hidden by its expression unblocks the gate lists, not isValid (#20)', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('minAge', 'int', {isOptional: true, options: {visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  typeInto(f, 'minAge', 'abc');
  assert.equal(validityOf(f, 'minAge'), 'Not a number');
  assert.deepEqual(f.invalidFields.value, ['Min Age'], 'the gate names the unparseable field');
  assert.deepEqual(f.missingRequired.value, []);

  f.getInput('type').value.value = 'Electric';
  assert.deepEqual(f.invalidFields.value, [], 'hidden: both gate lists empty — Run unblocks');
  assert.deepEqual(f.missingRequired.value, []);
  assert.equal(validityOf(f, 'minAge'), 'Not a number',
    'the constructor parse verdict is not active-gated');
  assert.equal(f.isValid, false, 'isValid stays visibility-blind by contract');

  f.getInput('type').value.value = 'ICE';
  assert.deepEqual(f.invalidFields.value, ['Min Age'], 'and blocks again the moment it shows');
  f.dispose();
});

// --- named validators (FP-W4-5, #21) ---

rules('a named verdict lands asynchronously and clears on a passing value', async () => {
  validator('isCode', (v) => String(v).startsWith('X') ? null : 'Code must start with X');
  const call = callOf(fp('code', 'string', {options: {validators: ['isCode']}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(validityOf(f, 'code'), null, 'empty is not validated');

  typeInto(f, 'code', 'abc');
  assert.equal(validityOf(f, 'code'), null, 'the verdict is never synchronous');
  await until(() => validityOf(f, 'code') === 'Code must start with X', 'the verdict to land');

  typeInto(f, 'code', 'X-1');
  await until(() => validityOf(f, 'code') === null, 'the pass to clear the verdict');
  f.dispose();
});

rules('false maps to the didn\'t-pass message; a full result object passes through', async () => {
  validator('flag', (v) => v !== 'bad');
  const call = callOf(fp('code', 'string', {options: {validators: ['flag']}}));
  call.setParamValue('code', 'bad');
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(validityOf(f, 'code'), 'Value didn\'t pass the flag check');
  f.dispose();
});

rules('a warning result (isError: false) rides the notice channel and clears on a clean landing', async () => {
  validator('legacy', (v) => v === 'old' ?
    {message: 'A **legacy** code form', isError: false, isHelper: false} : null);
  const call = callOf(fp('code', 'string', {options: {validators: ['legacy']}}));
  const f = mount(funcForm(call));
  await f.settled;

  typeInto(f, 'code', 'old');
  await until(() => stateOf(f, 'code') !== null && !stateOf(f, 'code').hidden, 'the notice to show');
  const state = stateOf(f, 'code');
  assert.ok(state.classList.contains('u2-param-source-notice'));
  assert.ok(state.textContent.includes('A legacy code form'), 'markup plained');
  assert.equal(validityOf(f, 'code'), null, 'a warning is not an error');

  typeInto(f, 'code', 'new');
  await until(() => stateOf(f, 'code').hidden, 'the clean landing to clear the notice');
  f.dispose();
});

rules('an unknown validator name lands inline as misconfigured — never a form-build throw (#21)', async () => {
  const call = callOf(fp('code', 'string', {options: {validators: ['nope']}}));
  call.setParamValue('code', 'abc');
  const warn = console.warn;
  let warns = 0;
  console.warn = () => warns++;
  let f;
  try {
    f = mount(funcForm(call));
    await f.settled;
  } finally {
    console.warn = warn;
  }
  assert.equal(validityOf(f, 'code'), 'Couldn\'t validate (validator \'nope\' is misconfigured)');
  const error = f.getInput('code').root.querySelector('.u2-input-error');
  assert.equal(error.title, 'Validator "nope" not found', 'the full message rides the title');
  assert.equal(warns, 1, 'the underlying message warned once');
  f.dispose();
});

rules('an async-registered validator func rejects as misconfigured, full reason on the title (#21)', async () => {
  validator('slowV', () => true, {sync: false});
  const call = callOf(fp('code', 'string', {options: {validators: ['slowV']}}));
  call.setParamValue('code', 'abc');
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(validityOf(f, 'code'), 'Couldn\'t validate (validator \'slowV\' is misconfigured)');
  assert.match(f.getInput('code').root.querySelector('.u2-input-error').title, /must be sync/);
  f.dispose();
});

rules('stale-drop: a slow first run loses to a fast second and never lands late', async () => {
  // the verdict encodes invocation ORDER: the slow first request evaluates LAST (run:2),
  // so a broken drop would overwrite the fast second's verdict after it landed
  let calls = 0;
  validator('orderV', () => {
    calls++;
    return 'run:' + calls;
  });
  const call = callOf(fp('code', 'string', {options: {validators: ['orderV']}}));
  const f = mount(funcForm(call));
  await f.settled;

  call.dart.evalDelayMs = 60;
  typeInto(f, 'code', 'a');
  call.dart.evalDelayMs = 1;
  typeInto(f, 'code', 'ab');
  await until(() => validityOf(f, 'code') === 'run:1', 'the fast second landing');
  await sleep(120);
  assert.equal(calls, 2, 'the slow first request still evaluated');
  assert.equal(validityOf(f, 'code'), 'run:1', 'but its landing was dropped, not applied late');
  f.dispose();
});

rules('emptying the value clears the named verdict synchronously, without an evaluation', async () => {
  let calls = 0;
  validator('countV', () => {
    calls++;
    return 'always wrong';
  });
  const call = callOf(fp('code', 'string', {options: {validators: ['countV']}}));
  const f = mount(funcForm(call));
  await f.settled;
  typeInto(f, 'code', 'abc');
  await until(() => validityOf(f, 'code') === 'always wrong', 'the verdict');
  assert.equal(calls, 1);

  typeInto(f, 'code', '');
  assert.equal(validityOf(f, 'code'), null, 'cleared in the same tick — empty is gated, not evaluated');
  await sleep(20);
  assert.equal(calls, 1, 'no evaluation for the empty value');
  f.dispose();
});

rules('settled includes the initial named runs (delay knob honored)', async () => {
  validator('isCode', (v) => String(v).startsWith('X') ? null : 'Code must start with X');
  const call = callOf(fp('code', 'string', {options: {validators: ['isCode']}}));
  call.setParamValue('code', 'abc');
  call.dart.evalDelayMs = 40;
  const f = mount(funcForm(call));
  assert.equal(validityOf(f, 'code'), null, 'not landed yet at construction');
  await f.settled;
  assert.equal(validityOf(f, 'code'), 'Code must start with X', 'landed within settled');
  f.dispose();
});

rules('a rebind re-runs the named validators once per field over the new call\'s values', async () => {
  let runs = 0;
  validator('countV', (v) => {
    runs++;
    return v === 'bad' ? 'no good' : null;
  });
  const call = callOf(fp('code', 'string', {options: {validators: ['countV']}}));
  call.setParamValue('code', 'bad');
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(runs, 1);
  assert.equal(validityOf(f, 'code'), 'no good');

  const next = callOf(fp('code', 'string', {options: {validators: ['countV']}}));
  next.setParamValue('code', 'fine');
  f.source = next;
  await f.settled;
  assert.equal(runs, 2, 'exactly one re-run per named field');
  assert.equal(validityOf(f, 'code'), null);
  f.dispose();
});

rules('an external two-way write re-runs the named validator (the _refreshField tail)', async () => {
  validator('isCode', (v) => String(v).startsWith('X') ? null : 'Code must start with X');
  const call = callOf(fp('code', 'string', {options: {validators: ['isCode']}}));
  const f = mount(funcForm(call));
  await f.settled;
  call.setParamValue('code', 'abc');
  await until(() => validityOf(f, 'code') === 'Code must start with X', 'the external write\'s verdict');
  f.dispose();
});

rules('an active flip clears the named verdict when hidden and re-runs when shown', async () => {
  validator('isCode', (v) => String(v).startsWith('X') ? null : 'Code must start with X');
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('code', 'string', {options: {validators: ['isCode'], visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  await f.settled;
  typeInto(f, 'code', 'abc');
  await until(() => validityOf(f, 'code') === 'Code must start with X', 'the verdict');

  f.getInput('type').value.value = 'Electric';
  assert.equal(validityOf(f, 'code'), null, 'hidden: exempt immediately, no landing awaited');
  f.getInput('type').value.value = 'ICE';
  await until(() => validityOf(f, 'code') === 'Code must start with X', 'the verdict to return');
  f.dispose();
});

// --- feature-detected inertness (W1-shape calls) ---

bare('a W1-shape call — none of the new members, no scriptSync — builds and warns once, inert', async () => {
  const visible = new FuncParam('gears', 'int', {options: {visible: 'type == "ICE"'}});
  const named = new FuncParam('code', 'string', {options: {validators: ['isCode']}});
  const params = [visible, named].map((property) => ({
    name: property.name, value: property.name === 'code' ? 'abc' : undefined, property,
    onChanged: {subscribe: () => ({unsubscribe() {}})},
  }));
  const call = {inputParams: {values: () => params.slice()}, setParamValue() {}};
  let f;
  const warns = countWarns(() => f = mount(funcForm(call)));
  await f.settled;
  assert.equal(hiddenOf(f, 'gears'), false, 'missing seam keeps the default visibility');
  assert.equal(validityOf(f, 'code'), null, 'missing evalParamValidators is inert');
  assert.equal(validityOf(f, 'gears'), null);
  assert.ok(warns >= 1, 'warned, not thrown');
  const again = countWarns(() => {
    f.getInput('gears').value.value = 3;
  });
  assert.equal(again, 0, 'warn-once per key: a later edit is silent');
  f.dispose();
});

test('the stub never leaks across tests: teardown deleted the global', () => {
  assert.equal('grok_ScriptSync' in globalThis, false);
  assert.equal(globalThis.grok_ScriptSync, undefined);
});

// --- categoryGroups (#18, #19, #22) ---

const GROUPS = '{"Power": ["Engine", "Battery"], "Details": ["Contact"]}';

function grouped() {
  return callWith({categoryGroups: GROUPS},
    fp('cylinders', 'int', {category: 'Engine'}),
    fp('capacity', 'double', {category: 'Battery'}),
    fp('email', 'string', {category: 'Contact'}),
    fp('notes', 'string', {category: 'Extra'}));
}

function levelOf(el) {
  return ['u2-form-category-l1', 'u2-form-category-l2', 'u2-form-category-l3']
    .find((c) => el.classList.contains(c)) ?? null;
}

rules('the plan renders headers in order with level classes; unlisted categories append (#18)', () => {
  const f = mount(funcForm(grouped()));
  const headers = headerEls(f);
  assert.deepEqual(headers.map((h) => h.textContent),
    ['Power', 'Engine', 'Battery', 'Details', 'Contact', 'Extra']);
  assert.deepEqual(headers.map(levelOf), [
    'u2-form-category-l1', 'u2-form-category-l2', 'u2-form-category-l2',
    'u2-form-category-l1', 'u2-form-category-l2', null],
  'group headers carry their depth; the appended category keeps the flat style');
  assert.deepEqual(f.inputs.map((i) => i.name), ['cylinders', 'capacity', 'email', 'notes'],
    'fields follow the plan order, the unlisted param present — Dart drops it (defect #12)');
  f.dispose();
});

rules('only field-owning categories collapse — a fieldless group header renders plain', () => {
  const f = mount(funcForm(grouped()));
  const power = headerOf(f, 'Power');
  assert.equal(power.getAttribute('role'), null, 'a group heading takes no clicks');
  assert.equal(power.querySelector('.u2-section-chevron'), null, 'and shows no chevron');
  const engine = headerOf(f, 'Engine');
  assert.equal(engine.getAttribute('role'), 'button');
  assert.notEqual(engine.querySelector('.u2-section-chevron'), null);
  f.dispose();

  const dup = mount(funcForm(callWith({categoryGroups: '{"Power": ["Engine", "Engine"]}'},
    fp('cylinders', 'int', {category: 'Engine'}))));
  const engines = headerEls(dup).filter((h) => h.textContent === 'Engine');
  assert.equal(engines[0].getAttribute('role'), 'button', 'the field-owning first render collapses');
  assert.equal(engines[1].getAttribute('role'), null, 'the fieldless duplicate is a plain heading');
  dup.dispose();
});

rules('a header item whose name matches a category renders that category\'s fields (fpe:849)', () => {
  const f = mount(funcForm(callWith({categoryGroups: '{"Engine": ["Sub"]}'},
    fp('cylinders', 'int', {category: 'Engine'}),
    fp('bore', 'int', {category: 'Sub'}))));
  assert.deepEqual(headerEls(f).map((h) => h.textContent), ['Engine', 'Sub']);
  assert.deepEqual(f.inputs.map((i) => i.name), ['cylinders', 'bore']);
  f.dispose();
});

rules('a category listed twice renders its fields once (#22b)', () => {
  const f = mount(funcForm(callWith({categoryGroups: '{"Power": ["Engine", "Engine"]}'},
    fp('cylinders', 'int', {category: 'Engine'}))));
  assert.deepEqual(f.inputs.map((i) => i.name), ['cylinders'], 'no duplicate bound inputs (Dart duplicates)');
  const engines = headerEls(f).filter((h) => h.textContent === 'Engine');
  assert.equal(engines.length, 2, 'the plan still emits both headers');
  assert.equal(engines[1].hidden, true, 'the fieldless duplicate auto-hides');
  f.dispose();
});

rules('malformed categoryGroups JSON falls back to the flat rendering with one warn (#22a)', () => {
  let f;
  const warns = countWarns(() => f = mount(funcForm(callWith({categoryGroups: '{oops'},
    fp('cylinders', 'int', {category: 'Engine'}),
    fp('email', 'string', {category: 'Contact'})))));
  assert.equal(warns, 1, 'warned, where Dart crashes the whole build (defect #14)');
  const headers = headerEls(f);
  assert.deepEqual(headers.map((h) => h.textContent), ['Engine', 'Contact'], 'W1 flat order');
  assert.deepEqual(headers.map(levelOf), [null, null], 'no level classes on the flat path');
  f.dispose();
});

rules('auto-hide: a header follows its fields, groups propagate up, and everything re-shows (#19)', () => {
  const call = callWith({categoryGroups: GROUPS},
    fp('type', 'string', {choices: ['ICE', 'Electric'], category: 'Contact'}),
    fp('cylinders', 'int', {category: 'Engine', options: {visible: 'type == "ICE"'}}),
    fp('capacity', 'double', {category: 'Battery', options: {visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  assert.equal(headerOf(f, 'Engine').hidden, false);
  assert.equal(headerOf(f, 'Power').hidden, false);

  f.getInput('type').value.value = 'Electric';
  assert.equal(hiddenOf(f, 'cylinders'), true);
  assert.equal(headerOf(f, 'Engine').hidden, true, 'the leaf header follows its only field');
  assert.equal(headerOf(f, 'Battery').hidden, true);
  assert.equal(headerOf(f, 'Power').hidden, true, 'the group propagates up');
  assert.equal(headerOf(f, 'Details').hidden, false, 'the visible branch stays');
  assert.equal(headerOf(f, 'Contact').hidden, false);

  f.getInput('type').value.value = 'ICE';
  assert.equal(headerOf(f, 'Engine').hidden, false, 're-shown');
  assert.equal(headerOf(f, 'Power').hidden, false);
  f.dispose();
});

rules('flat-mode headers auto-hide too when a category\'s only field hides', () => {
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric'], category: 'Main'}),
    fp('cylinders', 'int', {category: 'Engine', options: {visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  const f = mount(funcForm(call));
  assert.equal(headerOf(f, 'Engine').hidden, false);
  f.getInput('type').value.value = 'Electric';
  assert.equal(headerOf(f, 'Engine').hidden, true);
  assert.equal(headerOf(f, 'Main').hidden, false);
  f.dispose();
});

/* The dom-shim's MutationObserver never delivers attribute mutations, so the observer wiring
   itself is e2e-pinned (func-form-expressions/9); these pin the recompute over direct
   `root.hidden` writes by driving `_updateHeaders` directly. */

rules('a direct root.hidden write collapses the header on recompute, and a re-show restores it', () => {
  const f = mount(funcForm(callOf(
    fp('cylinders', 'int', {category: 'Engine'}),
    fp('email', 'string', {category: 'Contact'}))));
  assert.equal(headerOf(f, 'Engine').hidden, false);
  f.getInput('cylinders').root.hidden = true;
  f._updateHeaders();
  assert.equal(headerOf(f, 'Engine').hidden, true, 'the header follows a hide _applyVisible never saw');
  assert.equal(headerOf(f, 'Contact').hidden, false);
  f.getInput('cylinders').root.hidden = false;
  f._updateHeaders();
  assert.equal(headerOf(f, 'Engine').hidden, false);
  f.dispose();
});

rules('direct hides of every field under a group collapse the leaf and the group header', () => {
  const f = mount(funcForm(callWith({categoryGroups: '{"Power": ["Engine"]}'},
    fp('cylinders', 'int', {category: 'Engine'}),
    fp('capacity', 'double', {category: 'Engine'}))));
  f.getInput('cylinders').root.hidden = true;
  f.getInput('capacity').root.hidden = true;
  f._updateHeaders();
  assert.equal(headerOf(f, 'Engine').hidden, true);
  assert.equal(headerOf(f, 'Power').hidden, true, 'the fieldless group follows its emptied leaf');
  f.getInput('capacity').root.hidden = false;
  f._updateHeaders();
  assert.equal(headerOf(f, 'Engine').hidden, false);
  assert.equal(headerOf(f, 'Power').hidden, false);
  f.dispose();
});

// --- the verification-strategy traps ---

rules('echo audit: one keystroke → exactly N scriptSync evals, ZERO setParamValue from the rules layer', () => {
  let evals = 0;
  globalThis.grok_ScriptSync = (e, v) => {
    evals++;
    return scriptSyncStub(e, v);
  };
  const call = vehicle();
  const sets = [];
  const setParamValue = call.setParamValue.bind(call);
  call.setParamValue = (name, value) => {
    sets.push([name, value]);
    setParamValue(name, value);
  };
  const f = mount(funcForm(call));
  // vehicle rules: cylinders:visible, tankExtension:visible+enabled, batteryCapacity:visible
  const base = evals;
  sets.length = 0;
  f.getInput('type').value.value = 'Electric';
  assert.equal(evals - base, 4, 'one eval per expression, once');
  assert.deepEqual(sets, [['type', 'Electric']], 'the edit itself and nothing else');
  f.dispose();
});

rules('a revalidate storm never re-triggers expression evaluation', () => {
  let evals = 0;
  globalThis.grok_ScriptSync = (e, v) => {
    evals++;
    return scriptSyncStub(e, v);
  };
  const f = mount(funcForm(vehicle()));
  const base = evals;
  for (let i = 0; i < 10; i++)
    f.getInput('cylinders').revalidate();
  assert.equal(evals, base, 'revalidate writes _validity only');
  f.dispose();
});

rules('W2 debounce interplay: a choices landing after a hide never resurrects the hidden field\'s error', async () => {
  const provider = new Func('cities', {inputs: [], run: async () => ['Boston', 'Dallas']});
  Func.registry.push(provider);
  const call = callOf(
    fp('type', 'string', {choices: ['ICE', 'Electric']}),
    fp('city', 'string', {nullable: false,
      options: {choices: 'cities()', visible: 'type == "ICE"'}}));
  call.setParamValue('type', 'ICE');
  call.setParamValue('city', 'Paris');
  call.dart.evalDelayMs = 40;
  const f = mount(funcForm(call));
  // hide BEFORE the choices land; the landing prunes 'Paris' through the field effect
  f.getInput('type').value.value = 'Electric';
  assert.equal(hiddenOf(f, 'city'), true);
  await f.settled;
  await until(() => call.inputs.city === null, 'the prune to reach the call');
  assert.equal(validityOf(f, 'city'), null, 'the active guard held through the late landing');

  f.getInput('type').value.value = 'ICE';
  assert.equal(validityOf(f, 'city'), 'Value can\'t be empty', 'the verdict is back once shown');
  f.dispose();
});

// --- doubles honesty ---

rules('doubles: evalParamValidators never delivers synchronously and maps verdicts like Dart', async () => {
  validator('flag', (v) => v !== 'bad');
  const call = callOf(fp('code', 'string', {options: {validators: ['flag', 'notEmpty']}}));
  call.setParamValue('code', 'bad');
  let landed = null;
  const pending = call.evalParamValidators('code').then((r) => landed = r);
  assert.equal(landed, null, 'delivery is never synchronous');
  await pending;
  assert.deepEqual(landed, [{message: 'Value didn\'t pass the flag check', isError: true, isHelper: false}]);

  call.setParamValue('code', 'fine');
  assert.deepEqual(await call.evalParamValidators('code'), [], 'passing validators are omitted');

  call.setParamValue('code', '');
  assert.deepEqual(await call.evalParamValidators('code'),
    [{message: 'Value cannot be empty.', isError: true, isHelper: false}],
    'the notEmpty builtin answers without a registry entry');
});

rules('doubles: a missing validator rejects; {sync: false} rejects with the must-be-sync shape', async () => {
  validator('slowV', () => true, {sync: false});
  const call = callOf(
    fp('a', 'string', {options: {validators: ['ghost']}}),
    fp('b', 'string', {options: {validators: ['slowV']}}));
  call.setParamValue('a', 'x');
  call.setParamValue('b', 'x');
  await assert.rejects(() => call.evalParamValidators('a'), /Validator "ghost" not found/);
  await assert.rejects(() => call.evalParamValidators('b'), /must be sync/);
});
