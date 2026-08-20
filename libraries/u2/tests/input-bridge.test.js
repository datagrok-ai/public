/* Outbound bridge tests — the `dartInputFor` contract (plan WO-7.1). The platform base class and
   the `DG.TYPE` constants come from tests/dg-stub.mjs; everything under test is the shipped
   `DartInput`. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextInput} from '../src/components/text-input.js';
import {NumberInput} from '../src/components/number-input.js';
import {ChoiceInput} from '../src/components/choice-input.js';
import {inputForProperty} from '../src/dg/object-form.js';

register('./dg-stub.mjs', import.meta.url);
const {asDartInput, dartInputFor} = await import('../src/dg/input-bridge.js');
const {fileInput} = await import('../src/dg/file-input.js');

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function bridge(name, body) {
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

function type(input, text) {
  const el = input.root.querySelector('input');
  el.value = text;
  fire(el, 'input');
}

bridge('the editor is built on the first getInput(), with the property the platform bound', () => {
  const seen = [];
  const dart = dartInputFor((prop) => {
    seen.push(prop);
    return new TextInput({});
  });
  assert.equal(seen.length, 0, 'nothing is built before the platform asks');
  assert.equal(dart.u2, null);

  dart.getInput();
  assert.deepEqual(seen, [null], 'nothing bound — the builder gets a null property');
  assert.ok(dart.u2 instanceof TextInput);
  dart.detach();

  const prop = {name: 'dose', propertyType: 'double', min: 0, max: 10, units: 'mg',
    get: (x) => x.dose, set: (x, v) => x.dose = v};
  const bound = dartInputFor((p) => inputForProperty(p));
  bound.property = prop;
  bound.getInput();
  assert.ok(bound.u2 instanceof NumberInput, 'the bound property picks the editor');
  assert.ok(bound.u2.root.querySelector('.u2-number-slider'), 'and carries its metadata');
  bound.detach();
});

bridge('a value set before the build lands in the editor when it appears', () => {
  const dart = dartInputFor(() => new TextInput({}));
  dart.setValue('aspirin');
  assert.equal(dart.getValue(), 'aspirin', 'readable while the editor is still missing');

  dart.getInput();
  assert.equal(dart.u2.value.peek(), 'aspirin');
  assert.equal(dart.getValue(), 'aspirin');
  assert.equal(dart.u2.root.querySelector('input').value, 'aspirin');
  dart.detach();
});

bridge('getInput() twice returns the same editor and never rebuilds', () => {
  let builds = 0;
  const dart = dartInputFor(() => {
    builds++;
    return new TextInput({});
  });
  const first = dart.getInput();
  const second = dart.getInput();
  assert.equal(builds, 1);
  assert.equal(first, second);
  assert.equal(first, dart.u2.root.querySelector('.u2-input-box'));
  dart.detach();
});

bridge('getInput() hands over the editor WITH its options rail, so units survive the bridge', () => {
  const dart = dartInputFor(() => new NumberInput({label: 'Dose', postfix: 'mg', clicker: true,
    min: 0, max: 10}));
  const editor = dart.getInput();
  assert.equal(editor.classList.contains('u2-input-box'), true);
  assert.equal(editor.querySelector('.u2-input-postfix').textContent, 'mg');
  assert.equal(editor.querySelectorAll('.u2-number-click').length, 2);
  assert.equal(editor.querySelector('.u2-input-editor').classList.contains('u2-number'), true,
    'and the editor selector still points at the editor alone');
  assert.equal(editor.querySelector('.u2-input-label'), null, 'the caption stays the platform\'s');
  dart.detach();
});

/** …but not the units of a property the platform bound: the Dart host renders those itself into
 * the root it wraps this editor in (`input_base.dart:219-222`), so keeping them here would say
 * 'mg' twice. A form-generated input has no such host and keeps its own. */
bridge('units of a bound property are left to the Dart host', () => {
  const prop = {name: 'dose', propertyType: 'double', min: 0, max: 10, units: 'mg',
    showPlusMinus: true, get: (x) => x.dose, set: (x, v) => x.dose = v};
  const dart = dartInputFor((p) => inputForProperty(p, {assumeWritable: true}));
  dart.property = prop;
  const editor = dart.getInput();
  assert.equal(editor.querySelector('.u2-input-postfix'), null, 'no second postfix');
  assert.equal(editor.querySelectorAll('.u2-number-click').length, 2, 'the rest of the rail stays');
  dart.detach();

  const bare = dartInputFor((p) => inputForProperty(p, {assumeWritable: true}));
  bare.property = {name: 'dose', propertyType: 'double', min: 0, max: 10, units: 'mg'};
  assert.equal(bare.getInput().querySelector('.u2-input-options'), null,
    'a rail left with nothing on it goes with the postfix');
  bare.detach();

  const generated = inputForProperty(prop);
  assert.equal(generated.root.querySelector('.u2-input-postfix').textContent, 'mg',
    'the generated input keeps its postfix — nobody else renders it there');
  generated.dispose();
});

bridge('a dg input rides the same rail: the file picker crosses the bridge too', () => {
  const dart = dartInputFor(() => fileInput('Path'));
  const editor = dart.getInput();
  assert.equal(editor.classList.contains('u2-input-box'), true);
  const rail = editor.querySelector('.u2-input-options');
  assert.equal(rail.querySelectorAll('.u2-btn').length, 1, 'the folder-open icon');
  assert.equal(editor.querySelectorAll('input[type=file]').length, 1,
    'and the hidden picker, which crosses with the editor it sits in');
  dart.detach();
});

bridge('the deferred editor belongs to the bridge: detach() disposes it', () => {
  const dart = dartInputFor(() => new TextInput({}));
  dart.getInput();
  const input = dart.u2;
  dart.detach();
  assert.equal(input.scope.isDisposed, true, 'built under the bridge scope, released with it');
});

bridge('an input handed over is not owned: detach() leaves it alive', () => {
  const input = new TextInput({label: 'Name'});
  const dart = asDartInput(input);
  dart.detach();
  assert.equal(input.scope.isDisposed, false);
  input.dispose();
});

bridge('a user edit fires input and changed, a programmatic one only changed', async () => {
  const dart = dartInputFor(() => new TextInput({}));
  dart.getInput();
  type(dart.u2, 'ibuprofen');
  assert.deepEqual(dart.fired, {input: 0, changed: 0}, 'nothing fires from inside the value effect');
  await flush();
  assert.deepEqual(dart.fired, {input: 1, changed: 1});

  dart.setValue('paracetamol');
  await flush();
  assert.deepEqual(dart.fired, {input: 1, changed: 2});
  dart.detach();
});

bridge('one edit is one pair of events, however many writes it took', async () => {
  const dart = dartInputFor(() => new TextInput({}));
  dart.getInput();
  const input = dart.u2;
  input.value.value = 'a';
  input.value.value = 'ab';
  input.value.value = 'abc';
  await flush();
  assert.deepEqual(dart.fired, {input: 1, changed: 1});
  assert.equal(dart.getValue(), 'abc');
  dart.detach();
});

bridge('a subscriber that writes the signal back cannot trip the cycle guard', async () => {
  const dart = dartInputFor(() => new TextInput({}));
  dart.getInput();
  let hops = 0;
  // the shape a cross-synced form has: the platform subscriber writes a u2 signal back, and the
  // values never compare equal (a datetime crosses the bridge as a fresh object every time). Fired
  // from inside the value effect, the round trips would pile up in one batch until the signal
  // runtime aborted it — well before 100 of them
  dart.fireChanged = () => {
    dart.fired.changed++;
    if (hops < 150)
      dart.u2.value.value = `hop ${hops++}`;
  };
  type(dart.u2, 'start');
  await flush();
  assert.equal(hops, 150, 'every round trip ran');
  assert.equal(dart.u2.value.peek(), 'hop 149');
  assert.equal(dart.fired.changed, 151, 'one event per round trip, none of them lost');

  type(dart.u2, 'after');
  await flush();
  assert.equal(dart.fired.changed, 152, 'and the effect is still alive');
  dart.detach();
});

bridge('a value the input prunes on its own is not reported as a user edit', async () => {
  const items = new ChoiceInput({items: ['age', 'sex']});
  const dart = asDartInput(items);
  items.value.value = 'sex';
  await flush();
  assert.deepEqual(dart.fired, {input: 1, changed: 1});

  items.setItems(['age']);
  await flush();
  assert.equal(items.value.peek(), null, 'the item it stood for is gone');
  assert.deepEqual(dart.fired, {input: 1, changed: 2}, 'changed alone, as Dart does');
  dart.detach();
  items.dispose();
});

bridge('the validator follows the dart handle the platform re-homes the input on', async () => {
  const dart = dartInputFor(() => new NumberInput({min: 1, max: 10}));
  const first = dart.dart;
  // what `JsInputProxy.fromFunc` does before asking for the element
  dart.dart = {validators: []};
  dart.getInput();
  assert.equal(first.validators.length, 1, 'the handle built by JsInputBase kept its own');
  assert.equal(dart.dart.validators.length, 1, 'and the one the form binds got one too');

  type(dart.u2, '99');
  await flush();
  assert.match(dart.dart.validators[0](), /at most 10/, 'it reads the current validity');
  dart.getInput();
  assert.equal(dart.dart.validators.length, 1, 'and is added once per handle');
  dart.detach();
});

bridge('dataType comes from the component, not from the value', () => {
  const dart = dartInputFor(() => new NumberInput({}));
  assert.equal(dart.dataType, 'object', 'nothing built yet, nothing to infer from');
  dart.getInput();
  assert.equal(dart.dataType, 'double');
  dart.detach();

  const typed = dartInputFor(() => new TextInput({}), {dataType: 'string'});
  assert.equal(typed.dataType, 'string', 'the explicit type holds before the build');
  typed.detach();
});
