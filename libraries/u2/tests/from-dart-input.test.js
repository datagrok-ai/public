/* Inbound bridge tests. datagrok-api cannot load in node, so the platform input here is a minimal
   fake of the `DartInputLike` shape the bridge consumes — a real `DG.InputBase` satisfies it
   structurally (checked at compile time by src/dg/object-form.ts's auto mode) — and the core's
   validation surface is faked through the same globals the bridge feature-detects. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Form} from '../src/components/form.js';
import {TextInput} from '../src/components/text-input.js';
import {fromDartInput, PlatformInput} from '../src/dg/from-dart-input.js';
import {propertyForm} from '../src/dg/object-form.js';

/** Every test runs against a clean document, without leftover globals, and must leave the
 * live-scope count where it was. */
function bridge(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      delete globalThis.grok_InputBase_Get_ValidationMessages;
      delete globalThis.grok_InputBase_OnValidated;
      delete globalThis.DG;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** A platform input: its own root with the platform's label and editor, a value that fires
 * `onChanged` after it is stored (as `fireChanged` does), and a write counter. */
function dartInput(options = {}) {
  const root = document.createElement('div');
  root.className = 'ui-input-root';
  const label = document.createElement('label');
  label.className = 'ui-input-label';
  label.textContent = options.caption ?? 'Name';
  label.offsetWidth = options.labelWidth ?? 0;
  const editor = document.createElement('input');
  editor.className = 'ui-input-editor';
  root.append(label, editor);

  const subscribers = [];
  return {
    root, editor, subscribers, writes: 0,
    caption: options.caption ?? 'Name',
    dart: options.dart ?? null,
    _value: options.value ?? '',
    get value() {
      return this._value;
    },
    set value(x) {
      this.writes++;
      this._value = x;
      for (const fn of subscribers.slice())
        fn(x);
    },
    onChanged: {
      subscribe(fn) {
        subscribers.push(fn);
        return {unsubscribe: () => subscribers.splice(subscribers.indexOf(fn), 1)};
      },
    },
  };
}

function dart(messages = []) {
  return {messages, listeners: []};
}

/** The core surface the bridge feature-detects: `grok_InputBase_Get_ValidationMessages` and
 * `grok_InputBase_OnValidated`, added 2026-08-14 and absent from older platforms. */
function installValidation() {
  globalThis.grok_InputBase_Get_ValidationMessages = (x) => x.messages;
  globalThis.grok_InputBase_OnValidated = (x) => ({
    subscribe(fn) {
      x.listeners.push(fn);
      return {unsubscribe: () => x.listeners.splice(x.listeners.indexOf(fn), 1)};
    },
  });
}

function validate(x, messages) {
  x.messages = messages;
  for (const fn of x.listeners.slice())
    fn(messages);
}

/** What the bundler binds `datagrok-api/dg` to in a package build. */
function installDg(forProperty) {
  globalThis.DG = {InputBase: {forProperty}};
}

function mount(component) {
  document.body.append(component.root);
  return component;
}

function prop(name, options = {}) {
  return {
    name, propertyType: 'string',
    get: (x) => x[name],
    set: (x, v) => x[name] = v,
    ...options,
  };
}

bridge('mirrors the value both ways, with the platform value canonical', () => {
  const dg = dartInput({value: 'Aspirin', caption: 'Compound'});
  const input = mount(fromDartInput(dg));
  assert.equal(input instanceof PlatformInput, true);
  assert.equal(input.name, 'Compound', 'the caption is the u2 name');
  assert.equal(input.value.value, 'Aspirin');
  assert.equal(input.root.dataset.u2, 'dart-input');
  assert.equal(dg.root.parentNode, input.box, 'the platform root is the whole editor');
  assert.equal(input.root.querySelector('.u2-input-label'), null, 'no second caption');
  assert.equal(dg.root.classList.contains('u2-input-editor'), false);
  assert.equal(dg.editor.classList.contains('u2-input-editor'), true, 'the skin moved to the editor');

  dg.value = 'Ibuprofen';
  assert.equal(input.value.value, 'Ibuprofen');
  assert.equal(dg.writes, 1, 'a platform edit is not echoed back');

  input.value.value = 'Paracetamol';
  assert.equal(dg.value, 'Paracetamol');
  assert.equal(dg.writes, 2);
  input.value.value = 'Paracetamol';
  assert.equal(dg.writes, 2, 'an unchanged value writes nothing');
  input.dispose();
});

bridge('validity mirrors the platform validation messages and clears with them', () => {
  installValidation();
  const source = dart(['Not a valid SMILES']);
  const dg = dartInput({dart: source});
  const input = mount(fromDartInput(dg));
  assert.equal(input.validity.value, 'Not a valid SMILES', 'seeded from the current messages');
  assert.equal(input.root.querySelector('.u2-input-error'), null, 'the platform renders them itself');

  validate(source, ['Too long', 'Not a valid SMILES']);
  assert.equal(input.validity.value, 'Too long\nNot a valid SMILES');
  validate(source, []);
  assert.equal(input.validity.value, null);
  input.dispose();
  assert.equal(source.listeners.length, 0);
});

bridge('u2 validators keep working next to the platform ones, which come first', () => {
  installValidation();
  const source = dart();
  const dg = dartInput({dart: source, value: ''});
  const input = mount(fromDartInput(dg));
  input.addValidator((value) => value === '' ? 'Value can\'t be empty' : null);
  assert.equal(input.validity.value, 'Value can\'t be empty');

  validate(source, ['Server says no']);
  assert.equal(input.validity.value, 'Server says no');
  validate(source, []);
  assert.equal(input.validity.value, 'Value can\'t be empty', 'the u2 validator runs again');
  dg.value = 'CCO';
  assert.equal(input.validity.value, null);
  input.dispose();
});

bridge('without the core surface only u2 validators apply', () => {
  const dg = dartInput({dart: dart(['invisible']), value: ''});
  const input = mount(fromDartInput(dg));
  assert.equal(input.validity.value, null, 'an older platform is opaque to u2');
  input.addValidator((value) => value === '' ? 'Required' : null);
  assert.equal(input.validity.value, 'Required');
  assert.equal(input.root.querySelector('.u2-input-error').textContent, 'Required',
    'and u2 keeps its own error line, since the platform shows no u2 message');
  input.dispose();

  installValidation();
  const plain = mount(fromDartInput(dartInput({value: 'x'})));
  assert.equal(plain.validity.value, null, 'a dart-less input never reaches the globals');
  plain.dispose();
});

bridge('a platform input in a form flips form validity and joins the label column', async () => {
  installValidation();
  const source = dart();
  const dg = dartInput({dart: source, caption: 'Compound', labelWidth: 80});
  const platform = fromDartInput(dg);
  const plain = new TextInput({label: 'Batch'});
  plain.root.querySelector('.u2-input-label').offsetWidth = 40;
  const form = mount(new Form().addAll([platform, plain]));
  await flush();
  assert.equal(form.root.style.getPropertyValue('--u2-form-label-width'), '80px',
    'the platform label is measured too');

  assert.equal(form.validity.value, null);
  assert.equal(form.validate(), true);
  validate(source, ['Value can\'t be empty']);
  assert.equal(form.validity.value, 'Value can\'t be empty');
  assert.equal(form.validate(), false);
  assert.equal(document.activeElement, dg.editor, 'validate() focuses the platform editor');

  validate(source, []);
  assert.equal(form.validity.value, null);
  form.dispose();
  platform.dispose();
  plain.dispose();
});

bridge('an asDartInput product is unwrapped, never wrapped twice', () => {
  const original = new TextInput({label: 'Barcode'});
  // the brand is DartInput's `u2` reference; building a real one needs datagrok-api, which cannot
  // load in node, so this is the shape asDartInput hands back. The reverse guard (asDartInput of a
  // wrapper throws) lives in input-bridge.ts and needs the platform to be exercised.
  const asDart = {
    root: document.createElement('div'), value: '', caption: 'Barcode', dart: {}, u2: original,
    onChanged: {subscribe: () => ({unsubscribe() {}})},
  };
  assert.equal(fromDartInput(asDart), original, 'the original u2 input comes back');
  assert.equal(asDart.root.parentNode, null, 'and nothing was built around it');
  original.dispose();
});

bridge('dispose severs the bridge and leaves the platform input alone', () => {
  installValidation();
  const source = dart();
  const dg = dartInput({dart: source, value: 'a'});
  const input = mount(fromDartInput(dg));
  input.dispose();
  assert.equal(dg.subscribers.length, 0);
  assert.equal(source.listeners.length, 0);
  assert.equal(dg.root.parentNode, input.box, 'the platform input is never detached');

  dg.value = 'b';
  assert.equal(input.value.value, 'a', 'a severed bridge carries nothing');
  input.value.value = 'c';
  assert.equal(dg.value, 'b');
});

bridge('editors: auto asks the platform for the editor of a real property', () => {
  const dg = dartInput({caption: 'Structure', value: 'CCO'});
  const target = {smiles: 'CCO', name: 'Ethanol'};
  const asked = [];
  installDg((property, source) => {
    asked.push([property.name, source]);
    return dg;
  });
  const changes = [];
  const form = mount(propertyForm([prop('smiles', {dart: {}}), prop('name')], target,
    {editors: 'auto', onChanged: (name, value) => changes.push([name, value])}));
  assert.deepEqual(asked, [['smiles', target]], 'only a property with a dart handle is asked about');
  assert.deepEqual(form.inputs.map((i) => i.root.dataset.u2), ['dart-input', 'text-input']);
  assert.equal(form.inputs[0].name, 'smiles', 'keyed by property name, not the platform caption');
  assert.equal(form.input('smiles'), form.inputs[0]);
  assert.equal(form.getValues().smiles, 'CCO');
  assert.equal(dg.writes, 0, 'the platform editor owns its value, uncoerced');

  dg.value = 'CCC';
  assert.deepEqual(changes, [['smiles', 'CCC']], 'and its edits are reported');

  target.smiles = 'CCCC';
  form.refresh();
  assert.equal(dg.value, 'CCCC', 'refresh re-reads the raw property value');
  assert.deepEqual(changes, [['smiles', 'CCC']], 'refresh stays echo-suppressed');
  form.dispose();
});

bridge('auto mode falls back to the generated editor', () => {
  const target = {smiles: 'CCO'};
  installDg(() => {
    throw new Error('no editor for this property');
  });
  const thrown = mount(propertyForm([prop('smiles', {dart: {}})], target, {editors: 'auto'}));
  assert.deepEqual(thrown.inputs.map((i) => i.root.dataset.u2), ['text-input']);
  thrown.dispose();

  installDg(() => null);
  const none = mount(propertyForm([prop('smiles', {dart: {}})], target, {editors: 'auto'}));
  assert.deepEqual(none.inputs.map((i) => i.root.dataset.u2), ['text-input']);
  none.dispose();

  delete globalThis.DG;
  const offPlatform = mount(propertyForm([prop('smiles', {dart: {}})], target, {editors: 'auto'}));
  assert.deepEqual(offPlatform.inputs.map((i) => i.root.dataset.u2), ['text-input']);
  assert.equal(offPlatform.inputs[0].value.value, 'CCO', 'the generated editor is wired as always');
  offPlatform.dispose();
});

bridge('the default stays generated even for a real property', () => {
  installDg(() => dartInput());
  const form = mount(propertyForm([prop('smiles', {dart: {}})], {smiles: 'CCO'}));
  assert.deepEqual(form.inputs.map((i) => i.root.dataset.u2), ['text-input']);
  form.dispose();
});
