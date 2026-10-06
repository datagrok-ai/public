/* FunctionInput (dg): the tiny control over an nqName string, the FunctionsBrowser popup with
   click / activate commits, dismissal without a write, the `u2-function-input` platform meta, and
   the FunctionName semtype rule in the Editors registry. `DG`/`grok` come from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Registry} from '../src/spec/registry.js';
import {Func, Property} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {FunctionInput, functionInput} = await import('../src/dg/inputs/function-input.js');
const {propertyForm} = await import('../src/dg/forms/object-form.js');
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');

function ui(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    Func.registry = [
      new Func('MolWeight', {friendlyName: 'Molecular Weight', tags: ['Chem'],
        description: 'Molecular weight of a molecule',
        inputs: [new Property('smiles', 'string')], outputs: [new Property('mw', 'double')]}),
      new Func('Sketch', {namespace: 'Chem', tags: ['chem'],
        outputs: [new Property('view', 'view')]}),
      new Func('Sin', {tags: ['Math'],
        inputs: [new Property('x', 'double')], outputs: [new Property('y', 'double')]}),
    ];
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

function make(options = {}) {
  const input = new FunctionInput({label: 'Function', ...options});
  document.body.append(input.root);
  return input;
}

function control(input) {
  return input.root.querySelector('.u2-function-input');
}

function popup() {
  return document.body.querySelector('.u2-function-input-popup');
}

function open(input) {
  fire(control(input), 'click');
  const box = popup();
  const list = box.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  return box;
}

function rows(box) {
  return box.querySelectorAll('.u2-list-row[data-u2-func]');
}

ui('the control shows the friendly name of a resolvable value, the raw string otherwise', () => {
  const input = make({value: 'MolWeight'});
  assert.equal(input.root.dataset.u2, 'function-input');
  assert.equal(control(input).getAttribute('role'), 'button');
  assert.equal(control(input).getAttribute('aria-expanded'), 'false');
  const label = input.root.querySelector('.u2-function-input-name');
  assert.equal(label.textContent, 'Molecular Weight');
  assert.equal(label.title, 'MolWeight', 'the tooltip keeps the qualified name');
  input.value.value = 'Gone:Missing';
  assert.equal(input.root.querySelector('.u2-function-input-name').textContent, 'Gone:Missing');
  input.value.value = '';
  const name = input.root.querySelector('.u2-function-input-name');
  assert.equal(name.textContent, 'Pick a function…');
  assert.equal(name.classList.contains('u2-function-input-empty'), true);
  input.dispose();
});

ui('a row click commits the qualified name and closes; the popup hosts the compact browser', () => {
  const changes = [];
  const input = make({onChanged: (v) => changes.push(v)});
  const box = open(input);
  assert.equal(control(input).getAttribute('aria-expanded'), 'true');
  assert.notEqual(box.querySelector('[data-u2="functions-browser"]'), null);
  assert.equal(box.querySelector('[data-u2="fb-panes"]').hidden, true, 'no filter panes in the popup');
  const sketch = rows(box).find((r) => r.getAttribute('data-u2-func') === 'Chem:Sketch');
  fire(sketch, 'click');
  assert.equal(input.value.value, 'Chem:Sketch');
  assert.deepEqual(changes, ['Chem:Sketch']);
  assert.equal(popup(), null);
  assert.equal(control(input).getAttribute('aria-expanded'), 'false');
  input.dispose();
});

ui('double-click activates too; Esc and an outside click close without a write', () => {
  const input = make({value: 'Sin'});
  let box = open(input);
  const first = rows(box)[0];
  fire(first, 'dblclick');
  assert.equal(input.value.value, first.getAttribute('data-u2-func'));
  assert.equal(popup(), null);

  const committed = input.value.value;
  box = open(input);
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(popup(), null, 'Escape closes');
  box = open(input);
  fire(document.body, 'pointerdown');
  assert.equal(popup(), null, 'outside pointerdown closes');
  assert.equal(input.value.value, committed, 'neither dismissal wrote');
  input.dispose();
  assert.equal(popup(), null, 'disposing the input takes the popup with it');
});

ui('Backspace clears a nullable value; a required one keeps it; a disabled input does not open', () => {
  const input = make({value: 'Sin'});
  fire(control(input), 'keydown', {key: 'Backspace'});
  assert.equal(input.value.value, '');
  const required = make({value: 'Sin', nullable: false});
  fire(control(required), 'keydown', {key: 'Delete'});
  assert.equal(required.value.value, 'Sin');
  required.enabled = false;
  fire(control(required), 'click');
  assert.equal(popup(), null);
  input.dispose();
  required.dispose();
});

ui('a bound signal is adopted; the factory carries the label', () => {
  const bound = signal('MolWeight');
  const input = functionInput('Fn', {bind: bound});
  document.body.append(input.root);
  const box = open(input);
  fire(rows(box).find((r) => r.getAttribute('data-u2-func') === 'Sin'), 'click');
  assert.equal(bound.value, 'Sin');
  assert.equal(input.root.querySelector('.u2-input-label').textContent, 'Fn');
  input.dispose();
});

ui('the FunctionName semtype resolves to a FunctionInput in schema-driven forms', () => {
  const values = {handler: 'MolWeight'};
  const form = propertyForm([{name: 'handler', type: 'string', semType: 'FunctionName',
    get: (t) => t.handler, set: (t, v) => t.handler = v}], values);
  document.body.append(form.root);
  const editor = form.root.querySelector('[data-u2="function-input"]');
  assert.notEqual(editor, null);
  assert.equal(editor.querySelector('.u2-function-input-name').textContent, 'Molecular Weight');
  form.dispose();
});

ui('registerPlatformComponents carries the u2-function-input meta, and create() builds it', () => {
  const reg = new Registry();
  registerPlatformComponents(reg);
  const meta = reg.get('u2-function-input');
  assert.equal(meta.category, 'Inputs');
  const bound = signal('Sin');
  const input = meta.create({label: 'Fn', value: bound, scalarOnly: true});
  document.body.append(input.root);
  assert.equal(input.root.dataset.u2, 'function-input');
  assert.equal(input.root.querySelector('.u2-function-input-name').textContent, 'Sin');
  const box = open(input);
  assert.deepEqual(rows(box).map((r) => r.getAttribute('data-u2-func')).sort(),
    ['MolWeight', 'Sin'], 'scalarOnly leaves the view-typed output out');
  fire(document, 'keydown', {key: 'Escape'});
  input.dispose();
});

ui('opening and closing repeatedly leaves no scope behind', () => {
  const input = make();
  const live = Scope.liveCount;
  for (let i = 0; i < 5; i++) {
    open(input);
    fire(document.body, 'pointerdown');
  }
  assert.equal(Scope.liveCount, live);
  input.dispose();
});
