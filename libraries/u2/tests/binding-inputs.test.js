/* Input live chrome (UB-10): label, tooltipText, enabled, postfix and placeholder accept a bound
   signal and update in place — one mechanism on the Input base (`liveOption`), flags flipped once
   in the shared inputProps block. Strictly one-way; literals keep plain-setter ownership. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Signal, signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Input} from '../src/core/input-base.js';
import {TextInput} from '../src/components/inputs/text-input.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function spec(name, body) {
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

function registry() {
  const reg = new Registry();
  registerAll(reg);
  return reg;
}

const CHROME = ['label', 'postfix', 'tooltipText', 'enabled'];

spec('the five chrome props follow bound signals live, input element identity stable', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {
    l: 'Name', tt: 'A tip', en: true, px: 'mg', ph: 'Type here'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-text-input', name: 'driver',
      bind: {label: '$.l', tooltipText: '$.tt', enabled: '$.en', postfix: '$.px', placeholder: '$.ph'}},
  }, ctx, reg);
  const root = instance.node('driver').root;
  const editor = root.querySelector('input');
  assert.equal(root.querySelector('.u2-input-label').textContent, 'Name');
  assert.equal(root.title, 'A tip');
  assert.equal(editor.disabled, false);
  assert.equal(root.querySelector('.u2-input-postfix').textContent, 'mg');
  assert.equal(editor.placeholder, 'Type here');

  ctx.data.l.value = 'Renamed';
  ctx.data.tt.value = 'Another tip';
  ctx.data.en.value = false;
  ctx.data.px.value = 'kg';
  ctx.data.ph.value = 'Now here';
  assert.equal(root.querySelector('.u2-input-label').textContent, 'Renamed');
  assert.equal(root.title, 'Another tip');
  assert.equal(editor.disabled, true, 'enabled=false grays out live');
  assert.equal(root.classList.contains('u2-input-disabled'), true);
  assert.equal(root.querySelector('.u2-input-postfix').textContent, 'kg');
  assert.equal(editor.placeholder, 'Now here');
  assert.equal(root.querySelector('input'), editor, 'the editor was never re-created');
  assert.equal(instance.node('driver').root, root, 'nor was the input');

  ctx.data.en.value = true;
  assert.equal(editor.disabled, false, 'and un-grays live');
  instance.dispose();
});

spec('bound placeholder follows live on text-area, list and qnum editors too', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {ph: 'First'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-text-area', name: 'area', bind: {placeholder: '$.ph'}},
      {tag: 'u2-list-input', name: 'list', bind: {placeholder: '$.ph'}},
      {tag: 'u2-qnum-input', name: 'qnum', bind: {placeholder: '$.ph'}},
    ]},
  }, ctx, reg);
  const area = instance.node('area').root.querySelector('textarea');
  const listField = instance.node('list').root.querySelector('input');
  const listArea = instance.node('list').root.querySelector('textarea');
  const qnum = instance.node('qnum').root.querySelector('input');
  assert.equal(area.placeholder, 'First');
  assert.equal(listField.placeholder, 'First');
  assert.equal(listArea.placeholder, 'First', 'the expanded textarea follows too');
  assert.equal(qnum.placeholder, 'First');

  ctx.data.ph.value = 'Second';
  assert.equal(area.placeholder, 'Second');
  assert.equal(listField.placeholder, 'Second');
  assert.equal(listArea.placeholder, 'Second');
  assert.equal(qnum.placeholder, 'Second');
  instance.dispose();
});

spec('one-way by declaration: no chrome prop carries twoWay on any input meta', () => {
  const reg = registry();
  for (const meta of reg.metas()) {
    for (const name of [...CHROME, 'placeholder']) {
      const prop = meta.props.find((p) => p.name === name);
      if (prop !== undefined)
        assert.equal(prop.twoWay, undefined, `${meta.tag}.${name} must never be twoWay`);
    }
    // the deliberate exclusions: name and inline stay re-render tier
    for (const name of ['name', 'inline']) {
      const prop = meta.props.find((p) => p.name === name);
      if (prop !== undefined && meta.props.some((p) => p.name === 'value' && p.bindable))
        assert.equal(prop.bindable, undefined, `${meta.tag}.${name} stays re-render tier`);
    }
  }
  // the flags themselves, on a representative meta
  const text = reg.get('u2-text-input');
  for (const name of [...CHROME, 'placeholder'])
    assert.equal(text.props.find((p) => p.name === name)?.bindable, true, `${name} is live`);
});

spec('a deferred link wires chrome bound to a later node, through the _liveBinds bindStep', () => {
  const reg = registry();
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'div', children: [
      {tag: 'u2-text-input', name: 'driver', bind: {label: '$.later'}},
      {tag: 'u2-text-input', name: 'later', props: {value: 'world'}},
    ]},
  }, new SpecContext(), reg);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const label = () => instance.node('driver').root.querySelector('.u2-input-label').textContent;
  assert.equal(label(), 'world', 'linked after the render pass');
  const editor = instance.node('later').root.querySelector('input');
  editor.value = 'typed';
  fire(editor, 'input');
  assert.equal(label(), 'typed', 'and follows the source live afterwards');
  instance.dispose();
});

spec('bindProps advertises the bound chrome — real type, never writable; bindStep answers its signal', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {l: 'Name', en: true}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-text-input', name: 'driver', bind: {label: '$.l', enabled: '$.en'}},
  }, ctx, reg);
  const input = instance.node('driver');
  const props = input.bindProps();
  const names = props.map((p) => p.name);
  assert.ok(names.includes('label'), `label advertised: ${names}`);
  assert.ok(names.includes('value'), 'the value signal is still there');
  const label = props.find((p) => p.name === 'label');
  assert.equal(label.type, 'string', 'the type comes from the signal\'s current value');
  assert.equal(label.writable, false, 'chrome binds are strictly one-way (UB-10)');
  const enabled = props.find((p) => p.name === 'enabled');
  assert.equal(enabled.type, 'bool');
  assert.equal(enabled.writable, false);
  assert.equal(input.bindStep('label'), ctx.data.l, 'the bound signal itself answers the step');
  instance.dispose();
});

spec('disposal: a source write after dispose is inert', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {l: 'Name'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-text-input', name: 'driver', bind: {label: '$.l'}},
  }, ctx, reg);
  const labelEl = instance.node('driver').root.querySelector('.u2-input-label');
  instance.dispose();
  ctx.data.l.value = 'After';
  assert.equal(labelEl.textContent, 'Name', 'the chrome effect died with the input');
});

spec('literal chrome installs no effect: the plain setter stays the sole owner', async () => {
  const input = new TextInput({label: 'Lit', enabled: false, postfix: 'mg', tooltipText: 'tip'});
  assert.equal(input.enabled, false, 'the literal applied');
  assert.equal(input.bindStep('enabled'), null, 'no live bind registered for a literal');
  assert.equal(input.bindStep('label'), null);
  assert.equal(input.bindStep('postfix'), null);
  assert.equal(input.bindStep('tooltipText'), null);

  input.enabled = true;
  input.label = 'Renamed';
  await flush();
  assert.equal(input.enabled, true, 'the imperative write sticks — nothing fights it');
  assert.equal(input.root.querySelector('input').disabled, false);
  assert.equal(input.label, 'Renamed');
  input.dispose();
});

spec('literal postfix \'\' never materializes the span or the rail', () => {
  const input = new TextInput({postfix: ''});
  assert.equal(input.root.querySelector('.u2-input-postfix'), null);
  assert.equal(input.root.querySelector('.u2-input-options'), null, 'the rail never materialized');
  input.dispose();
});

spec('a signal postfix materializes on the first non-empty value only', () => {
  const px = signal('');
  const input = new TextInput({postfix: px});
  assert.equal(input.root.querySelector('.u2-input-postfix'), null);
  px.value = 'mg';
  assert.equal(input.root.querySelector('.u2-input-postfix').textContent, 'mg');
  px.value = '';
  assert.equal(input.root.querySelector('.u2-input-postfix').textContent, '', 'the span stays, emptied');
  input.dispose();
});

spec('literal enabled: true runs no sweep — createEditor-disabled internals survive', () => {
  class Probe extends Input {
    createEditor() {
      const el = document.createElement('div');
      this.btn = document.createElement('button');
      this.btn.disabled = true;
      el.append(this.btn);
      return el;
    }
  }
  const input = new Probe({enabled: true}, null);
  assert.equal(input.btn.disabled, true, 'the internal control stays as the editor built it');
  assert.equal(input.enabled, true);
  input.dispose();
});

spec('literal chrome advertises nothing: bindProps carries no chrome entries', () => {
  const input = new TextInput({label: 'Lit', enabled: false, postfix: 'mg', tooltipText: 'tip'});
  const names = input.bindProps().map((p) => p.name);
  for (const name of [...CHROME, 'placeholder'])
    assert.equal(names.includes(name), false, `${name} is not advertised for a literal`);
  input.dispose();
});

spec('a spec-rendered input with bound label: dump keeps the bind, name never becomes a Signal', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {l: 'Caption'}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-form', name: 'f', children: [
      {tag: 'u2-text-input', name: 'driver', bind: {label: '$.l'}},
    ]},
  }, ctx, reg);
  const input = instance.node('driver');
  assert.equal(input.name instanceof Signal, false, 'a Signal never leaks into the form key');
  assert.ok(input.name === undefined || typeof input.name === 'string');
  const dumped = JSON.stringify(instance.dump());
  assert.ok(dumped.includes('"label":"$.l"'), `dump keeps the bind path: ${dumped}`);
  assert.equal(dumped.includes('Caption'), false, 'dump never carries the resolved value');
  instance.dispose();
});

spec('a signal option outside the spec engine works the same, one-way', () => {
  const live = signal('Caption');
  const input = new TextInput({label: live});
  assert.equal(input.label, 'Caption');
  live.value = 'Changed';
  assert.equal(input.label, 'Changed');
  assert.equal(input.bindStep('label'), live, 'answers the signal for deferred links');
  assert.equal(input.name, undefined, 'a Signal never leaks into the form key');
  input.dispose();
  live.value = 'Dead';
  assert.equal(input.root.querySelector('.u2-input-label').textContent, 'Changed',
    'the effect died with the input');
});
