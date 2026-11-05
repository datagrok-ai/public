/* What a dg-ui/1 spec says about actions and spacing: a `u2-button` off with the reason it is off (the README's first
   rule, sayable in a spec), buttons inside a `u2-form` in the form's button row, and a container's `gap` on the token
   scale. Every prop is live. Every test leaves the live-scope count where it found it. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';

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

spec('u2-button disabled: the reason it is off, live; a press on it runs nothing', async () => {
  const ran = [];
  const ctx = new SpecContext({data: {why: 'Fill in Name'}, commands: {save: () => ran.push('save')}});
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'div', children: [
    {tag: 'u2-button', name: 'save', props: {text: 'Save', primary: true}, bind: {disabled: '$.why'},
      on: {click: 'cmd:save'}},
    {tag: 'u2-button', name: 'later', props: {text: 'Later', disabled: 'Comes in the next release'}},
  ]}}, ctx, registry());
  document.body.append(instance.root);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const save = instance.node('save').root;
  assert.equal(save.getAttribute('aria-disabled'), 'true');
  assert.equal(save.getAttribute('aria-description'), 'Fill in Name', 'the reason is said on focus');
  assert.equal(instance.node('save').disabled, 'Fill in Name');
  fire(save, 'click');
  assert.deepEqual(ran, []);
  ctx.data.why.value = null;
  await flush();
  assert.equal(save.getAttribute('aria-disabled'), null);
  assert.equal(save.getAttribute('aria-description'), null);
  fire(save, 'click');
  assert.deepEqual(ran, ['save']);
  assert.equal(instance.node('later').root.getAttribute('aria-description'), 'Comes in the next release');
  instance.dispose();
});

spec('u2-form-buttons in a u2-form is its button row, under the rows; align places it', () => {
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-form', name: 'form', children: [
    {tag: 'u2-text-input', name: 'name', props: {label: 'Name'}},
    {tag: 'u2-form-buttons', name: 'buttons', props: {align: 'editors'}, children: [
      {tag: 'u2-button', name: 'cancel', props: {text: 'Cancel'}},
      {tag: 'u2-button', name: 'save', props: {text: 'Save', primary: true}},
    ]},
    {tag: 'u2-number-input', name: 'amount', props: {label: 'Amount'}},
  ]}}, new SpecContext(), registry());
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const form = instance.node('form').root;
  const row = form.querySelector('.u2-form-buttons');
  assert.ok(row, 'the form has a button row');
  assert.ok(row.parentElement === form, 'under the rows, not one of them');
  assert.ok(row.children[0] === instance.node('buttons').root);
  assert.deepEqual([...instance.node('buttons').root.children].map((b) => b.textContent), ['Cancel', 'Save']);
  assert.equal(instance.node('buttons').root.classList.contains('u2-gap-m'), true, 'the buttons 8px apart');
  assert.equal(row.classList.contains('u2-form-buttons-editors'), true);
  assert.equal(instance.node('form').inputs.length, 2);
  instance.dispose();
});

spec('a u2 tag in a registry without u2\'s components says to call registerAll()', async () => {
  const warn = console.warn;
  console.warn = () => {};
  try {
    const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-text-input', props: {label: 'Name'}}},
      new SpecContext(), new Registry());
    assert.match(instance.root.textContent, /registerAll\(\)/);
    instance.dispose();
  } finally {
    console.warn = warn;
  }
});

spec('a container\'s gap is a token-scale class, live', async () => {
  const ctx = new SpecContext({data: {g: 'm'}});
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-div-h', name: 'row', bind: {gap: '$.g'},
    children: [{tag: 'u2-button', props: {text: 'A'}}, {tag: 'u2-button', props: {text: 'B'}}]}}, ctx, registry());
  const row = instance.node('row').root;
  assert.equal(row.classList.contains('u2-gap-m'), true);
  ctx.data.g.value = 'xl';
  await flush();
  assert.equal(row.classList.contains('u2-gap-m'), false);
  assert.equal(row.classList.contains('u2-gap-xl'), true);
  instance.dispose();
});
