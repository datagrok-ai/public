/* u2 on the u2core engine: one signal engine per page (P1), element handlers owned, batched and
   contained (R17, P83), empty text for empty values (P67), and content built once inside an effect
   run kept alive by an explicit owner (R5 — a lazy tab). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal, span, link, button, bindValue, TabStrip, TextInput} from '../src/index.js';

function smoke(name, body) {
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

test('a second signal engine on the page fails at import, naming the fix', async () => {
  const saved = globalThis.DG;
  try {
    globalThis.DG = {U2: {Signal: class {}}};
    await assert.rejects(import('../src/core/signals.js?second-engine'),
      /a second copy of the signal engine is bundled — build with @datagrok\/build-config/);
    globalThis.DG = {U2: await import('datagrok-api/u2core')};
    await import('../src/core/signals.js?one-engine');
  } finally {
    globalThis.DG = saved;
  }
});

smoke('link() and button() hand their handler the DOM event', async () => {
  const scope = new Scope();
  const got = [];
  const [a, b] = Scope.runWith(scope, () => [
    link('Open', (e) => got.push(`link ${e.type}`)),
    button('Run', (e) => got.push(`button ${e.type}`)),
  ]);
  fire(a, 'click');
  fire(a, 'keydown', {key: 'Enter'});
  fire(a, 'keydown', {key: 'a'});
  fire(b, 'click');
  assert.deepEqual(got, ['link click', 'link keydown', 'button click']);
  scope.dispose();
  fire(b, 'click');
  assert.equal(got.length, 3, 'the owner disposed, the handler is gone');
});

smoke('an element handler is batched and a throw is reported to its owner, not thrown', async () => {
  const scope = new Scope();
  const errors = [];
  scope.onError = (e) => errors.push(e.message);
  const a = signal(0);
  const b = signal(0);
  const seen = [];
  scope.effect(() => seen.push(`${a.value}/${b.value}`));
  const el = Scope.runWith(scope, () => button('Both', () => {
    a.value = 1;
    b.value = 1;
  }));
  const bad = Scope.runWith(scope, () => button('Bad', () => {
    throw new Error('boom');
  }));
  fire(el, 'click');
  assert.deepEqual(seen, ['0/0', '1/1'], 'one effect run for two writes');
  fire(bad, 'click');
  assert.deepEqual(errors, ['boom']);
  scope.dispose();
});

smoke('a handler on an element built with no owner still runs, and its throw is reported', async () => {
  const saved = Scope.onUnhandledError;
  const errors = [];
  Scope.onUnhandledError = (e) => errors.push(e.message);
  try {
    let clicks = 0;
    const ok = button('Ok', () => clicks++);
    const bad = button('Bad', () => {
      throw new Error('unowned');
    });
    fire(ok, 'click');
    fire(bad, 'click');
    assert.equal(clicks, 1);
    assert.deepEqual(errors, ['unowned']);
  } finally {
    Scope.onUnhandledError = saved;
  }
});

smoke('span(signal) shows null and undefined as empty', async () => {
  const scope = new Scope();
  const value = signal(null);
  const el = Scope.runWith(scope, () => span(value));
  assert.equal(el.textContent, '');
  value.value = undefined;
  assert.equal(el.textContent, '');
  value.value = 0;
  assert.equal(el.textContent, '0');
  scope.dispose();
});

smoke('bindValue: a commit writes the signal and is not echoed back; an outside write is shown', async () => {
  const scope = new Scope();
  const input = document.createElement('input');
  let text = '';
  let fieldWrites = 0;
  Object.defineProperty(input, 'value', {get: () => text, set: (x) => {
    fieldWrites++;
    text = x;
  }});
  const type = (x) => {
    text = x;
    fire(input, 'input');
  };
  const source = signal('a');
  bindValue(scope, input, source);
  assert.equal(text, 'a');
  fieldWrites = 0;
  type('ab ');
  assert.equal(source.value, 'ab ');
  assert.equal(fieldWrites, 0, 'the commit is not written back into the field');
  source.value = 'xyz';
  assert.equal(text, 'xyz', 'an outside write reaches the field');
  type('xyz');
  source.value = 'next';
  assert.equal(text, 'next', 'a commit that changed nothing leaves no echo behind');
  scope.dispose();
});

smoke('a lazy tab built once survives switching away and back, and stays live', async () => {
  const count = signal(1);
  let builds = 0;
  let input;
  const tabs = new TabStrip({tabs: [
    {id: 'a', label: 'A', content: () => {
      builds++;
      input = new TextInput({label: 'Name'});
      const el = document.createElement('div');
      el.append(input.root, span(count, 'count'));
      return el;
    }},
    {id: 'b', label: 'B', content: document.createElement('div')},
  ]});
  document.body.append(tabs.root);
  assert.equal(builds, 1);
  for (const id of ['b', 'a', 'b', 'a'])
    tabs.activeTab.value = id;
  assert.equal(builds, 1, 'built once');
  count.value = 2;
  assert.equal(tabs.root.querySelector('.count').textContent, '2',
    'the binding built on first show still follows its signal');
  assert.equal(input.scope.isDisposed, false, 'the control built on first show is alive');
  tabs.dispose();
  assert.equal(input.scope.isDisposed, true, 'and dies with the strip');
});
