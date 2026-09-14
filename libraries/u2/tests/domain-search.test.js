/* `domains.search` (WO 2-7) over the memory backend: the box writes the source's `search` after a
   pause, on Enter at once, Escape clears, and a search set in code shows in the box. `DG` and
   `grok` come from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';
import {backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {registerDomainComponents} = await import('../src/dg/domain/registrations.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);

async function issues() {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  await flush();
  return src;
}

function type(input, value) {
  input.value = value;
  fire(input, 'input');
}

scoped('typing searches after the pause; Enter at once; Escape clears; the total follows', async () => {
  const src = await issues();
  const search = domains.search(src, {debounceMs: 20});
  const input = search.root.querySelector('input');
  assert.equal(search.root.dataset.u2, 'domain-search');
  assert.equal(input.placeholder, 'Search issues…');
  type(input, 'asp');
  assert.equal(src.search.value, '', 'not yet');
  await wait(30);
  assert.equal(src.search.value, 'asp');
  await flush();
  assert.equal(src.total.value, 1, 'the memory backend searched the title');
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Aspirin']);
  type(input, 'ibu');
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(src.search.value, 'ibu', 'Enter does not wait');
  fire(input, 'keydown', {key: 'Escape'});
  assert.equal(input.value, '');
  assert.equal(src.search.value, '');
  await flush();
  assert.equal(src.total.value, 3);
  search.dispose();
  src.dispose();
});

scoped('two-way: a search set in code shows in the box and cancels a pending keystroke', async () => {
  const src = await issues();
  const search = domains.search(src, {debounceMs: 20});
  const input = search.root.querySelector('input');
  type(input, 'pending');
  src.search.value = 'nap';
  assert.equal(input.value, 'nap');
  await wait(30);
  assert.equal(src.search.value, 'nap', 'the keystroke that was pending never lands');
  search.dispose();
  src.dispose();
});

scoped('spec: u2-domain-search is registered with usage', () => {
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const meta = reg.get('u2-domain-search');
  assert.equal(meta.usage.length > 0, true);
  assert.deepEqual(meta.props.slice(0, 3).map((p) => p.name), ['source', 'placeholder', 'debounceMs']);
});

scoped('a search while the session is dirty asks: cancel puts the text back, discard searches', async () => {
  const src = await issues();
  const search = domains.search(src, {debounceMs: 20});
  const input = search.root.querySelector('input');
  src.rows.byKey('i1').title = 'Edited';
  assert.equal(src.isDirty.value, true);
  type(input, 'nap');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.notEqual(document.body.querySelector('.u2-dialog'), null, 'the gate asks');
  assert.equal(src.search.value, '', 'nothing written yet');
  fire(buttonNamed('CANCEL'), 'click');
  await flush();
  assert.equal(input.value, '', 'the box shows the search in force again');
  assert.equal(src.isDirty.value, true, 'the edit is kept');
  type(input, 'nap');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(src.search.value, 'nap');
  assert.equal(src.isDirty.value, false);
  await flush();
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Naproxen'], 're-queried');
  search.dispose();
  src.dispose();
});
