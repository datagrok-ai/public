/* AsyncView: the state → content projection over a hand-resolved AsyncSource — loading, ready,
   empty, error with Retry — and the render callback running untracked inside the state effect,
   so a render that writes a signal (a builder normalizing its bound tree) neither throws nor
   re-triggers the view. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, signal} from '../src/index.js';
import {AsyncSource} from '../src/core/async-source.js';
import {AsyncView} from '../src/components/display/async-view.js';

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

function tracker() {
  const calls = [];
  const fetch = (query, abort) => {
    const call = {query, abort};
    const promise = new Promise((resolve, reject) => Object.assign(call, {resolve, reject}));
    calls.push(call);
    return promise;
  };
  return {calls, fetch};
}

/** Past the debounce `AsyncView.owned` keeps by default. */
const settle = () => new Promise((resolve) => setTimeout(resolve, 200));

function list(items) {
  const el = document.createElement('ul');
  for (const item of items) {
    const li = document.createElement('li');
    li.textContent = String(item);
    el.append(li);
  }
  return el;
}

smoke('states: loading → ready renders the items, empty and error show their rows', async () => {
  const {calls, fetch} = tracker();
  const source = new AsyncSource(fetch, {debounceMs: 0});
  const view = new AsyncView(source, list);
  assert.equal(view.root.children.length, 0, 'idle renders nothing');

  view.refresh();
  assert.equal(view.root.getAttribute('aria-busy'), 'true');
  assert.ok(view.root.querySelector('.u2-loader'));
  await flush();
  calls[0].resolve([1, 2]);
  await flush();
  assert.equal(view.root.getAttribute('aria-busy'), 'false');
  assert.equal(view.root.querySelectorAll('li').length, 2);

  view.refresh();
  await flush();
  calls[1].resolve([]);
  await flush();
  assert.equal(view.root.querySelector('.u2-async-empty')?.textContent, 'No data');

  view.refresh();
  await flush();
  calls[2].reject(new Error('boom'));
  await flush();
  assert.equal(view.root.querySelector('.u2-async-error-message')?.textContent, 'boom');
  view.root.querySelector('.u2-async-error button').click();
  await flush();
  assert.equal(calls.length, 4, 'Retry re-runs the fetch');
  view.dispose();
  source.dispose();
});

smoke('render runs untracked: writing a signal it read neither throws nor re-triggers the view', async () => {
  const {calls, fetch} = tracker();
  const source = new AsyncSource(fetch, {debounceMs: 0});
  const normalized = signal(0);
  let renders = 0;
  const view = new AsyncView(source, (items) => {
    renders++;
    normalized.value = normalized.value + items.length;
    return list(items);
  });
  view.refresh();
  await flush();
  calls[0].resolve(['a', 'b', 'c']);
  await flush();
  assert.equal(renders, 1);
  assert.equal(normalized.value, 3);

  normalized.value = 10;
  await flush();
  assert.equal(renders, 1, 'the signal the render wrote is not a dependency of the view');
  assert.equal(view.root.querySelectorAll('li').length, 3);

  view.refresh();
  await flush();
  calls[1].resolve(['d']);
  await flush();
  assert.equal(renders, 2, 'a state change still re-renders');
  assert.equal(normalized.value, 11);
  view.dispose();
  source.dispose();
});

smoke('dispose: the content scope goes with the view', async () => {
  const {calls, fetch} = tracker();
  let inner;
  const view = AsyncView.owned(fetch, (items) => {
    inner = Scope.ambient;
    return list(items);
  });
  view.refresh();
  await settle();
  calls[0].resolve([1]);
  await flush();
  assert.ok(inner && !inner.isDisposed, 'the render ran inside a live content scope');
  view.dispose();
  assert.equal(inner.isDisposed, true);
});
