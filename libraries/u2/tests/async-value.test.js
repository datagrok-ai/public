/* AsyncValue: the single-value async engine sources run on — auto/manual kick, deps-driven
   re-runs cancelling the in-flight call, stale-response drop, error/refresh recovery, debounce
   coalescing and disposal, all over hand-resolved promises (the AsyncSource suite's technique). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, signal} from '../src/index.js';
import {AsyncValue, asyncValue} from '../src/core/async-value.js';

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
  const fetch = (abort) => {
    const call = {abort};
    const promise = new Promise((resolve, reject) => Object.assign(call, {resolve, reject}));
    calls.push(call);
    return promise;
  };
  return {calls, fetch};
}

smoke('auto: runs at construction, projects ready into value', async () => {
  const {calls, fetch} = tracker();
  const av = new AsyncValue(fetch, {debounceMs: 0});
  assert.equal(av.state.value.kind, 'loading');
  assert.equal(av.value.value, undefined);
  assert.equal(calls.length, 0, 'the fetch waits for the debounce');
  await flush();
  assert.equal(calls.length, 1);
  calls[0].resolve(42);
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 42});
  assert.equal(av.value.value, 42);
  av.dispose();
});

smoke('manual: idle until refresh, then the same lifecycle', async () => {
  const {calls, fetch} = tracker();
  const av = asyncValue(fetch, {auto: false, debounceMs: 0});
  await flush();
  assert.equal(av.state.value.kind, 'idle');
  assert.equal(calls.length, 0);
  av.refresh();
  assert.equal(av.state.value.kind, 'loading');
  await flush();
  assert.equal(calls.length, 1);
  calls[0].resolve('done');
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 'done'});
  av.dispose();
});

smoke('deps: a dep change re-runs and aborts the in-flight call', async () => {
  const {calls, fetch} = tracker();
  const dep = signal(1);
  const av = new AsyncValue(fetch, {deps: () => dep.value, debounceMs: 0});
  await flush();
  assert.equal(calls.length, 1);

  dep.value = 2;
  assert.equal(av.state.value.kind, 'loading');
  await flush();
  assert.equal(calls.length, 2);
  assert.equal(calls[0].abort.aborted, true, 'the superseded call is aborted');

  calls[0].resolve('stale');
  await flush();
  assert.equal(av.state.value.kind, 'loading', 'the stale response is dropped');
  calls[1].resolve('fresh');
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 'fresh'});
  av.dispose();
});

smoke('stale drop: a slow first response never overwrites the fast second', async () => {
  const {calls, fetch} = tracker();
  const av = new AsyncValue(fetch, {debounceMs: 0});
  await flush();
  av.refresh();
  await flush();
  assert.equal(calls.length, 2);
  calls[1].resolve('second');
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 'second'});
  calls[0].resolve('first');
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 'second'}, 'the late first run lost');
  av.dispose();
});

smoke('error: reported with the message, refresh recovers', async () => {
  const {calls, fetch} = tracker();
  const av = new AsyncValue(fetch, {debounceMs: 0});
  await flush();
  calls[0].reject(new Error('boom'));
  await flush();
  assert.deepEqual(av.state.value, {kind: 'error', message: 'boom'});
  assert.equal(av.value.value, undefined);

  av.refresh();
  assert.equal(av.state.value.kind, 'loading');
  await flush();
  calls[1].resolve(7);
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 7});
  av.dispose();
});

smoke('dispose: mid-flight aborts and drops the response, mid-debounce never fetches', async () => {
  const {calls, fetch} = tracker();
  const av = new AsyncValue(fetch, {debounceMs: 0});
  await flush();
  assert.equal(calls.length, 1);
  av.dispose();
  assert.equal(calls[0].abort.aborted, true);
  calls[0].resolve('too late');
  await flush();
  assert.equal(av.state.value.kind, 'loading', 'a disposed value never settles');

  const second = tracker();
  const av2 = new AsyncValue(second.fetch, {debounceMs: 0});
  av2.dispose();
  await flush();
  assert.equal(second.calls.length, 0, 'disposal during the debounce cancels the timer');
  av2.refresh();
  await flush();
  assert.equal(second.calls.length, 0, 'refresh after dispose is inert');
});

smoke('debounce: rapid dep changes coalesce into one run', async () => {
  const {calls, fetch} = tracker();
  const dep = signal(0);
  const av = new AsyncValue(fetch, {deps: () => dep.value, debounceMs: 20});
  dep.value = 1;
  dep.value = 2;
  dep.value = 3;
  assert.equal(av.state.value.kind, 'loading');
  assert.equal(calls.length, 0);
  await new Promise((resolve) => setTimeout(resolve, 50));
  assert.equal(calls.length, 1, 'one run for the whole burst');
  calls[0].resolve('coalesced');
  await flush();
  assert.deepEqual(av.state.value, {kind: 'ready', value: 'coalesced'});
  av.dispose();
});

smoke('ambient adoption: an AsyncValue built under a scope dies with it', async () => {
  const {calls, fetch} = tracker();
  const owner = new Scope();
  let av;
  Scope.runWith(owner, () => {
    av = new AsyncValue(fetch, {auto: false, debounceMs: 0});
  });
  owner.dispose();
  av.refresh();
  await flush();
  assert.equal(calls.length, 0, 'the adopted value was disposed with its owner');
  assert.equal(av.state.value.kind, 'idle');
});
