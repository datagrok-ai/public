/* subscribe(): a platform or rxjs stream subscription owned by a scope or a control — unsubscribed with its
   owner, its handler batched and contained like every u2 handler, and visible to the leak audit as a
   `subscription` scope. */

import {test as nodeTest} from 'node:test';
import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {debug} from 'datagrok-api/u2core';
import {Scope, Control, signal, subscribe} from '../src/index.js';
import {Stream} from './platform-doubles.mjs';

function owned(name, body, run = test) {
  run(name, async () => {
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

/** A test of what happens with no owner: no test scope around its body. */
function unowned(name, body) {
  owned(name, body, nodeTest);
}

owned('a subscription made while a control builds ends with the control', () => {
  const stream = new Stream();
  const seen = [];
  const control = Control.build(() => {
    subscribe(stream, (x) => seen.push(x));
    return [];
  });
  stream.fire(1);
  assert.deepEqual(seen, [1]);
  assert.equal(stream.count, 1);
  control.dispose();
  assert.equal(stream.count, 0, 'unsubscribed with its owner');
  stream.fire(2);
  assert.deepEqual(seen, [1]);
});

owned('an explicit owner — a scope or a component — is what it ends with, after an await too', async () => {
  const stream = new Stream();
  const scope = new Scope();
  const control = new Control();
  await Promise.resolve();
  assert.equal(Scope.ambient, undefined, 'nothing is building after an await');
  scope.subscribe(stream, () => {});
  control.subscribe(stream, () => {});
  assert.equal(stream.count, 2);
  scope.dispose();
  assert.equal(stream.count, 1);
  control.dispose();
  assert.equal(stream.count, 0);
});

unowned('without an owner it throws, and under a disposed one it throws and never subscribes', () => {
  const stream = new Stream();
  assert.throws(() => subscribe(stream, () => {}), /subscribe\(\) needs an owner/);
  const gone = new Scope();
  gone.dispose();
  assert.throws(() => gone.subscribe(stream, () => {}), {name: 'U2UsageError', code: 'owner-disposed'});
  assert.equal(stream.count, 0);
});

owned('the answer unsubscribes and leaves the owner alive', () => {
  const stream = new Stream();
  const scope = new Scope();
  const off = scope.subscribe(stream, () => {});
  off();
  assert.equal(stream.count, 0);
  assert.equal(scope.isDisposed, false);
  scope.dispose();
});

owned('a throwing handler is reported to the owner, not thrown at the stream', () => {
  const stream = new Stream();
  const scope = new Scope();
  const errors = [];
  scope.onError = (e) => errors.push(e.message);
  scope.subscribe(stream, () => {
    throw new Error('boom');
  });
  assert.doesNotThrow(() => stream.fire(1));
  assert.deepEqual(errors, ['boom']);
  scope.dispose();
});

owned('the handler runs batched: two writes reach an effect once', () => {
  const stream = new Stream();
  const a = signal(0);
  const b = signal(0);
  const scope = new Scope();
  let runs = 0;
  scope.effect(() => {
    a.value;
    b.value;
    runs++;
  });
  scope.subscribe(stream, (x) => {
    a.value = x;
    b.value = x;
  });
  stream.fire(5);
  assert.equal(runs, 2, 'the initial run and one for the batch');
  scope.dispose();
});

owned('with u2 debug on, a live subscription is a `subscription` scope with its stack', () => {
  const was = debug.enabled;
  debug.enable();
  try {
    const stream = new Stream();
    const scope = new Scope();
    const before = debug.snapshot();
    scope.subscribe(stream, () => {});
    const leaked = debug.diff(before).scopes;
    assert.equal(leaked.length, 1);
    assert.equal(leaked[0].label, 'subscription');
    assert.ok(leaked[0].stack.includes('subscribe.test.js'), 'the stack names who subscribed');
    scope.dispose();
    assert.equal(debug.diff(before).total, 0);
  } finally {
    if (!was)
      debug.disable();
  }
});
