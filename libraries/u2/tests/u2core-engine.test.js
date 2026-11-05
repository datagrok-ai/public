/* The u2core engine (js-api `u2core/scope.ts`, `debug.ts`, `Component`): owner capture per effect run,
   error containment, disposal that survives a throwing cleanup, the scope helpers and `DG.U2.debug`. */

import './dom-shim.js';
import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {Scope, Control, Component, signal, debug} from 'datagrok-api/u2core';

function capturing(fn) {
  const saved = Scope.onUnhandledError;
  const errors = [];
  Scope.onUnhandledError = (e) => errors.push(e);
  try {
    fn(errors);
  }
  finally {
    Scope.onUnhandledError = saved;
  }
  return errors;
}

function tracking(fn) {
  debug.enable({orphanSeconds: 0});
  let result;
  try {
    result = fn();
  }
  finally {
    if (!(result instanceof Promise))
      debug.disable();
  }
  return result instanceof Promise ? result.finally(() => debug.disable()) : result;
}

test('a throwing disposer does not stop the next; the error reaches onUnhandledError', () => {
  const ran = [];
  const errors = capturing(() => {
    const scope = new Scope();
    scope.own(() => ran.push(1));
    scope.own(() => {
      throw new Error('boom');
    });
    scope.own(() => ran.push(3));
    scope.dispose();
  });
  assert.deepEqual(ran, [1, 3]);
  assert.equal(errors.length, 1);
  assert.equal(errors[0].message, 'boom');
});

test('an effect that throws on a write: the writer does not throw, onError gets it, it re-runs', () => {
  const owner = new Scope();
  const seen = [];
  owner.onError = (e) => seen.push(e.message);
  const s = signal(0);
  const runs = [];
  owner.effect(() => {
    runs.push(s.value);
    if (s.value === 1)
      throw new Error('bad ' + s.value);
  });
  s.value = 1;
  assert.deepEqual(seen, ['bad 1']);
  s.value = 2;
  assert.deepEqual(runs, [0, 1, 2]);
  owner.dispose();
});

test('a control built inside an effect run is disposed on the next run and on owner dispose', () => {
  const owner = new Scope();
  const s = signal(0);
  const built = [];
  owner.effect(() => {
    s.value;
    built.push(new Control());
  });
  assert.equal(built.length, 1);
  assert.equal(built[0].scope.isDisposed, false);
  s.value = 1;
  assert.equal(built.length, 2);
  assert.equal(built[0].scope.isDisposed, true);
  assert.equal(built[1].scope.isDisposed, false);
  owner.dispose();
  assert.equal(built[1].scope.isDisposed, true);
  s.value = 2;
  assert.equal(built.length, 2);
});

test('an effect run that owns nothing is not counted; one that does is, until the next run', () => {
  const base = Scope.liveCount;
  const owner = new Scope();
  const s = signal(0);
  owner.effect(() => s.value);
  assert.equal(Scope.liveCount, base + 1);
  const cleanups = [];
  owner.effect(() => {
    const v = s.value;
    return () => cleanups.push(v);
  });
  assert.equal(Scope.liveCount, base + 2);
  s.value = 1;
  assert.deepEqual(cleanups, [0]);
  assert.equal(Scope.liveCount, base + 2);
  owner.dispose();
  assert.equal(Scope.liveCount, base);
});

test('a child disposed early leaves its parent; no growth over 1,000 cycles', () => tracking(() => {
  const parent = new Scope();
  let disposed = 0;
  const first = parent.child('row');
  first.own(() => disposed++);
  first.dispose();
  parent.dispose();
  assert.equal(disposed, 1);

  const host = new Scope();
  const before = debug.liveScopes().length;
  for (let i = 0; i < 1000; i++) {
    const c = host.child('row');
    c.own(() => {});
    c.dispose();
  }
  assert.equal(debug.liveScopes().length, before);
  const child = host.child();
  const records = debug.liveScopes();
  assert.equal(records[records.length - 1].parent, records[records.length - 2].id);
  child.dispose();
  host.dispose();
}));

test('bind after dispose does nothing; bind batches two writes into one effect run', () => {
  const scope = new Scope();
  const a = signal(0);
  const b = signal(0);
  const runs = [];
  scope.effect(() => runs.push(a.value + b.value));
  const both = scope.bind(() => {
    a.value = 1;
    b.value = 2;
    return 'done';
  });
  assert.equal(both(), 'done');
  assert.deepEqual(runs, [0, 3]);
  let called = false;
  const late = scope.bind(() => called = true);
  scope.dispose();
  assert.equal(late(), undefined);
  assert.equal(called, false);
});

test('bind reports a throw and answers undefined', () => {
  const scope = new Scope();
  const seen = [];
  scope.onError = (e, s) => seen.push([e.message, s === scope]);
  assert.equal(scope.bind(() => {
    throw new Error('handler');
  })(), undefined);
  assert.deepEqual(seen, [['handler', true]]);
  scope.dispose();
});

const rowLike = () => new Proxy({}, {get: (_, name) => {
  throw new Error(`No such column: ${String(name)}`);
}});

test('Scope.thenable: a promise, or an object or function with a callable then; never a throw on reading then', () => {
  const promise = Promise.reject(new Error('unread'));
  promise.catch(() => {});
  const custom = {then: () => {}};
  const callable = Object.assign(() => {}, {then: () => {}});
  assert.equal(Scope.thenable(promise), promise);
  assert.equal(Scope.thenable(custom), custom);
  assert.equal(Scope.thenable(callable), callable);
  for (const x of [rowLike(), {then: 1}, {}, () => {}, null, undefined, 0, 'then'])
    assert.equal(Scope.thenable(x), undefined);
});

test('a handler answering a proxy that throws on then reports nothing; a rejecting thenable is reported', () => {
  const scope = new Control().scope;
  const seen = [];
  scope.onError = (e) => seen.push(e.message);
  const row = rowLike();
  assert.equal(scope.bind(() => row)(), row);
  assert.deepEqual(seen, []);
  let fail;
  const pending = {then: (_ok, reject) => fail = reject};
  assert.equal(scope.bind(() => pending)(), pending);
  fail(new Error('rejected'));
  assert.deepEqual(seen, ['rejected']);
  scope.dispose();
});

test('a static control with text only leaves Scope.liveCount unchanged, alone or under an owner', () => {
  class Caption extends Control {
    constructor(text) {
      super();
      this.root.dataset.u2 = 'caption';
      this.scope.run(() => this.root.textContent = text);
    }
  }
  const base = Scope.liveCount;
  let owned;
  const owner = Control.build(() => [owned = new Caption('Campaigns')]);
  const alone = new Caption('Hit Design');
  alone.name = 'title';
  assert.equal(alone.getState().kind, 'caption');
  assert.ok(alone.rootId.startsWith('u2-'));
  assert.equal(Scope.liveCount, base + 1, 'the owner only');
  owned.own(() => {});
  assert.equal(Scope.liveCount, base + 2, 'counted once it owns something');
  alone.dispose();
  owner.dispose();
  assert.equal(owned.scope.isDisposed, true, 'disposed with its owner');
  assert.equal(Scope.liveCount, base);
});

test('a control marked singleton before it owns anything keeps the mark, and owns the listener it adds', () =>
  tracking(() => {
    const snap = debug.snapshot();
    const control = new Control();
    debug.markSingleton(control);
    assert.equal(debug.diff(snap).scopes.length, 0, 'not tracked while it owns nothing');
    control.scope.listen(document, 'u2-late', () => {});
    const report = debug.diff(snap);
    assert.equal(report.scopes.length, 1);
    assert.equal(report.listeners.length, 1);
    assert.equal(report.listeners[0].owner, report.scopes[0].id);
    assert.equal(report.listeners[0].label, 'singleton');
    assert.equal(report.total, 0);
    control.dispose();
  }));

test('bind reports the rejection of an async handler, and still answers its promise', async () => {
  const scope = new Scope();
  const seen = [];
  scope.onError = (e) => seen.push(e.message);
  const answer = scope.bind(async () => {
    throw new Error('later');
  })();
  assert.ok(answer instanceof Promise);
  await assert.rejects(answer, /later/);
  assert.deepEqual(seen, ['later']);
  assert.equal(await scope.bind(async () => 'fine')(), 'fine');
  assert.deepEqual(seen, ['later']);
  scope.dispose();
});

test('a raw global listener the browser drops itself — once, or on its signal — is not reported', () => tracking(() => {
  const shim = globalThis.document;
  const doc = new EventTarget();
  globalThis.document = doc;
  const scope = new Scope();
  try {
    const handler = () => {};
    const abort = new AbortController();
    scope.build(() => {
      doc.addEventListener('u2-once', handler, {once: true});
      doc.addEventListener('u2-signal', handler, {signal: abort.signal});
      doc.addEventListener('u2-plain', handler);
    });
    assert.deepEqual(debug.listeners().map((l) => l.type).sort(), ['u2-plain', 'u2-signal']);
    abort.abort();
    assert.deepEqual(debug.listeners().map((l) => l.type), ['u2-plain']);
    doc.removeEventListener('u2-plain', handler);
    assert.equal(debug.listeners().length, 0);
  }
  finally {
    globalThis.document = shim;
    scope.dispose();
  }
}));

test('listen and setInterval stop on dispose and show in debug while alive', async () => tracking(async () => {
  const scope = new Scope();
  const hits = [];
  scope.listen(document, 'u2-ping', () => hits.push('ping'));
  let ticks = 0;
  scope.setInterval(() => ticks++, 1);
  let timeouts = 0;
  scope.setTimeout(() => timeouts++, 1);
  assert.equal(debug.listeners().length, 1);
  assert.equal(debug.listeners()[0].target, 'document');
  assert.deepEqual(debug.timers().map((t) => t.kind).sort(), ['interval', 'timeout']);
  document.dispatchEvent(new Event('u2-ping'));
  await new Promise((r) => setTimeout(r, 20));
  assert.ok(ticks > 0);
  assert.equal(timeouts, 1);
  assert.deepEqual(debug.timers().map((t) => t.kind), ['interval']);
  scope.dispose();
  const after = ticks;
  document.dispatchEvent(new Event('u2-ping'));
  await new Promise((r) => setTimeout(r, 20));
  assert.equal(ticks, after);
  assert.deepEqual(hits, ['ping']);
  assert.equal(debug.listeners().length, 0);
  assert.equal(debug.timers().length, 0);
}));

test('abortSignal aborts on dispose', () => {
  const scope = new Scope();
  const signal = scope.abortSignal;
  assert.equal(signal.aborted, false);
  scope.dispose();
  assert.equal(signal.aborted, true);
});

test('snapshot/diff: three scopes after the snapshot, each with a stack; back to 0 after dispose', () => tracking(() => {
  const snap = debug.snapshot();
  const scopes = [new Scope(), new Scope(), new Scope()];
  const report = debug.diff(snap);
  assert.equal(report.scopes.length, 3);
  assert.ok(report.scopes.every((r) => typeof r.stack === 'string' && r.stack.length > 0));
  for (const s of scopes)
    s.dispose();
  assert.equal(debug.diff(snap).total, 0);
}));

test('a control scope is described by its kind and name; singletons do not count', () => tracking(() => {
  const snap = debug.snapshot();
  const control = new Control();
  control.root.dataset.u2 = 'probe';
  control.name = 'p1';
  control.own(() => {});
  const [record] = debug.diff(snap).scopes;
  assert.equal(record.label, 'probe');
  assert.equal(record.name, 'p1');
  assert.equal(record.connected, false);
  const host = new Scope();
  host.listen(document, 'x', () => {});
  debug.markSingleton(host);
  const report = debug.diff(snap);
  assert.equal(report.listeners[0].label, 'singleton');
  assert.equal(report.total, 1);
  control.dispose();
  host.dispose();
}));

test('liveCount is exact with tracking off, and nothing is listed', () => {
  assert.equal(debug.enabled, false);
  const base = debug.liveCount();
  const scopes = Array.from({length: 5}, () => new Scope());
  assert.equal(debug.liveCount(), base + 5);
  assert.equal(Scope.liveCount, base + 5);
  assert.deepEqual(debug.liveScopes(), []);
  for (const s of scopes)
    s.dispose();
  assert.equal(debug.liveCount(), base);
});

test('the effect loop warning fires once at 25 re-runs', () => tracking(() => {
  const warn = console.warn;
  const warnings = [];
  console.warn = (...args) => warnings.push(args[0]);
  try {
    const scope = new Scope();
    const s = signal(0);
    scope.effect(() => {
      if (s.value < 40)
        s.value++;
    });
    assert.equal(s.peek(), 40);
    const loops = warnings.filter((w) => String(w).startsWith('u2: effect re-ran 25 times'));
    assert.equal(loops.length, 1);
    scope.dispose();
  }
  finally {
    console.warn = warn;
  }
}));

test('Component.on passes (arg, source) and stops once its owner is disposed', () => {
  class Emitter extends Component {
    emit(x) {
      this.fireEvent('ping', x);
    }
  }
  const emitter = new Emitter();
  const owner = new Scope();
  const got = [];
  owner.build(() => emitter.on('ping', (x, source) => got.push([x, source === emitter])));
  emitter.emit(1);
  owner.dispose();
  emitter.emit(2);
  assert.deepEqual(got, [[1, true]]);
  emitter.dispose();
});

test('fireEvent: a throwing listener is reported and the rest still run', () => {
  class Emitter extends Component {
    emit() {
      this.fireEvent('go');
    }
  }
  const emitter = new Emitter();
  const seen = [];
  emitter.scope.onError = (e) => seen.push(e.message);
  const ran = [];
  emitter.onEvent('go').subscribe(() => {
    throw new Error('first');
  });
  emitter.onEvent('go').subscribe(() => ran.push('second'));
  emitter.emit();
  assert.deepEqual(seen, ['first']);
  assert.deepEqual(ran, ['second']);
  emitter.dispose();
});

test('an adopted component disowns itself from its adopter on dispose', () => {
  const host = new Scope();
  const base = Scope.liveCount;
  for (let i = 0; i < 100; i++) {
    const c = host.build(() => new Control());
    c.dispose();
  }
  assert.equal(Scope.liveCount, base);
  assert.equal(host._disposers?.size ?? 0, 0);
  let calls = 0;
  const c = host.build(() => new Control());
  c.own(() => calls++);
  host.dispose();
  assert.equal(calls, 1);
});

test('runWith(undefined) builds with no owner', () => {
  const host = new Scope();
  const free = host.build(() => Scope.runWith(undefined, () => new Control()));
  host.dispose();
  assert.equal(free.scope.isDisposed, false);
  free.dispose();
});
