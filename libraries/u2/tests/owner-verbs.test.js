import {test as nodeTest} from 'node:test';
import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {debug} from '../src/core/scope.js';
import {Scope, Control, Dialog, Form, RecordsInput, TextInput, Signal, signal, computed, batch, lens, afterFlush,
  component, div, span} from '../src/index.js';
import {isWritableSignal} from '../src/core/signals.js';
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

const later = () => new Promise((resolve) => setTimeout(resolve, 1));

unowned('build: runs now with the owner ambient, answers what it returns, and a throw propagates', () => {
  const owner = Scope.root('page');
  const input = owner.build(() => {
    assert.equal(Scope.ambient, owner);
    return new TextInput({label: 'Name'});
  });
  assert.equal(Scope.ambient, undefined);
  assert.throws(() => owner.build(() => {
    throw new Error('broken');
  }), /broken/);
  assert.equal(Scope.ambient, undefined, 'the owner left the stack');
  owner.dispose();
  assert.equal(input.scope.isDisposed, true, 'what build built goes with the owner');
});

owned('every owner a consumer holds has the verbs: dialog.build, form.task, control.subscribe go with it', async () => {
  const dialog = new Dialog('Edit');
  const name = dialog.build(() => new TextInput({label: 'Name'}));
  const form = new Form();
  let late;
  await form.task(async (t) => {
    await later();
    late = t.build(() => new TextInput({label: 'Late'}));
  }).done;
  const stream = new Stream();
  const card = component(() => []);
  card.subscribe(stream, () => {});
  assert.equal(stream.count, 1);
  dialog.dispose();
  form.dispose();
  card.dispose();
  assert.equal(name.scope.isDisposed, true);
  assert.equal(late.scope.isDisposed, true);
  assert.equal(stream.count, 0);
});

const OLD_ENGINE = !Object.hasOwn(Control.prototype, 'dispose') &&
  'the js-api build predates Control.dispose removing its root';

test('dispose() by hand takes a control off the page; its owner going leaves the DOM to whoever holds it',
  {skip: OLD_ENGINE}, async () => {
  const live = Scope.liveCount;
  try {
    const host = div([]);
    document.body.append(host);
    const owner = Scope.root('page');
    const card = owner.build(() => component(() => [span('x')]));
    const name = owner.build(() => new TextInput({label: 'Name'}));
    host.append(card.root, name.root);
    name.dispose();
    assert.equal(name.root.isConnected, false, 'disposed by hand: off the page');
    assert.equal(card.root.parentNode, host);
    owner.dispose();
    assert.equal(card.scope.isDisposed, true);
    assert.equal(card.root.parentNode, host, 'a region swapping its content or a spec re-rendering removes it');
  } finally {
    resetDom();
    await flush();
  }
  assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
});

test('the engine\'s own types declare no owner verb on Component: a platform class with its own build() compiles', () => {
  const dts = readFileSync(new URL('component.d.ts', import.meta.resolve('datagrok-api/u2core')), 'utf8');
  const declared = dts.slice(dts.indexOf('export declare class Component'), dts.indexOf('export declare class Control'));
  assert.ok(declared.length > 0);
  assert.doesNotMatch(declared, /^\s+(build|task|subscribe)</m);
  assert.match(readFileSync(new URL('../src/core/component.ts', import.meta.url), 'utf8'),
    /interface Component extends OwnerVerbs/, 'u2 declares them for its own consumers');
});

owned('task: after an await, t.build is owned by the owner; the run settles done', async () => {
  const owner = Scope.root('page');
  let built;
  const run = owner.task(async (t) => {
    assert.equal(t.owner, owner);
    await later();
    assert.equal(Scope.ambient, undefined, 'nothing is building after an await');
    built = t.build(() => new TextInput({label: 'Loaded'}));
  });
  await run.done;
  assert.ok(built instanceof TextInput);
  owner.dispose();
  assert.equal(built.scope.isDisposed, true);
});

owned('task: disposing the owner aborts the signal and the build after the await is skipped', async () => {
  const owner = Scope.root('page');
  let signal;
  let after = 'not reached';
  const run = owner.task(async (t) => {
    signal = t.signal;
    await later();
    after = t.build(() => new TextInput({label: 'Late'}));
  });
  owner.dispose();
  assert.equal(signal.aborted, true);
  await run.done;
  assert.equal(after, undefined, 'skipped once the owner is gone');
});

owned('task: cancel aborts one run and skips its build; the owner lives on', async () => {
  const owner = Scope.root('page');
  let built = 'not reached';
  const run = owner.task(async (t) => {
    await later();
    built = t.build(() => 'built');
  });
  run.cancel();
  await run.done;
  assert.equal(built, undefined);
  assert.equal(owner.isDisposed, false);
  owner.dispose();
});

owned('task: serial — three runs asked for while one runs make two, the first and the last', async () => {
  const owner = Scope.root('page');
  const started = [];
  let release;
  const gate = new Promise((resolve) => release = resolve);
  const save = (n) => owner.task(async () => {
    started.push(n);
    if (n === 1)
      await gate;
  }, {serial: 'save'});
  const runs = [save(1), save(2), save(3)];
  assert.deepEqual(started, [1]);
  release();
  await Promise.all(runs.map((r) => r.done));
  assert.deepEqual(started, [1, 3], 'the second was replaced by the third');
  const other = owner.task(async () => started.push('other'), {serial: 'load'});
  await other.done;
  assert.deepEqual(started, [1, 3, 'other'], 'another queue does not wait');
  owner.dispose();
});

owned('task: a throw or a rejection goes to the owner\'s onError, and done still settles', async () => {
  const owner = Scope.root('page');
  const seen = [];
  owner.onError = (e) => seen.push(e.message);
  await owner.task(() => {
    throw new Error('sync');
  }).done;
  await owner.task(async () => {
    throw new Error('async');
  }).done;
  assert.deepEqual(seen, ['sync', 'async']);
  owner.dispose();
});

owned('subscribe: ends with its owner, and shows in DG.U2.debug as a subscription', () => {
  const was = debug.enabled;
  debug.enable({orphanSeconds: 0});
  try {
    const stream = new Stream();
    const owner = Scope.root('page');
    const before = debug.snapshot();
    const seen = [];
    owner.subscribe(stream, (x) => seen.push(x));
    assert.deepEqual(debug.diff(before).scopes.map((s) => s.label), ['subscription']);
    stream.fire(1);
    owner.dispose();
    stream.fire(2);
    assert.deepEqual(seen, [1]);
    assert.equal(stream.count, 0);
    assert.equal(debug.diff(before).total, 0);
  } finally {
    if (!was)
      debug.disable();
  }
});

owned('Scope.root: counted at once and named in leak reports', () => {
  const was = debug.enabled;
  debug.enable({orphanSeconds: 0});
  try {
    const before = debug.snapshot();
    const live = Scope.liveCount;
    const root = Scope.root('campaign session');
    assert.equal(Scope.liveCount, live + 1);
    assert.deepEqual(debug.diff(before).scopes.map((s) => s.label), ['campaign session']);
    root.dispose();
  } finally {
    if (!was)
      debug.disable();
  }
});

owned('lens: reads and writes a part of its source, and is a writable signal', () => {
  const person = signal({name: 'Ann', age: 41});
  const name = lens(person, (p) => p.name, (p, value) => ({...p, name: value}));
  assert.ok(name instanceof Signal);
  assert.equal(isWritableSignal(name), true);
  assert.equal(name.value, 'Ann');
  name.value = 'Bo';
  assert.deepEqual(person.peek(), {name: 'Bo', age: 41});
  person.value = {name: 'Cy', age: 7};
  assert.equal(name.peek(), 'Cy');
  const owner = component(() => [span(name)]);
  assert.equal(owner.root.textContent, 'Cy');
  name.value = 'Di';
  assert.equal(owner.root.textContent, 'Di', 'it binds like any signal');
  owner.dispose();
});

owned('lens: writing what it already reads writes nothing, so a two-way mirror settles', () => {
  const source = signal({n: 1});
  const n = lens(source, (s) => s.n, (s, value) => ({...s, n: value}));
  const mirror = signal(0);
  let writes = 0;
  const owner = Scope.root('page');
  owner.effect(() => {
    source.value;
    writes++;
  });
  owner.effect(() => mirror.value = n.value);
  owner.effect(() => n.value = mirror.value);
  const before = source.peek();
  n.value = 1;
  assert.equal(source.peek(), before, 'the same value leaves the source as it is');
  mirror.value = 5;
  assert.equal(source.peek().n, 5);
  assert.equal(writes, 2, 'one run at start, one for the change');
  owner.dispose();
});

owned('afterFlush: inside a batch it waits for the effects; outside a batch it runs now', () => {
  const items = signal(['a']);
  const owner = Scope.root('page');
  const shown = [];
  owner.effect(() => shown.push(items.value.length));
  const seen = [];
  batch(() => {
    items.value = ['a', 'b'];
    afterFlush(() => seen.push(shown.at(-1)));
    seen.push('in batch');
  });
  afterFlush(() => seen.push('now'));
  assert.deepEqual(seen, ['in batch', 2, 'now']);
  owner.dispose();
});

owned('RecordsInput.add from a handler focuses the new row\'s first cell before the handler returns', () => {
  const input = new RecordsInput({label: 'Fields',
    columns: [{name: 'name', header: 'Name', input: (value) => new TextInput({inline: true, value: value ?? ''})}],
    newItem: () => ({name: ''}), addText: 'Add field', value: [{name: 'Chemist'}]});
  document.body.append(input.root);
  fire(input.root.querySelector('[data-u2-part="add"]'), 'click');
  const rows = [...input.root.querySelector('.u2-records-body').children];
  assert.equal(rows.length, 2);
  assert.ok(document.activeElement === rows[1].querySelector('input'), 'no microtask needed');
  input.dispose();
});

owned('component: a handle makes the control carry its API, getters live, and name stamps the id', () => {
  const count = signal(0);
  const counter = component((self) => {
    const text = span(computed(() => String(count.value)));
    self.effect(() => text.dataset.count = String(count.value));
    return {root: text, get doubled() {
      return count.value * 2;
    }, increment: () => count.value++};
  }, {name: 'counter'});
  assert.equal(counter.root.dataset.u2Name, 'counter');
  assert.equal(counter.root.tagName, 'SPAN');
  counter.increment();
  counter.increment();
  assert.equal(counter.doubled, 4, 'a getter reads the live value');
  assert.equal(counter.root.textContent, '2');
  assert.equal(counter.root.dataset.count, '2');
  counter.dispose();
  count.value = 9;
  assert.equal(counter.root.dataset.count, '2', 'its effects went with it');
});

owned('component: self is the control\'s scope; a control or a text root goes into a div', () => {
  let self;
  const inner = component(() => [new TextInput({label: 'Inner'})]);
  const wrapped = component((s) => {
    self = s;
    return inner;
  });
  assert.equal(self, wrapped.scope);
  assert.equal(wrapped.root.tagName, 'DIV');
  assert.equal(wrapped.root.firstChild, inner.root);
  assert.equal(Control.forElement(inner.root), inner, 'the inner control keeps its root');
  const text = component(() => 'Nothing yet');
  assert.equal(text.root.textContent, 'Nothing yet');
  text.dispose();
  wrapped.dispose();
  inner.dispose();
});

owned('component: a static build costs no live scope until it owns something; a throwing builder leaks none', () => {
  const live = Scope.liveCount;
  const plain = component(() => div(['static']));
  assert.equal(Scope.liveCount, live);
  plain.own(() => {});
  assert.equal(Scope.liveCount, live + 1);
  plain.dispose();
  const alias = Control.build(() => div(['deprecated alias']), {name: 'alias'});
  assert.equal(alias.root.dataset.u2Name, 'alias', 'Control.build, which published bundles call, is the same');
  alias.dispose();
  assert.throws(() => component(() => {
    new TextInput({label: 'Half'});
    throw new Error('builder broke');
  }), /builder broke/);
  assert.equal(Scope.liveCount, live);
});

owned('getState: a control disabled with a reason says it; enabled or reasonless, no disabledReason', () => {
  const owner = Scope.root('page');
  const input = owner.build(() => new TextInput({label: 'Name'}));
  assert.equal(input.getState().disabledReason, undefined);
  input.disabled = 'Pick a campaign first';
  assert.equal(input.getState().enabled, false);
  assert.equal(input.getState().disabledReason, 'Pick a campaign first');
  input.root.removeAttribute('aria-description');
  assert.equal(input.getState().disabledReason, undefined, 'off with no reason (one taken off by hand)');
  input.disabled = null;
  input.root.setAttribute('aria-description', 'Pick a campaign first');
  assert.equal(input.getState().disabledReason, undefined, 'a description is no reason while enabled');
  owner.dispose();
});
