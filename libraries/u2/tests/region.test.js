/* each / show (core/region.ts): a kept key keeps its element, its scope and the focus inside it; a
   removed key's scope is disposed; a reorder moves nodes without building any; a changed item is
   patched or rebuilt; a repeated key warns once; show swaps two builds, each in its own scope; both
   throw without an owner and leave no scope behind. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, signal, each, show} from '../src/index.js';

function ui(name, body) {
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

/** A list region over `items` (`{id, text}`), each row an input; counts the builds. */
function rows(items, options = {}) {
  const scope = new Scope();
  const host = document.createElement('div');
  document.body.append(host);
  const built = [];
  Scope.runWith(scope, () => each(host, items, (item) => item.id, (item) => {
    built.push(item.id);
    const row = document.createElement('div');
    const input = document.createElement('input');
    input.value = item.text;
    row.append(input);
    return row;
  }, options));
  const ids = () => host.children.map((c) => c.dataset.u2Item);
  return {scope, host, built, ids};
}

ui('a kept key keeps its element and the focus inside it while another item changes', () => {
  const a = {id: 'a', text: 'A'};
  const items = signal([a, {id: 'b', text: 'B'}]);
  const {scope, host, built, ids} = rows(items, {update: (el, item) => el.querySelector('input').value = item.text});
  assert.deepEqual(ids(), ['a', 'b'], 'data-u2-item is the key');
  const input = host.querySelector('input');
  input.focus();
  const first = host.children[0];
  items.value = [a, {id: 'b', text: 'B2'}, {id: 'c', text: 'C'}];
  assert.equal(host.children[0] === first, true);
  assert.equal(document.activeElement === input, true, 'focus stays where it was');
  assert.equal(host.children[1].querySelector('input').value, 'B2', 'the changed item was patched');
  assert.deepEqual(built, ['a', 'b', 'c'], 'only the new key was built');
  scope.dispose();
});

ui('a removed key\'s scope is disposed: what its build started stops', () => {
  const tick = signal(0);
  const runs = {a: 0, b: 0};
  const items = signal(['a', 'b']);
  const scope = new Scope();
  const host = document.createElement('div');
  Scope.runWith(scope, () => each(host, items, (id) => id, (id) => {
    Scope.ambient.effect(() => {
      tick.value;
      runs[id]++;
    });
    return document.createElement('span');
  }));
  items.value = ['b'];
  tick.value++;
  assert.deepEqual(runs, {a: 1, b: 2}, 'a stopped with its key, b runs on');
  assert.equal(host.children.length, 1);
  scope.dispose();
  tick.value++;
  assert.deepEqual(runs, {a: 1, b: 2}, 'the owner disposes the kept ones');
});

ui('a reorder moves nodes without building any; a node already in place is not moved', () => {
  const items = signal(['a', 'b', 'c'].map((id) => ({id, text: id})));
  const {scope, host, built, ids} = rows(items);
  const [a, b, c] = host.children;
  const moves = [];
  const insertBefore = host.insertBefore.bind(host);
  host.insertBefore = (node, ref) => {
    moves.push(node.dataset.u2Item);
    return insertBefore(node, ref);
  };
  items.value = [items.value[2], items.value[0], items.value[1]];
  assert.deepEqual(ids(), ['c', 'a', 'b']);
  assert.deepEqual(moves, ['c'], 'one node moved');
  assert.deepEqual([host.children[0] === c, host.children[1] === a, host.children[2] === b], [true, true, true]);
  assert.deepEqual(built, ['a', 'b', 'c']);
  scope.dispose();
});

ui('a changed item without update, or with update answering false, is rebuilt; a function source works too', () => {
  let list = [{id: 'a', text: 'A'}];
  const version = signal(0);
  const source = () => {
    version.value;
    return list;
  };
  const {scope, host, built} = rows(source, {update: (el, item) => item.text !== 'rebuild'});
  const first = host.children[0];
  list = [{id: 'a', text: 'patched'}];
  version.value++;
  assert.equal(host.children[0] === first, true, 'patched in place');
  list = [{id: 'a', text: 'rebuild'}];
  version.value++;
  assert.equal(host.children[0] === first, false);
  assert.equal(host.children[0].querySelector('input').value, 'rebuild');
  assert.deepEqual(built, ['a', 'a']);
  scope.dispose();
});

ui('an empty list shows the empty build; items replace it; an own data-u2-item is kept', () => {
  const items = signal([]);
  const scope = new Scope();
  const host = document.createElement('div');
  Scope.runWith(scope, () => each(host, items, (x) => x, (x) => {
    const el = document.createElement('span');
    if (x === 'own')
      el.dataset.u2Item = 'mine';
    return el;
  }, {empty: () => {
    const el = document.createElement('em');
    el.textContent = 'Nothing';
    return el;
  }}));
  assert.equal(host.textContent, 'Nothing');
  items.value = ['x', 'own'];
  assert.deepEqual(host.children.map((c) => [c.tagName, c.dataset.u2Item]), [['SPAN', 'x'], ['SPAN', 'mine']]);
  items.value = [];
  assert.deepEqual(host.children.map((c) => c.tagName), ['EM']);
  scope.dispose();
});

ui('a repeated key is warned about once and suffixed', () => {
  const warn = console.warn;
  const warnings = [];
  console.warn = (...args) => warnings.push(args.join(' '));
  try {
    const items = signal(['a', 'a', 'a']);
    const scope = new Scope();
    const host = document.createElement('div');
    Scope.runWith(scope, () => each(host, items, (x) => x, () => document.createElement('span')));
    items.value = ['a', 'a'];
    assert.deepEqual(host.children.map((c) => c.dataset.u2Item), ['a', 'a#2']);
    assert.deepEqual(warnings, ['u2: each(): duplicate key "a"']);
    scope.dispose();
  } finally {
    console.warn = warn;
  }
});

ui('show builds a branch while its condition holds, each in its own scope, disposed when left', () => {
  const on = signal(true);
  const tick = signal(0);
  let runs = 0;
  const scope = new Scope();
  const host = document.createElement('div');
  Scope.runWith(scope, () => show(host, on, () => {
    Scope.ambient.effect(() => {
      tick.value;
      runs++;
    });
    const el = document.createElement('b');
    el.textContent = 'yes';
    return el;
  }, () => {
    const el = document.createElement('i');
    el.textContent = 'no';
    return el;
  }));
  assert.equal(host.textContent, 'yes');
  on.value = false;
  assert.equal(host.textContent, 'no');
  tick.value++;
  assert.equal(runs, 1, 'the left branch stopped');
  on.value = true;
  assert.equal(host.textContent, 'yes');
  assert.equal(runs, 2, 'built again on entry');
  scope.dispose();
});

ui('show without an otherwise empties the host; the same verdict builds nothing new', () => {
  const on = signal(false);
  let builds = 0;
  const scope = new Scope();
  const host = document.createElement('div');
  Scope.runWith(scope, () => show(host, () => on.value, () => {
    builds++;
    return document.createElement('b');
  }));
  assert.equal(host.children.length, 0);
  on.value = true;
  on.value = true;
  assert.equal(builds, 1);
  on.value = false;
  assert.equal(host.children.length, 0);
  scope.dispose();
});

ui('no owner: each and show throw', () => {
  const host = document.createElement('div');
  assert.throws(() => each(host, () => [], (x) => x, () => document.createElement('span')), /each\(\) needs an owner/);
  assert.throws(() => show(host, () => true, () => document.createElement('span')), /show\(\) needs an owner/);
});
