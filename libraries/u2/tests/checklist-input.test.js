/* ChecklistInput: checkboxes the user may also put in order — the value is the checked items in the order
   shown, `order` every item's; orderable rows move by a handle or Alt+↑/↓. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, ChecklistInput, signal} from '../src/index.js';

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

function mount(input) {
  document.body.append(input.root);
  return input;
}

const labels = (input) => [...input.root.querySelectorAll('.u2-checklist-item')].map((row) => row.textContent);
const box = (input, value) => input.root.querySelector(`[data-u2-item="${value}"] input`);

function toggle(input, value) {
  const el = box(input, value);
  el.checked = !el.checked;
  fire(el, 'change');
}

smoke('a checkbox per item; the value is the checked ones in the order shown', () => {
  const input = mount(new ChecklistInput({label: 'Run', items: ['Descriptors', 'Toxicity', 'Solubility']}));
  assert.equal(input.root.dataset.u2, 'checklist-input');
  assert.deepEqual(labels(input), ['Descriptors', 'Toxicity', 'Solubility']);
  toggle(input, 'Solubility');
  toggle(input, 'Descriptors');
  assert.deepEqual(input.value, ['Descriptors', 'Solubility']);
  assert.deepEqual(input.getItems().list.map((i) => i.checked), [true, false, true]);
  assert.deepEqual(input.getState().selection, ['Descriptors', 'Solubility']);
  input.dispose();
});

smoke('orderable: the checked items given come first, in the value\'s order; move() reorders the value', () => {
  const input = mount(new ChecklistInput({items: [{value: 'd', label: 'Descriptors'}, {value: 't', label: 'Toxicity'},
    {value: 's', label: 'Solubility'}], value: ['s', 'd'], orderable: true}));
  assert.deepEqual(labels(input), ['Solubility', 'Descriptors', 'Toxicity']);
  assert.deepEqual(input.order.value, ['s', 'd', 't']);
  assert.equal(input.root.querySelectorAll('[data-u2-part="handle"]').length, 3);
  input.move(2, 0);
  assert.deepEqual(input.order.value, ['t', 's', 'd']);
  assert.deepEqual(input.value, ['s', 'd'], 'an unchecked item moved changes no pick');
  input.move(2, 0);
  assert.deepEqual(input.value, ['d', 's']);
  input.dispose();
});

smoke('Alt+↑/↓ moves the row under the focus', () => {
  const input = mount(new ChecklistInput({items: ['a', 'b', 'c'], value: ['a', 'b', 'c'], orderable: true}));
  const b = box(input, 'b');
  b.focus();
  fire(b, 'keydown', {key: 'ArrowUp', altKey: true});
  assert.deepEqual(input.value, ['b', 'a', 'c']);
  fire(b, 'keydown', {key: 'ArrowUp', altKey: true});
  assert.deepEqual(input.value, ['b', 'a', 'c'], 'the first row stays first');
  input.dispose();
});

smoke('a value from elsewhere is checked and, orderable, put first in its order', () => {
  const input = mount(new ChecklistInput({items: ['a', 'b', 'c'], orderable: true}));
  input.value = ['c', 'a'];
  assert.deepEqual(input.order.value, ['c', 'a', 'b']);
  assert.equal(box(input, 'c').checked, true);
  assert.equal(box(input, 'b').checked, false);
  input.dispose();
});

smoke('items that change keep the order and checks of those still there; a vanished pick leaves the value', () => {
  const items = signal(['a', 'b', 'c']);
  const changes = [];
  const input = mount(new ChecklistInput({items, value: ['c', 'b'], orderable: true,
    onInput: (v) => changes.push(v)}));
  items.value = ['b', 'c', 'd'];
  assert.deepEqual(input.order.value, ['c', 'b', 'd']);
  assert.deepEqual(input.value, ['c', 'b']);
  items.value = ['b', 'd'];
  assert.deepEqual(input.value, ['b']);
  assert.deepEqual(changes, [], 'the input\'s own pruning is no user edit');
  input.dispose();
});

smoke('no items: a row of "No items" (or emptyText), gone once items come', () => {
  const items = signal([]);
  const input = mount(new ChecklistInput({items, emptyText: 'No functions'}));
  assert.equal(input.root.querySelector('.u2-checklist-empty').textContent, 'No functions');
  assert.equal(input.getItems().count, 0);
  items.value = ['a'];
  assert.equal(input.root.querySelector('.u2-checklist-empty'), null);
  assert.equal(input.getItems().count, 1);
  input.dispose();
});

smoke('not orderable: the items\' order, no handles; required refuses nothing checked', () => {
  const input = mount(new ChecklistInput({items: ['a', 'b'], value: ['b', 'a'], required: true}));
  assert.deepEqual(labels(input), ['a', 'b']);
  assert.equal(input.root.querySelector('[data-u2-part="handle"]'), null);
  toggle(input, 'a');
  toggle(input, 'b');
  assert.deepEqual(input.value, []);
  assert.equal(input.validity.value, 'Value can\'t be empty');
  input.dispose();
});

smoke('not orderable: new items in a new order show in that order, as MultiChoiceInput\'s; each label is a part', () => {
  const items = signal(['a', 'b', 'c']);
  const input = mount(new ChecklistInput({items, value: ['b']}));
  items.value = ['c', 'b', 'a', 'd'];
  assert.deepEqual(labels(input), ['c', 'b', 'a', 'd']);
  assert.deepEqual(input.value, ['b']);
  assert.deepEqual([...input.root.querySelectorAll('[data-u2-part="item-label"]')].map((el) => el.textContent),
    ['c', 'b', 'a', 'd']);
  const ordered = mount(new ChecklistInput({items: ['a', 'b'], orderable: true}));
  ordered.move(1, 0);
  ordered.items = ['a', 'b', 'c'];
  assert.deepEqual(labels(ordered), ['b', 'a', 'c'], 'orderable keeps the user\'s order, the new item last');
  assert.equal(ordered.root.querySelectorAll('.u2-reorder-grip[data-u2-part="handle"]').length, 3);
  input.dispose();
  ordered.dispose();
});
