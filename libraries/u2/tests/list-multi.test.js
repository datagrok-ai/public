/* F5 — multi-selection on the virtualized list and its tree projection: Ctrl toggles, Shift
   ranges from the anchor (the last plain click), the contextmenu keeps a set it lands in, and
   every single-select path — plain click, keyboard, setItems, a programmatic lead write —
   collapses the set to the lead alone, so single-select consumers never see a ghost. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {batch} from '../src/core/signals.js';
import {VirtualList} from '../src/components/list.js';
import {VirtualTree} from '../src/components/tree.js';

function scoped(name, body) {
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

function list(options = {}) {
  const l = new VirtualList({
    itemHeight: 20,
    render: (item) => {
      const el = document.createElement('span');
      el.textContent = String(item);
      return el;
    },
    ...options,
  });
  document.body.append(l.root);
  l.root.clientHeight = 220;
  l.setItems(Array.from({length: 10}, (_, i) => `item${i}`));
  return l;
}

const row = (l, i) => l.root.querySelector(`.u2-list-row[data-index="${i}"]`);
const selected = (l) => [...l.selectedIndices.value];
const click = (l, i, init = {}) => fire(row(l, i), 'click', init);

scoped('Ctrl-click toggles membership and the lead follows the toggle', () => {
  const l = list();
  click(l, 2);
  assert.deepEqual(selected(l), [2]);
  assert.equal(l.selectedIndex.value, 2);

  click(l, 5, {ctrlKey: true});
  assert.deepEqual(selected(l), [2, 5], 'selection order, not index order');
  assert.equal(l.selectedIndex.value, 5, 'the lead follows the added member');
  assert.equal(row(l, 2).getAttribute('aria-selected'), 'true');
  assert.equal(row(l, 5).getAttribute('aria-selected'), 'true');
  assert.equal(row(l, 3).getAttribute('aria-selected'), 'false');

  click(l, 7, {metaKey: true});
  assert.deepEqual(selected(l), [2, 5, 7], 'Cmd toggles too');

  click(l, 5, {ctrlKey: true});
  assert.deepEqual(selected(l), [2, 7]);
  assert.equal(l.selectedIndex.value, 7, 'removing a non-lead leaves the lead alone');

  click(l, 7, {ctrlKey: true});
  assert.deepEqual(selected(l), [2]);
  assert.equal(l.selectedIndex.value, 2, 'removing the lead hands it to a remaining member');

  click(l, 2, {ctrlKey: true});
  assert.deepEqual(selected(l), [2], 'toggling the sole member off is a no-op — a selection must exist');
  assert.equal(l.selectedIndex.value, 2);
  l.dispose();
});

scoped('a gesture wrapped in an outer batch keeps its set past the deferred flush', () => {
  const l = list();
  click(l, 2);
  batch(() => click(l, 5, {ctrlKey: true}));
  assert.deepEqual(selected(l), [2, 5], 'the collapse effect consumed the gesture, not the set');
  assert.equal(l.selectedIndex.value, 5);

  batch(() => {
    click(l, 7, {ctrlKey: true});
    click(l, 8, {ctrlKey: true});
  });
  assert.deepEqual(selected(l), [2, 5, 7, 8], 'two gestures in one batch both survive');

  l.selectedIndex.value = 1;
  assert.deepEqual(selected(l), [1], 'a programmatic lead write still collapses afterwards');
  l.dispose();
});

scoped('Shift-click ranges from the last plain click, replacing the set', () => {
  const l = list();
  click(l, 3);
  click(l, 6, {shiftKey: true});
  assert.deepEqual(selected(l), [3, 4, 5, 6]);
  assert.equal(l.selectedIndex.value, 6, 'the lead is the clicked end');

  click(l, 1, {shiftKey: true});
  assert.deepEqual(selected(l), [3, 2, 1], 'the range runs anchor-to-click, replacing the set');
  assert.equal(l.selectedIndex.value, 1);

  click(l, 8, {ctrlKey: true});
  assert.deepEqual(selected(l), [3, 2, 1, 8]);
  click(l, 5, {shiftKey: true});
  assert.deepEqual(selected(l), [3, 4, 5], 'a Ctrl toggle never moves the anchor');
  l.dispose();
});

scoped('keyboard, programmatic lead writes and setItems collapse to the lead alone', () => {
  const l = list({keyOf: (item) => item});
  click(l, 2);
  click(l, 5, {ctrlKey: true});
  assert.deepEqual(selected(l), [2, 5]);

  fire(l.root, 'keydown', {key: 'ArrowDown'});
  assert.deepEqual(selected(l), [6], 'an arrow move is a single selection');

  click(l, 8, {ctrlKey: true});
  l.selectedIndex.value = 1;
  assert.deepEqual(selected(l), [1], 'a programmatic lead write collapses');
  assert.equal(l.selectedIndex.value, 1);

  click(l, 4, {ctrlKey: true});
  assert.deepEqual(selected(l), [1, 4]);
  l.setItems(['item4', 'other', 'another']);
  assert.equal(l.selectedIndex.value, 0, 'the key preserved the lead across setItems');
  assert.deepEqual(selected(l), [0], 'and the lead alone');
  l.dispose();
});

scoped('the contextmenu keeps the set it lands in, and collapses to a row outside it', () => {
  const l = list({contextActions: (item) => [{name: `Act on ${item}`, run: () => {}}]});
  click(l, 2);
  click(l, 5, {ctrlKey: true});

  fire(row(l, 2).firstChild, 'contextmenu');
  assert.deepEqual(selected(l), [2, 5], 'a right-click inside the selection keeps it');
  assert.equal(l.selectedIndex.value, 2, 'with the clicked row as the lead');
  fire(document.querySelector('[role="menuitem"]'), 'click');

  fire(row(l, 7).firstChild, 'contextmenu');
  assert.deepEqual(selected(l), [7], 'outside the selection it collapses to the row');
  assert.equal(l.selectedIndex.value, 7);
  fire(document.querySelector('[role="menuitem"]'), 'click');
  l.dispose();
});

scoped('tree: selectedNodes projects the set in selection order, and clearSelection empties both', () => {
  const tree = new VirtualTree();
  document.body.append(tree.root);
  tree.root.querySelector('.u2-list').clientHeight = 220;
  tree.expanded.value = new Set(['b']);
  tree.setRoots([
    {id: 'a', label: 'a', data: 'A'},
    {id: 'b', label: 'b', data: 'B', children: [
      {id: 'c', label: 'c', data: 'C'},
      {id: 'd', label: 'd', data: 'D'},
    ]},
    {id: 'e', label: 'e', data: 'E'},
  ]);
  const rowAt = (i) => tree.root.querySelector(`.u2-list-row[data-index="${i}"]`);

  fire(rowAt(2).querySelector('.u2-tree-label'), 'click');
  fire(rowAt(4).querySelector('.u2-tree-label'), 'click', {ctrlKey: true});
  assert.deepEqual(tree.selectedNodes.value.map((n) => n.data), ['C', 'E']);
  assert.equal(tree.selectedNode.value.data, 'E', 'selectedNode stays the lead');

  fire(rowAt(1).querySelector('.u2-tree-label'), 'click', {shiftKey: true});
  assert.deepEqual(tree.selectedNodes.value.map((n) => n.data), ['C', 'B'],
    'the range runs from the anchor toward the click');

  tree.clearSelection();
  assert.deepEqual(tree.selectedNodes.value, []);
  assert.equal(tree.selectedNode.value, null);
  tree.dispose();
});
