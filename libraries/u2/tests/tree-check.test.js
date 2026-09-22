/* VirtualTree checkboxes: a node with `checked` carries a box, `locked` disables it, `disabled`
   greys the row; a click on the box reports the flipped state through `onCheck` without moving
   the selection, and Space toggles the selected row. Every test leaves the live-scope count
   where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {VirtualTree} from '../src/components/collections/tree.js';

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

function mount(checks) {
  const tree = new VirtualTree({onCheck: (node, checked) => checks.push([node.id, checked])});
  document.body.append(tree.root);
  tree.root.querySelector('.u2-list').clientHeight = 220;
  tree.setRoots([
    {id: 'a', label: 'Alpha', checked: true, children: [
      {id: 'a1', label: 'Alpha one', checked: true, locked: true},
      {id: 'a2', label: 'Alpha two', checked: false, disabled: true},
    ]},
    {id: 'b', label: 'Beta'},
  ]);
  tree.expanded.value = new Set(['a']);
  return tree;
}

const row = (tree, index) => tree.root.querySelector(`.u2-list-row[data-index="${index}"]`);
const box = (tree, index) => row(tree, index).querySelector('.u2-tree-check');

ui('a node with `checked` carries a box; locked disables it; disabled greys the row; no state, no box', () => {
  const tree = mount([]);
  assert.equal(box(tree, 0).checked, true);
  assert.equal(box(tree, 0).disabled, false);
  assert.equal(box(tree, 1).checked, true);
  assert.equal(box(tree, 1).disabled, true, 'locked');
  assert.equal(box(tree, 2).checked, false);
  assert.equal(row(tree, 2).querySelector('.u2-tree-row').classList.contains('u2-tree-row-disabled'), true);
  assert.equal(row(tree, 2).getAttribute('aria-disabled'), 'true');
  assert.equal(box(tree, 3), null, 'Beta has no checkbox');
  assert.equal(box(tree, 0).getAttribute('aria-label'), 'Include Alpha');
  tree.dispose();
});

ui('a click on the box reports the flipped state and leaves the selection alone; a locked box reports nothing', () => {
  const checks = [];
  const tree = mount(checks);
  fire(row(tree, 3).querySelector('.u2-tree-row'), 'click');
  assert.equal(tree.selectedNode.value.id, 'b');
  fire(box(tree, 0), 'click');
  assert.deepEqual(checks, [['a', false]]);
  assert.equal(tree.selectedNode.value.id, 'b', 'the selection did not move to Alpha');
  fire(box(tree, 1), 'click');
  assert.deepEqual(checks, [['a', false]], 'locked');
  tree.dispose();
});

ui('Space toggles the selected row; a row without a box ignores it', () => {
  const checks = [];
  const tree = mount(checks);
  const list = tree.root.querySelector('.u2-list');
  fire(row(tree, 2).querySelector('.u2-tree-row'), 'click');
  assert.equal(tree.selectedNode.value.id, 'a2');
  fire(list, 'keydown', {key: ' '});
  assert.deepEqual(checks, [['a2', true]]);
  fire(row(tree, 3).querySelector('.u2-tree-row'), 'click');
  fire(list, 'keydown', {key: ' '});
  assert.deepEqual(checks, [['a2', true]]);
  tree.dispose();
});
