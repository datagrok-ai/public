/* Keyed virtual rows (R6, P20, P2, P104): an items change keeps the rows whose key stays in the
   window and renders only the items that changed; every row owns a scope its render binds in; the
   selection follows its keys; items are reported and revealed whether rendered or not; a tree box
   shows its node, never a click the owner refused. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {span} from '../src/core/elements.js';
import {VirtualList} from '../src/components/collections/list.js';
import {VirtualTree} from '../src/components/collections/tree.js';

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

/** 10,000 items, 20 px rows in a 400 px scroller: 20 visible + 3 overscan = 23 rendered. */
function bigList(options = {}, count = 10000) {
  const calls = {render: 0};
  const items = Array.from({length: count}, (_, i) => ({id: `r${i}`, name: `Item ${i}`}));
  const list = new VirtualList({
    itemHeight: 20,
    keyOf: (item) => item.id,
    render: (item) => {
      calls.render++;
      return span(item.name);
    },
    ...options,
  });
  document.body.append(list.root);
  layout(list.root, {clientHeight: 400});
  list.setItems(items);
  return {list, items, calls};
}

const rowAt = (list, i) => list.root.querySelector(`.u2-list-row[data-index="${i}"]`);
const texts = (list) => [...list.root.querySelectorAll('.u2-list-row')]
  .sort((a, b) => Number(a.dataset.index) - Number(b.dataset.index)).map((r) => r.textContent);

function scrollTo(list, top) {
  list.root.scrollTop = top;
  fire(list.root, 'scroll');
}

scoped('one immutable rename in a 10,000-row list renders exactly one row (was 23 of 23)', () => {
  const {list, items, calls} = bigList();
  assert.equal(list.renderedCount, 23);
  const kept = rowAt(list, 6);
  calls.render = 0;
  const next = items.slice();
  next[5] = {...next[5], name: 'Renamed'};
  list.setItems(next);
  assert.equal(calls.render, 1);
  assert.equal(rowAt(list, 5).textContent, 'Renamed');
  assert.equal(rowAt(list, 6), kept, 'the other rows keep their elements');
  list.dispose();
});

scoped('an insert at the top moves the kept rows and renders only the new item', () => {
  const {list, items, calls} = bigList();
  const first = rowAt(list, 0);
  calls.render = 0;
  list.setItems([{id: 'new', name: 'New'}, ...items]);
  assert.equal(calls.render, 1, 'the new item alone');
  assert.equal(rowAt(list, 1), first, 'the old first row moved down');
  assert.equal(first.getAttribute('aria-posinset'), '2');
  assert.equal(first.getAttribute('aria-setsize'), '10001');
  assert.equal(first.style.top, '20px');
  assert.deepEqual(texts(list).slice(0, 3), ['New', 'Item 0', 'Item 1']);
  list.dispose();
});

scoped('an item mutated in place shows its change through refresh(key), one row only', () => {
  const {list, items, calls} = bigList();
  items[3].name = 'Changed';
  list.setItems(items.slice());
  assert.equal(rowAt(list, 3).textContent, 'Item 3', 'the list cannot see a mutation');
  calls.render = 0;
  list.refresh('r3');
  assert.equal(calls.render, 1);
  assert.equal(rowAt(list, 3).textContent, 'Changed');
  list.dispose();
});

scoped('update patches a replaced item in place, in the row scope; false renders it again', () => {
  const patched = [];
  const {list, items, calls} = bigList({
    update: (item, index, row, el) => {
      if (item.name === 'refuse')
        return false;
      patched.push(index);
      el.textContent = item.name;
    },
  });
  calls.render = 0;
  const next = items.slice();
  next[2] = {...next[2], name: 'Patched'};
  next[4] = {...next[4], name: 'refuse'};
  list.setItems(next);
  assert.deepEqual(patched, [2]);
  assert.equal(calls.render, 1, 'only the refused patch rendered again');
  assert.deepEqual(texts(list).slice(2, 5), ['Patched', 'Item 3', 'refuse']);
  list.dispose();
});

scoped('liveItems renders the window again on every items change', () => {
  const {list, items, calls} = bigList({liveItems: true});
  calls.render = 0;
  list.setItems(items.slice());
  assert.equal(calls.render, 23);
  list.dispose();
});

scoped('a signal bound in render follows its item and stops with the row', async () => {
  const names = Array.from({length: 100}, (_, i) => signal(`n${i}`));
  const list = new VirtualList({itemHeight: 20, keyOf: (s) => names.indexOf(s).toString(), render: (s) => span(s)});
  document.body.append(list.root);
  layout(list.root, {clientHeight: 200});
  const before = Scope.liveCount;
  list.setItems(names);
  assert.equal(Scope.liveCount - before, list.renderedCount, 'one live scope per row that binds');
  names[2].value = 'live';
  assert.equal(rowAt(list, 2).textContent, 'live', 'the binding follows without a render');

  const element = rowAt(list, 2);
  scrollTo(list, 60 * 20);
  assert.equal(rowAt(list, 2), null, 'row 2 left the window');
  assert.equal(Scope.liveCount - before, list.renderedCount, 'its scope went with it');
  names[2].value = 'stale';
  assert.notEqual(element.textContent, 'stale', 'the released row no longer follows item 2');
  list.dispose();
});

scoped('scrolling back and forth ten times: every row shows its own item, bindings included', () => {
  const names = Array.from({length: 1000}, (_, i) => signal(`n${i}`));
  const list = new VirtualList({itemHeight: 20, keyOf: (s) => s.peek(), render: (s) => span(s)});
  document.body.append(list.root);
  layout(list.root, {clientHeight: 200});
  list.setItems(names);
  for (let round = 0; round < 10; round++) {
    scrollTo(list, (round % 2 === 0 ? 400 + round * 13 : round * 7) * 20);
    for (const row of list.root.querySelectorAll('.u2-list-row'))
      assert.equal(row.textContent, names[Number(row.dataset.index)].peek());
  }
  for (const [i, s] of names.entries())
    s.value = `m${i}`;
  for (const row of list.root.querySelectorAll('.u2-list-row'))
    assert.equal(row.textContent, `m${row.dataset.index}`, 'each binding writes into its own row');
  list.dispose();
});

scoped('a Ctrl-selection of three keys survives setItems of a reordered list', () => {
  const {list, items} = bigList();
  fire(rowAt(list, 2), 'click');
  fire(rowAt(list, 4), 'click', {ctrlKey: true});
  fire(rowAt(list, 6), 'click', {ctrlKey: true});
  assert.deepEqual(list.getState().selection, ['r2', 'r4', 'r6']);
  list.setItems([items[6], items[0], items[4], items[9], items[2], ...items.slice(10)]);
  assert.deepEqual([...list.selectedIndices.value].sort(), [0, 2, 4]);
  assert.equal(list.selectedIndex.value, 0, 'the lead (r6) kept its key');
  assert.deepEqual(list.getState().selection.sort(), ['r2', 'r4', 'r6']);
  assert.equal(rowAt(list, 4).getAttribute('aria-selected'), 'true');
  assert.equal(rowAt(list, 1).getAttribute('aria-selected'), 'false');
  list.dispose();
});

scoped('items report every item, not the window; rows carry data-u2-item', () => {
  const {list} = bigList();
  layout(list.root, {clientHeight: 340});
  scrollTo(list, 0);
  assert.equal(list.renderedCount, 20);
  const items = list.getItems();
  assert.equal(items.count, 10000);
  assert.equal(items.list.length, 500);
  assert.equal(items.truncated, true);
  assert.deepEqual(items.list[3], {key: 'r3', label: 'Item 3', index: 3, selected: false, disabled: false, rendered: true});
  assert.deepEqual(items.list[499], {key: 'r499', label: 'r499', index: 499, selected: false, disabled: false,
    rendered: false}, 'an item off screen is labelled by its key');
  assert.equal(list.getItems({all: true}).list.length, 10000);
  assert.equal(rowAt(list, 7).dataset.u2Item, 'r7');
  list.dispose();
});

scoped('keys default to the string item, then name or id, then the position', () => {
  const list = new VirtualList({render: (x) => span(String(x?.label ?? x))});
  document.body.append(list.root);
  layout(list.root, {clientHeight: 200});
  list.setItems(['a', 1, {name: 'n'}, {id: 'i'}, {label: 'plain'}, null]);
  assert.deepEqual(list.getItems().list.map((i) => i.key), ['a', '1', 'n', 'i', '#4', '#5']);
  list.dispose();
});

scoped('revealItem scrolls item 900 into the window and resolves to its row', async () => {
  const {list} = bigList({labelOf: (item) => item.name});
  assert.equal(rowAt(list, 900), null);
  const row = await list.revealItem('r900');
  assert.equal(row, rowAt(list, 900));
  assert.equal(row.textContent, 'Item 900');
  assert.equal(await list.revealItem('item 901'), rowAt(list, 901), 'by label, any case');
  assert.equal(await list.revealItem('nothing'), null);
  list.dispose();
});

scoped('an ambiguous label rejects naming the candidates', async () => {
  const list = new VirtualList({keyOf: (x) => x.id, labelOf: (x) => x.name, render: (x) => span(x.name)});
  document.body.append(list.root);
  layout(list.root, {clientHeight: 200});
  list.setItems([{id: 'a', name: 'Same'}, {id: 'b', name: 'Other'}, {id: 'c', name: 'same'}]);
  await assert.rejects(list.revealItem('Same'), {message: 'u2: 2 items match "Same": a, c'});
  assert.equal(await list.revealItem('c'), rowAt(list, 2), 'the key wins over a label');
  list.dispose();
});

function checkTree(onCheck) {
  const tree = new VirtualTree({onCheck});
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 220});
  tree.setRoots([{id: 'a', label: 'Alpha', checked: false}, {id: 'b', label: 'Beta', checked: true}]);
  return tree;
}

const treeRow = (tree, i) => tree.root.querySelector(`.u2-list-row[data-index="${i}"]`);

scoped('tree: a refused onCheck leaves the box and aria-checked on the node state (P104)', () => {
  const asked = [];
  const tree = checkTree((node, checked) => asked.push([node.id, checked]));
  const box = treeRow(tree, 0).querySelector('.u2-tree-check');
  box.checked = true; // what the browser does before the handlers run
  const proceeded = fire(box, 'click');
  assert.equal(proceeded, false, 'the click is prevented, so the browser puts the box back');
  assert.deepEqual(asked, [['a', true]]);
  assert.equal(treeRow(tree, 0).getAttribute('aria-checked'), 'false', 'the owner refused: still unchecked');
  assert.equal(tree.getItems().list[0].checked, false);
  tree.dispose();
});

scoped('tree: the owner\'s answer re-renders that row only; no onCheck disables the boxes', () => {
  let roots;
  const tree = checkTree((node, checked) => {
    roots = roots.map((n) => n.id === node.id ? {...n, checked} : n);
    tree.setRoots(roots);
  });
  roots = [{id: 'a', label: 'Alpha', checked: false}, {id: 'b', label: 'Beta', checked: true}];
  tree.setRoots(roots);
  const beta = treeRow(tree, 1);
  fire(treeRow(tree, 0).querySelector('.u2-tree-check'), 'click');
  assert.equal(treeRow(tree, 0).getAttribute('aria-checked'), 'true');
  assert.equal(treeRow(tree, 0).querySelector('.u2-tree-check').checked, true);
  assert.equal(treeRow(tree, 1), beta, 'Beta kept its row');
  tree.dispose();

  const passive = new VirtualTree();
  document.body.append(passive.root);
  layout(passive.root.querySelector('.u2-list'), {clientHeight: 220});
  passive.setRoots([{id: 'a', label: 'Alpha', checked: true}]);
  assert.equal(treeRow(passive, 0).querySelector('.u2-tree-check').disabled, true, 'nothing would answer it');
  passive.dispose();
});

scoped('tree: a node changed in place shows once the roots are set again, and only that row renders', () => {
  const renders = [];
  const roots = [{id: 'a', label: 'Alpha'}, {id: 'b', label: 'Beta'}, {id: 'c', label: 'Gamma'}];
  const tree = new VirtualTree({render: (node) => {
    renders.push(node.id);
    return span(node.label);
  }});
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 220});
  tree.setRoots(roots);
  renders.length = 0;
  roots[1].label = 'Bravo';
  tree.setRoots(roots);
  assert.deepEqual(renders, ['b']);
  assert.equal(treeRow(tree, 1).querySelector('.u2-tree-label').textContent, 'Bravo');
  tree.dispose();
});

scoped('tree: items are the visible nodes with level, expansion, checks and selection; reveal opens ancestors', async () => {
  const tree = new VirtualTree();
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 220});
  tree.setRoots([
    {id: 'p', label: 'Project', children: [
      {id: 'p/t', label: 'Tables', children: [{id: 'p/t/demog', label: 'demog'}]},
    ]},
    {id: 'q', label: 'Other', checked: true},
  ]);
  assert.deepEqual(tree.getItems().list.map((i) => [i.key, i.level, i.expanded ?? null, i.checked ?? null]),
    [['p', 1, false, null], ['q', 1, null, true]]);
  const row = await tree.revealItem('demog');
  assert.equal(row.dataset.u2Item, 'p/t/demog');
  assert.deepEqual([...tree.expanded.value].sort(), ['p', 'p/t']);
  assert.equal(tree.selectedNode.value, null, 'revealed, not selected');
  fire(row.querySelector('.u2-tree-label'), 'click');
  const state = tree.getState();
  assert.deepEqual(state.selection, ['p/t/demog']);
  assert.deepEqual(state.items.list.map((i) => [i.key, i.index, i.selected]),
    [['p', 0, false], ['p/t', 1, false], ['p/t/demog', 2, true], ['q', 3, false]]);
  const inner = tree.root.querySelector('.u2-list');
  assert.deepEqual((await import('../src/core/component.js')).Control.forElement(inner).getItems().count, 4,
    'the scroller answers for the tree');
  assert.equal(row.querySelector('[data-u2-part="item-label"]').textContent, 'demog');
  tree.dispose();
});

scoped('tree: a lazy read is owned by the tree, not by the rebuild that started it; a failure is read again', async () => {
  const signals = [];
  let fail = true;
  const tree = new VirtualTree();
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 220});
  tree.setRoots([
    {id: 'lazy', label: 'Lazy', children: () => {
      signals.push(Scope.ambient.abortSignal);
      return new Promise((resolve, reject) => setTimeout(() => fail ?
        reject(new Error('offline')) : resolve([{id: 'child', label: 'Child'}]), 5));
    }},
    {id: 'other', label: 'Other', children: [{id: 'o1', label: 'o1'}]},
  ]);
  tree.expanded.value = new Set(['lazy']);
  tree.expanded.value = new Set(['lazy', 'other']);
  assert.equal(signals.length, 1, 'read once');
  assert.equal(signals[0].aborted, false, 'a rebuild does not abort the read in flight');
  await new Promise((resolve) => setTimeout(resolve, 20));
  assert.equal(tree.root.querySelector('.u2-tree-error').textContent, 'offline');
  fire(treeRow(tree, 0).querySelector('.u2-tree-twistie'), 'click');
  fail = false;
  fire(treeRow(tree, 0).querySelector('.u2-tree-twistie'), 'click');
  assert.equal(signals.length, 2, 'collapsing and opening a failed branch reads it again');
  await new Promise((resolve) => setTimeout(resolve, 20));
  assert.equal(treeRow(tree, 1).dataset.u2Item, 'child');
  tree.dispose();
  assert.equal(signals[1].aborted, true, 'the tree aborts what it still owns');
});

scoped('tree: expandPath keeps the branches opened or closed while it awaited a lazy read', async () => {
  const tree = new VirtualTree();
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 220});
  let answer;
  tree.setRoots([
    {id: 'lazy', label: 'Lazy', children: () => new Promise((resolve) => answer = resolve)},
    {id: 'a', label: 'A', children: [{id: 'a1', label: 'a1'}]},
    {id: 'b', label: 'B', children: [{id: 'b1', label: 'b1'}]},
  ]);
  tree.expanded.value = new Set(['b']);
  const done = tree.expandPath(['lazy', 'child']);
  tree.expanded.value = new Set(['a']);
  answer([{id: 'child', label: 'Child'}]);
  await done;
  assert.deepEqual([...tree.expanded.value].sort(), ['a', 'lazy']);
  assert.equal(tree.selectedNode.value?.id, 'child');
  tree.dispose();
});

scoped('tree: a lazy read is dropped once its node leaves the tree; a collapsed branch keeps its own', async () => {
  const reads = new Map();
  let count = 0;
  const lazy = (id, children = []) => ({id, label: id, children: () => {
    count++;
    reads.set(id, Scope.ambient.abortSignal);
    return Promise.resolve(children);
  }});
  const x = lazy('x', [lazy('x/sub')]);
  const y = lazy('y');
  const tree = new VirtualTree();
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 220});
  tree.setRoots([x, y]);
  tree.expanded.value = new Set(['x', 'y']);
  await flush();
  tree.expanded.value = new Set(['x', 'x/sub', 'y']);
  await flush();
  assert.deepEqual([...reads.keys()].sort(), ['x', 'x/sub', 'y']);
  tree.expanded.value = new Set(['y']);
  tree.setRoots([x, y]);
  assert.equal(reads.get('x/sub').aborted, false, 'under a collapsed branch, still in the tree');
  tree.setRoots([y]);
  assert.deepEqual(['x', 'x/sub', 'y'].map((id) => reads.get(id).aborted), [true, true, false]);
  assert.equal(count, 3, 'nothing read again');
  tree.dispose();
  assert.equal(reads.get('y').aborted, true);
});
