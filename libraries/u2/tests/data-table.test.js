/* DataTable (WO 3-6) on the DOM shim: the shim lays nothing out, so each test gives the root a
   viewport (clientHeight) and fires `scroll` to re-render. What it must get right: only the
   visible window plus overscan in the DOM, rows AND cells recycled through the pool, the
   selection kept by ROW KEY across setItems, the cell state read by key (never by index) into
   the `u2-cell-*` classes and the title, the row context menu, and a dispose that leaves
   nothing listening. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {DataTable} from '../src/components/collections/data-table.js';
import {arrayRows} from '../src/sources/rows-like.js';

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

const COLUMNS = [{name: 'name', header: 'Name'}, {name: 'amount', align: 'right', width: '80px'}];

const items = (n) => Array.from({length: n}, (_, i) => ({id: `r${i}`, name: `Row ${i}`, amount: i}));

/** 20px rows in a 220px viewport: the header eats one row, so 10 rows are visible. */
function table(options = {}, count = 100) {
  const t = new DataTable({columns: COLUMNS, rowHeight: 20, keyOf: (item) => item.id, ...options});
  document.body.append(t.root);
  t.root.clientHeight = 220;
  if (options.items === undefined)
    t.setItems(items(count));
  return t;
}

const row = (t, i) => t.root.querySelector(`.u2-data-table-row[data-index="${i}"]`);
const rendered = (t) => t.root.querySelectorAll('.u2-data-table-row')
  .map((r) => Number(r.dataset.index)).sort((a, b) => a - b);
const cells = (t, i) => row(t, i).children;
const scrollTo = (t, top) => {
  t.root.scrollTop = top;
  fire(t.root, 'scroll');
};

scoped('the header is a sticky row of columnheaders; the spacer carries the full row count', () => {
  const t = table();
  const header = t.root.querySelector('.u2-data-table-header');
  assert.equal(header.getAttribute('role'), 'row');
  assert.deepEqual(header.children.map((c) => c.textContent), ['Name', 'amount'],
    'a column without a header falls back to its name');
  assert.equal(header.children[0].getAttribute('role'), 'columnheader');
  assert.equal(header.style.gridTemplateColumns, 'minmax(0, 1fr) 80px');
  assert.equal(t.root.getAttribute('role'), 'grid');
  assert.equal(t.root.getAttribute('data-u2'), 'data-table');
  assert.equal(t.root.querySelector('.u2-data-table-content').style.height, '2000px');
  t.dispose();
});

scoped('only the visible window plus overscan is rendered, and the range follows the scroll', () => {
  const t = table();
  // 10 rows fit under the header, 3 overscan below (none above at the top)
  assert.equal(t.renderedCount, 13);
  assert.deepEqual(rendered(t).slice(0, 3), [0, 1, 2]);
  assert.equal(row(t, 0).style.top, '0px');
  assert.equal(row(t, 5).style.top, '100px');
  assert.equal(row(t, 0).getAttribute('aria-rowindex'), '2', 'the header is row 1');
  assert.deepEqual(cells(t, 3).map((c) => c.textContent), ['Row 3', '3']);
  assert.ok(cells(t, 3)[1].className.includes('u2-data-table-align-right'));

  scrollTo(t, 400);
  assert.deepEqual(rendered(t), Array.from({length: 16}, (_, i) => 17 + i));
  assert.equal(row(t, 0), null, 'the rows that left the window are gone');
  assert.equal(t.renderedCount, 16);
  t.dispose();
});

scoped('rows and their cells are recycled: scrolling creates no new elements', () => {
  const t = table();
  scrollTo(t, 400);
  // every element the window has needed so far, and the cells they own
  const pooled = new Set(t.root.querySelectorAll('.u2-data-table-row'));
  const ownCells = new Set([...pooled].flatMap((r) => r.children));
  assert.equal(pooled.size, 16);

  for (const top of [0, 1000, 240, 1980]) {
    scrollTo(t, top);
    for (const r of t.root.querySelectorAll('.u2-data-table-row')) {
      assert.ok(pooled.has(r), `row ${r.dataset.index} came from the pool`);
      for (const cell of r.children)
        assert.ok(ownCells.has(cell), 'and its cells are the row\'s own, not rebuilt per render');
    }
    assert.equal(t.root.querySelectorAll('.u2-data-table-row').length, t.renderedCount,
      'nothing is left behind in the DOM');
  }
  assert.equal(row(t, 99).children[0].textContent, 'Row 99', 'a recycled row shows its own item');
  t.dispose();
});

scoped('scrollToIndex brings a row under the header, and the keyboard moves the lead', () => {
  const t = table();
  t.scrollToIndex(50);
  assert.equal(t.root.scrollTop, 50 * 20 + 20 - 200);
  fire(t.root, 'keydown', {key: 'Home'});
  assert.equal(t.selectedIndex.value, 0);
  assert.equal(t.root.scrollTop, 0);
  fire(t.root, 'keydown', {key: 'ArrowDown'});
  fire(t.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(t.selectedIndex.value, 2);
  assert.equal(row(t, 2).getAttribute('aria-selected'), 'true');
  assert.equal(t.root.getAttribute('aria-activedescendant'), row(t, 2).id);
  fire(t.root, 'keydown', {key: 'End'});
  assert.equal(t.selectedIndex.value, 99);
  t.dispose();
});

scoped('Enter and a double-click activate the row; a plain click selects, Ctrl and Shift build a set', () => {
  const activated = [];
  const t = table({onActivate: (item, index) => activated.push([item.id, index])});
  fire(row(t, 3), 'click');
  assert.deepEqual([...t.selectedIndices.value], [3]);
  fire(t.root, 'keydown', {key: 'Enter', target: t.root});
  assert.deepEqual(activated, [['r3', 3]]);
  fire(row(t, 4), 'dblclick');
  assert.deepEqual(activated[1], ['r4', 4]);

  fire(row(t, 1), 'click');
  fire(row(t, 5), 'click', {ctrlKey: true});
  assert.deepEqual([...t.selectedIndices.value], [1, 5]);
  assert.equal(t.selectedIndex.value, 5);
  fire(row(t, 3), 'click', {shiftKey: true});
  assert.deepEqual([...t.selectedIndices.value], [1, 2, 3], 'a range from the anchor');
  fire(row(t, 7), 'click');
  assert.deepEqual([...t.selectedIndices.value], [7], 'a plain click collapses the set');
  t.dispose();
});

scoped('the selection follows the ROW KEY across setItems, and a RowsLike brings its own keyOf', () => {
  const t = table({keyOf: undefined, items: undefined});
  const rows = signal(items(10));
  t.setItems(arrayRows(rows, (item) => item.id));
  fire(row(t, 4), 'click');
  assert.equal(t.selectedIndex.value, 4);

  // the same rows with two inserted in front: the selected row is now at 6
  rows.value = [{id: 'x1', name: 'X1', amount: 0}, {id: 'x2', name: 'X2', amount: 0}, ...items(10)];
  assert.equal(t.selectedIndex.value, 6, 'the key, not the index, is what the selection keeps');
  assert.deepEqual([...t.selectedIndices.value], [6], 'and the set collapses to the lead');

  rows.value = items(10).filter((item) => item.id !== 'r4');
  assert.equal(t.selectedIndex.value, -1, 'a row that left takes the selection with it');

  // a second RowsLike REPLACES the first one's identity — the keyed re-select must use the new one
  const uppercase = signal(items(6).map((item) => ({...item, id: item.id.toUpperCase()})));
  t.setItems(arrayRows(uppercase, (item) => item.id));
  fire(row(t, 2), 'click');
  uppercase.value = [{id: 'X', name: 'X', amount: 0}, ...uppercase.value];
  assert.equal(t.selectedIndex.value, 3, 'keyed by the second RowsLike, not the first');
  t.dispose();
});

scoped('an explicit keyOf outranks whatever a RowsLike brings, across setItems', () => {
  const t = table({keyOf: (item) => item.name, items: undefined});
  t.setItems(arrayRows(signal(items(5)), (item) => item.id));
  fire(row(t, 1), 'click');
  // the same names under different ids: the option's keyOf is what the re-select follows
  t.setItems(arrayRows(signal(items(5).map((item) => ({...item, id: `other-${item.id}`}))),
    (item) => item.id));
  assert.equal(t.selectedIndex.value, 1, 'matched by name, which the ids no longer agree with');
  t.dispose();
});

scoped('cell state is read BY KEY: changed, error (with its message as the title) and readonly', () => {
  const changed = new Set(['r1.name']);
  const errors = new Map([['r2.amount', {message: 'Must be at least 1', kind: 'validation'}]]);
  const t = table({cellState: {
    isChanged: (key, column) => changed.has(`${key}.${column}`),
    canEdit: (key, column) => !(key === 'r0' && column === 'amount'),
    errorOf: (key, column) => errors.get(`${key}.${column}`) ?? null,
  }});
  assert.ok(cells(t, 1)[0].className.includes('u2-cell-changed'));
  assert.ok(!cells(t, 1)[1].className.includes('u2-cell-changed'), 'the column is part of the key');
  assert.ok(!cells(t, 0)[0].className.includes('u2-cell-changed'));
  assert.ok(cells(t, 2)[1].className.includes('u2-cell-error'));
  assert.equal(cells(t, 2)[1].getAttribute('title'), 'Must be at least 1');
  assert.equal(cells(t, 2)[0].getAttribute('title'), null);
  assert.ok(cells(t, 0)[1].className.includes('u2-cell-readonly'));

  // a verdict raised over the window already on screen: refresh re-reads it, and the stale
  // classes of a recycled cell are gone
  changed.delete('r1.name');
  changed.add('r3.amount');
  errors.clear();
  t.refresh();
  assert.ok(!cells(t, 1)[0].className.includes('u2-cell-changed'));
  assert.ok(cells(t, 3)[1].className.includes('u2-cell-changed'));
  assert.ok(!cells(t, 2)[1].className.includes('u2-cell-error'));
  assert.equal(cells(t, 2)[1].getAttribute('title'), null);
  t.dispose();
});

scoped('a custom render owns the cell; scrolling never carries a cell\'s classes to another row', () => {
  const t = table({columns: [{name: 'name', render: (item, index, cell) => {
    cell.classList.add(`u2-demo-${item.id}`);
    const el = document.createElement('b');
    el.textContent = item.name.toUpperCase();
    return el;
  }}]});
  assert.equal(cells(t, 0)[0].textContent, 'ROW 0');
  assert.ok(cells(t, 0)[0].className.includes('u2-demo-r0'));
  scrollTo(t, 400);
  const recycled = [...t.root.querySelectorAll('.u2-data-table-row')]
    .filter((r) => r.children[0].className.includes('u2-demo-r0'));
  assert.deepEqual(recycled.map((r) => r.dataset.index), [], 'the class set is reset before each render');
  t.dispose();
});

scoped('the row context menu: inside the selection it keeps the set, outside it collapses', () => {
  const t = table({contextActions: (item) => [{name: `Act on ${item.id}`, run: () => {}}]});
  fire(row(t, 2), 'click');
  fire(row(t, 5), 'click', {ctrlKey: true});
  fire(cells(t, 2)[0], 'contextmenu');
  assert.deepEqual([...t.selectedIndices.value], [2, 5]);
  assert.equal(t.selectedIndex.value, 2, 'with the clicked row as the lead');
  assert.equal(document.querySelector('[role="menuitem"]').textContent, 'Act on r2');
  fire(document.querySelector('[role="menuitem"]'), 'click');

  fire(cells(t, 7)[0], 'contextmenu');
  assert.deepEqual([...t.selectedIndices.value], [7]);
  fire(document.querySelector('[role="menuitem"]'), 'click');
  t.dispose();
});

scoped('a key pressed inside a row belongs to the row: the scroller takes only its own', () => {
  const t = table({columns: [{name: 'name', render: (item, index, cell) => {
    const input = document.createElement('input');
    input.value = item.name;
    cell.append(input);
    return input;
  }}]});
  fire(row(t, 3), 'click');
  fire(row(t, 3).children[0].children[0], 'keydown', {key: 'End'});
  assert.equal(t.selectedIndex.value, 3, 'End in a cell editor keeps its caret, and the lead');
  fire(t.root, 'keydown', {key: 'End'});
  assert.equal(t.selectedIndex.value, 99, 'the same key on the scroller still moves it');
  t.dispose();
});

scoped('an empty table renders no rows; dispose leaves nothing listening', () => {
  const t = table({}, 0);
  assert.equal(t.renderedCount, 0);
  assert.equal(t.root.querySelector('.u2-data-table-content').style.height, '0px');

  t.setItems(items(20));
  assert.equal(t.renderedCount, 13);
  fire(row(t, 2), 'click');
  const before = t.selectedIndex.value;
  t.dispose();
  scrollTo(t, 400);
  fire(t.root, 'keydown', {key: 'End'});
  assert.equal(t.selectedIndex.value, before, 'the keyboard and the scroll handler are gone');
  assert.equal(t.renderedCount, 13, 'and nothing re-rendered');
});
