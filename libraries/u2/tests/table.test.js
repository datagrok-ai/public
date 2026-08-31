/* BasicTable and tableFromMap on the node DOM shim. Same contract as smoke.test.js — every test
   leaves the live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, computed, Scope, Control} from '../src/index.js';
import {span} from '../src/core/elements.js';
import {BasicTable, tableFromMap} from '../src/components/collections/table.js';

function table(name, body) {
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

const ITEMS = [
  {name: 'Aspirin', mw: 180.16},
  {name: 'Caffeine', mw: 194.19},
  {name: 'Ibuprofen', mw: 206.28},
];

function columns() {
  return [
    {header: 'Name', render: (c) => c.name, width: '60%'},
    {header: 'MW', render: (c) => c.mw.toFixed(2), align: 'right'},
  ];
}

function rows(t) {
  return t.root.querySelectorAll('.u2-table-row');
}

function cells(row) {
  return row.querySelectorAll('td').map((td) => td.textContent);
}

table('renders a header, one row per item and the declared widths', () => {
  const t = new BasicTable({columns: columns(), items: ITEMS});
  document.body.append(t.root);

  assert.equal(t.root.tagName, 'TABLE');
  assert.equal(t.root.dataset.u2, 'table');
  assert.deepEqual(t.root.querySelectorAll('th').map((th) => th.textContent), ['Name', 'MW']);
  assert.equal(t.root.querySelectorAll('col')[0].style.width, '60%');
  assert.equal(rows(t).length, 3);
  assert.deepEqual(cells(rows(t)[0]), ['Aspirin', '180.16']);
  assert.deepEqual(rows(t).map((r) => r.dataset.index), ['0', '1', '2']);
  t.dispose();
});

table('headerVisible: false drops the header row entirely', () => {
  const t = new BasicTable({columns: columns(), items: ITEMS, headerVisible: false});
  assert.equal(t.root.querySelectorAll('thead').length, 0);
  assert.equal(t.root.querySelectorAll('th').length, 0);
  assert.equal(rows(t).length, 3);
  t.dispose();
});

table('string cells become text, element cells are appended, alignment is a class', () => {
  const t = new BasicTable({
    columns: [
      {header: 'Name', render: (c) => c.name},
      {header: 'MW', render: (c) => span(String(c.mw)), align: 'right'},
      {header: 'Mid', render: () => 'x', align: 'center'},
    ],
    items: ITEMS,
  });
  document.body.append(t.root);

  const cell = rows(t)[0].querySelectorAll('td');
  assert.equal(cell[0].children.length, 0, 'a string cell has no element children');
  assert.equal(cell[0].textContent, 'Aspirin');
  assert.equal(cell[1].children[0].tagName, 'SPAN');
  assert.equal(cell[1].classList.contains('u2-table-align-right'), true);
  assert.equal(cell[2].classList.contains('u2-table-align-center'), true);
  assert.equal(cell[0].classList.contains('u2-table-align-left'), false, 'no class without align');
  assert.equal(t.root.querySelectorAll('th')[1].classList.contains('u2-table-align-right'), true);
  t.dispose();
});

table('a signal source re-renders the rows live', () => {
  const items = signal(ITEMS.slice(0, 2));
  const t = new BasicTable({columns: columns(), items});
  document.body.append(t.root);
  assert.equal(rows(t).length, 2);

  items.value = ITEMS;
  assert.equal(rows(t).length, 3);
  assert.deepEqual(cells(rows(t)[2]), ['Ibuprofen', '206.28']);

  t.setItems([ITEMS[0]]);
  assert.equal(rows(t).length, 1);
  items.value = [];
  assert.equal(rows(t).length, 1, 'setItems switched the source away from the signal');
  t.dispose();
});

table('selection: click, arrows, Home/End, aria-selected and the clamp', () => {
  const t = new BasicTable({columns: columns(), items: ITEMS, selectable: true});
  document.body.append(t.root);
  assert.equal(t.root.tabIndex, 0);
  assert.equal(t.selectedIndex.value, -1);

  fire(rows(t)[1], 'click');
  assert.equal(t.selectedIndex.value, 1);
  assert.equal(rows(t)[1].getAttribute('aria-selected'), 'true');
  assert.equal(rows(t)[1].classList.contains('u2-table-row-selected'), true);
  assert.equal(rows(t)[0].hasAttribute('aria-selected'), false);

  fire(t.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(t.selectedIndex.value, 2);
  fire(t.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(t.selectedIndex.value, 2, 'stops at the last row');
  fire(t.root, 'keydown', {key: 'Home'});
  assert.equal(t.selectedIndex.value, 0);
  fire(t.root, 'keydown', {key: 'End'});
  assert.equal(t.selectedIndex.value, 2);
  assert.equal(rows(t)[2].getAttribute('aria-selected'), 'true');
  fire(t.root, 'keydown', {key: 'ArrowUp'});
  assert.equal(t.selectedIndex.value, 1);

  t.setItems([ITEMS[0]]);
  assert.equal(t.selectedIndex.value, -1, 'a shorter list clears a selection past its end');
  t.dispose();
});

table('onRowClick fires on click and on Enter, and is not wired to keys without a selection', () => {
  const clicks = [];
  const t = new BasicTable({
    columns: columns(),
    items: ITEMS,
    selectable: true,
    onRowClick: (item, index) => clicks.push(`${item.name}:${index}`),
  });
  document.body.append(t.root);

  fire(t.root, 'keydown', {key: 'Enter'});
  assert.deepEqual(clicks, [], 'nothing selected, nothing fired');

  fire(rows(t)[2].querySelectorAll('td')[0], 'click');
  assert.deepEqual(clicks, ['Ibuprofen:2'], 'a click on a cell finds its row');
  fire(t.root, 'keydown', {key: 'Enter'});
  assert.deepEqual(clicks, ['Ibuprofen:2', 'Ibuprofen:2']);
  t.dispose();
});

table('cell renderers run under a per-render scope that is released on re-render', () => {
  const ticks = signal(0);
  const items = signal(ITEMS);
  const t = new BasicTable({
    columns: [
      {header: 'Name', render: (c) => c.name},
      {header: 'Ticks', render: () => span(computed(() => String(ticks.value)))},
    ],
    items,
  });
  document.body.append(t.root);

  const bound = () => rows(t).map((r) => r.querySelectorAll('td')[1].textContent);
  assert.deepEqual(bound(), ['0', '0', '0']);
  ticks.value = 1;
  assert.deepEqual(bound(), ['1', '1', '1']);

  const live = Scope.liveCount;
  for (let i = 0; i < 5; i++)
    items.value = ITEMS.slice();
  assert.equal(Scope.liveCount, live, 'each re-render disposes the previous row scope');

  ticks.value = 2;
  assert.deepEqual(bound(), ['2', '2', '2'], 'the rebuilt rows are bound to the same signal');
  t.dispose();
  ticks.value = 3;
  assert.deepEqual(bound(), ['2', '2', '2'], 'dispose releases the row bindings');
});

table('tableFromMap: two columns, muted keys, element values', () => {
  const value = span('Approved');
  const map = tableFromMap({Name: 'Aspirin', MW: 180.16, Status: value, Missing: null});
  document.body.append(map);

  assert.equal(map.tagName, 'TABLE');
  assert.equal(map.classList.contains('u2-table-map'), true);
  assert.equal(map.querySelectorAll('thead').length, 0);
  const all = map.querySelectorAll('.u2-table-row');
  assert.equal(all.length, 4);
  assert.deepEqual(all.map((r) => r.querySelectorAll('td')[0].textContent), ['Name', 'MW', 'Status', 'Missing']);
  assert.equal(all[0].querySelectorAll('td')[0].classList.contains('u2-table-key'), true);
  assert.equal(all[1].querySelectorAll('td')[1].textContent, '180.16');
  assert.equal(all[2].querySelectorAll('td')[1].children[0], value);
  assert.equal(all[3].querySelectorAll('td')[1].textContent, '');
});

table('dispose: listeners and effects die with the component scope', () => {
  const items = signal(ITEMS);
  let clicks = 0;
  const t = Control.build(() => [new BasicTable({
    columns: columns(),
    items,
    selectable: true,
    onRowClick: () => clicks++,
  })]);
  document.body.append(t.root);

  const root = t.root.querySelector('.u2-table');
  fire(root.querySelectorAll('.u2-table-row')[0], 'click');
  assert.equal(clicks, 1);

  t.dispose();
  fire(root.querySelectorAll('.u2-table-row')[0], 'click');
  assert.equal(clicks, 1, 'the click listener is gone');
  items.value = ITEMS.slice(0, 1);
  assert.equal(root.querySelectorAll('.u2-table-row').length, 3, 'the items effect is gone');
});
