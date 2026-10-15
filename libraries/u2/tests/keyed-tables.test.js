/* Keyed tables and the property grid (R6, P25, P82): BasicTable rows are built once per key (row
   hook, section rows, selection by key), DataTable keeps unchanged rows of plain items and re-renders
   live ones, VirtualGrid cells are keyed and addressable, a PropertyGrid keystroke wakes one row, and
   each of them answers the state block's items. Every test leaves the live-scope count where it
   found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal, computed} from '../src/core/signals.js';
import {span} from '../src/core/elements.js';
import {BasicTable} from '../src/components/collections/table.js';
import {DataTable} from '../src/components/collections/data-table.js';
import {VirtualGrid} from '../src/components/collections/grid.js';
import {FunctionsBrowser} from '../src/components/collections/functions-browser.js';
import {PropertyGrid} from '../src/components/forms/property-grid.js';
import {propertyEditor} from '../src/dg/forms/property-editor.js';

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

const DRUGS = [
  {id: 'a', name: 'Aspirin', mw: 180},
  {id: 'c', name: 'Caffeine', mw: 194},
  {id: 'i', name: 'Ibuprofen', mw: 206},
];

/** A table whose name column counts its renders. */
function drugTable(options = {}) {
  const calls = {render: 0};
  const t = new BasicTable({
    keyOf: (x) => x.id,
    columns: [
      {header: 'Name', render: (x) => {
        calls.render++;
        return x.name;
      }},
      {header: 'MW', render: (x) => String(x.mw), align: 'right'},
    ],
    ...options,
  });
  document.body.append(t.root);
  return {t, calls};
}

const tableRows = (t) => t.root.querySelectorAll('tbody tr');
const rowOf = (t, key) => tableRows(t).find((r) => r.dataset.u2Item === key);

scoped('BasicTable: one item replaced re-renders that row only, the others keep their elements', () => {
  const items = signal(DRUGS);
  const {t, calls} = drugTable({items, selectable: true});
  assert.equal(calls.render, 3);
  assert.deepEqual(tableRows(t).map((r) => r.dataset.u2Item), ['a', 'c', 'i']);
  assert.deepEqual(tableRows(t).map((r) => r.getAttribute('role')), ['row', 'row', 'row']);
  const aspirin = rowOf(t, 'a');
  const ibuprofen = rowOf(t, 'i');

  calls.render = 0;
  items.value = [DRUGS[0], {...DRUGS[1], name: 'Caffeine citrate'}, DRUGS[2]];
  assert.equal(calls.render, 1, 'only the replaced item rendered');
  assert.equal(rowOf(t, 'a'), aspirin);
  assert.equal(rowOf(t, 'i'), ibuprofen);
  assert.equal(rowOf(t, 'c').querySelector('td').textContent, 'Caffeine citrate');

  fire(rowOf(t, 'i'), 'click');
  assert.equal(t.selectedIndex.value, 2);
  calls.render = 0;
  items.value = [DRUGS[2], DRUGS[0]];
  assert.equal(calls.render, 0, 'a reorder renders nothing');
  assert.equal(t.selectedIndex.value, 0, 'the selection follows its key');
  assert.equal(rowOf(t, 'i'), ibuprofen);
  assert.equal(ibuprofen.getAttribute('aria-selected'), 'true');
  assert.deepEqual(tableRows(t).map((r) => r.dataset.index), ['0', '1']);
  assert.deepEqual(t.getState().selection, ['i']);
  t.dispose();
});

scoped('BasicTable: {section} items are spanning header rows that selection passes by', async () => {
  const clicks = [];
  const {t} = drugTable({selectable: true, onRowClick: (x, i) => clicks.push(`${x.name}:${i}`),
    items: [{section: 'Analgesics'}, DRUGS[0], DRUGS[2], {section: 'Stimulants'}, DRUGS[1]]});
  const section = rowOf(t, 'section:Analgesics');
  assert.equal(section.dataset.u2Part, 'section');
  assert.equal(section.classList.contains('u2-table-row'), false);
  const th = section.querySelector('th');
  assert.equal(th.colSpan, 2);
  assert.equal(th.getAttribute('role'), 'rowheader');
  assert.equal(th.textContent, 'Analgesics');

  fire(th, 'click');
  assert.equal(t.selectedIndex.value, -1, 'a section is not selectable');
  assert.deepEqual(clicks, []);

  fire(t.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(t.selectedIndex.value, 1, 'the first row, past the section');
  fire(t.root, 'keydown', {key: 'ArrowDown'});
  fire(t.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(t.selectedIndex.value, 4, 'the second section is skipped');
  fire(t.root, 'keydown', {key: 'Enter'});
  assert.deepEqual(clicks, ['Caffeine:4']);
  fire(t.root, 'keydown', {key: 'Home'});
  assert.equal(t.selectedIndex.value, 1);

  const items = t.getItems();
  assert.equal(items.count, 3, 'sections are not items');
  assert.deepEqual(items.list.map((x) => [x.key, x.label, x.index, x.selected]),
    [['a', 'Aspirin', 0, true], ['i', 'Ibuprofen', 1, false], ['c', 'Caffeine', 2, false]]);
  assert.equal(await t.revealItem('Analgesics'), null);
  assert.equal(await t.revealItem('ibuprofen'), rowOf(t, 'i'));
  t.dispose();
});

scoped('BasicTable: the row hook names, classes and titles a row, in the row\'s scope', () => {
  const flagged = signal(new Set(['c']));
  const items = signal(DRUGS);
  const {t} = drugTable({items, row: (x, i, tr) => {
    tr.dataset.u2Name = `drug-${x.id}`;
    tr.title = `${x.name} (#${i})`;
    Scope.ambient.effect(() => tr.classList.toggle('flagged', flagged.value.has(x.id)));
  }});
  const caffeine = rowOf(t, 'c');
  assert.equal(caffeine.dataset.u2Name, 'drug-c');
  assert.equal(caffeine.title, 'Caffeine (#1)');
  assert.equal(caffeine.classList.contains('flagged'), true);
  flagged.value = new Set(['a']);
  assert.equal(caffeine.classList.contains('flagged'), false);
  assert.equal(rowOf(t, 'a').classList.contains('flagged'), true);

  const live = Scope.liveCount;
  items.value = [DRUGS[0], DRUGS[2]];
  assert.equal(Scope.liveCount, live - 1, 'the dropped row\'s scope is disposed');
  flagged.value = new Set(['c']);
  assert.equal(caffeine.classList.contains('flagged'), false, 'and its binding with it');
  t.dispose();
});

scoped('BasicTable: items without an identity are keyed by position (the import wizard\'s issues)', () => {
  const issues = [{row: 3, message: 'bad'}, {row: 7, message: 'worse'}];
  const t = new BasicTable({items: issues, columns: [{header: 'Row', render: (x) => String(x.row)},
    {header: 'Message', render: (x) => x.message}]});
  document.body.append(t.root);
  assert.deepEqual(tableRows(t).map((r) => r.dataset.u2Item), ['#0', '#1']);
  assert.deepEqual(t.getItems().list.map((x) => x.label), ['3', '7']);
  t.setItems([issues[1]]);
  assert.deepEqual(tableRows(t).map((r) => r.textContent), ['7worse']);
  t.dispose();
});

/** 10,000 items, 20 px rows under a 20 px header in a 420 px scroller: 20 visible + 3 = 23. */
function bigTable(options = {}, source = undefined) {
  const calls = {render: 0};
  const items = Array.from({length: 10000}, (_, i) => ({id: `r${i}`, name: `Item ${i}`, mw: i}));
  const t = new DataTable({rowHeight: 20, keyOf: (x) => x.id, columns: [
    {name: 'name', render: (x) => {
      calls.render++;
      return x.name;
    }},
    {name: 'mw'},
  ], ...options});
  document.body.append(t.root);
  t.root.clientHeight = 420;
  t.setItems(source ?? items);
  return {t, items, calls};
}

const dataRow = (t, i) => t.root.querySelector(`.u2-data-table-row[data-index="${i}"]`);

scoped('DataTable: an array edit re-renders one row (was the window); kept rows follow their index', () => {
  const {t, items, calls} = bigTable();
  assert.equal(t.renderedCount, 23);
  const row3 = dataRow(t, 3);
  calls.render = 0;
  const edited = items.slice();
  edited[5] = {...edited[5], name: 'Renamed'};
  t.setItems(edited);
  assert.equal(calls.render, 1, 'one row of 23');
  assert.equal(dataRow(t, 5).children[0].textContent, 'Renamed');

  calls.render = 0;
  t.setItems([{id: 'new', name: 'First', mw: -1}, ...edited]);
  assert.equal(calls.render, 1, 'an insert at the top renders the new row only');
  assert.equal(dataRow(t, 4), row3, 'the kept row moved down');
  assert.equal(row3.getAttribute('aria-rowindex'), '6', 'and its aria-rowindex with it (header is 1)');
  assert.equal(row3.dataset.u2Item, 'r3');
  assert.equal(row3.children[0].dataset.u2Part, 'item-label');
  t.dispose();
});

scoped('DataTable: a RowsLike source is live by default; liveItems: false keeps unchanged rows', () => {
  const live = Array.from({length: 50}, (_, i) => ({id: `r${i}`, name: `Item ${i}`, mw: i}));
  const items = signal(live);
  const rowsLike = {items, keyOf: (x) => x.id, byKey: (k) => live.find((x) => x.id === k)};
  const {t, calls} = bigTable({}, rowsLike);
  calls.render = 0;
  live[2].name = 'Edited in place';
  items.value = live.slice();
  assert.equal(calls.render, t.renderedCount, 'live: the window renders again');
  assert.equal(dataRow(t, 2).children[0].textContent, 'Edited in place');
  t.dispose();

  const {t: still, calls: stillCalls} = bigTable({liveItems: false}, rowsLike);
  stillCalls.render = 0;
  items.value = live.slice();
  assert.equal(stillCalls.render, 0, 'not live: the same objects render nothing');
  still.refresh('r2');
  assert.equal(stillCalls.render, 1, 'refresh(key) renders that row');
  still.dispose();
});

scoped('DataTable: an item\'s label is its first column\'s value, rendered or not', async () => {
  const {t} = bigTable();
  const items = t.getItems();
  assert.equal(items.count, 10000);
  assert.equal(items.list[3].label, 'Item 3');
  assert.equal(items.list[400].label, 'Item 400');
  assert.equal(items.list[400].rendered, false);
  const row = await t.revealItem('item 900');
  assert.equal(row, dataRow(t, 900));
  t.dispose();
});

/** 10 columns × 5 rows of 20 px cells over string items. */
function grid(items, options = {}) {
  const calls = {render: 0};
  const g = new VirtualGrid({cellWidth: 20, cellHeight: 20, render: (x, _i, cell) => {
    calls.render++;
    cell.title = x;
    return span(x);
  }, ...options});
  document.body.append(g.root);
  g.root.clientWidth = 200;
  g.root.clientHeight = 100;
  fire(g.root, 'scroll');
  g.setItems(items);
  return {g, calls};
}

const ICONS = Array.from({length: 500}, (_, i) => `icon-${i}`);
const cellOf = (g, key) => g.root.querySelectorAll('.u2-grid-cell').find((c) => c.dataset.u2Item === key);

scoped('VirtualGrid: cells are keyed and addressable; narrowing keeps the cells that stay', async () => {
  const {g, calls} = grid(ICONS);
  assert.equal(g.renderedCount, 70);
  const kept = cellOf(g, 'icon-5');
  assert.equal(kept.getAttribute('role'), 'option');
  fire(kept, 'click');
  assert.equal(g.selectedIndex.value, 5);

  calls.render = 0;
  g.setItems(ICONS.filter((x) => x.endsWith('5')));
  assert.equal(cellOf(g, 'icon-5'), kept, 'the same element');
  assert.equal(kept.dataset.index, '0');
  assert.equal(g.selectedIndex.value, 0, 'the selection follows its key');
  assert.equal(calls.render, g.renderedCount - 7, 'only the cells new to the window rendered');

  g.setItems(ICONS);
  const items = g.getItems();
  assert.equal(items.count, 500);
  assert.equal(items.truncated, false);
  assert.deepEqual([items.list[5].key, items.list[5].label, items.list[5].selected], ['icon-5', 'icon-5', true]);
  assert.equal(items.list[400].rendered, false);
  assert.deepEqual(g.getState().selection, ['icon-5']);
  const far = await g.revealItem('icon-450');
  assert.equal(far, cellOf(g, 'icon-450'));
  assert.equal(far.title, 'icon-450');
  g.dispose();
});

scoped('VirtualGrid: a render binds in the cell\'s scope, which dies with the cell', () => {
  const tick = signal(0);
  const {g} = grid(ICONS.slice(0, 20), {render: (x) => span(computed(() => `${x}:${tick.value}`))});
  const cell = cellOf(g, 'icon-1');
  tick.value = 1;
  assert.equal(cell.textContent, 'icon-1:1');
  const live = Scope.liveCount;
  g.setItems(ICONS.slice(2, 20));
  assert.equal(Scope.liveCount, live - 2, 'the two dropped cells released their scopes');
  g.dispose();
});

function propertyGrid(count, options) {
  const props = Array.from({length: count}, (_, i) => ({name: `p${i}`, type: i % 2 ? 'string' : 'int'}));
  const values = Object.fromEntries(props.map((p, i) => [p.name, i % 2 ? `v${i}` : i]));
  const pg = new PropertyGrid(options);
  document.body.append(pg.root);
  pg.setProperties(props, values);
  return {pg, props, values};
}

const pgRow = (pg, name) => pg.root.querySelectorAll('.u2-propgrid-row').find((r) => r.dataset.u2Item === name);

scoped('PropertyGrid: one keystroke over 200 rows runs at most 2 row effects (was 201)', () => {
  const {pg} = propertyGrid(200);
  let runs = 0;
  const coerce = PropertyGrid._coerce;
  PropertyGrid._coerce = function (...args) {
    runs++;
    return coerce.apply(this, args);
  };
  try {
    const field = pgRow(pg, 'p1').querySelector('input');
    field.value = 'typed';
    fire(field, 'input');
    assert.equal(pg.values.peek().p1, 'typed');
    assert.ok(runs <= 2, `row effects run: ${runs}`);
  } finally {
    PropertyGrid._coerce = coerce;
  }
  pg.dispose();
});

scoped('PropertyGrid: rows are keyed by name — kept, rebuilt when the descriptor changes, items reported', () => {
  const {pg, props, values} = propertyGrid(6);
  const row1 = pgRow(pg, 'p1');
  const row2 = pgRow(pg, 'p2');
  assert.equal(row1.getAttribute('role'), 'row');
  assert.equal(pg.root.getAttribute('role'), 'grid');

  pg.setProperties(props.map((p) => p.name === 'p2' ? {...p, readonly: true} : {...p}), {...values, p1: 'outside'});
  assert.equal(pgRow(pg, 'p1'), row1, 'an equal descriptor keeps its row');
  assert.equal(row1.querySelector('input').value, 'outside', 'and shows the new value');
  assert.notEqual(pgRow(pg, 'p2'), row2, 'a changed descriptor rebuilds its row');
  assert.equal(pgRow(pg, 'p2').querySelector('.u2-propgrid-value').textContent, '2');

  const items = pg.getItems();
  assert.equal(items.count, 6);
  assert.deepEqual(items.list[2], {key: 'p2', label: 'P2', index: 2, disabled: true, rendered: true},
    'the label is the caption the row shows');
  pg.setProperties([{name: 'only', type: 'bool', category: 'Flags'}], {only: true});
  assert.deepEqual(pg.getItems().list.map((x) => x.key), ['only']);
  assert.equal(pg.root.querySelector('[data-u2-part="category"]').textContent, '▸Flags');
  pg.dispose();
});

scoped('PropertyGrid: an edit reaches the onChanged option and the change event, never a write from code', () => {
  const seen = [];
  const {pg} = propertyGrid(4, {onChanged: (change, grid) => seen.push(['option', change.name, change.value, grid])});
  pg.on('change', (change) => seen.push(['event', change.name, change.value]));
  const field = pgRow(pg, 'p3').querySelector('input');
  field.value = 'x';
  fire(field, 'input');
  field.value = 'y';
  fire(field, 'input');
  assert.deepEqual(seen, [['option', 'p3', 'x', pg], ['event', 'p3', 'x'], ['option', 'p3', 'y', pg], ['event', 'p3', 'y']]);
  pg.values.value = {...pg.values.peek(), p3: 'from code'};
  assert.equal(field.value, 'from code');
  assert.equal(seen.length, 4, 'a write from code is not an edit');
  pg.dispose();
});

scoped('PropertyEditor: its items are the fields of both forms, keyed by name', () => {
  const editor = propertyEditor({name: 'dose', type: 'double'});
  document.body.append(editor.root);
  const items = editor.getItems();
  const keys = items.list.map((x) => x.key);
  assert.equal(keys[0], 'name');
  assert.ok(keys.includes('min') && keys.includes('showSlider'), keys.join());
  assert.equal(items.list.find((x) => x.key === 'semType').label, 'Semantic type');
  assert.deepEqual(items.list.map((x) => x.index), keys.map((_, i) => i));
  editor.dispose();
});

function funcs(count) {
  return Array.from({length: count}, (_, i) => {
    const name = `Pkg:f${String(i).padStart(3, '0')}`;
    const label = `Function ${i}`;
    return {name, label, tags: [], roles: [], signature: '(x) : int', search: `${name} ${label}`.toLowerCase()};
  });
}

function browser(options) {
  const fb = new FunctionsBrowser(options);
  document.body.append(fb.root);
  const list = fb.root.querySelector('.u2-list');
  list.clientHeight = 280;
  fire(list, 'scroll');
  return fb;
}

scoped('FunctionsBrowser: items are the shown functions, labelled by their friendly name', async () => {
  const fb = browser({items: funcs(300)});
  fb.query.value = 'pkg:f1';
  const items = fb.getItems();
  assert.equal(items.count, 100);
  assert.deepEqual([items.list[0].key, items.list[0].label], ['Pkg:f100', 'Function 100']);
  const row = await fb.revealItem('Function 199');
  assert.equal(row.dataset.u2Item, 'Pkg:f199');
  fb.dispose();
});

scoped('FunctionsBrowser: rows the selection effect scrolls into view outlive that effect\'s next run', () => {
  const tick = signal(0);
  const fb = browser({items: funcs(100),
    render: (item) => span(computed(() => `${item.label}/${tick.value}`))});
  const list = fb.root.querySelector('.u2-list');
  const row = (key) => list.querySelectorAll('.u2-list-row').find((r) => r.dataset.u2Item === key);
  // select the last function, then change the query so the list keeps it and the browser's
  // selection effect scrolls to it: its row is rendered INSIDE that effect run
  fb._list.selectedIndex.value = 99;
  assert.equal(row('Pkg:f099'), undefined, 'not rendered yet');
  fb.query.value = 'pkg';
  const rendered = row('Pkg:f099');
  assert.ok(rendered !== undefined, 'scrolled into view by the selection effect');
  assert.equal(rendered.textContent, 'Function 99/0');

  fb.query.value = 'pkg:';
  assert.equal(row('Pkg:f099'), rendered, 'the effect ran again; the row is kept');
  tick.value = 1;
  assert.equal(rendered.textContent, 'Function 99/1', 'its binding is alive');
  fb.dispose();
});
