// where a UTC midnight read locally is the day before
process.env.TZ = 'America/New_York';

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {Scope, Format, backends} from '../src/index.js';
import {CatalogTable} from '../src/components/collections/catalog-table.js';
import {DataTable} from '../src/components/collections/data-table.js';
import {tableFromMap} from '../src/components/collections/table.js';

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      delete backends.format;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

scoped('value: empty for null, undefined and NaN; booleans true/false; numbers and strings as they are', () => {
  assert.deepEqual([null, undefined, NaN, ''].map((v) => Format.value(v)), ['', '', '', '']);
  assert.deepEqual([true, false].map((v) => Format.value(v)), ['true', 'false'], 'never Yes/No');
  assert.deepEqual([0, 46, -1.5, 1234567].map((v) => Format.value(v)), ['0', '46', '-1.5', '1234567']);
  assert.equal(Format.value('KHL-12'), 'KHL-12');
  assert.equal(Format.value(12n), '12');
});

scoped('value: a list item by item; a record of no known shape is empty, never [object Object]', () => {
  assert.equal(Format.value(['a', null, 2, true]), 'a, 2, true');
  assert.equal(Format.value({login: 'ann'}), '');
  assert.equal(Format.value(Object.create(null)), '');
  assert.equal(Format.value({toString: () => 'x'}), '', 'a record is a record, whatever it carries');
  class Entity {
    toString() {
      return 'Ann Lee';
    }
  }
  assert.equal(Format.value(new Entity()), 'Ann Lee', 'an object of a class says what it is');
});

scoped('date: the platform default M/d/yyyy, with the finest time part the value carries', () => {
  assert.equal(Format.value(new Date(2026, 5, 29)), '6/29/2026');
  assert.equal(Format.value(new Date(2026, 5, 29, 14, 5)), '6/29/2026 14:05');
  assert.equal(Format.value(new Date(2026, 5, 29, 0, 0, 7)), '6/29/2026 00:00:07');
  assert.equal(Format.value(new Date(2026, 5, 29, 9, 3, 7, 42)), '6/29/2026 09:03:07.042');
  assert.equal(Format.value({toDate: () => new Date(2026, 0, 2)}), '1/2/2026', 'a dayjs reads as its date');
  assert.equal(Format.value(new Date('not a date')), '');
});

scoped('date: every platform token; any other run of letters prints as it is', () => {
  const d = new Date(2026, 2, 5, 13, 7, 9, 345);
  assert.equal(Format.date(d, {format: 'yyyy-MM-dd HH:mm:ss.fff'}), '2026-03-05 13:07:09.345');
  assert.equal(Format.date(d, {format: 'yy M/d H:m:s f ff'}), '26 3/5 13:7:9 3 34');
  assert.equal(Format.date(d, {format: 'dddd, MMMM d, yyyy'}), 'Thursday, March 5, 2026');
  assert.equal(Format.date(d, {format: 'ddd MMM d h:mm tt'}), 'Thu Mar 5 1:07 PM');
  assert.equal(Format.date(new Date(2026, 2, 5, 0, 30), {format: 'hh:mm tt'}), '12:30 AM');
  assert.equal(Format.date(new Date(2026, 2, 5, 12, 30), {format: 'h:mm tt'}), '12:30 PM');
  assert.equal(Format.date(d, {format: 'yyyy-MM-ddTHH:mm'}), '2026-03-05T13:07');
});

scoped('date: utcDates reads a midnight-UTC date in UTC; every other value reads in the local zone', () => {
  const day = new Date('2026-10-02T00:00:00Z');
  assert.equal(Format.value(day), '10/1/2026 20:00', 'locally it is the evening before');
  assert.equal(Format.value(day, {utcDates: true}), '10/2/2026');
  const moment = new Date('2026-10-02T15:30:00Z');
  assert.equal(Format.value(moment, {utcDates: true}), '10/2/2026 11:30', 'a moment stays local');
});

scoped('the platform backend: the user\'s date format and the platform\'s number formats', () => {
  const asked = [];
  backends.format = {dateFormat: () => 'yyyy/M/d', number: (x, format) => {
    asked.push([x, format]);
    return `${format}:${x}`;
  }};
  assert.equal(Format.value(new Date(2026, 5, 29)), '2026/6/29');
  assert.equal(Format.value(new Date(2026, 5, 29, 8, 15)), '2026/6/29 08:15');
  assert.equal(Format.value(new Date(2026, 5, 29), {format: 'dd.MM.yyyy'}), '29.06.2026', 'a column\'s own wins');
  assert.equal(Format.value(3.14159, {format: '0.00'}), '0.00:3.14159');
  assert.equal(Format.value(3.14159), '3.14159', 'no format, no platform call');
  assert.deepEqual(asked, [[3.14159, '0.00']]);
  backends.format = {dateFormat: () => null, number: (x) => String(x)};
  assert.equal(Format.value(new Date(2026, 5, 29)), '6/29/2026', 'no date format chosen: the platform default');
  delete backends.format;
  assert.equal(Format.value(3.14159, {format: '0.00'}), '3.14159', 'headless, a number format has no one to apply it');
});

scoped('count groups the thousands; spelled is the locale\'s short date, with the time on request', () => {
  assert.equal(Format.count(1234567), (1234567).toLocaleString());
  const d = new Date(2026, 9, 2, 14, 30);
  const day = {month: 'short', day: 'numeric', year: 'numeric'};
  assert.equal(Format.spelled(d), d.toLocaleDateString(undefined, day));
  assert.equal(Format.spelled(d, {time: true}), d.toLocaleString(undefined, {...day, hour: '2-digit', minute: '2-digit'}));
  const utc = new Date('2026-10-02T00:00:00Z');
  assert.equal(Format.spelled(utc, {utc: true}), utc.toLocaleDateString(undefined, {...day, timeZone: 'UTC'}));
});

// ---- the catalog contract ---------------------------------------------------------------------------

const ITEMS = [
  {code: 'KHL-1', live: true, day: new Date(2026, 5, 29), moment: new Date(2026, 5, 29, 9, 30), count: 8, kind: 'Kinase'},
  {code: 'KHL-2', live: false, day: new Date(2026, 5, 29), moment: new Date(2026, 5, 29, 9, 30), count: 12, kind: 'GPCR'},
  {code: 'KHL-3', live: true, day: new Date(2026, 6, 1), moment: new Date(2026, 6, 1, 17, 0), count: 8, kind: 'Kinase'},
  {code: 'KHL-4', live: false, day: new Date(2026, 9, 2), moment: new Date(2026, 9, 2, 0, 5), count: 1.5,
    kind: 'GPCR'},
];
const FIELDS = ['live', 'day', 'moment', 'count', 'kind'];

function catalog() {
  const t = new CatalogTable({name: 'campaigns', items: ITEMS, keyOf: (c) => c.code,
    columns: [{name: 'code', value: (c) => c.code, pinned: true}, ...FIELDS.map((f) => ({name: f, value: (c) => c[f]}))],
    groupings: FIELDS.map((f) => ({name: f, key: (c) => c[f]}))});
  document.body.append(t.root);
  return t;
}

/** Each group header's label with the texts its rows show in the column of the same name. */
async function groups(t, field) {
  t.groupBy(field);
  await flush();
  const out = [];
  for (const item of t.getItems({all: true}).list) {
    if (item.key.startsWith('section:'))
      out.push({label: item.label, cells: []});
    else {
      const tr = [...t.root.querySelectorAll('tr')].find((r) => r.dataset.u2Item === item.key);
      out.at(-1).cells.push([...tr.children].find((td) => td.dataset.column === field).textContent);
    }
  }
  return out;
}

async function headersReadAsCells(t, field, count) {
  const found = await groups(t, field);
  assert.equal(found.length, count, `${field}: ${found.map((g) => g.label).join(', ')}`);
  for (const g of found) {
    assert.ok(g.cells.length > 0);
    for (const cell of g.cells)
      assert.equal(cell, g.label, `the ${field} header reads "${g.label}" over a row that reads "${cell}"`);
  }
  return found.map((g) => g.label);
}

scoped('catalog: for bool, date, number and string a group header reads as the cells of its rows', async () => {
  const t = catalog();
  assert.deepEqual(await headersReadAsCells(t, 'live', 2), ['true', 'false'], 'booleans read true/false');
  assert.deepEqual(await headersReadAsCells(t, 'day', 3), ['6/29/2026', '7/1/2026', '10/2/2026']);
  assert.deepEqual(await headersReadAsCells(t, 'moment', 3), ['6/29/2026 09:30', '7/1/2026 17:00', '10/2/2026 00:05']);
  assert.deepEqual(await headersReadAsCells(t, 'count', 3), ['8', '12', '1.5']);
  assert.deepEqual(await headersReadAsCells(t, 'kind', 2), ['Kinase', 'GPCR']);
  t.dispose();
});

scoped('catalog: with the platform\'s date format, the date headers and cells both read in it', async () => {
  backends.format = {dateFormat: () => 'yyyy/M/d', number: (x) => String(x)};
  const t = catalog();
  assert.deepEqual(await headersReadAsCells(t, 'day', 3), ['2026/6/29', '2026/7/1', '2026/10/2']);
  t.dispose();
});

scoped('a DataTable cell and a tableFromMap value with no render read as a CatalogTable cell does', async () => {
  const values = {moment: new Date(2026, 5, 29, 9, 30), live: true, none: null, owner: {login: 'ann'}};
  const names = Object.keys(values);
  const item = {code: 'KHL-1', ...values};
  const catalog = new CatalogTable({items: [item], keyOf: (c) => c.code,
    columns: names.map((name) => ({name, value: (c) => c[name]}))});
  const data = new DataTable({columns: names.map((name) => ({name})), rowHeight: 20, keyOf: (c) => c.code});
  document.body.append(catalog.root, data.root);
  layout(data.root, {clientHeight: 220});
  data.setItems([item]);
  await flush();
  const tr = [...catalog.root.querySelectorAll('tr')].find((r) => r.dataset.u2Item === 'KHL-1');
  const cells = names.map((name) => [...tr.children].find((td) => td.dataset.column === name).textContent);
  assert.deepEqual(cells, ['6/29/2026 09:30', 'true', '', '']);
  assert.deepEqual([...data.root.querySelector('.u2-data-table-row[data-index="0"]').children].map((c) => c.textContent),
    cells, 'DataTable');
  assert.deepEqual([...tableFromMap(values).querySelectorAll('tr')].map((r) => r.children[1].textContent), cells,
    'tableFromMap');
  catalog.dispose();
  data.dispose();
});

scoped('catalog: a key with no text — a record — goes to the empty group, never "[object Object]"', async () => {
  const t = new CatalogTable({items: ITEMS, keyOf: (c) => c.code, columns: [{name: 'code', value: (c) => c.code}],
    groupings: [{name: 'owner', key: (c) => ({login: c.kind})}], singleGroupHeader: true});
  document.body.append(t.root);
  t.groupBy('owner');
  await flush();
  assert.deepEqual(t.getItems({all: true}).list.filter((i) => i.key.startsWith('section:')).map((i) => i.label),
    ['No information']);
  t.dispose();
});
