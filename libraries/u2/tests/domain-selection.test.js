/* `DomainSelection.bind` (3-4) — the one binding between a collection control and its source,
   here over a `DataTable`, which never had the multi-selection mirror: the lead and `currentRow`
   are one value in either direction, the mark is armed by a gesture and mirrored into the frame's
   selection bitset, a collection read again keeps it BY KEY, and Escape clears both. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {DataTable} from '../src/components/collections/data-table.js';
import {backends} from '../src/sources/backends.js';
import {backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {DomainSelection} = await import('../src/dg/domain/selection.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** A bare data table over the source's rows, bound and mounted with a viewport. */
async function bound(options = {}) {
  backends.domain = backend();
  const handle = await domains.table('grit.issue');
  const source = handle.source(options);
  await flush();
  const owner = new Control();
  const table = owner.runInScope(() => new DataTable({rowHeight: 20, items: source.rows,
    columns: [{name: 'title'}]}));
  DomainSelection.bind(owner, source, table);
  owner.root.append(table.root);
  document.body.append(owner.root);
  table.root.clientHeight = 300;
  fire(table.root, 'scroll');
  return {source, owner, table};
}

const ids = (source) => source.selection.value.map((r) => r.id);
const rowAt = (table, i) => table.root.querySelector(`.u2-data-table-row[data-index="${i}"]`);

scoped('the lead and the current row are one value, and a load\'s lead is not a selection', async () => {
  const {source, owner, table} = await bound();
  assert.equal(table.selectedIndex.value, -1);
  source.currentRow.value = source.rows.byKey('i1');
  await flush();
  assert.equal(table.selectedIndex.value, 0, 'the current row leads');
  assert.deepEqual(ids(source), [], 'a lead nobody picked is not a selection');

  table.selectedIndex.value = 2;
  await flush();
  assert.equal(source.currentRow.value.id, 'i3', 'and the other way round');
  owner.dispose();
  source.dispose();
});

scoped('a gesture arms the mirror: the frame\'s bitset is what the user picked', async () => {
  const {source, owner, table} = await bound();
  fire(rowAt(table, 0), 'click');
  await flush();
  assert.deepEqual(ids(source), ['i1'], 'a click is a selection');
  assert.equal(source.currentRow.value.id, 'i1');
  fire(rowAt(table, 2), 'click', {ctrlKey: true});
  await flush();
  assert.deepEqual(ids(source), ['i1', 'i3']);

  fire(table.root, 'keydown', {key: 'Escape'});
  await flush();
  assert.deepEqual(ids(source), [], 'Escape clears the mark');
  assert.equal(table.selectedIndex.value, -1, 'and the lead');
  assert.equal(source.currentRow.value, null);
  owner.dispose();
  source.dispose();
});

scoped('a collection read again keeps the selection by key; one holding none of the keys clears it', async () => {
  const {source, owner, table} = await bound();
  fire(rowAt(table, 0), 'click');
  fire(rowAt(table, 2), 'click', {ctrlKey: true});
  await flush();
  assert.deepEqual(ids(source), ['i1', 'i3']);

  await source.refresh();
  await flush();
  assert.deepEqual(ids(source), ['i1', 'i3'], 'the same rows, a new frame: the keys are the identity');

  source.query.value = 'title = "Ibuprofen"';
  await flush();
  assert.deepEqual(ids(source), [], 'a collection answering none of the kept keys selects nothing');
  fire(rowAt(table, 0), 'click');
  await flush();
  assert.deepEqual(ids(source), ['i2'], 'and the next gesture arms it again');
  owner.dispose();
  source.dispose();
});

scoped('the binding lives as long as the owner: disposing it releases the listeners', async () => {
  const {source, owner, table} = await bound();
  fire(rowAt(table, 0), 'click');
  await flush();
  assert.deepEqual(ids(source), ['i1']);
  const row = rowAt(table, 1);
  owner.dispose();
  fire(row, 'click');
  await flush();
  assert.deepEqual(ids(source), ['i1'], 'nothing follows a disposed binding');
  source.dispose();
});
