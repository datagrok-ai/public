/* The two shell helpers lifted out of EMS (`collectionSelection`, `under`), driven by u2's own
   primitives — a `FrameRows` over a `MemoryFrame` — so their contract holds with no domain source. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush} from './dom-shim.js';
import {layout} from './layout.mjs';
import {signal} from '../src/core/signals.js';
import {Control} from '../src/core/component.js';
import {DataTable} from '../src/components/collections/data-table.js';
import {FrameRows} from '../src/sources/df-rows.js';
import {MemoryFrame} from '../src/sources/memory-frame.js';
import {collectionSelection} from '../src/dg/shell/collection-selection.js';
import {u2test, mount, user} from '../src/testing/index.js';

register('./dg-stub.mjs', import.meta.url);
const {under} = await import('../src/dg/shell/view-sync.js');

const COLUMNS = [{name: 'id', type: 'string'}, {name: 'title', type: 'string'}];
const frameOf = (ids) => new MemoryFrame(COLUMNS, ids.map((id) => ({id, title: id.toUpperCase()})));

async function bound() {
  const owner = mount(() => new Control());
  const df = signal(frameOf(['a', 'b', 'c']));
  const source = {currentRow: signal(null), rows: new FrameRows(df, {owner: owner.scope}), df};
  const table = owner.build(() => new DataTable({rowHeight: 20, items: source.rows, columns: [{name: 'title'}]}));
  collectionSelection(source, table, {owner: owner});
  owner.root.append(table.root);
  layout(table.root, {clientHeight: 300});
  fire(table.root, 'scroll');
  await flush();
  return {source, owner, table};
}

const marked = (df) => df.rows.map((r, i) => df.selection.get(i) ? r.id : null).filter((id) => id !== null);
const rowAt = (table, i) => table.root.querySelector(`.u2-data-table-row[data-index="${i}"]`);

test('collectionSelection: the lead and currentRow are one value; only a gesture marks the frame', u2test(async () => {
  const {source, owner, table} = await bound();
  source.currentRow.value = source.rows.byKey('b');
  await flush();
  assert.equal(table.selectedIndex.value, 1);
  assert.deepEqual(marked(source.df.value), [], 'a lead nobody picked is not a selection');

  user.click(rowAt(table, 0));
  user.click(rowAt(table, 2), {ctrlKey: true});
  await flush();
  assert.deepEqual(marked(source.df.value), ['a', 'c']);

  user.press(table.root, 'Escape');
  await flush();
  assert.deepEqual(marked(source.df.value), []);
  assert.equal(source.currentRow.value, null);
  owner.dispose();
}));

test('collectionSelection: a new frame keeps the selection by key, one without the keys clears it', u2test(async () => {
  const {source, owner, table} = await bound();
  user.click(rowAt(table, 0));
  user.click(rowAt(table, 2), {ctrlKey: true});
  await flush();

  source.df.value = frameOf(['c', 'b', 'a']);
  await flush();
  assert.deepEqual(marked(source.df.value), ['c', 'a']);

  source.df.value = frameOf(['x', 'y']);
  await flush();
  assert.deepEqual(marked(source.df.value), []);
  owner.dispose();
}));

test('under: a base continues at a segment or a query, never mid-segment', () => {
  assert.equal(under('/apps/grit', '/apps/grit'), true);
  assert.equal(under('/apps/grit?entity=i1', '/apps/grit'), true);
  assert.equal(under('/apps/grit/i1', '/apps/grit'), true);
  assert.equal(under('/apps/grit_labs', '/apps/grit'), false);
});
