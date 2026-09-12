/* RowsLike (WO-5): `arrayRows` over an array or a signal, and `frameRows` over the DataFrame
   double — id-keyed live proxies, items following values / filter / rows-added events minus the
   rows marked deleted, writes to the frame or to the `onWrite` hook, repoint without leaks. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {BitArray} from 'datagrok-api/u2core';
import {arrayRows} from '../src/sources/rows-like.js';
import {frameRows, FrameRows} from '../src/sources/df-rows.js';
import {Rows} from '../src/sources/rows-like.js';
import {BitSet, DataFrame} from './platform-doubles.mjs';

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

function issues() {
  return new DataFrame(
    [{name: 'id', type: 'string'}, {name: 'title', type: 'string'}, {name: '~state', type: 'string'}],
    [{id: 'a', title: 'Aspirin', '~state': ''}, {id: 'b', title: 'Ibuprofen', '~state': ''},
      {id: 'c', title: 'Naproxen', '~state': ''}]);
}

scoped('arrayRows: an array or a signal, keyed lookups follow the signal', () => {
  const plain = arrayRows([{id: 1, n: 'a'}, {id: 2, n: 'b'}], (x) => String(x.id));
  assert.deepEqual(plain.items.value.map((x) => x.n), ['a', 'b']);
  assert.equal(plain.keyOf(plain.items.value[1]), '2');
  assert.equal(plain.byKey('2').n, 'b');
  assert.equal(plain.byKey('9'), undefined);

  const live = signal([{id: 1, n: 'a'}]);
  const rows = arrayRows(live, (x) => String(x.id));
  assert.equal(rows.items, live, 'the signal is the items');
  assert.equal(rows.byKey('1').n, 'a');
  live.value = [{id: 3, n: 'c'}];
  assert.equal(rows.byKey('1'), undefined);
  assert.equal(rows.byKey('3').n, 'c');
});

scoped('frameRows: one live proxy per id — reads follow the frame, writes reach it', () => {
  const scope = new Scope();
  const df = issues();
  const rows = frameRows(signal(df), scope);
  assert.deepEqual(rows.items.value.map((r) => r.id), ['a', 'b', 'c']);
  const b = rows.byKey('b');
  assert.equal(b.title, 'Ibuprofen');
  assert.equal(rows.items.value[1], b, 'the same proxy everywhere');
  assert.equal(rows.keyOf(b), 'b');
  assert.equal('title' in b, true);
  assert.deepEqual(Object.keys(b), ['id', 'title'], 'service columns are never enumerated (H7)');
  assert.equal('~state' in b, false);
  assert.equal(b['~state'], '', 'but read by name');
  assert.equal(b.nope, undefined, 'an unknown column reads undefined');

  df.set('title', 1, 'Ibuprofen 400');
  assert.equal(b.title, 'Ibuprofen 400', 'live: no copy to refresh');
  b.title = 'Ibu';
  assert.deepEqual(df.dart.writes.at(-1), ['title', 1, 'Ibu']);
  assert.equal(rows.byKey('z'), undefined);
  assert.equal(rows.indexOf('c'), 2);
  scope.dispose();
  assert.equal(df.liveSubscriptions(), 0, 'nothing left on the frame');
});

scoped('frameRows: items skip filtered-out rows, keep deleted ones, and follow rows added', () => {
  const scope = new Scope();
  const df = issues();
  const rows = frameRows(signal(df), scope);
  let versions = 0;
  scope.effect(() => {
    rows.items.value;
    versions++;
  });

  const a = rows.byKey('a');
  df.set('~state', 0, 'deleted');
  assert.deepEqual(rows.items.value.map((r) => r.id), ['a', 'b', 'c'],
    'a row marked deleted stays an item until the save, as the platform grid keeps it');
  assert.equal(rows.byKey('a')['~state'], 'deleted', 'its state says so');
  assert.equal(a.title, 'Aspirin', 'a proxy already held still reads');

  df.dart.filter = BitSet.fromBitArray(BitArray.create(3, (i) => i === 2));
  df.onFilterChanged.fire();
  assert.deepEqual(rows.items.value.map((r) => r.id), ['a', 'c'],
    'the frame\'s filter applies — except to the deleted row the editor masked out of it');
  df.dart.filter = BitSet.fromBitArray(BitArray.create(3, (i) => i !== 1));
  df.onFilterChanged.fire();

  df.dart.rows.push({id: 'd', title: 'Diclofenac', '~state': 'new'});
  df.dart.filter = BitSet.fromBitArray(BitArray.create(4, (i) => i !== 1));
  df.onRowsAdded.fire();
  assert.deepEqual(rows.items.value.map((r) => r.id), ['a', 'c', 'd']);
  assert.equal(rows.byKey('d').title, 'Diclofenac');
  assert.equal(versions, 5, 'one rebuild per event');
  scope.dispose();
});

scoped('frameRows: onWrite takes the writes instead of the frame; an id-less row is keyed by index', () => {
  const scope = new Scope();
  const df = new DataFrame([{name: 'id', type: 'string'}, {name: 'title', type: 'string'}],
    [{id: 'a', title: 'x'}, {id: null, title: 'draft'}]);
  const writes = [];
  const rows = frameRows(signal(df), scope, {onWrite: (id, column, value) => writes.push([id, column, value])});
  assert.deepEqual(rows.items.value.map((r) => r.id), ['a', '~row:1']);
  rows.byKey('a').title = 'y';
  assert.deepEqual(writes, [['a', 'title', 'y']]);
  assert.deepEqual(df.dart.writes, [], 'the frame was not written');
  assert.equal(rows.byKey('~row:1').title, 'draft');
  scope.dispose();
});

scoped('frameRows: a repoint resubscribes and rebuilds; no frame is no rows', () => {
  const scope = new Scope();
  const first = issues();
  const source = signal(first);
  const rows = frameRows(source, scope);
  assert.equal(first.liveSubscriptions(), 5);

  const second = new DataFrame([{name: 'id', type: 'string'}, {name: 'title', type: 'string'}],
    [{id: 'z', title: 'Zinc'}]);
  source.value = second;
  assert.equal(first.liveSubscriptions(), 0, 'let go of the old frame');
  assert.equal(second.liveSubscriptions(), 5);
  assert.deepEqual(rows.items.value.map((r) => r.id), ['z']);
  assert.equal(rows.byKey('a'), undefined);

  source.value = undefined;
  assert.deepEqual(rows.items.value, []);
  assert.equal(second.liveSubscriptions(), 0);
  scope.dispose();
});

scoped('Rows: the one holder of the ~ conventions', () => {
  assert.equal(Rows.STATE, '~state');
  assert.equal(Rows.isService('~can_edit'), true);
  assert.equal(Rows.isService('title'), false);
  assert.equal(Rows.draftKey(3), '~row:3');
  assert.equal(Rows.draftIndex('~row:3'), 3);
  assert.equal(Rows.isDraft('~row:3'), true);
  assert.equal(Rows.isDraft({id: '~row:0'}), true);
  assert.equal(Rows.isDraft({id: 'i1'}), false);
});

scoped('FrameRows.proxy: the row shape over any reader', () => {
  const store = {title: 'x', count: 1};
  const row = FrameRows.proxy(() => 'k', {get: (c) => store[c], set: (c, v) => store[c] = v,
    keys: () => Object.keys(store)});
  assert.equal(row.id, 'k');
  assert.equal(row.title, 'x');
  row.count = 2;
  assert.equal(store.count, 2);
  row.id = 'other';
  assert.equal(row.id, 'k', 'the key is not writable through the row');
  assert.deepEqual({...row}, {id: 'k', title: 'x', count: 2}, 'a spread copies the columns');
});
