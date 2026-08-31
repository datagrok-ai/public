/* dfBindings over a fake DataFrameLike: fixed steps, lazy cached two-way column signals with
   echo suppression, live-object selection/filter signals that re-fire without changing identity,
   version-driven rowCount/columns, and the repoint contract (resubscribe + refresh + no leaks). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope, signal} from '../src/index.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {Column, DataFrame} from './platform-doubles.mjs';

function smoke(name, body) {
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

function orders() {
  return new DataFrame(
    [{name: 'name', type: 'string', semType: 'Text'},
      {name: 'total', type: 'double'},
      {name: 'Mol Weight', type: 'double', semType: 'Molecule Weight'}],
    [{'name': 'aspirin', 'total': 10, 'Mol Weight': 180.16},
      {'name': 'ibuprofen', 'total': 20, 'Mol Weight': 206.29}]);
}

smoke('fixed steps: df identity, rowCount, columns, and "" defaults to df', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);

  const dfSig = source.bindStep('df');
  assert.equal(dfSig.value, df);
  assert.equal(source.bindStep(''), dfSig, 'the default binding is the frame itself');
  assert.equal(source.bindStep('rowCount').value, 2);
  assert.deepEqual(source.bindStep('columns').value, ['name', 'total', 'Mol Weight']);
  assert.equal(source.bindStep('nope'), null);

  df.dart.rows.push({'name': 'naproxen', 'total': 30, 'Mol Weight': 230.26});
  df.onValuesChanged.fire({});
  assert.equal(source.bindStep('rowCount').value, 3);

  df.columns.add(new Column('flag', 'bool'));
  assert.deepEqual(source.bindStep('columns').value, ['name', 'total', 'Mol Weight', 'flag']);
  scope.dispose();
});

smoke('currentRowIdx: two-way with echo suppression on both sides', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);
  const idx = source.bindStep('currentRowIdx');
  assert.equal(idx.value, 0);

  const before = df.dart.idxWrites;
  idx.value = 1;
  assert.equal(df.currentRowIdx, 1);
  assert.equal(df.dart.idxWrites, before + 1, 'exactly one write reaches the frame');

  df.currentRowIdx = 0;
  assert.equal(idx.value, 0);
  assert.equal(df.dart.idxWrites, before + 2, 'a frame-side change does not echo a write back');
  scope.dispose();
});

smoke('currentRow.<col>: cached lazy signal, read/write/echo with two consumers', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);
  const row = source.bindStep('currentRow');

  const name = row.bindStep('name');
  assert.equal(row.bindStep('name'), name, 'the column signal is cached');
  assert.equal(name.value, 'aspirin');

  let firstRuns = 0;
  let secondRuns = 0;
  scope.effect(() => {
    name.value;
    firstRuns++;
  });
  scope.effect(() => {
    name.value;
    secondRuns++;
  });

  name.value = 'acetylsalicylic acid';
  await flush();
  assert.deepEqual(df.dart.writes, [['name', 0, 'acetylsalicylic acid']], 'one write, no oscillation');
  assert.equal(df.dart.rows[0].name, 'acetylsalicylic acid');
  assert.equal(firstRuns, 2);
  assert.equal(secondRuns, 2);

  df.dart.rows[0].name = 'external';
  df.onValuesChanged.fire({});
  assert.equal(name.value, 'external');
  assert.equal(df.dart.writes.length, 1, 'a frame-side edit does not echo a write back');

  df.currentRowIdx = 1;
  assert.equal(name.value, 'ibuprofen', 'a row change re-reads the cached signals');
  assert.equal(df.dart.writes.length, 1);
  scope.dispose();
});

smoke('bracket column names: a name with spaces walks and writes', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);
  const weight = source.bindStep('currentRow').bindStep('Mol Weight');
  assert.equal(weight.value, 180.16);
  weight.value = 181;
  assert.deepEqual(df.dart.writes, [['Mol Weight', 0, 181]]);
  scope.dispose();
});

smoke('currentRow: an unknown or frameless column reads undefined and refuses writes', async () => {
  const scope = new Scope();
  const df = orders();
  const dfSig = signal(undefined);
  const source = dfBindings(dfSig, scope);
  const row = source.bindStep('currentRow');
  assert.deepEqual(row.bindProps(), [], 'no frame, nothing to enumerate');

  const name = row.bindStep('name');
  assert.equal(name.value, undefined);
  name.value = 'ghost';
  assert.equal(source.bindStep('rowCount').value, 0);

  dfSig.value = df;
  assert.equal(name.value, 'aspirin', 'the frame arriving fills the pre-bound signal');
  assert.equal(df.dart.writes.length, 0, 'the frameless write never lands');

  const missing = row.bindStep('missing');
  assert.equal(missing.value, undefined);
  missing.value = 'nope';
  assert.equal(df.dart.writes.length, 0, 'an unknown column refuses writes');
  scope.dispose();
});

smoke('selection/filter: identity stable across version bumps while effects re-fire', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);
  const sel = source.bindStep('selection');
  const flt = source.bindStep('filter');
  assert.equal(sel.value, df.selection);
  assert.equal(flt.value, df.filter);

  let selRuns = 0;
  const seen = [];
  scope.effect(() => {
    seen.push(sel.value);
    selRuns++;
  });
  assert.equal(selRuns, 1);

  df.onSelectionChanged.fire({});
  assert.equal(selRuns, 2, 'the effect re-fires on a version bump');
  df.onSelectionChanged.fire({});
  assert.equal(selRuns, 3);
  assert.equal(seen.every((x) => x === df.selection), true, 'the identity never changed');

  let fltRuns = 0;
  scope.effect(() => {
    flt.value;
    fltRuns++;
  });
  df.onFilterChanged.fire({});
  assert.equal(fltRuns, 2);

  // the prototype-descriptor writability walk (core isWritableSignal), local so the fake
  // classification is pinned even while WO-11's helper is in flight
  const writable = (sig) => {
    for (let p = Object.getPrototypeOf(sig); p !== null; p = Object.getPrototypeOf(p)) {
      const d = Object.getOwnPropertyDescriptor(p, 'value');
      if (d)
        return d.set !== undefined;
    }
    return false;
  };
  assert.equal(writable(sel), false, 'the live-object signal classifies read-only');
  assert.equal(writable(flt), false);
  assert.equal(writable(source.bindStep('currentRowIdx')), true);
  scope.dispose();
});

smoke('repoint: resubscribes, refreshes every cached signal, drops the old subscriptions', async () => {
  const scope = new Scope();
  const df1 = orders();
  const df2 = new DataFrame(
    [{name: 'name', type: 'string'}, {name: 'city', type: 'string'}],
    [{name: 'zoe', city: 'Kyiv'}, {name: 'ada', city: 'Lviv'}, {name: 'bob', city: 'Riga'}]);
  df2.dart.currentRowIdx = 2;
  const dfSig = signal(df1);
  const source = dfBindings(dfSig, scope);

  const name = source.bindStep('currentRow').bindStep('name');
  const total = source.bindStep('currentRow').bindStep('total');
  const idx = source.bindStep('currentRowIdx');
  const sel = source.bindStep('selection');
  assert.equal(df1.liveSubscriptions(), 5);
  assert.equal(df2.liveSubscriptions(), 0);

  dfSig.value = df2;
  assert.equal(df1.liveSubscriptions(), 0, 'the old frame is fully unsubscribed');
  assert.equal(df2.liveSubscriptions(), 5);
  assert.equal(idx.value, 2);
  assert.equal(name.value, 'bob');
  assert.equal(total.value, undefined, 'a column the new frame lacks reads undefined');
  assert.equal(sel.value, df2.selection);
  assert.equal(source.bindStep('rowCount').value, 3);
  assert.deepEqual(source.bindStep('columns').value, ['name', 'city']);

  df1.onValuesChanged.fire({});
  df1.currentRowIdx = 1;
  assert.equal(name.value, 'bob', 'the abandoned frame no longer feeds the signals');

  name.value = 'boris';
  assert.deepEqual(df2.dart.writes, [['name', 2, 'boris']]);
  assert.equal(df1.dart.writes.length, 0);

  scope.dispose();
  assert.equal(df2.liveSubscriptions(), 0, 'disposal drops the live subscriptions');
});

smoke('laziness: bindProps enumerates columns without allocating signals or reading cells', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);

  const props = source.bindProps();
  assert.deepEqual(props.map((p) => p.name),
    ['df', 'currentRowIdx', 'currentRow', 'selection', 'filter', 'rowCount', 'columns']);
  assert.equal(props.find((p) => p.name === 'currentRow').walkable, true);
  assert.equal(props.find((p) => p.name === 'currentRowIdx').writable, true);

  const rowProps = source.bindStep('currentRow').bindProps();
  assert.deepEqual(rowProps, [
    {name: 'name', type: 'string', semType: 'Text', writable: true},
    {name: 'total', type: 'double', semType: undefined, writable: true},
    {name: 'Mol Weight', type: 'double', semType: 'Molecule Weight', writable: true},
  ]);
  assert.equal(df.dart.reads, 0, 'enumeration never read a cell — no column signal was allocated');

  source.bindStep('currentRow').bindStep('name');
  assert.ok(df.dart.reads > 0, 'the first bindStep allocates and reads');
  scope.dispose();
});

smoke('an index past the last row reads empty — and comes back when it is a row again', async () => {
  const scope = new Scope();
  const df = orders();
  const source = dfBindings(signal(df), scope);
  const idx = source.bindStep('currentRowIdx');
  const name = source.bindStep('currentRow').bindStep('name');
  assert.equal(name.value, 'aspirin');

  idx.value = 2;
  assert.equal(name.value, undefined,
    'the row does not exist: the previous row\'s text sitting there reads as real data');
  idx.value = 5;
  assert.equal(name.value, undefined, 'and the same answer whichever out-of-range index it is');
  assert.deepEqual(df.dart.writes, [], 'nothing is written back to a row that is not there');
  assert.equal(df.currentRowIdx, 0,
    'nor is the frame moved there — the platform throws on an out-of-bounds row');

  idx.value = 1;
  assert.equal(name.value, 'ibuprofen', 'back in range, the fields fill again');
  idx.value = -1;
  assert.equal(name.value, undefined, 'no current row is no value either');
  scope.dispose();
});
