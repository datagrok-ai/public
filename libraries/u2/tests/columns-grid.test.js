/* The ColumnsInput platform popup (src/dg/inputs/columns.ts `_openPlatform`): the REAL ColumnsGrid
   in checkbox mode hosted with DIALOG semantics — checks buffer in the grid, OK commits ONE write
   through the existing value path, Cancel/Esc/outside-dismiss write nothing. Driven against the
   shared FakeColumnGrid through the dg stub's live binding; the native popup stays the
   feature-absent backend (columns.test.js pins it under the ColumnGrid-less platform stub). */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Input} from '../src/core/input-base.js';
import {DataFrame} from './platform-doubles.mjs';
import {FakeColumnGrid} from './fake-column-grid.mjs';

register('./dg-stub.mjs', import.meta.url);
const DG = await import('datagrok-api/dg');
const {ColumnsInput} = await import('../src/dg/inputs/columns.js');

function selector(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      DG._installColumnGrid(FakeColumnGrid);
      await body();
    } finally {
      DG._installColumnGrid(undefined);
      FakeColumnGrid.calls.length = 0;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const DEMO = () => new DataFrame([
  {name: 'age', type: 'int'}, {name: 'height', type: 'double'},
  {name: 'sex', type: 'string'}, {name: 'started', type: 'datetime'},
]);

/** The popup is anchored to the control, and an overlay whose anchor left the document closes
 * itself — so every input that opens one is mounted. */
function make(options = {}) {
  const input = new ColumnsInput({label: 'Columns', table: DEMO(), ...options});
  document.body.append(input.root);
  return input;
}

function control(input) {
  return input.root.querySelector('.u2-columns');
}

function host() {
  return document.body.querySelector('.u2-overlay .u2-columns-popup');
}

function open(input) {
  fire(control(input), 'click');
  const cg = FakeColumnGrid.calls.at(-1);
  assert.ok(cg, 'the feature-detected platform path is taken');
  return cg;
}

function press(text) {
  fire(host().querySelectorAll('button').find((b) => b.textContent === text), 'click');
}

/** Every value write and who spoke it — `[value, isSystemWrite]`, the pickers convention. */
function track(input) {
  const writes = [];
  input.effect(() => writes.push([input.value.value, Input.isSystemWrite]));
  writes.length = 0;
  return writes;
}

selector('the platform grid is the popup: hosted with OK/CANCEL, no native list, filter intact', () => {
  const table = DEMO();
  const input = new ColumnsInput({label: 'Columns', table, filter: (c) => c.type !== 'datetime'});
  document.body.append(input.root);
  const cg = open(input);

  assert.equal(cg.table, table, 'the grid is built over the bound table');
  assert.equal(cg.options.checkAll, undefined, 'nothing pre-checked but the seed');
  assert.equal(typeof cg.options.filter, 'function');
  assert.equal(cg.options.filter(table.columns.byName('age')), true);
  assert.equal(cg.options.filter(table.columns.byName('started')), false, 'the predicate crosses intact');
  assert.deepEqual(cg.columns.map((c) => c.name), ['age', 'height', 'sex']);
  assert.ok(cg.root.isConnected, 'the grid root is hosted');
  assert.equal(cg.showSearch, true, 'search visible from open');
  assert.equal(document.activeElement, cg.search, 'and focused');
  assert.equal(host().querySelector('.u2-columns-list'), null, 'no native list rides along');
  assert.deepEqual(host().querySelectorAll('.u2-columns-buttons button').map((b) => b.textContent),
    ['OK', 'CANCEL'], 'the dialog buttons are the u2 row');
  input.dispose();
  assert.ok(cg.closed > 0, 'disposal closes the grid (the scope-owned teardown)');
});

selector('opening seeds the checks from the current value — through setSelectedColumns, not the ctor', () => {
  const input = make({value: ['sex', 'age']});
  const cg = open(input);
  assert.deepEqual(cg.getCheckedColumnNames(), ['age', 'sex'],
    'the current picks arrive checked (grid order)');
  input.dispose();
});

selector('toggle + OK is exactly one write: kept picks in value order, new ones appended', () => {
  const input = make({value: ['sex']});
  const writes = track(input);
  const cg = open(input);

  cg.toggle('age');
  assert.deepEqual(writes, [], 'checks buffer in the grid — no write before OK');
  press('OK');
  assert.equal(host(), null, 'OK closes');
  assert.deepEqual(writes, [[['sex', 'age'], false]],
    'one write, the user speaking, selected-on-top order');
  input.dispose();
});

selector('unchecking to empty commits the empty array; OK with nothing touched writes nothing', () => {
  const input = make({value: ['sex']});
  const writes = track(input);
  let cg = open(input);
  press('OK');
  assert.deepEqual(writes, [], 'an untouched OK is not a write');

  cg = open(input);
  cg.toggle('sex');
  press('OK');
  assert.deepEqual(writes, [[[], false]], 'clearing every check commits []');
  input.dispose();
});

selector('CANCEL discards the toggles — zero writes', () => {
  const input = make({value: ['sex']});
  const writes = track(input);
  const cg = open(input);
  cg.toggle('age');
  cg.toggle('sex');
  press('CANCEL');
  assert.equal(host(), null, 'CANCEL closes');
  assert.deepEqual(writes, [], 'and commits nothing');
  assert.deepEqual(input.value.value, ['sex']);
  input.dispose();
});

selector('Esc from the search box cancels even though the Dart box swallows it', () => {
  const input = make({value: ['sex']});
  const writes = track(input);
  const cg = open(input);
  cg.toggle('age');

  fire(cg.search, 'keydown', {key: 'Escape'});
  assert.equal(host(), null, 'the capture-phase listener closes before the swallow');
  assert.ok(cg.closed > 0);
  assert.deepEqual(writes, [], 'no write');
  assert.deepEqual(input.value.value, ['sex']);
  assert.equal(document.activeElement, control(input), 'focus returns to the editor');
  input.dispose();
});

selector('an outside pointerdown and a second control click both cancel', () => {
  const input = make({value: ['sex']});
  const writes = track(input);
  let cg = open(input);
  cg.toggle('age');
  fire(document.body, 'pointerdown');
  assert.equal(host(), null, 'an outside pointerdown closes the popup');
  assert.ok(cg.closed > 0);

  cg = open(input);
  cg.toggle('age');
  fire(control(input), 'click');
  assert.equal(host(), null, 'the second click closes it');
  assert.deepEqual(writes, [], 'both dismissals discard the pending checks');
  assert.deepEqual(input.value.value, ['sex']);
  input.dispose();
});

selector('Enter in the search box commits like OK', () => {
  const input = make();
  const writes = track(input);
  const cg = open(input);
  cg.toggle('height');
  fire(cg.search, 'keydown', {key: 'Enter'});
  assert.equal(host(), null, 'Enter in the search commits and closes');
  assert.deepEqual(writes, [[['height'], false]]);
  input.dispose();
});

selector('changeTable closes the platform popup; the next open builds a fresh grid over the new table', () => {
  const input = make({value: ['sex']});
  const writes = track(input);
  const cg = open(input);
  cg.toggle('age');

  const other = new DataFrame([{name: 'x', type: 'int'}, {name: 'y', type: 'int'}], [], 'other');
  input.changeTable(other);
  assert.ok(cg.closed > 0, 'the open popup goes with the old table');
  assert.equal(host(), null);
  assert.deepEqual(writes.splice(0), [[[], false]], 'only changeTable\'s own clear is written');

  const cg2 = open(input);
  assert.notEqual(cg2, cg, 'a fresh grid per open');
  assert.equal(cg2.table, other);
  press('CANCEL');
  input.dispose();
});

selector('reopening the popup leaves nothing behind', () => {
  const input = make();
  const disposers = input.scope._disposers.length;
  const live = Scope.liveCount;
  for (let i = 0; i < 10; i++) {
    open(input);
    fire(document.body, 'pointerdown');
  }
  assert.equal(host(), null);
  assert.equal(FakeColumnGrid.calls.length, 10, 'a fresh grid per open');
  assert.ok(FakeColumnGrid.calls.every((cg) => cg.closed > 0), 'each one closed');
  assert.equal(Scope.liveCount, live, 'each open owns its own scope, and closing releases it');
  assert.equal(input.scope._disposers.length, disposers, 'no per-open cleanup piles up on the input');
  input.dispose();
});

selector('without the platform grid the native popup serves, byte-unchanged', () => {
  DG._installColumnGrid(undefined);
  const input = make();
  fire(control(input), 'click');
  assert.equal(FakeColumnGrid.calls.length, 0);
  assert.ok(host().querySelector('.u2-columns-list'), 'the u2 check-list is what opens');
  press('CANCEL');
  input.dispose();
});
