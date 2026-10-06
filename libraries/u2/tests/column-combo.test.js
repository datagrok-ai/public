/* The pickers-parity ColumnInput (src/dg/inputs/column-combo.ts) to the plan-pickers FP-P-4
   contract, blind to the implementation (the zero-reconciliation method): the name-valued editor,
   the native popup backend end-to-end (open, filter, pick, keyboard, dismiss), changeTable and the
   internal frame following, the closed-combo keys, and the backend selection — the dg stub ships
   no ColumnGrid so the native path serves headlessly, and a fake installed through the stub's
   live binding proves the platform path. Tables are the getter-backed frames of
   tests/platform-doubles.mjs. */

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
const {ColumnInput} = await import('../src/dg/inputs/column-combo.js');
const {columnInput} = await import('../src/dg/inputs/pickers.js');

function combo(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
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
  {name: 'sex', type: 'string', semType: 'Sex'}, {name: 'started', type: 'datetime'},
]);

/** The popup is anchored to the control, and an overlay whose anchor left the document closes
 * itself — so every combo that opens one is mounted. */
function make(options = {}) {
  const input = new ColumnInput({label: 'Column', table: DEMO(), ...options});
  document.body.append(input.root);
  return input;
}

function control(input) {
  return input.root.querySelector('.u2-columns');
}

function summary(input) {
  return input.root.querySelector('.u2-columns-summary').textContent;
}

function popup() {
  return document.body.querySelector('.u2-overlay .u2-list');
}

function open(input) {
  fire(control(input), 'click');
  const list = popup();
  assert.ok(list, 'the native popup opens');
  // the shim lays nothing out: give the list a viewport — and rewind the highlight scroll,
  // which the layoutless shim never clamps — so every row is rendered
  list.clientHeight = 400;
  list.scrollTop = 0;
  fire(list, 'scroll');
  return list;
}

/** In index order — recycling appends rendered rows in render order, not list order. */
function rows(list) {
  return list.querySelectorAll('.u2-list-row').sort((a, b) => Number(a.dataset.index) - Number(b.dataset.index));
}

function rowName(row) {
  return row.querySelector('.u2-column-name')?.textContent ?? '';
}

function rowOf(list, name) {
  return rows(list).find((row) => rowName(row) === name);
}

function searchBox() {
  return document.body.querySelector('.u2-overlay .u2-input-search input');
}

/** Every value write and who spoke it — `[value, isSystemWrite]`, the pickers convention. */
function track(input) {
  const writes = [];
  input.effect(() => writes.push([input.value.value, Input.isSystemWrite]));
  writes.length = 0;
  return writes;
}

// --- the editor and the value contract ---

combo('the value is the column name, shown on the editor; the placeholder covers null', () => {
  const input = make({value: 'sex', placeholder: 'Pick a column'});
  assert.equal(input.root.dataset.u2, 'column-combo');
  const ctl = control(input);
  assert.ok(ctl, 'the shared u2-columns editor serves the combo');
  assert.equal(ctl.getAttribute('role'), 'button');
  assert.equal(ctl.getAttribute('aria-haspopup'), 'listbox');
  assert.equal(summary(input), 'sex');

  input.value.value = 'age';
  assert.equal(summary(input), 'age');
  input.value.value = null;
  assert.equal(summary(input), 'Pick a column');
  input.dispose();

  const bare = make();
  assert.equal(summary(bare), '', 'no placeholder is an empty box');
  bare.dispose();
});

combo('names(): the filter decides the offer; a null table offers nothing', () => {
  const table = DEMO();
  const input = new ColumnInput({label: 'Column', table, filter: (c) => c.isNumerical});
  assert.equal(input.table, table);
  assert.deepEqual(input.names(), ['age', 'height', 'started'],
    'datetime is numerical to the platform');
  input.dispose();

  const none = new ColumnInput({label: 'Column', table: null});
  assert.deepEqual(none.names(), []);
  none.dispose();
});

combo('changeTable closes the popup, retargets the offer and clears the value', () => {
  const table = DEMO();
  const input = make({table, value: 'sex'});
  open(input);

  const other = new DataFrame([{name: 'x', type: 'int'}, {name: 'y', type: 'string'}], [], 'other');
  input.changeTable(other, (c) => c.isNumerical);
  assert.equal(popup(), null, 'the open popup goes with the old table');
  assert.equal(input.value.value, null, 'the selection belonged to the old table');
  assert.deepEqual(input.names(), ['x'], 'the new filter rides along');
  assert.equal(table.liveSubscriptions(), 0, 'the old table is let go');
  input.dispose();
  assert.equal(other.liveSubscriptions(), 0);
});

combo('a null table reads inert — aria-disabled, tooltip, closed popup — and flips with changeTable', () => {
  const input = make({table: null});
  const ctl = control(input);
  assert.equal(ctl.getAttribute('aria-disabled'), 'true');
  assert.equal(input.box.getAttribute('title'), 'Select a table first');
  fire(ctl, 'click');
  assert.equal(popup(), null, 'the popup refuses');

  input.changeTable(DEMO());
  assert.equal(ctl.getAttribute('aria-disabled'), 'false');
  assert.equal(input.box.hasAttribute('title'), false);
  fire(ctl, 'click');
  assert.ok(popup(), 'a table brings it back');
  input.dispose();
});

// --- the native popup backend ---

combo('native popup: renderer rows with the empty row first, the current value highlighted', () => {
  const input = make({value: 'sex'});
  const list = open(input);
  assert.deepEqual(rows(list).map(rowName), ['', 'age', 'height', 'sex', 'started'],
    'a nullable combo leads with the empty row (FP-P-4)');
  assert.ok(rowOf(list, 'sex').querySelector('.u2-column-semtype'),
    'columnRenderer rows — the semType hint rides along');
  const selected = rows(list).filter((row) => row.getAttribute('aria-selected') === 'true');
  assert.deepEqual(selected.map(rowName), ['sex'], 'the held value is the highlighted row');
  input.dispose();
});

combo('nullable: false drops the empty row', () => {
  const input = make({nullable: false});
  const list = open(input);
  assert.deepEqual(rows(list).map(rowName), ['age', 'height', 'sex', 'started']);
  input.dispose();
});

combo('the search box opens focused and filters as typed', () => {
  const input = make();
  const list = open(input);
  const search = searchBox();
  assert.ok(search, 'the search sits on top (divergence #15 — visible from open)');
  assert.equal(document.activeElement, search, 'the popup opens ready to filter');

  search.value = 'he';
  fire(search, 'input');
  assert.deepEqual(rows(list).map(rowName).filter((name) => name !== ''), ['height']);
  input.dispose();
});

combo('a row click commits as a user edit and closes', () => {
  const input = make();
  const writes = track(input);
  const list = open(input);
  fire(rowOf(list, 'height'), 'click');
  assert.equal(popup(), null, 'a pick closes the popup');
  assert.deepEqual(writes, [['height', false]], 'one write, the user speaking');
  assert.equal(summary(input), 'height');
  input.dispose();
});

combo('the empty row picks null', () => {
  const input = make({value: 'sex'});
  const writes = track(input);
  const list = open(input);
  fire(rowOf(list, ''), 'click');
  assert.equal(popup(), null);
  assert.deepEqual(writes, [[null, false]], 'the empty row answers null, as a user edit');
  assert.equal(summary(input), '');
  input.dispose();
});

combo('arrows move the highlight from the search box; Enter picks the highlighted row', () => {
  const input = make();
  const writes = track(input);
  const list = open(input);
  const search = searchBox();
  fire(search, 'keydown', {key: 'ArrowDown'});
  fire(search, 'keydown', {key: 'ArrowDown'});
  const selected = rows(list).find((row) => row.getAttribute('aria-selected') === 'true');
  assert.ok(selected, 'the arrows move a visible highlight');
  const picked = rowName(selected);
  assert.notEqual(picked, '', 'two steps down are past the empty row');

  fire(search, 'keydown', {key: 'Enter'});
  assert.equal(popup(), null, 'Enter picks and closes');
  assert.deepEqual(writes, [[picked, false]], 'the highlighted row is what lands, a user edit');
  input.dispose();
});

combo('narrowing the filter below the held highlight revives Enter — the unique match lands', () => {
  const input = make({value: 'started'});
  const writes = track(input);
  open(input);
  const search = searchBox();
  // the held value highlights row 4 (past the empty row); 'he' leaves two rows
  search.value = 'he';
  fire(search, 'input');
  fire(search, 'keydown', {key: 'Enter'});
  assert.equal(popup(), null, 'Enter picks and closes');
  assert.deepEqual(writes, [['height', false]], 'the single passing match is what lands');
  input.dispose();
});

combo('Esc and an outside click both close without a write', () => {
  const input = make({value: 'sex'});
  const writes = track(input);
  open(input);
  fire(searchBox(), 'keydown', {key: 'ArrowDown'});
  fire(searchBox(), 'keydown', {key: 'Escape'});
  assert.equal(popup(), null, 'Esc closes');
  assert.equal(input.value.value, 'sex');

  open(input);
  fire(document.body, 'pointerdown');
  assert.equal(popup(), null, 'an outside pointerdown closes');
  assert.deepEqual(writes, [], 'neither dismissal wrote anything');
  input.dispose();
});

combo('reopening the popup leaves nothing behind', () => {
  const input = make();
  const disposers = input.scope._disposers.length;
  const live = Scope.liveCount;
  for (let i = 0; i < 10; i++) {
    open(input);
    fire(document.body, 'pointerdown');
  }
  assert.equal(popup(), null);
  assert.equal(Scope.liveCount, live, 'each open owns its own scope, and closing releases it');
  assert.equal(input.scope._disposers.length, disposers, 'no per-open cleanup piles up on the input');
  input.dispose();
});

// --- frame following (the pickers convention, moved into the control per FP-P-4) ---

combo('frame following: a rename remaps, a drop prunes — both the input speaking for itself', () => {
  const table = DEMO();
  const input = make({table, value: 'sex'});
  const writes = track(input);

  table.columns.byName('sex').name = 'gender';
  assert.deepEqual(writes.splice(0), [['gender', true]], 'the remap is a system write');
  assert.deepEqual(input.names(), ['age', 'height', 'gender', 'started']);

  table.columns.byName('age').name = 'years';
  assert.deepEqual(writes, [], 'renaming another column leaves the value alone');

  table.columns.remove('gender');
  assert.deepEqual(writes.splice(0), [[null, true]], 'a dropped column prunes to null');
  assert.deepEqual(input.names(), ['years', 'height', 'started']);

  table.columns.remove('years');
  assert.deepEqual(writes, [], 'dropping an unheld column changes nothing');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0, 'disposal takes the follow with it');
});

combo('dispose closes an open popup and releases the frame', () => {
  const table = DEMO();
  const input = make({table});
  open(input);
  input.dispose();
  assert.equal(popup(), null, 'the popup goes with the input');
  assert.equal(table.liveSubscriptions(), 0);
});

// --- the closed-combo keys ---

combo('closed-combo arrows step through the passing columns as user edits', () => {
  const input = make({value: 'age', filter: (c) => c.isNumerical});
  const writes = track(input);
  const ctl = control(input);

  fire(ctl, 'keydown', {key: 'ArrowDown'});
  assert.deepEqual(writes.splice(0), [['height', false]],
    'ArrowDown is the next passing column, a real user edit');
  fire(ctl, 'keydown', {key: 'ArrowRight'});
  assert.equal(input.value.value, 'started', 'ArrowRight steps too, skipping what fails the filter');
  fire(ctl, 'keydown', {key: 'ArrowUp'});
  assert.equal(input.value.value, 'height', 'ArrowUp steps back');
  fire(ctl, 'keydown', {key: 'ArrowLeft'});
  assert.equal(input.value.value, 'age');
  assert.equal(popup(), null, 'stepping never opens the popup');
  input.dispose();
});

combo('Space and Enter open the popup; a printable key opens it with the search seeded', () => {
  const input = make();
  fire(control(input), 'keydown', {key: 'Enter'});
  assert.ok(popup(), 'Enter opens');
  fire(document.body, 'pointerdown');
  assert.equal(popup(), null);

  fire(control(input), 'keydown', {key: ' '});
  assert.ok(popup(), 'Space opens');
  fire(document.body, 'pointerdown');

  fire(control(input), 'keydown', {key: 'h'});
  const list = popup();
  assert.ok(list, 'a printable key opens');
  assert.equal(searchBox().value, 'h', 'seeded with what was typed');
  list.clientHeight = 400;
  fire(list, 'scroll');
  assert.deepEqual(rows(list).map(rowName).filter((name) => name !== ''), ['height'],
    'and already filtering');
  input.dispose();
});

// --- the factory overload ---

combo('columnInput {grid: true} answers the ColumnInput; the plain default is untouched', () => {
  const input = columnInput('Column', DEMO(), {grid: true, filter: (c) => c.isNumerical});
  assert.ok(input instanceof ColumnInput);
  assert.equal(input.root.dataset.u2, 'column-combo');
  assert.deepEqual(input.names(), ['age', 'height', 'started']);
  input.dispose();

  const plain = columnInput('Column', DEMO());
  assert.ok(!(plain instanceof ColumnInput), 'the existing call shape keeps the plain choice');
  plain.dispose();
});

// --- backend selection (FP-P-4: chosen once at open, structurally) ---

combo('no DG.ColumnGrid: the native backend serves the popup', () => {
  assert.equal(DG.ColumnGrid, undefined, 'the stub ships no ColumnGrid');
  const input = make();
  assert.ok(open(input), 'so the u2 list is what opens');
  input.dispose();
});

combo('DG.ColumnGrid present: the platform popup is the backend', () => {
  DG._installColumnGrid(FakeColumnGrid);
  const table = DEMO();
  const input = make({table, value: 'sex', filter: (c) => c.type !== 'datetime'});
  fire(control(input), 'click');

  assert.equal(FakeColumnGrid.calls.length, 1, 'the feature-detected path is taken');
  const cg = FakeColumnGrid.calls[0];
  assert.equal(cg.table, table, 'over the bound table');
  assert.equal(cg.options.addEmpty, true, 'nullable crosses as addEmpty');
  assert.equal(typeof cg.options.filter, 'function');
  assert.equal(cg.options.filter(table.columns.byName('age')), true);
  assert.equal(cg.options.filter(table.columns.byName('started')), false, 'the predicate crosses intact');
  assert.ok(cg.root.isConnected, 'the grid root is hosted');
  assert.equal(popup(), null, 'and no native list rides along');
  input.dispose();
  assert.ok(cg.closed > 0, 'disposal closes the grid (the scope-owned teardown)');
});

combo('platform backend: nullable: false crosses as addEmpty: false', () => {
  DG._installColumnGrid(FakeColumnGrid);
  const input = make({nullable: false});
  fire(control(input), 'click');
  assert.equal(FakeColumnGrid.calls[0].options.addEmpty, false);
  input.dispose();
});

combo('platform backend: onCurrentRowChanged commits the current column and closes', () => {
  DG._installColumnGrid(FakeColumnGrid);
  const input = make();
  const writes = track(input);
  fire(control(input), 'click');
  const cg = FakeColumnGrid.calls[0];

  cg.pick('height');
  assert.deepEqual(writes.splice(0), [['height', false]], 'the pick lands as a user edit');
  assert.ok(cg.closed > 0, 'and the popup closes');
  assert.ok(!cg.root.isConnected);

  fire(control(input), 'click');
  const cg2 = FakeColumnGrid.calls[1];
  assert.notEqual(cg2, cg, 'a fresh grid per open (the CCB lifecycle)');
  cg2.pick('');
  assert.deepEqual(writes.splice(0), [[null, false]], 'the addEmpty row answers null');
  input.dispose();
});

combo('platform backend: Enter with no walked row commits the grid\'s verdict — a friendlyName match', () => {
  DG._installColumnGrid(FakeColumnGrid);
  const table = DEMO();
  table.columns.byName('height').friendlyName = 'stature';
  const input = new ColumnInput({label: 'Column', table});
  document.body.append(input.root);
  const writes = track(input);
  fire(control(input), 'click');
  const cg = FakeColumnGrid.calls[0];

  // 'stat' matches no column NAME — only the friendlyName the grid itself would show
  cg.search.value = 'stat';
  fire(cg.search, 'keydown', {key: 'Enter'});
  assert.deepEqual(writes, [['height', false]], 'the grid\'s unique passing column lands');
  assert.ok(cg.closed > 0, 'and the popup closes');
  input.dispose();
});

combo('platform backend: Esc closes the grid without a write', () => {
  DG._installColumnGrid(FakeColumnGrid);
  const input = make({value: 'sex'});
  const writes = track(input);
  fire(control(input), 'click');
  const cg = FakeColumnGrid.calls[0];

  fire(cg.root, 'keydown', {key: 'Escape'});
  assert.ok(cg.closed > 0, 'Esc tears the grid down');
  assert.ok(!cg.root.isConnected);
  assert.deepEqual(writes, [], 'and writes nothing');
  assert.equal(input.value.value, 'sex');
  input.dispose();
});

combo('platform backend: Esc from the search box closes even though the Dart box swallows it', () => {
  DG._installColumnGrid(FakeColumnGrid);
  const input = make({value: 'sex'});
  const writes = track(input);
  fire(control(input), 'click');
  const cg = FakeColumnGrid.calls[0];
  assert.equal(document.activeElement, cg.search, 'focus sits in the grid search box');

  fire(cg.search, 'keydown', {key: 'Escape'});
  assert.ok(cg.closed > 0, 'the capture-phase listener closes before the swallow');
  assert.ok(!cg.root.isConnected);
  assert.deepEqual(writes, [], 'no write');
  assert.equal(input.value.value, 'sex');
  assert.equal(document.activeElement, control(input), 'focus returns to the combo editor');
  input.dispose();
});
