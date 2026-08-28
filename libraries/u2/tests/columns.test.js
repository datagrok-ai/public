/* The columns family (src/dg/inputs/columns.ts) over the stub in tests/platform-stub.mjs: the ported
   aggregation rules, the summary label and the anchored selection popup of columnsInput, the
   per-key pickers of columnsMapInput, and the aggregation rows. Tables are the getter-backed
   frames of tests/platform-doubles.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Input} from '../src/core/input-base.js';
import {DataFrame} from './platform-doubles.mjs';

register('./platform-stub.mjs', import.meta.url);
const {ColumnsInput, columnsInput, columnsMapInput, aggregatedColumnsInput, aggregationsFor,
  defaultAggregation} = await import('../src/dg/inputs/columns.js');

function columns(name, body) {
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

const DEMO = () => new DataFrame([
  {name: 'age', type: 'int'}, {name: 'height', type: 'double'},
  {name: 'sex', type: 'string'}, {name: 'started', type: 'datetime'},
]);

/** The popup is anchored to the control, and an overlay whose anchor left the document closes
 * itself — so every columnsInput that opens one is mounted. */
function picker(table = DEMO(), options = {}) {
  const input = columnsInput('Columns', table, options);
  document.body.append(input.root);
  return input;
}

function summary(input) {
  return input.root.querySelector('.u2-columns-summary').textContent;
}

function control(input) {
  return input.root.querySelector('.u2-columns');
}

function popup() {
  return document.body.querySelector('.u2-columns-popup');
}

function open(input) {
  fire(control(input), 'click');
  const box = popup();
  const list = box.querySelector('.u2-columns-list');
  // the shim lays nothing out: give the list a viewport so every row is rendered
  list.clientHeight = 400;
  fire(list, 'scroll');
  return box;
}

function rowNames(box) {
  return box.querySelectorAll('.u2-column-name').map((el) => el.textContent);
}

function toggle(box, name) {
  const row = box.querySelectorAll('.u2-list-row')
    .find((el) => el.querySelector('.u2-column-name').textContent === name);
  fire(row.querySelector('input[type=checkbox]'), 'click');
}

function press(box, text) {
  fire(box.querySelectorAll('button').find((b) => b.textContent === text), 'click');
}

/** What the platform's outside-mousedown dismissal does: the overlay listens on the document. */
function clickOutside() {
  fire(document.body, 'pointerdown');
}

columns('aggregationsFor: the per-type exclusions of the Dart editor are honoured', () => {
  assert.deepEqual(aggregationsFor('string'), ['count', 'values', 'unique', 'nulls'],
    'the string functions and first are excluded');
  assert.deepEqual(aggregationsFor('bool'), ['count', 'values', 'unique', 'nulls']);
  assert.deepEqual(aggregationsFor('datetime'), ['count', 'values', 'unique', 'nulls', 'range'],
    'datetime keeps range, loses first/min/max/avg');
  const numeric = aggregationsFor('double');
  assert.deepEqual(numeric, aggregationsFor('int'));
  assert.deepEqual(numeric, aggregationsFor('qnum'));
  assert.ok(numeric.includes('first') && numeric.includes('avg') && numeric.includes('geomean'));
  assert.deepEqual(aggregationsFor('byte_array'), ['count', 'values', 'unique', 'nulls'],
    'an unregistered type counts rows and nothing else');
});

columns('defaultAggregation: avg where the type offers it, the first entry otherwise', () => {
  const table = DEMO();
  const column = (name) => table.columns.byName(name);
  assert.equal(defaultAggregation(column('age')), 'avg');
  assert.equal(defaultAggregation(column('sex')), 'count');
  assert.equal(defaultAggregation(column('started')), 'count',
    'datetime is numerical to the platform but has no avg left');
});

columns('columnsInput: the summary reads as Dart writes it', () => {
  const table = DEMO();
  const input = columnsInput('Columns', table);
  assert.equal(summary(input), '(0)');

  input.value.value = ['age', 'sex'];
  assert.equal(summary(input), '(2) age, sex');

  input.value.value = ['age', 'height', 'sex', 'started'];
  assert.equal(summary(input), '(4) All');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

columns('columnsInput: the popup commits on OK, keeps the value on cancel', () => {
  const input = picker();
  input.value.value = ['sex'];

  let box = open(input);
  toggle(box, 'age');
  press(box, 'CANCEL');
  assert.deepEqual(input.value.value, ['sex'], 'cancel commits nothing');
  assert.equal(popup(), null);

  box = open(input);
  toggle(box, 'age');
  press(box, 'OK');
  assert.deepEqual(input.value.value, ['sex', 'age'], 'selected-on-top order is the value order');
  assert.equal(summary(input), '(2) sex, age');
  assert.equal(popup(), null);
  input.dispose();
});

/** Dart parity (`dialog.dart:486-494`): the selector is dismissed by anything outside it, and
 * every such dismissal discards what was checked. */
columns('columnsInput: an outside click and Esc both cancel', () => {
  const input = picker();
  input.value.value = ['sex'];

  let box = open(input);
  toggle(box, 'age');
  clickOutside();
  assert.equal(popup(), null, 'an outside pointerdown closes the popup');
  assert.deepEqual(input.value.value, ['sex'], 'and discards the pending checks');

  box = open(input);
  toggle(box, 'age');
  fire(box, 'keydown', {key: 'Escape'});
  assert.equal(popup(), null);
  assert.deepEqual(input.value.value, ['sex'], 'Esc cancels too');
  input.dispose();
});

/** `columns_input.dart:84-86`: a second click on the control closes the selector it opened. */
columns('columnsInput: clicking the control toggles the popup, and cancels what is open', () => {
  const input = picker();
  const box = open(input);
  toggle(box, 'age');

  fire(control(input), 'click');
  assert.equal(popup(), null, 'the second click closes it');
  assert.deepEqual(input.value.value, [], 'as a cancel');
  assert.equal(control(input).getAttribute('aria-expanded'), 'false');

  fire(control(input), 'click');
  assert.ok(popup(), 'and the third opens it again');
  assert.equal(control(input).getAttribute('aria-expanded'), 'true');
  input.dispose();
});

columns('columnsInput: reopening the popup leaves nothing behind', () => {
  const input = picker();
  const disposers = input.scope._disposers.length;
  const live = Scope.liveCount;
  for (let i = 0; i < 10; i++) {
    open(input);
    press(popup(), 'CANCEL');
  }
  assert.equal(popup(), null);
  assert.equal(Scope.liveCount, live, 'each open owns its own scope, and closing releases it');
  assert.equal(input.scope._disposers.length, disposers, 'no per-open cleanup piles up on the input');
  input.dispose();
});

columns('columnsInput: a disabled input does not open the popup', () => {
  const input = picker();
  input.enabled = false;
  fire(control(input), 'click');
  assert.equal(popup(), null);

  input.enabled = true;
  fire(control(input), 'click');
  assert.ok(popup());
  input.dispose();
});

columns('columnsInput: selected columns come first, the search box filters', () => {
  const input = picker(DEMO(), {filter: (c) => c.name !== 'started'});
  input.value.value = ['sex'];
  const box = open(input);
  assert.deepEqual(rowNames(box), ['sex', 'age', 'height'], 'picked first, filter applied');

  const search = box.querySelector('.u2-input-search input');
  search.value = 'he';
  fire(search, 'input');
  assert.deepEqual(rowNames(box), ['height']);

  toggle(box, 'height');
  press(box, 'OK');
  assert.deepEqual(input.value.value, ['sex', 'height']);
  input.dispose();
});

columns('columnsInput: arrows and Space pick without a mouse, Enter is OK', () => {
  const input = picker();
  const box = open(input);
  const list = box.querySelector('.u2-columns-list');
  fire(list, 'keydown', {key: 'ArrowDown'});
  fire(list, 'keydown', {key: 'ArrowDown'});
  fire(list, 'keydown', {key: ' '});
  fire(list, 'keydown', {key: 'Enter'});
  assert.equal(popup(), null, 'Enter over the list commits and closes');
  assert.deepEqual(input.value.value, ['height']);
  input.dispose();
});

/** The focus opens in the search box, so the popup carries the row keys for it — arrow, Space and
 * Enter all work without ever leaving it. */
columns('columnsInput: the row keys work from the search box, Enter in it commits', () => {
  const input = picker();
  const box = open(input);
  const search = box.querySelector('.u2-input-search input');
  assert.equal(document.activeElement, search, 'the popup opens ready to filter');

  fire(box.querySelectorAll('button').find((b) => b.textContent === 'CANCEL'), 'keydown', {key: 'Enter'});
  assert.ok(popup(), 'Enter on a button is left to the button');

  fire(search, 'keydown', {key: 'ArrowDown'});
  fire(search, 'keydown', {key: ' '});
  assert.equal(box.querySelector('.u2-columns-count').textContent, '1 selected');
  fire(search, 'keydown', {key: 'Enter'});
  assert.equal(popup(), null, 'Enter in the search box commits and closes');
  assert.deepEqual(input.value.value, ['age']);
  input.dispose();
});

columns('columnsInput: a dropped column leaves the value; changeTable resets it', () => {
  const table = DEMO();
  const input = columnsInput('Columns', table);
  input.value.value = ['age', 'sex'];
  table.columns.remove('sex');
  assert.deepEqual(input.value.value, ['age']);

  const other = new DataFrame([{name: 'x', type: 'int'}, {name: 'y', type: 'int'}]);
  input.changeTable(other);
  assert.deepEqual(input.value.value, []);
  assert.equal(summary(input), '(0)');
  assert.equal(table.liveSubscriptions(), 0, 'the old table is let go');

  input.value.value = ['x', 'y'];
  assert.equal(summary(input), '(2) All');
  input.dispose();
  assert.equal(other.liveSubscriptions(), 0);
});

columns('columnsInput: a null table reads inert — aria-disabled, tooltip, closed popup — and flips with changeTable', () => {
  const input = new ColumnsInput({label: 'Columns', table: null});
  document.body.append(input.root);
  const ctl = control(input);
  assert.equal(ctl.getAttribute('aria-disabled'), 'true');
  assert.equal(input.box.getAttribute('title'), 'Select a table first');
  fire(ctl, 'click');
  assert.equal(popup(), null, 'the popup stays closed');

  input.changeTable(DEMO());
  assert.equal(ctl.getAttribute('aria-disabled'), 'false');
  assert.equal(input.box.hasAttribute('title'), false);
  fire(ctl, 'click');
  assert.ok(popup(), 'a table brings the popup back');

  input.changeTable(null);
  assert.equal(popup(), null, 'changeTable closes it');
  assert.equal(ctl.getAttribute('aria-disabled'), 'true', 'and the affordance flips back');
  assert.equal(input.box.getAttribute('title'), 'Select a table first');
  input.dispose();
});

columns('columnsInput: a renamed column is remapped, and as a system write', () => {
  const table = DEMO();
  let systemic = null;
  const input = columnsInput('Columns', table,
    {onChanged: () => systemic = Input.isSystemWrite});
  input.value.value = ['age', 'sex'];
  assert.equal(systemic, false, 'a plain write is a user edit');

  table.columns.byName('sex').name = 'gender';
  assert.deepEqual(input.value.value, ['age', 'gender'], 'remapped, not pruned');
  assert.equal(summary(input), '(2) age, gender');
  assert.equal(systemic, true, 'the remap is the input speaking for itself');

  table.columns.byName('height').name = 'cm';
  assert.deepEqual(input.value.value, ['age', 'gender'], 'a column outside the value changes nothing');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

columns('columnsInput: All and None act on what the search left, the count follows', () => {
  const input = picker();
  const box = open(input);
  const count = () => box.querySelector('.u2-columns-count').textContent;
  const act = (text) => fire(box.querySelectorAll('a').find((a) => a.textContent === text), 'click');
  assert.equal(count(), '0 selected');

  act('All');
  assert.equal(count(), '4 selected');
  act('None');
  assert.equal(count(), '0 selected');

  const search = box.querySelector('.u2-input-search input');
  search.value = 'e';
  fire(search, 'input');
  assert.deepEqual(rowNames(box), ['age', 'height', 'sex', 'started']);
  search.value = 'he';
  fire(search, 'input');
  act('All');
  assert.equal(count(), '1 selected', 'only the visible rows are picked');
  press(box, 'OK');
  assert.deepEqual(input.value.value, ['height']);
  input.dispose();
});

columns('columnsInput: every option in the popup names itself for a screen reader', () => {
  const input = picker();
  const box = open(input);
  assert.deepEqual(box.querySelectorAll('.u2-columns-option input[type=checkbox]')
    .map((b) => b.getAttribute('aria-label')), ['age', 'height', 'sex', 'started']);
  press(box, 'CANCEL');
  input.dispose();
});

columns('columnsMapInput: a picker per key, values both ways, type filter per key', () => {
  const input = columnsMapInput('Map', ['x', {name: 'when', type: 'datetime'}], DEMO());
  const pickers = input.root.querySelectorAll('select');
  assert.equal(pickers.length, 2);
  assert.deepEqual(pickers[1].querySelectorAll('option').map((o) => o.value).filter((v) => v !== ''),
    ['started'], 'a typed key only offers columns of that type');

  pickers[0].value = 'age';
  fire(pickers[0], 'change');
  assert.deepEqual(input.value.value, {x: 'age'}, 'unmapped keys stay out of the value');

  input.value.value = {x: 'height', when: 'started'};
  assert.equal(pickers[0].value, 'height');
  assert.equal(pickers[1].value, 'started');
  assert.equal(input.input('when').value.value, 'started');
  input.dispose();
});

columns('columnsMapInput: int and double are one type, as in Dart', () => {
  const input = columnsMapInput('Map', [{name: 'value', type: 'int'}], DEMO());
  const items = input.root.querySelectorAll('option').map((o) => o.value).filter((v) => v !== '');
  assert.deepEqual(items, ['age', 'height']);
  input.dispose();
});

columns('columnsMapInput: a renamed column is remapped in the record and in the picker', () => {
  const table = DEMO();
  const input = columnsMapInput('Map', ['x', 'y'], table);
  input.value.value = {x: 'age', y: 'sex'};

  table.columns.byName('age').name = 'ageYears';
  assert.deepEqual(input.value.value, {x: 'ageYears', y: 'sex'});
  assert.equal(input.root.querySelectorAll('select')[0].value, 'ageYears');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

columns('aggregatedColumnsInput: a renamed column is remapped in the rows and the choices', () => {
  const table = DEMO();
  const input = aggregatedColumnsInput('Aggregations', table);
  input.value.value = [{column: 'age', aggregation: 'avg'}];

  table.columns.byName('age').name = 'ageYears';
  assert.deepEqual(input.value.value, [{column: 'ageYears', aggregation: 'avg'}]);
  const column = input.root.querySelectorAll('.u2-aggregations-row select')[0];
  assert.equal(column.value, 'ageYears');
  assert.ok(column.querySelectorAll('option').map((o) => o.value).includes('ageYears'),
    'the choices carry the new name too');
  input.dispose();
  assert.equal(table.liveSubscriptions(), 0);
});

columns('aggregatedColumnsInput: rows add and remove, an empty value shows one add button', () => {
  const input = aggregatedColumnsInput('Aggregations', DEMO());
  const buttons = () => input.root.querySelectorAll('.u2-aggregations-row [data-u2="icon-button"]');
  assert.equal(input.root.querySelectorAll('select').length, 0);
  assert.equal(buttons().length, 1, 'nothing selected — just the add button');

  fire(buttons()[0], 'click');
  assert.deepEqual(input.value.value, [{column: 'age', aggregation: 'avg'}],
    'the first numerical column at its default aggregation');

  fire(buttons()[0], 'click');
  assert.equal(input.value.value.length, 2);
  fire(buttons()[1], 'click');
  assert.equal(input.value.value.length, 1, 'the second button of the row removes it');
  input.dispose();
});

columns('aggregatedColumnsInput: the aggregation list follows the column type', () => {
  const input = aggregatedColumnsInput('Aggregations', DEMO());
  input.value.value = [{column: 'age', aggregation: 'avg'}];
  const column = () => input.root.querySelectorAll('.u2-aggregations-row select')[0];
  const aggregation = () => input.root.querySelectorAll('.u2-aggregations-row select')[1];
  assert.ok(aggregation().querySelectorAll('option').map((o) => o.value).includes('geomean'));

  column().value = 'sex';
  fire(column(), 'change');
  assert.deepEqual(aggregation().querySelectorAll('option').map((o) => o.value),
    ['count', 'values', 'unique', 'nulls']);
  assert.deepEqual(input.value.value, [{column: 'sex', aggregation: 'count'}],
    'an aggregation the new type does not offer falls back to its default');

  aggregation().value = 'unique';
  fire(aggregation(), 'change');
  assert.deepEqual(input.value.value, [{column: 'sex', aggregation: 'unique'}]);

  column().value = 'height';
  fire(column(), 'change');
  assert.deepEqual(input.value.value, [{column: 'height', aggregation: 'unique'}],
    'an aggregation both types offer survives the column change');
  input.dispose();
});
