/* RecordsInput: an editable list of typed records — keyed rows of inline inputs under a header, add, remove
   and (orderable) move; the value is the records, a new array per edit; a failing cell names its row. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, RecordsInput, TextInput, ChoiceInput, Form} from '../src/index.js';

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

const text = (value) => new TextInput({inline: true, value: value ?? '', required: true});
const fields = (options = {}) => new RecordsInput({
  label: 'Fields',
  columns: [{name: 'name', header: 'Name', input: text},
    {name: 'type', header: 'Type', width: '120px',
      input: (value) => new ChoiceInput({inline: true, items: ['String', 'Number'], value: value ?? 'String'})}],
  newItem: () => ({name: '', type: 'String'}),
  addText: 'Add field',
  ...options,
});

function mount(input) {
  document.body.append(input.root);
  return input;
}

const rows = (input) => [...input.root.querySelector('.u2-records-body').children];
const cell = (input, row, column) => rows(input)[row].querySelectorAll('.u2-records-cell')[column];

function type(el, value) {
  el.value = value;
  fire(el, 'input');
}

smoke('a row per record, an inline input per field, under a header that hides with no rows', () => {
  const input = mount(fields({value: [{name: 'Chemist', type: 'String'}, {name: 'Dose', type: 'Number'}]}));
  assert.equal(input.root.dataset.u2, 'records-input');
  assert.equal(rows(input).length, 2);
  assert.equal(cell(input, 1, 0).querySelector('input').value, 'Dose');
  assert.equal(cell(input, 1, 1).querySelector('select').value, 'Number');
  const head = input.root.querySelector('.u2-records-head');
  assert.deepEqual([...head.querySelectorAll('[role="columnheader"]')].map((h) => h.textContent), ['Name', 'Type']);
  assert.equal(head.hidden, false);
  input.value = [];
  assert.equal(rows(input).length, 0);
  assert.equal(head.hidden, true, 'no header over no rows');
  input.dispose();
});

smoke('typing in a cell rewrites that record — a new array, the other records as they were', () => {
  const first = {name: 'Chemist', type: 'String'};
  const second = {name: 'Dose', type: 'Number'};
  const input = mount(fields({value: [first, second]}));
  const kept = rows(input)[1];
  type(cell(input, 0, 0).querySelector('input'), 'Biologist');
  const value = input.value;
  assert.deepEqual(value, [{name: 'Biologist', type: 'String'}, second]);
  assert.ok(value[1] === second, 'an untouched record keeps its identity');
  assert.ok(value[0] !== first);
  assert.ok(rows(input)[1] === kept, 'an edit never rebuilds another row');
  input.dispose();
});

smoke('add appends a new record and focuses its first cell; remove takes a row out', async () => {
  const input = mount(fields({value: [{name: 'Chemist', type: 'String'}]}));
  const add = input.root.querySelector('[data-u2-part="add"]');
  assert.equal(add.textContent, 'Add field');
  fire(add, 'click');
  assert.equal(rows(input).length, 2);
  assert.deepEqual(input.value, [{name: 'Chemist', type: 'String'}, {name: '', type: 'String'}]);
  await flush();
  assert.ok(document.activeElement === cell(input, 1, 0).querySelector('input'), 'the new row\'s first cell');
  fire(rows(input)[0].querySelector('.u2-records-remove'), 'click');
  assert.deepEqual(input.value, [{name: '', type: 'String'}]);
  input.dispose();
});

smoke('the verdict names the row and the field; a whole-list rule goes through addValidator', () => {
  const input = mount(fields({value: [{name: 'Chemist', type: 'String'}, {name: '', type: 'String'}],
    rowName: (_, i) => `Field ${i + 1}`}));
  assert.equal(input.validity.value, 'Field 2: Name: Value can\'t be empty');
  assert.equal(input.getItems().list[1].invalid, true);
  type(cell(input, 1, 0).querySelector('input'), 'Chemist');
  assert.equal(input.validity.value, null);
  input.addValidator((records) => new Set(records.map((r) => r.name)).size < records.length ?
    'Names must be unique' : null);
  assert.equal(input.validity.value, 'Names must be unique');
  input.dispose();
});

smoke('one column: the verdict is the row\'s; required refuses an empty list', () => {
  const input = mount(new RecordsInput({label: 'Stages', required: true,
    columns: [{name: 'name', input: text}], value: [{name: ''}], rowName: (_, i) => `Stage ${i + 1}`}));
  assert.equal(input.validity.value, 'Stage 1: Value can\'t be empty');
  input.remove(0);
  assert.equal(input.validity.value, 'Value can\'t be empty');
  assert.equal(input.root.classList.contains('u2-input-required'), true);
  input.dispose();
});

smoke('orderable: a handle per row, move() and Alt+↑/↓ reorder the records, the rows kept', async () => {
  const input = mount(new RecordsInput({columns: [{name: 'name', input: text}], orderable: true,
    value: [{name: 'Design'}, {name: 'Synthesis'}, {name: 'Assay'}]}));
  assert.equal(rows(input)[0].querySelectorAll('[data-u2-part="handle"]').length, 1);
  const assay = rows(input)[2];
  input.move(2, 0);
  assert.deepEqual(input.value.map((r) => r.name), ['Assay', 'Design', 'Synthesis']);
  assert.ok(rows(input)[0] === assay, 'the row moved, not rebuilt');
  const editor = cell(input, 0, 0).querySelector('input');
  editor.focus();
  fire(editor, 'keydown', {key: 'ArrowDown', altKey: true});
  assert.deepEqual(input.value.map((r) => r.name), ['Design', 'Assay', 'Synthesis']);
  await flush();
  assert.ok(document.activeElement === editor, 'the moved row keeps the focus');
  input.dispose();
});

smoke('a value written from elsewhere reloads the rows; the items are keyed by the first field', () => {
  const input = mount(fields({value: [{name: 'Chemist', type: 'String'}]}));
  input.value = [{name: 'Target', type: 'String'}, {name: 'Dose', type: 'Number'}];
  assert.deepEqual(input.getItems().list.map((i) => i.key), ['Target', 'Dose']);
  assert.deepEqual(rows(input).map((r) => r.dataset.u2Item), ['Target', 'Dose'], 'the DOM keys are the items\' keys');
  type(cell(input, 1, 0).querySelector('input'), 'Volume');
  assert.equal(rows(input)[1].dataset.u2Item, 'Volume');
  assert.equal(input.getState().kind, 'records-input');
  assert.equal(input.getState().items.count, 2);
  input.reset([{name: 'Only', type: 'String'}]);
  assert.equal(rows(input).length, 1);
  input.dispose();
});

smoke('read-only: every cell refuses edits, rows added later too', () => {
  const input = mount(fields({value: [{name: 'Chemist', type: 'String'}], readOnly: true}));
  assert.equal(cell(input, 0, 0).querySelector('input').readOnly, true);
  input.add();
  assert.equal(cell(input, 1, 0).querySelector('input').readOnly, true);
  input.dispose();
});

smoke('in a form: the form goes by the records\' verdict', () => {
  const input = fields({value: [{name: '', type: 'String'}]});
  const form = new Form().add(input);
  mount(form);
  assert.match(form.validity.value, /Row 1: Name: Value can't be empty/);
  form.dispose();
  input.dispose();
});
