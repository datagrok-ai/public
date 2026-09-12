/* `EditorEditState` (WO-8) over a fake js-api editor and the DataFrame double: keys map onto row
   indices (the id cell, `~row:<index>` for a draft), every write reaches the editor, the signals
   follow the editor's observables, validity is the first blocking error on a live row, and
   dispose detaches. The real editor runs only in the platform — the U2Demo `U2: domain source`
   category covers that. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {DataFrame, Stream} from './platform-doubles.mjs';
import {EditorEditState} from '../src/dg/domain/editor-state.js';
import {FrameRows} from '../src/sources/df-rows.js';
import {Rows} from '../src/sources/rows-like.js';

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

/** The editor's surface the host reads and writes, over a frame double. */
function fakeEditor(rows) {
  const df = new DataFrame([{name: 'id', type: 'string'}, {name: 'title', type: 'string'}], rows);
  const states = rows.map(() => '');
  const errors = rows.map(() => ({}));
  const editor = {
    table: 'grit.issue',
    isDirty: false, changeCount: 0, isSaving: false, detached: 0, calls: [],
    onChanged: new Stream(), onDirtyChanged: new Stream(), onSavingChanged: new Stream(),
    get dataFrame() { return df; },
    stateOf: (row) => states[row] ?? '',
    errorsOf: (row) => errors[row] ?? {},
    isChanged(row, column) { return states[row] === 'new' || this.calls.some((c) => c[0] === 'set' && c[1] === row && c[2] === column); },
    errorOf(row, column) { return errors[row]?.[column] ?? null; },
    setValue(row, column, value) {
      this.calls.push(['set', row, column, value]);
      df.set(column, row, value);
      this.changeCount++;
      this.isDirty = true;
      this.onChanged.fire(this);
    },
    addRow(values, options) {
      if (this.isSaving)
        return -1;
      df.dart.rows.push({id: null, ...values});
      df.onRowsAdded.fire();
      states.push('new');
      errors.push({});
      if (!options?.pristine)
        this.changeCount++;
      this.onChanged.fire(this);
      return df.rowCount - 1;
    },
    markDeleted(row) {
      this.calls.push(['delete', row]);
      states[row] = 'deleted';
      this.onChanged.fire(this);
    },
    discard() { this.calls.push(['discard']); },
    save: async () => true,
    detach() { this.detached++; },
    fail(row, column, message, kind = 'error') {
      errors[row][column] = {message, kind};
      this.onChanged.fire(this);
    },
  };
  return {df, editor};
}

const ROWS = [{id: 'i1', title: 'Aspirin'}, {id: 'i2', title: 'Ibuprofen'}];

scoped('keys are the id cells; writes, deletes and state reads go to the editor by row index', async () => {
  const {editor} = fakeEditor(ROWS);
  const state = new EditorEditState(editor);
  assert.equal(state.indexOf('i2'), 1);
  assert.equal(state.indexOf('nope'), -1);
  state.setValue('i2', 'title', 'Ibu');
  assert.deepEqual(editor.calls, [['set', 1, 'title', 'Ibu']]);
  assert.equal(state.isChanged('i2', 'title'), true);
  assert.equal(state.isChanged('i1', 'title'), false);
  state.setValue('nope', 'title', 'x');
  assert.equal(editor.calls.length, 1, 'an unknown key writes nothing');
  state.markDeleted('i1');
  assert.deepEqual(editor.calls[1], ['delete', 0]);
  state.discard();
  assert.deepEqual(editor.calls[2], ['discard']);
  assert.equal(await state.save(), true);
  state.dispose();
  assert.equal(editor.detached, 1);
  assert.equal(editor.onChanged.count + editor.onDirtyChanged.count + editor.onSavingChanged.count, 0,
    'dispose releases the editor subscriptions');
  assert.equal(editor.dataFrame.liveSubscriptions(), 0, 'and the frame ones');
});

scoped('a draft is keyed by its index until it has an id; the key is what `newRow` answers', async () => {
  const {df, editor} = fakeEditor(ROWS);
  const state = new EditorEditState(editor);
  const key = state.newRow({title: 'Draft'}, {pristine: true});
  assert.equal(key, '~row:2');
  assert.equal(state.indexOf(key), 2);
  assert.equal(state.changeCount.value, 0, 'pristine');
  assert.equal(state.isChanged(key, 'title'), true, 'every cell of a new row is a change');
  state.setValue(key, 'title', 'Draft 2');
  assert.equal(df.get('title', 2), 'Draft 2');
  assert.equal(FrameRows.keyOf(df, 0), 'i1');
  assert.equal(FrameRows.keyOf(df, 2), Rows.draftKey(2));
  df.set('id', 2, 'i9');
  assert.equal(state.indexOf('i9'), 2, 'the index follows the frame\'s cells');
  assert.equal(state.indexOf('~row:9'), -1, 'past the frame');
  editor.isSaving = true;
  assert.throws(() => state.newRow({}), /being saved/);
  state.dispose();
});

scoped('the signals follow the editor: change count, dirty, saving; listeners hear every change', async () => {
  const {editor} = fakeEditor(ROWS);
  const state = new EditorEditState(editor);
  const heard = [];
  const sub = state.onChanged.subscribe((key) => heard.push(key));
  assert.equal(state.isDirty.value, false);
  assert.equal(state.changeCount.value, 0);
  state.setValue('i1', 'title', 'A');
  assert.equal(state.changeCount.value, 1);
  assert.equal(state.isDirty.value, true, 'read off the editor on its onChanged');
  assert.deepEqual(heard, [null], 'the editor does not say which row: null');
  editor.isDirty = false;
  editor.onDirtyChanged.fire(false);
  assert.equal(state.isDirty.value, false);
  assert.equal(state.isSaving.value, false);
  editor.onSavingChanged.fire(true);
  assert.equal(state.isSaving.value, true);
  sub.unsubscribe();
  state.setValue('i1', 'title', 'B');
  assert.equal(heard.length, 1, 'unsubscribed');
  state.dispose();
});

scoped('validity is the first blocking error on a row that is not deleted; a conflict does not block', async () => {
  const {editor} = fakeEditor(ROWS);
  const walked = [];
  editor.errorsOf = ((inner) => (row) => walked.push(row) && inner(row))(editor.errorsOf);
  const state = new EditorEditState(editor);
  assert.equal(walked.length, 0, 'lazy: no row walk until validity is read');
  assert.equal(state.validity.value, null);
  assert.equal(walked.length, 2);
  editor.setValue(0, 'title', 'A');
  editor.setValue(0, 'title', 'B');
  assert.equal(walked.length, 2, 'two keystrokes, no walk');
  assert.equal(state.validity.value, null);
  assert.equal(walked.length, 4, 'one walk on the next read');
  editor.fail(1, 'title', 'Too long');
  assert.equal(state.validity.value, 'Too long');
  assert.equal(state.errorOf('i2', 'title'), 'Too long');
  assert.equal(state.errorOf('i1', 'title'), null);
  assert.equal(state.errorOf('nope', 'title'), null);
  editor.markDeleted(1);
  assert.equal(state.validity.value, null, 'a deleted row\'s errors do not block');
  editor.fail(0, 'title', 'Someone else changed it', 'conflict');
  assert.equal(state.validity.value, null, 'a dismissed conflict only marks the cell');
  assert.equal(state.errorOf('i1', 'title'), 'Someone else changed it');
  state.dispose();
});
