import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {Scope, FrameRows, Rows} from '@datagrok-libraries/u2';
import {DataFrame, Stream} from '@datagrok-libraries/u2/tests/platform-doubles.mjs';
import {register} from 'node:module';
import {EditorEditState} from '../src/dg/editor-state.js';
import {DomainSource} from '../src/core/domain-source.js';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const DG = await import('datagrok-api/dg');
const {DgDomainBackend} = await import('../src/dg/backend.js');
const {MemoryDomainBackend} = await import('../src/core/memory-domain.js');

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
    onChanged: new Stream(), onDirtyChanged: new Stream(), onSavingChanged: new Stream(), onSaved: new Stream(),
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
      df.dart.rows.push({...values, id: Rows.draftId()});
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
    markRestored(row) {
      this.calls.push(['restore', row]);
      states[row] = 'restored';
      this.onChanged.fire(this);
    },
    unmarkRestored(row) {
      this.calls.push(['unrestore', row]);
      states[row] = '';
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
  state.markRestored('i2');
  assert.deepEqual(editor.calls[2], ['restore', 1], 'the adapter reaches the editor by index');
  state.unmarkRestored('i2');
  assert.deepEqual(editor.calls[3], ['unrestore', 1]);
  state.markRestored('nope');
  assert.equal(editor.calls.length, 4, 'an unknown key reaches nothing');
  state.discard();
  assert.deepEqual(editor.calls[4], ['discard']);
  assert.equal(await state.save(), true);
  state.dispose();
  assert.equal(editor.detached, 1);
  assert.equal(editor.onChanged.count + editor.onDirtyChanged.count + editor.onSavingChanged.count +
    editor.onSaved.count, 0, 'dispose releases the editor subscriptions');
  assert.equal(editor.dataFrame.liveSubscriptions(), 0, 'and the frame ones');
});

scoped('a draft is keyed by the ~new: id the editor stamps; onSaved carries the assigned ids', async () => {
  const {df, editor} = fakeEditor(ROWS);
  const state = new EditorEditState(editor);
  const key = state.newRow({title: 'Draft'}, {pristine: true});
  assert.equal(Rows.isDraft(key), true);
  assert.equal(state.indexOf(key), 2);
  assert.equal(state.changeCount.value, 0, 'pristine');
  assert.equal(state.isChanged(key, 'title'), true, 'every cell of a new row is a change');
  state.setValue(key, 'title', 'Draft 2');
  assert.equal(df.get('title', 2), 'Draft 2');
  assert.equal(FrameRows.keyOf(df, 0), 'i1');
  assert.equal(FrameRows.keyOf(df, 2), key);
  const heard = [];
  const sub = state.onSaved.subscribe((e) => heard.push(e));
  df.set('id', 2, 'i9');
  editor.onSaved.fire({assigned: {[key]: 'i9'}});
  assert.deepEqual(heard, [{assigned: {[key]: 'i9'}}]);
  sub.unsubscribe();
  assert.equal(state.indexOf('i9'), 2, 'the index follows the frame\'s cells');
  assert.equal(state.indexOf(key), -1, 'the draft id left with the cell');
  assert.equal(state.indexOf('~row:9'), -1, 'no index-keyed path');
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
  editor.markRestored(1);
  assert.equal(state.validity.value, null, 'nor a restored one\'s');
  editor.fail(0, 'title', 'Someone else changed it', 'conflict');
  assert.equal(state.validity.value, null, 'a dismissed conflict only marks the cell');
  assert.equal(state.errorOf('i1', 'title'), 'Someone else changed it');
  state.dispose();
});

// EQ-033: the platform editor marks the server's field error on the cell as kind 'error', like its
// own checks; saveAll keeps that mark apart, so it stays shown, holds nothing, and the next save is sent
scoped('EQ-033: a cell the server refused stays shown, does not hold Save, and the next save is sent', async () => {
  const {editor} = fakeEditor(ROWS);
  let validated = 0;
  editor.commitEdit = function(row, column) { // the editor re-runs the cell's own checks: the server's mark goes
    validated++;
    delete this.errorsOf(row)[column];
    this.onChanged.fire(this);
  };
  const state = new EditorEditState(editor);
  state.setValue('i1', 'title', 'Aspirin');
  const real = DG.DomainSession.prototype.save;
  let sent = 0;
  DG.DomainSession.prototype.save = async function() {
    this.lastRefusal = ++sent > 1 ? null : {message: 'Title Aspirin already belongs to i7', editor: this.editors[0]};
    if (this.lastRefusal) // the platform maps the server's field error onto the cell it names
      this.editors[0].fail(0, 'title', this.lastRefusal.message);
    return this.lastRefusal === null;
  };
  try {
    const dg = new DgDomainBackend();
    assert.equal(await dg.saveAll([state]), false);
    assert.equal(state.validity.value, 'Title Aspirin already belongs to i7', 'the refusal stays shown');
    assert.equal((state.blocking ?? state.validity).value, null, 'a server refusal does not hold Save');
    assert.equal(await dg.saveAll([state]), true);
    assert.equal(validated, 1, 'the next save took the server\'s mark back before the editor built its batch');
    assert.equal(sent, 2);
  } finally {
    DG.DomainSession.prototype.save = real;
  }
  state.dispose();
});

scoped('EQ-033: a client error the retry finds is the client\'s, and holds', async () => {
  const {editor} = fakeEditor(ROWS);
  editor.commitEdit = function(row, column) { // the server's mark goes; the client finds the title empty
    this.errorsOf(row)[column] = {message: 'Title is required', kind: 'error'};
    this.onChanged.fire(this);
  };
  const state = new EditorEditState(editor);
  state.setValue('i1', 'title', '');
  const real = DG.DomainSession.prototype.save;
  let sent = 0;
  DG.DomainSession.prototype.save = async function() {
    if (++sent === 1)
      this.editors[0].fail(0, 'title', 'Taken');
    // then the editor refuses its own cells and sends nothing
    this.lastRefusal = sent === 1 ? {message: 'Taken', editor: this.editors[0]} :
      {message: 'Fix the cell errors first', editor: null};
    return false;
  };
  try {
    const dg = new DgDomainBackend();
    assert.equal(await dg.saveAll([state]), false);
    assert.equal(await dg.saveAll([state]), false);
    assert.equal(state.refused('i1', 'title'), false, 'never taken for the server\'s');
    assert.equal(state.blocking.value, 'Title is required', 'the client\'s error holds the batch');
  } finally {
    DG.DomainSession.prototype.save = real;
  }
  state.dispose();
});

scoped('a save after a server refusal is not held by a stored value of the row that the user never touched',
  async () => {
    const {editor} = fakeEditor(ROWS);
    // as the platform's: every rule over every cell of the pending rows, the stored empty status included
    editor.validate = function() {
      delete this.errorsOf(0).title;
      this.errorsOf(0).status = {message: 'Value can\'t be empty', kind: 'error'};
      this.onChanged.fire(this);
      return 1;
    };
    editor.commitEdit = function(row, column) { // the rules of that one cell: the typed title passes
      delete this.errorsOf(row)[column];
      this.onChanged.fire(this);
    };
    const state = new EditorEditState(editor);
    state.setValue('i1', 'title', 'Aspirin 2');
    const real = DG.DomainSession.prototype.save;
    let sent = 0;
    DG.DomainSession.prototype.save = async function() {
      if (Object.values(editor.errorsOf(0)).some((e) => e.kind === 'error'))
        return (this.lastRefusal = {message: DG.DomainSession.CELL_ERRORS, editor: null}, false);
      this.lastRefusal = ++sent > 1 ? null : {message: 'Title Aspirin 2 already belongs to i7', editor: this.editors[0]};
      if (this.lastRefusal)
        this.editors[0].fail(0, 'title', this.lastRefusal.message);
      return this.lastRefusal === null;
    };
    try {
      const dg = new DgDomainBackend();
      assert.equal(await dg.saveAll([state]), false);
      assert.equal(await dg.saveAll([state]), true, 'the second Save is sent');
      assert.equal(state.errorOf('i1', 'status'), null, 'the untouched status is not marked');
    } finally {
      DG.DomainSession.prototype.save = real;
    }
    state.dispose();
  });

scoped('JSON-looking text retyped in a plain string column is written as typed, by both writers', async () => {
  const {df, editor} = fakeEditor([{id: 'i1', title: '{"a":1}'}]);
  editor.properties = [{name: 'title', propertyType: 'string'}];
  editor.changesOf = () => ({});
  const platform = new EditorEditState(editor);
  platform.setValue('i1', 'title', '{"a": 1}');
  assert.equal(df.get('title', 0), '{"a": 1}', 'the platform writer writes the text');
  platform.dispose();
  const memory = new MemoryDomainBackend({name: 'qa', tables: {item: {columns: {
    notes: {type: 'string'}, meta: {type: 'json'}}}}}, {rows: {item: [{id: 'r1', notes: '{"a":1}', meta: {a: 1}}]}});
  const frame = await memory.tableSync('qa.item').frame({});
  frame.edit.setValue('r1', 'notes', '{"a": 1}');
  assert.equal(frame.df.get('notes', 0), '{"a": 1}', 'the memory writer writes the text');
  frame.edit.setValue('r1', 'meta', '{"a": 1}');
  assert.equal(frame.edit.isChanged('r1', 'meta'), false, 'a json cell retyped as it was is no change');
  assert.equal(frame.edit.changeCount.value, 1, 'only the string column counts');
  frame.dispose();
});

/* EQ-057: a json cell written through the platform editor. The double behaves as js-api does where it
   matters: `toDart` makes a Dart map of a PLAIN object only (anything else reaches the frame as it is),
   a Dart map is cyclic through its type information, the editor's change compare is `JSON.stringify`
   (`wireEquals`), and `validateCellValue` is the server's js-api — before 1.27.14 it refuses a json
   cell's text, after it leaves the text to the server. */
function jsApiEditor(acceptsText) {
  const df = new DataFrame([{name: 'id', type: 'string'}, {name: 'meta', type: 'string'}],
    [{id: 'r1', meta: '{"k":"v"}'}]);
  const changes = [{}];
  const toDart = (x) => {
    if (typeof x !== 'object' || x === null || Object.getPrototypeOf(x) !== Object.prototype)
      return x;
    const map = Object.assign(Object.create({constructor: function LinkedHashMapCell() {}}), x);
    map.$ti = {owner: map};
    return map;
  };
  const wireEquals = (a, b) => a === b || JSON.stringify(a) === JSON.stringify(b);
  class Editor {
    static validateCellValue(p, value) {
      if (value === null || value === undefined || (acceptsText && typeof value === 'string'))
        return null;
      return typeof value === 'object' && !Array.isArray(value) ? null :
        `Types differ. Expected: map, passed: ${Array.isArray(value) ? 'list' : 'string'}`;
    }
    constructor() {
      this.table = 'qa.item';
      this.properties = [{name: 'meta', propertyType: 'map', nullable: true}];
      this.isDirty = false; this.changeCount = 0; this.isSaving = false;
      this.onChanged = new Stream(); this.onDirtyChanged = new Stream(); this.onSavingChanged = new Stream();
      this.onSaved = new Stream();
    }
    get dataFrame() { return df; }
    stateOf(row) { return Object.keys(changes[row]).length > 0 ? 'modified' : ''; }
    changesOf(row) { return changes[row]; }
    errorsOf() { return {}; }
    errorOf(row, column) {
      const message = Object.hasOwn(changes[row], column) ? Editor.validateCellValue(this.properties[0], df.get(column, row)) : null;
      return message === null ? null : {message, kind: 'error'};
    }
    isChanged(row, column) { return Object.hasOwn(changes[row], column); }
    setValue(row, column, value) {
      const original = Object.hasOwn(changes[row], column) ? undefined : df.get(column, row);
      df.set(column, row, toDart(value));
      if (original !== undefined)
        changes[row][column] = original;
      if (wireEquals(changes[row][column], df.get(column, row)))
        delete changes[row][column];
      this.onChanged.fire(this);
    }
    revertCell(row, column) {
      df.set(column, row, changes[row][column]);
      delete changes[row][column];
      this.onChanged.fire(this);
    }
    detach() {}
  }
  return {df, editor: new Editor()};
}

for (const acceptsText of [true, false]) {
  const server = acceptsText ? 'a js-api that takes json text' : 'a js-api that refuses json text';
  scoped(`EQ-057: a typed json object reaches the platform frame as something it can compare and send (${server})`, async () => {
    const {df, editor} = jsApiEditor(acceptsText);
    const state = new EditorEditState(editor);
    const write = (value) => state.setValue('r1', 'meta', DomainSource.frameValue(df, 'meta', value));
    assert.doesNotThrow(() => write({k: 'v', k2: 'v2'}),
      'writing a typed json object into the frame does not throw (no Dart map in the cell)');
    const cell = df.get('meta', 0);
    assert.equal(typeof cell === 'string', acceptsText, acceptsText ? 'as its JSON text' : 'as an object, the fallback');
    assert.deepEqual(JSON.parse(acceptsText ? cell : JSON.stringify(cell)), {k: 'v', k2: 'v2'});
    assert.equal(state.errorOf('r1', 'meta'), null, 'and the editor takes it');
    write({k: 'v'});
    assert.equal(state.isChanged('r1', 'meta'), false, 'retyping the original pair is no change');
    state.dispose();
  });
}
