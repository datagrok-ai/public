// EQ-060: the built bundle's DomainFrameEditor.buildOps and resolveConflict for a stale delete and restore, run on
// an editor whose frame, client and access are stubs (the Dart DataFrame is not available in Node).
const test = require('node:test');
const assert = require('node:assert/strict');
const {loadBundle} = require('../load-bundle.cjs');

const {DG, sandbox} = loadBundle(process.env.JS_API_BUNDLE || undefined);
const E = DG.DomainFrameEditor;
// the bundle runs in its own realm: compare its objects by value
const plain = (v) => JSON.parse(JSON.stringify(v ?? null));

function frame(cols) {
  const column = (name) => ({
    get: (r) => cols[name][r],
    isNone: (r) => cols[name][r] == null || cols[name][r] === '',
    set: (r, v) => { cols[name][r] = v; },
  });
  // requestFilter recomputes the filter at once, letting [onRowsFiltering] AND in its exclusions
  const bits = cols.id.map(() => true);
  const df = {
    get rowCount() { return cols.id.length; },
    columns: {byName: (n) => n in cols ? column(n) : null, contains: (n) => n in cols},
    set: (n, r, v) => { cols[n][r] = v; },
    filter: {get: (r) => bits[r], set: (r, v) => { bits[r] = v; }},
    rows: {requestFilter: () => { bits.fill(true); df.onRowsFiltering?.(); }},
    onRowsFiltering: null,
    cols,
  };
  return df;
}

/** An editor over one persisted row in [state], its client answering get() with [fresh] and recording the options. */
function editor(state, fresh, {concurrency = 'version'} = {}) {
  const ed = Object.create(E.prototype);
  const reads = [];
  Object.assign(ed, {
    _df: frame({id: ['r1'], name: ['Mine'], version: [1], [DG.DOMAIN_DELETED_COLUMN]: [state === 'restored'],
      [E.STATE_COLUMN]: [state], [E.CHANGES_COLUMN]: [''], [E.ERRORS_COLUMN]: ['']}),
    _propByName: new Map(), _parsed: new Map(), _snapshots: new Map(), _unguarded: new Set(), _pristine: new Set(),
    _info: {nameColumn: 'name'},
    access: {fields: {name: 'editable'}, support: {concurrency}},
    client: {table: 'item', get: async (id, options) => { reads.push({id, options}); return fresh; }},
    _onConflict: {next() {}}, _onChanged: {next() {}}, _onDirtyChanged: {next() {}},
    _changeCount: 1, _dirty: true,
  });
  ed._recount = () => {};
  ed._df.onRowsFiltering = () => ed._maskDeleted();
  return {ed, reads};
}

const conflict = (current) => new DG.DomainVersionConflictError(`Version conflict: expected 1, current ${current}`, 409,
  {error: 'version-conflict', id: 'r1', currentVersion: current, expectedVersion: 1, opIndex: 0});

/** The dialog the platform shows, answered with [decision]; records what it was asked. */
function dialog(decision, {showFor = true} = {}) {
  const asked = [];
  sandbox.grok_DomainConflictDialog_Show = (subject) => { asked.push({subject}); return Promise.resolve(decision); };
  sandbox.grok_DomainConflictDialog_ShowFor = showFor ?
    (subject, op) => { asked.push({subject, op}); return Promise.resolve(decision); } : undefined;
  return asked;
}

for (const op of ['delete', 'restore']) {
  const state = op === 'delete' ? 'deleted' : 'restored';

  test(`a ${op} op carries the row's expectedVersion on a version-guarded table, none elsewhere`, () => {
    assert.deepEqual(plain(editor(state, null).ed.buildOps().map((p) => p.op)),
      [{op, table: 'item', id: 'r1', expectedVersion: 1}]);
    assert.deepEqual(plain(editor(state, null, {concurrency: 'none'}).ed.buildOps().map((p) => p.op)),
      [{op, table: 'item', id: 'r1'}]);
  });

  test(`a stale ${op}, RELOAD: the dialog is worded for a ${op}, the row is re-read and the op dropped`, async () => {
    const asked = dialog('reload');
    const {ed, reads} = editor(state, {id: 'r1', name: 'Theirs', version: 2});
    const failing = ed.buildOps()[0];
    assert.equal(await ed.resolveConflict(conflict(2), failing), true);
    assert.deepEqual(plain(asked), [{subject: 'Mine', op}]);
    assert.deepEqual(plain(reads), [op === 'restore' ? {id: 'r1', options: {deleted: 'include'}} : {id: 'r1'}],
      'a restore reads the trashed row with deleted: include');
    assert.equal(ed.stateOf(0), '', `the ${op} is no longer staged`);
    assert.equal(ed.buildOps().length, 0, 'no op is left');
    assert.equal(ed._df.cols.name[0], 'Theirs');
    assert.equal(ed._df.cols.version[0], 2);
  });

  test(`a stale ${op}, OVERWRITE: the ${op} is resent at the current version`, async () => {
    dialog('overwrite');
    const {ed, reads} = editor(state, null);
    assert.equal(await ed.resolveConflict(conflict(5), ed.buildOps()[0]), true);
    assert.deepEqual(reads, [], 'overwrite reads nothing');
    assert.deepEqual(plain(ed.buildOps().map((p) => p.op)), [{op, table: 'item', id: 'r1', expectedVersion: 5}]);
  });
}

test('a server without the delete wording shows the edit dialog', async () => {
  const asked = dialog('reload', {showFor: false});
  const {ed} = editor('deleted', {id: 'r1', name: 'Theirs', version: 2});
  assert.equal(await ed.resolveConflict(conflict(2), ed.buildOps()[0]), true);
  assert.deepEqual(plain(asked), [{subject: 'Mine'}]);
  assert.equal(ed.stateOf(0), '');
});

test('a stale delete, RELOAD: the row is shown again before the retry', async () => {
  dialog('reload');
  const {ed} = editor('deleted', {id: 'r1', name: 'Theirs', version: 2});
  ed._df.rows.requestFilter();
  assert.equal(ed._df.filter.get(0), false, 'the staged delete hides the row');
  assert.equal(await ed.resolveConflict(conflict(2), ed.buildOps()[0]), true);
  assert.equal(ed._df.filter.get(0), true, 'the dropped delete no longer hides the row, whatever the retry does');
});

test('a stale update, RELOAD: json and string_list values land as the JSON text the frame holds', async () => {
  dialog('reload');
  const {ed} = editor('modified', {id: 'r1', name: 'Theirs', version: 2, spec: {a: 1}, tags: ['x', 'y']});
  Object.assign(ed._df.cols, {spec: ['{}'], tags: ['[]']});
  assert.equal(await ed.resolveConflict(conflict(2), {op: {op: 'update', table: 'item', id: 'r1'}, row: 0}), true);
  assert.equal(ed._df.cols.spec[0], '{"a":1}');
  assert.equal(ed._df.cols.tags[0], '["x","y"]');
});
