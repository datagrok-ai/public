// The built bundle's DomainFrameEditor.buildOps and resolveConflict, run on an editor whose frame, client and
// access are stubs (the Dart DataFrame is not available in Node).
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
  const bits = cols.id.map(() => true);
  const df = {
    get rowCount() { return cols.id.length; },
    columns: {byName: (n) => n in cols ? column(n) : null, contains: (n) => n in cols},
    set: (n, r, v) => { cols[n][r] = v; },
    filter: {get: (r) => bits[r], set: (r, v) => { bits[r] = v; }},
    rows: {
      requestFilter: () => { bits.fill(true); df.onRowsFiltering?.(); },
      removeAt: (r, n) => { for (const c of Object.values(cols)) c.splice(r, n); bits.splice(r, n); },
    },
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
    _propByName: new Map(), _parsed: new Map(), _contributions: new Map(), _snapshots: new Map(), _unguarded: new Set(), _pristine: new Set(),
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

/** The dialog the platform shows, answered with [decision]. */
function dialog(decision) {
  sandbox.grok_DomainConflictDialog_Show = sandbox.grok_DomainConflictDialog_ShowFor = () => Promise.resolve(decision);
}

for (const op of ['delete', 'restore']) {
  const state = op === 'delete' ? 'deleted' : 'restored';

  test(`a ${op} op carries the row's expectedVersion on a version-guarded table, none elsewhere`, () => {
    assert.deepEqual(plain(editor(state, null).ed.buildOps().map((p) => p.op)),
      [{op, table: 'item', id: 'r1', expectedVersion: 1}]);
    assert.deepEqual(plain(editor(state, null, {concurrency: 'none'}).ed.buildOps().map((p) => p.op)),
      [{op, table: 'item', id: 'r1'}]);
  });

  test(`a stale ${op}, RELOAD: the row is re-read and the op dropped`, async () => {
    dialog('reload');
    const {ed, reads} = editor(state, {id: 'r1', name: 'Theirs', version: 2});
    assert.equal(await ed.resolveConflict(conflict(2), ed.buildOps()[0]), true);
    assert.deepEqual(plain(reads), [op === 'restore' ? {id: 'r1', options: {deleted: 'include'}} : {id: 'r1'}],
      'a restore reads the trashed row with deleted: include');
    assert.equal(ed.buildOps().length, 0);
    assert.deepEqual([ed._df.cols.name[0], ed._df.cols.version[0]], ['Theirs', 2]);
  });

  test(`a stale ${op}, OVERWRITE: the ${op} is resent at the current version`, async () => {
    dialog('overwrite');
    const {ed, reads} = editor(state, null);
    assert.equal(await ed.resolveConflict(conflict(5), ed.buildOps()[0]), true);
    assert.deepEqual(reads, [], 'overwrite reads nothing');
    assert.deepEqual(plain(ed.buildOps().map((p) => p.op)), [{op, table: 'item', id: 'r1', expectedVersion: 5}]);
  });
}

test('a stale delete, RELOAD: the row is shown again before the retry', async () => {
  dialog('reload');
  const {ed} = editor('deleted', {id: 'r1', name: 'Theirs', version: 2});
  ed._df.rows.requestFilter();
  assert.equal(ed._df.filter.get(0), false, 'the staged delete hides the row');
  assert.equal(await ed.resolveConflict(conflict(2), ed.buildOps()[0]), true);
  assert.equal(ed._df.filter.get(0), true);
});

test('a stale delete, RELOAD of a row the server no longer has: the row leaves the frame, as a saved delete', async () => {
  dialog('reload');
  const balloons = [];
  Object.assign(sandbox, {grok_JS_To_Map: (o) => o, grok_Balloon: (s) => { balloons.push(s); }});
  const {ed} = editor('deleted', null);
  assert.equal(await ed.resolveConflict(conflict(2), ed.buildOps()[0]), true);
  assert.deepEqual(balloons, ['Mine no longer exists.']);
  assert.equal(ed._df.rowCount, 0, 'the row is not shown again as a live one');
  assert.equal(ed.buildOps().length, 0);
});

test('a stale update, RELOAD: json and string_list values land as the JSON text the frame holds', async () => {
  dialog('reload');
  const {ed} = editor('modified', {id: 'r1', name: 'Theirs', version: 2, spec: {a: 1}, tags: ['x', 'y']});
  Object.assign(ed._df.cols, {spec: ['{}'], tags: ['[]']});
  assert.equal(await ed.resolveConflict(conflict(2), {op: {op: 'update', table: 'item', id: 'r1'}, row: 0}), true);
  assert.equal(ed._df.cols.spec[0], '{"a":1}');
  assert.equal(ed._df.cols.tags[0], '["x","y"]');
});

test('a text starting with ~new: is sent as typed, and only a draft id of the batch as a $ref', () => {
  const draft = '~new:7b1f8c3e-5d2a-4e6b-9c1d-0a2b3c4d5e6f';
  const {ed} = editor('modified', null, {concurrency: 'none'});
  Object.assign(ed.access.fields, {parent_id: 'editable'});
  Object.assign(ed._df.cols, {id: ['r1', 'r2', draft, '~new:9a8b7c6d-5e4f-4a3b-8c2d-1e0f9a8b7c6d'],
    name: ['~new:release notes', '$x', 'parent', 'child'], parent_id: [null, null, null, draft], version: [1, 1, null, null],
    [DG.DOMAIN_DELETED_COLUMN]: [false, false, false, false], [E.STATE_COLUMN]: ['modified', 'modified', 'new', 'new'],
    [E.CHANGES_COLUMN]: ['{"name":"old"}', '{"name":"old"}', '', ''], [E.ERRORS_COLUMN]: ['', '', '', '']});
  const values = ed.buildOps().map((p) => p.op.values);
  assert.deepEqual(plain([values[0].name, values[1].name, values[3].parent_id]), ['~new:release notes', '$$x', `$${draft}`]);
});

test('a draft added and then deleted is not a draft of the batch', () => {
  const {ed} = editor('deleted', null);
  ed._df.cols.id[0] = '~new:7b1f8c3e-5d2a-4e6b-9c1d-0a2b3c4d5e6f';
  assert.equal(ed.draftIds().size, 0);
});

test('a reference to a deleted parent is refused on its cell, by the column caption and without the id', async () => {
  const {ed} = editor('restored', null);
  ed._df.cols.location_id = ['p1'];
  ed._propByName.set('location_id', {friendlyName: 'Location'});
  const refusal = 'Column "location_id" references a deleted row in "location" ("p1")';
  const message = await ed.refusalFor(new DG.DomainRestrictError(refusal, 409,
    {error: 'restrict', id: 'r1', column: 'location_id', table: 'location', message: refusal}),
    {op: {op: 'restore'}, row: 0});
  assert.equal(message, 'Location refers to a deleted location: restore that location first');
  assert.equal(ed.errorsOf(0).location_id.message, message);
  ed.unmarkRestored(0);
  assert.equal(ed.errorOf(0, 'location_id'), null, 'an unstaged restore keeps no mark');
  sandbox.grok_DomainRegistry_TableInfo = () => Promise.resolve({singularName: 'Location', nameColumn: 'name'});
  sandbox.grok_Dapi_Domains = () => ({});
  sandbox.grok_Dapi_Domains_GetRow = () => Promise.resolve({id: 'p1', name: 'Shelf 2'});
  const peer = await ed.refusalFor(new DG.DomainRestrictError(refusal, 409, {error: 'restrict', id: 'r1',
    column: 'location_id', table: 'location', parentId: 'p1', message: refusal}), {op: {op: 'update'}, row: 0});
  assert.equal(peer, 'Location "Shelf 2" of "Mine" was deleted: pick another location');
  assert.equal(ed.errorsOf(0).location_id.kind, 'conflict');
});

test('a delete of a row someone else already deleted lands, and the rest of the batch is sent', async () => {
  const {ed} = editor('deleted', null);
  for (const [name, value] of [['id', 'r2'], ['name', 'Other'], ['version', 1], [DG.DOMAIN_DELETED_COLUMN, false],
    [E.STATE_COLUMN, 'deleted'], [E.CHANGES_COLUMN, ''], [E.ERRORS_COLUMN, '']])
    ed._df.cols[name].push(value);
  Object.assign(ed, {_onSavingChanged: {next() {}}, _onSaved: {next() {}}, writeBack: async () => {}, rebind() {}});
  ed.client.schema = 's';
  const sent = [];
  sandbox.grok_Dapi_Domains = () => ({});
  sandbox.grok_Dapi_Domains_Transaction = (_, schema, ops) => {
    sent.push(plain(ops).map((op) => op.id));
    return sent.length > 1 ? Promise.resolve([{id: 'r2'}]) : Promise.reject(new DG.DomainNotFoundError(
      'Row "r1" not found in s.item', 404, {error: 'not-found', id: 'r1', opIndex: 0}));
  };
  const session = Object.assign(Object.create(DG.DomainSession.prototype),
    {editors: [ed], _quiet: true, _onSaved: {next() {}}, _lastRefusal: null});
  assert.equal(await session.save(), true);
  assert.deepEqual(sent, [['r1', 'r2'], ['r2']]);
  assert.equal(ed._df.rowCount, 0, 'both rows left the frame');
});
