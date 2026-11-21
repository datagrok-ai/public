// EQ-006: the built bundle's DomainFrameEditor.buildOps sends a text that only starts with `~new:` as typed: only
// the id of a draft row of the batch becomes a `$` reference; a literal leading `$` keeps its `$$` escape. Run on an
// editor whose frame, client and access are stubs (the Dart DataFrame is not available in Node).
const test = require('node:test');
const assert = require('node:assert/strict');
const {loadBundle} = require('../load-bundle.cjs');

const {DG} = loadBundle(process.env.JS_API_BUNDLE || undefined);
const E = DG.DomainFrameEditor;
// the bundle runs in its own realm: compare its objects by value
const plain = (v) => JSON.parse(JSON.stringify(v ?? null));

function frame(cols) {
  const column = (name) => ({
    get: (r) => cols[name][r],
    isNone: (r) => cols[name][r] == null || cols[name][r] === '',
  });
  return {
    get rowCount() { return cols.id.length; },
    columns: {byName: (n) => n in cols ? column(n) : null, contains: (n) => n in cols},
  };
}

/** An editor over [rows] ({id, state, note, parent_id, changes}) of a table with writable `note` and `parent_id`. */
function editor(rows) {
  const ed = Object.create(E.prototype);
  const col = (name) => rows.map((r) => r[name] ?? null);
  Object.assign(ed, {
    _df: frame({id: col('id'), note: col('note'), parent_id: col('parent_id'),
      [E.STATE_COLUMN]: col('state'), [E.CHANGES_COLUMN]: rows.map((r) => JSON.stringify(r.changes ?? {})),
      [E.ERRORS_COLUMN]: rows.map(() => '')}),
    _propByName: new Map(), _parsed: new Map(), _snapshots: new Map(), _unguarded: new Set(), _pristine: new Set(),
    _info: {nameColumn: 'note'},
    access: {fields: {note: 'editable', parent_id: 'editable'}, support: {concurrency: 'none'}},
    client: {table: 'item'},
  });
  return ed;
}

const DRAFT = '~new:7b1f8c3e-5d2a-4e6b-9c1d-0a2b3c4d5e6f';
const OTHER = '~new:0f9e8d7c-6b5a-4c3d-8e2f-1a0b9c8d7e6f';

test('a text that only starts with ~new: is sent as typed; a draft id of the batch is a $ref', () => {
  const ed = editor([
    {id: 'r1', state: 'modified', note: '~new:release notes', changes: {note: 'old'}},
    {id: 'r2', state: 'modified', note: '$x', changes: {note: 'old'}},
    {id: DRAFT, state: 'new', note: '~new:draft title'},
    {id: '~new:9a8b7c6d-5e4f-4a3b-8c2d-1e0f9a8b7c6d', state: 'new', note: 'child', parent_id: DRAFT},
  ]);
  const values = plain(ed.buildOps().map((p) => p.op.values));
  assert.equal(values[0].note, '~new:release notes', 'a ~new: text is sent as typed, not as a $ref');
  assert.equal(values[1].note, '$$x', 'a literal $ keeps its $$ escape');
  assert.equal(values[2].note, '~new:draft title', 'and in an insert');
  assert.equal(values[3].parent_id, `$${DRAFT}`, 'a draft id of the batch is the server\'s $ref');
});

test('a draft id of another participant of the batch is a $ref', () => {
  const child = editor([{id: '~new:9a8b7c6d-5e4f-4a3b-8c2d-1e0f9a8b7c6d', state: 'new', note: 'child', parent_id: OTHER}]);
  assert.equal(plain(child.buildOps(new Set([OTHER]))[0].op.values).parent_id, `$${OTHER}`);
});
