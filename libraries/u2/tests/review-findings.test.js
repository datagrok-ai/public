/* The seven findings of the phase-1 external review, reproduced against the phase-2 code. A `todo`
   test is one that FAILS today (the finding stands) — it asserts the behaviour the fix would give,
   and node reports it without failing the suite; a plain test is one that PASSES, recording what
   the code actually does where the review asked a question or got the shape wrong.

   Findings 3 and 5 run the REAL compiled js-api (`domains-editor.js`, `domains-session.js`) through
   a resolve hook: its ESM is webpack-shaped (extensionless relative imports, CSS imports), which
   node cannot load on its own. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {DomainSource} from '../src/sources/domain-source.js';
import {SharedSession} from '../src/sources/session.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {EditorEditState} from '../src/dg/domain/editor-state.js';
import {DataFrame, Stream} from './platform-doubles.mjs';
import {Rows} from '../src/sources/rows-like.js';
import {notify} from '../src/components/display/notify.js';
import {backend} from './domain-fixtures.mjs';

/** js-api's compiled output as node can read it: `./x` → `./x.js`, `rxjs/operators` →
 * `.../index.js`, a `.css` import → an empty module. */
const HOOK = 'data:text/javascript,' + encodeURIComponent(`
import {existsSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
export async function resolve(specifier, context, next) {
  try { return await next(specifier, context); }
  catch (e) {
    if (!specifier.startsWith('.') && !specifier.startsWith('/'))
      return next(specifier + '/index.js', context);
    const base = new URL(specifier, context.parentURL);
    for (const candidate of [base.href + '.js', base.href + '/index.js'])
      if (existsSync(fileURLToPath(candidate)))
        return {url: candidate, shortCircuit: true, format: 'module'};
    throw e;
  }
}
export async function load(url, context, next) {
  return url.endsWith('.css') ? {source: '', format: 'module', shortCircuit: true} : next(url, context);
}
`);
register(HOOK);

function scoped(name, options, body) {
  test(name, options, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const env = {designTime: false, subBinds: {}, resolve: () => null};

async function issues(options = {}) {
  const src = new DomainSource({table: 'grit.issue', ...options}, env);
  src.start();
  await flush();
  return src;
}

// ─────────────────────────── 1. a refresh drops edits made while it loaded ──────────────────────

/** The memory backend with a gate on `frame`: the test lets the response land when it chooses. */
function gated(inner) {
  const gates = [];
  return {backend: {table: async (address) => {
    const t = await inner.table(address);
    return {...t, access: () => t.access(), count: (f, s) => t.count(f, s), query: (s) => t.query(s),
      transaction: (ops) => t.transaction(ops),
      frame: async (spec) => {
        if (gates.length > 0)
          await gates.shift();
        return t.frame(spec);
      }};
  }}, hold: () => {
    let open;
    gates.push(new Promise((r) => open = r));
    return open;
  }};
}

scoped('finding 1: a query change while the response is in flight keeps the edits made meanwhile',
  {}, async () => {
    const g = gated(backend());
    backends.domain = g.backend;
    const src = await issues();
    src.rows.byKey('i1').title = 'Edited';
    assert.equal(src.isDirty.value, true);
    src.query.value = 'done = false';
    await flush();
    assert.equal(src.isDirty.value, true, 'a dirty source does not re-query at all (the start() gate)');
    src.edit.value.discard();
    assert.equal(src.isDirty.value, false);
    const open = g.hold();
    src.query.value = 'weight > 0';
    await flush();
    // the old frame is still the live one, and still editable, while the response is awaited
    src.rows.byKey('i1').title = 'Edited while loading';
    assert.equal(src.isDirty.value, true, 'the edit landed on the frame the refresh is about to drop');
    open();
    await flush();
    await flush();
    assert.equal(src.isDirty.value, true, 'the edit made during the load is still pending');
    src.dispose();
  });

// ───────────────────────── 2. the editor's row index after a saved deletion ─────────────────────

/** The js-api editor's surface as `applyResults` leaves it: a saved delete removes the row with
 * `removeAt(row, 1, false)` — no `onRowsRemoved` (row_list.dart:263, `_endUpdate(notify, …)`) —
 * and announces the batch through `onChanged` / `onSaved` only (domains-editor.ts:893-900). */
function fakeEditor(rows) {
  const df = new DataFrame([{name: 'id', type: 'string'}, {name: 'title', type: 'string'}], rows);
  const states = rows.map(() => '');
  return {df, editor: {
    table: 'grit.issue', isDirty: false, changeCount: 0, isSaving: false,
    onChanged: new Stream(), onDirtyChanged: new Stream(), onSavingChanged: new Stream(), onSaved: new Stream(),
    get dataFrame() { return df; },
    stateOf: (row) => states[row] ?? '',
    errorsOf: () => ({}),
    isChanged: () => false,
    errorOf: () => null,
    setValue(row, column, value) { df.set(column, row, value); },
    markDeleted(row) {
      states[row] = 'deleted';
      this.changeCount++;
      this.isDirty = true;
      this.onChanged.fire(this);
    },
    discard() {},
    detach() {},
    /** applyResults: the deleted rows leave the frame WITHOUT a row event, then onChanged/onSaved. */
    async save() {
      for (let row = df.rowCount - 1; row >= 0; row--)
        if (states[row] === 'deleted') {
          df.dart.rows.splice(row, 1);
          states.splice(row, 1);
        }
      this.changeCount = 0;
      this.isDirty = false;
      this.onChanged.fire(this);
      this.onSaved.fire({inserted: 0, updated: 0, deleted: 1, assigned: {}});
      return true;
    },
  }};
}

scoped('finding 2: the cached row index is rebuilt after a saved deletion', {}, async () => {
    const {editor} = fakeEditor([{id: 'i1', title: 'Aspirin'}, {id: 'i2', title: 'Ibuprofen'}]);
    const state = new EditorEditState(editor);
    assert.equal(state.indexOf('i1'), 0);
    assert.equal(state.indexOf('i2'), 1, 'the index is built and cached');
    state.markDeleted('i1');
    assert.equal(await state.save(), true);
    assert.equal(editor.dataFrame.rowCount, 1, 'the row left the frame');
    assert.equal(state.indexOf('i2'), 0, 'the survivor is row 0 now');
    assert.equal(state.indexOf('i1'), -1, 'the deleted key is gone');
    state.dispose();
  });

// ───────────────────── 3. the grid hook cannot express row-level permissions ────────────────────

scoped('finding 3: editability answers per row, not per table', {}, async () => {
  const {DomainFrameEditor} = await import('../../../js-api/src/ui/domains/domains-editor.js');
  const get = Object.getOwnPropertyDescriptor(DomainFrameEditor.prototype, 'writableColumns').get;
  // a row-mode table: the table-level right is false, the rows carry their own `~can_edit`
  const rowMode = {securityMode: 'row', can: {view: true, insert: true, edit: false, delete: false},
    fields: {id: 'readonly', title: 'editable', number: 'editable'}};
  assert.deepEqual(get.call({access: rowMode}), ['title', 'number'],
    'a row-mode table with editable rows is not a read-only grid');
  const editor = Object.create(DomainFrameEditor.prototype);
  editor.access = rowMode;
  editor._df = {columns: {byName: (n) => n === '~can_edit' ? {get: (row) => row === 0} : null}};
  editor.stateOf = (row) => row === 2 ? 'new' : '';
  assert.equal(editor.canEdit(0, 'title'), true, 'the row says it may be edited');
  assert.equal(editor.canEdit(1, 'title'), false, '`~can_edit` false locks that row only');
  assert.equal(editor.canEdit(0, 'id'), false, 'a readonly column stays readonly');
  assert.equal(editor.canEdit(2, 'title'), true, 'a draft answers to `insert`');
});

scoped('finding 3: a table-mode frame carries `~can_edit` too, and a row added to it is editable',
  {}, async () => {
    const {DomainFrameEditor} = await import('../../../js-api/src/ui/domains/domains-editor.js');
    // the server computes `~can_edit` in BOTH security modes (`repository.dart` _accessColumns);
    // a row appended locally has no answer in it — a bool column has no null slot, so it reads
    // false — and the right that governs it is `insert`
    const tableMode = {securityMode: 'table', can: {view: true, insert: true, edit: true, delete: false},
      fields: {id: 'readonly', label: 'readonly', title: 'editable'}};
    const editor = Object.create(DomainFrameEditor.prototype);
    editor.access = tableMode;
    editor._df = {columns: {byName: (n) => n === '~can_edit' ? {get: (row) => row === 0} : null}};
    editor.stateOf = (row) => row === 2 ? 'new' : '';
    editor._info = {singularName: 'Container'};
    editor._propByName = new Map([['label', {friendlyName: 'Label'}]]);
    editor._displayOf = () => null;
    assert.equal(editor.canEdit(0, 'title'), true, 'the row the server cleared is editable');
    assert.equal(editor.canEdit(1, 'title'), false, '`~can_edit` false locks that row in table mode too');
    assert.equal(editor.canEdit(2, 'title'), true,
      'a locally added row must not read its false-defaulted `~can_edit` cell');
    assert.equal(editor.refusalOf(2, 'title'), null, 'no refusal for a cell that may be edited');
    assert.equal(editor.refusalOf(2, 'label'), 'Label is read-only', 'the column is named, not the row');
    assert.equal(editor.refusalOf(1, 'title'), 'This container is read-only for you',
      'a locked row is refused in the table\'s words');
  });

// ───────────────────── 4. saving invalidates the offset the source keeps paging from ────────────

scoped('finding 4: a saved deletion re-bases the paging offset', {}, async () => {
    backends.domain = backend();
    const src = await issues({pageSize: 2});
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Aspirin', 'Ibuprofen']);
    src.edit.value.markDeleted('i1');
    assert.equal(await src.save(), true);
    assert.equal(src.total.value, 2, 'the count followed the save');
    await src.loadMore();
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Ibuprofen', 'Naproxen'],
      'the third row is not skipped by a stale offset');
    src.dispose();
  });

scoped('finding 4: a saved insert re-bases it too — no row skipped, none twice', {}, async () => {
    backends.domain = backend();
    const src = await issues({pageSize: 2});
    src.newRow({project_id: 'p1', title: 'Fourth'});
    assert.equal(await src.save(), true);
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen'],
      'the window is what the server answers at those offsets, not the frame the batch left');
    await src.loadMore();
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen', 'Fourth'],
      'the saved row comes back with the next page, once');
    assert.equal(src.total.value, 4);
    src.dispose();
  });

scoped('finding 4: an update that takes the row out of the query re-bases it as well', {}, async () => {
    backends.domain = backend();
    const src = await issues({pageSize: 1, query: 'done = false'});
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Ibuprofen']);
    src.rows.byKey('i2').done = true;
    assert.equal(await src.save(), true);
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Naproxen'],
      'the row that left the query is gone and the one behind it is not skipped');
    await src.loadMore();
    assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Naproxen'], 'and nothing is loaded twice');
    src.dispose();
  });

// ───────────────────────── 5. every standalone save leaks a session ─────────────────────────────

/** The participant surface `DomainSession` subscribes to and drives; an empty batch is enough —
 * the constructor has already subscribed by the time `save` returns. */
async function sessionFake() {
  const rxjs = await import('../../../js-api/node_modules/rxjs/index.js');
  return {
    table: 'grit.issue', quiet: true, isDirty: false, changeCount: 0, isSaving: false,
    client: {schema: 'grit', table: 'issue'},
    onChanged: new rxjs.Subject(), onDirtyChanged: new rxjs.Subject(),
    onSavingChanged: new rxjs.Subject(), onSaved: new rxjs.Subject(),
    prepareSave: () => [],
    applyResults: () => ({inserted: 0, updated: 0, deleted: 0, assigned: {}}),
    writeBack: async () => {},
    setSaving() {}, discard() {},
  };
}

scoped('finding 5: a standalone save leaves no subscription on the editor', {}, async () => {
  const {DomainFrameEditor} = await import('../../../js-api/src/ui/domains/domains-editor.js');
  const {DomainSession} = await import('../../../js-api/src/ui/domains/domains-session.js');
  const editor = await sessionFake();
  // the standalone path: `DomainFrameEditor.save` builds a session of one and disposes it
  for (let i = 0; i < 5; i++)
    assert.equal(await DomainFrameEditor.prototype.save.call(editor), true);
  assert.equal(editor.onChanged.observers.length, 0, 'onChanged listeners released');
  assert.equal(editor.onSavingChanged.observers.length, 0, 'onSavingChanged listeners released');
  const live = new DomainSession([editor], {quiet: true});
  assert.equal(editor.onChanged.observers.length, 1, 'a session that is kept stays subscribed');
  live.dispose();
  assert.equal(editor.onChanged.observers.length, 0);
});

// ─────────────────────────── 6. the delete-ordering rule ────────────────────────────────────────

const SELF = {name: 'stock', tables: {location: {friendlyName: 'Locations', columns: {
  name: {type: 'string', isName: true},
  parent_id: {type: 'ref', ref: 'location'},
}}}};

const rooms = () => new MemoryDomainBackend(SELF,
  {rows: {location: [{id: 'l1', name: 'Room'}, {id: 'l2', name: 'Shelf', parent_id: 'l1'}]}});

scoped('finding 6: a self-referencing table orders its own parent/child deletes', {}, async () => {
  const table = await rooms().table('stock.location');
  const results = await table.transaction([
    {op: 'delete', table: 'stock.location', id: 'l1'},
    {op: 'delete', table: 'stock.location', id: 'l2'}]);
  assert.deepEqual(results.map((r) => r.id), ['l1', 'l2'], 'the child delete was run first, results stay at request index');
});

scoped('finding 6 (recorded): child-first request order lands too; an unbatched child still vetoes', {}, async () => {
  const ok = await (await rooms().table('stock.location')).transaction([
    {op: 'delete', table: 'stock.location', id: 'l2'},
    {op: 'delete', table: 'stock.location', id: 'l1'}]);
  assert.deepEqual(ok.map((r) => r.id), ['l2', 'l1']);
  const table = await rooms().table('stock.location');
  await assert.rejects(() => table.transaction([{op: 'delete', table: 'stock.location', id: 'l1'}]),
    /Operation 0: row "l1" is referenced by location.parent_id/,
    'the FK veto stands for a child the batch does not delete');
});

scoped('finding 6 (recorded): a three-table delete cycle is refused, not silently mis-ordered', {}, async () => {
  const cycle = new MemoryDomainBackend({name: 'cyc', tables: {
    a: {columns: {name: {type: 'string', isName: true}, b_id: {type: 'ref', ref: 'b'}}},
    b: {columns: {name: {type: 'string', isName: true}, c_id: {type: 'ref', ref: 'c'}}},
    c: {columns: {name: {type: 'string', isName: true}, a_id: {type: 'ref', ref: 'a'}}},
  }}, {rows: {a: [{id: 'a1', name: 'A', b_id: 'b1'}], b: [{id: 'b1', name: 'B', c_id: 'c1'}],
    c: [{id: 'c1', name: 'C', a_id: 'a1'}]}});
  const table = await cycle.table('cyc.a');
  await assert.rejects(() => table.transaction([
    {op: 'delete', table: 'cyc.a', id: 'a1'}, {op: 'delete', table: 'cyc.b', id: 'b1'},
    {op: 'delete', table: 'cyc.c', id: 'c1'}]),
  /circular reference among refs/,
  'three non-reciprocal edges close a cycle: the sort fails — with a message that names no ref');
});

// ───────────────────── 7. the same persisted row in two sources of one session ──────────────────

async function twoSources() {
  backends.domain = backend();
  const session = new SharedSession();
  const one = new DomainSource({table: 'grit.issue', session}, env);
  const two = new DomainSource({table: 'grit.issue', session}, env);
  one.start();
  two.start();
  await flush();
  return {session, one, two};
}

const serverRow = async (id) => (await (await backends.domain.table('grit.issue'))
  .query({filter: {property: 'id', operator: '=', value: id}}))[0];

scoped('finding 12: the same persisted row pending in two sources is refused before anything is sent', {},
  async () => {
    const {session, one, two} = await twoSources();
    let sent = 0;
    const real = backends.domain.saveAll.bind(backends.domain);
    backends.domain.saveAll = (edits) => {
      sent++;
      return real(edits);
    };
    one.rows.byKey('i1').weight = 42;
    two.rows.byKey('i1').title = 'From two';
    assert.equal(session.changeCount.value, 2, 'the session counts both — it has no row identity');
    assert.equal(await session.save(), false);
    assert.equal(sent, 0, 'refused before the transaction, not by it');
    const message = String(one.error.value);
    assert.match(message, /Cannot save: Issue "Aspirin" is edited in two places/, 'the refusal names the row');
    assert.equal(String(two.error.value), message, 'every source holding it takes the refusal');
    assert.equal(one.problemRow.value, 'i1');
    assert.equal(two.problemRow.value, 'i1');
    assert.equal((await serverRow('i1')).version, 1, 'nothing landed');
    assert.equal(session.isDirty.value, true, 'both edits stay pending');

    // the columns being disjoint changes nothing: rejection is the whole policy, there is no merge
    two.revert();
    two.rows.byKey('i1').weight = 7;
    assert.equal(await session.save(), false, 'the same column in both is refused just the same');
    assert.equal(sent, 0);
    one.revert();
    assert.equal(await session.save(), true, 'the row pending in one source only saves');
    assert.equal((await serverRow('i1')).weight, 7);
    one.dispose();
    two.dispose();
  });

scoped('finding 7b: a sibling source of the same table reloads after another one saves', {}, async () => {
  const {session, one, two} = await twoSources();
  one.rows.byKey('i1').title = 'Only one';
  assert.equal(await session.save(), true);
  assert.equal((await serverRow('i1')).version, 2);
  await flush();
  assert.equal(two.rows.byKey('i1').title, 'Only one', 'the other source reloaded the row it holds');
  assert.equal(two.rows.byKey('i1').version, 2, 'and it is at the landed version');
  two.rows.byKey('i1').title = 'Now two';
  assert.equal(await session.save(), true, 'its next save lands on the version it now holds');
  await flush();
  assert.equal((await serverRow('i1')).title, 'Now two');
  one.dispose();
  two.dispose();
});
