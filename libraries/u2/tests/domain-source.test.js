/* DomainSource (WO-7) over the memory backend's frame: load → rows/total/access, a query change
   resets, loadMore appends into the same frame, drafts are rows (newRow → dirty → save → id
   assigned, clean), discard, the `currentRow.<col>` binding writes through the edit state, the
   spec round trip, the session and the guards a form pairs through, the frame's current row
   mirrored both ways, H8 disposal, and the containment of a missing table or backend. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal, Signal} from '../src/core/signals.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {backends} from '../src/sources/backends.js';
import {DomainSource} from '../src/sources/domain-source.js';
import {SingleSession} from '../src/sources/session.js';
import {MemoryFrame} from '../src/sources/memory-frame.js';
import {Rows} from '../src/sources/rows-like.js';
import {Access} from '../src/core/access.js';
import {notify} from '../src/components/display/notify.js';
import {SCHEMA, ROWS, backend} from './domain-fixtures.mjs';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';

/** Every test runs over its own backend and must leave the live-scope count where it was. */
function source(name, body) {
  test(name, async () => {
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

function env(options = {}) {
  return {designTime: false, subBinds: options.subBinds ?? {}, resolve: options.resolve ?? (() => null)};
}

async function issues(options = {}) {
  const src = new DomainSource({table: 'grit.issue', ...options}, env());
  src.start();
  await flush();
  return src;
}

const titles = (src) => src.rows.items.value.map((r) => r.title);
const fn = (src, name) => src.getFunctions().find((f) => f.name === name);
const UUID = /^[0-9a-f-]{36}$/;

/** The memory backend behind a table whose frames count their disposals (H8). */
function counting(memory, frames) {
  return {table: async (address) => {
    const t = await memory.table(address);
    return {
      ...t, access: () => t.access(), count: (f) => t.count(f), query: (s) => t.query(s),
      transaction: (ops) => t.transaction(ops),
      frame: async (spec) => {
        const frame = await t.frame(spec);
        const wrapped = {...frame, disposed: 0, dispose() {
          this.disposed++;
          frame.dispose();
        }};
        frames.push(wrapped);
        return wrapped;
      },
    };
  }};
}

source('load: rows by id, total, access, schema and the summary', async () => {
  backends.domain = backend();
  const src = new DomainSource({table: 'grit.issue', pageSize: 2}, env());
  assert.equal(src.state.value, 'idle');
  assert.equal(src.access.value, Access.readOnly, 'nothing is editable until the server says so');
  src.start();
  assert.equal(src.state.value, 'loading');
  assert.equal(src.summary.value, 'Loading…');
  await flush();
  assert.equal(src.state.value, 'ready');
  assert.deepEqual(titles(src), ['Aspirin', 'Ibuprofen']);
  assert.equal(src.total.value, 3);
  assert.equal(src.summary.value, '2 of 3');
  assert.equal(src.rows.byKey('i2').title, 'Ibuprofen');
  assert.equal(src.rows.keyOf(src.rows.items.value[0]), 'i1');
  assert.equal(src.access.value.field('title'), 'editable');
  assert.equal(src.access.value.field('id'), 'readonly');
  assert.equal(src.access.value.can('delete'), true);
  assert.equal(src.rows.items.value[0]['~can_edit'], true, 'withAccess: the rows carry their own');
  assert.equal('~can_edit' in src.rows.items.value[0], false, 'read by name, never enumerated (H7)');
  assert.equal(src.edit.value.isChanged('i1', 'title'), false, 'the writer is reachable for the controls');
  assert.equal(src.schema.info.pluralName, 'Issues');
  assert.equal(src.schema.properties.find((p) => p.name === 'project_id').ref, 'grit.project');
  assert.ok(src.df.value instanceof MemoryFrame, 'the memory backend is a frame host too');
  assert.ok(src.session instanceof SingleSession, 'phase 1: a source is a session of one');
  assert.equal(src.isDirty.value, false);
  assert.equal(src.isDraft, false);
  src.dispose();
});

source('query: a change reloads from the first page; loadMore appends until the last short page', async () => {
  backends.domain = backend();
  const src = await issues({pageSize: 2});
  const df = src.df.value;
  await src.loadMore();
  assert.deepEqual(titles(src), ['Aspirin', 'Ibuprofen', 'Naproxen']);
  assert.equal(src.df.value, df, 'appended into the same frame');
  assert.equal(src.summary.value, '3 issues');
  await src.loadMore();
  assert.equal(src.rows.items.value.length, 3, 'done: nothing more to load');

  src.query.value = 'done = false';
  await flush();
  assert.deepEqual(titles(src), ['Ibuprofen', 'Naproxen']);
  assert.notEqual(src.df.value, df, 'a new query is a new frame');
  assert.equal(src.total.value, 2);
  assert.equal(src.summary.value, '2 issues');
  assert.equal(src.bindStep('total').value, 2);
  src.dispose();
});

source('newRow → dirty → save: the draft is a row keyed by its index, gets its id, and the source is clean', async () => {
  backends.domain = backend();
  const src = await issues({defaults: {project_id: 'p1'}});
  const draft = src.newRow({title: ''});
  assert.equal(draft.id, Rows.draftKey(3), 'keyed by its index until saved');
  assert.equal(Rows.isDraft(draft), true);
  assert.equal(draft.project_id, 'p1', 'defaults apply');
  assert.equal(src.currentRow.value, draft, 'the draft is current: a form binds to it');
  assert.equal(src.rows.items.value.at(-1), draft, 'drafts are rows');
  assert.equal(src.df.value.currentRowIdx, 3, 'and the frame\'s current row');
  assert.equal(src.isDirty.value, true);
  assert.equal(src.changeCount.value, 1);
  assert.equal(src.summary.value, '1 unsaved change');
  assert.equal(src.validity.value, 'Value can\'t be empty', 'a required column');
  assert.equal(await src.save(), false, 'refused while invalid');
  assert.equal(src.error.value.message, 'Cannot save: Value can\'t be empty', 'and says why');
  assert.equal(src.isDirty.value, true);

  draft.title = 'New issue';
  assert.equal(src.validity.value, null);
  assert.equal(await src.save(), true);
  const saved = src.currentRow.value;
  assert.match(saved.id, UUID, 'the same row, now under the server\'s id, is current');
  assert.equal(saved.title, 'New issue');
  assert.equal(saved.version, 1);
  assert.equal(src.rows.byKey(draft.id), undefined, 'the draft key is gone');
  assert.equal(src.isDirty.value, false);
  assert.equal(src.error.value, undefined);
  assert.equal(src.total.value, 4);
  assert.equal(src.summary.value, '4 issues');
  const store = backends.domain.tableSync('grit.issue').rows;
  assert.equal(store.at(-1).title, 'New issue');
  assert.equal(store.at(-1).author_id, 'me', 'stamped on insert');
  src.dispose();
});

source('edit → save updates with the expected version; a conflict is an error, the change stays pending', async () => {
  backends.domain = backend();
  const src = await issues();
  const row = src.rows.byKey('i2');
  row.title = 'Ibuprofen 400';
  assert.equal(src.isDirty.value, true);
  assert.equal(src.changeCount.value, 1);
  assert.equal(row['~state'], 'modified');
  row.title = 'Ibuprofen';
  assert.equal(src.isDirty.value, false, 'back to the original is no change');
  row.title = 'Ibuprofen 400';
  assert.equal(await src.save(), true);
  assert.equal(row.version, 2);
  assert.equal(backends.domain.tableSync('grit.issue').rows[1].title, 'Ibuprofen 400');

  backends.domain.tableSync('grit.issue').rows[1].version = 9;
  row.title = 'Ibuprofen 600';
  assert.equal(await src.save(), false);
  assert.equal(src.error.value.code, 'version-conflict', 'the error object, code and all');
  assert.match(src.bindStep('error').value, /expected 2/);
  assert.match(src.summary.value, /expected 2/, 'the status bar says why');
  assert.equal(src.isDirty.value, true, 'nothing was lost');
  assert.equal(row.title, 'Ibuprofen 600');

  backends.domain.tableSync('grit.issue').rows[1].version = 2;
  assert.equal(await src.save(), true);
  assert.equal(src.error.value, undefined, 'a later success clears it');
  assert.equal(src.summary.value, '3 issues');
  src.dispose();
});

source('discard: modified values come back, drafts go, deletes are undone', async () => {
  backends.domain = backend();
  const src = await issues();
  const row = src.rows.byKey('i1');
  row.title = 'x';
  const draft = src.newRow({title: 'draft'});
  src.edit.value.markDeleted('i2');
  assert.equal(src.changeCount.value, 3);
  assert.deepEqual(titles(src), ['x', 'Ibuprofen', 'Naproxen', 'draft'], 'a deleted row stays an item until the save');
  assert.equal(src.rows.byKey('i2')['~state'], 'deleted', 'saying so');
  src.edit.value.unmarkDeleted('i2');
  assert.equal(src.rows.byKey('i2')['~state'], '');
  assert.equal(src.changeCount.value, 2);
  src.edit.value.markDeleted('i2');
  src.rows.byKey('i2').title = 'written while deleted';
  assert.equal(src.rows.byKey('i2')['~state'], 'deleted', 'a write does not undelete');
  assert.equal(src.edit.value.isChanged('i2', 'title'), true, 'but is tracked');
  src.edit.value.unmarkDeleted('i2');
  assert.equal(src.rows.byKey('i2')['~state'], 'modified', 'restored as edited');
  src.edit.value.markDeleted('i2');
  src.discard();
  assert.equal(row.title, 'Aspirin');
  assert.equal(src.rows.byKey('i2').title, 'Ibuprofen', 'the write on the deleted row is undone too');
  assert.deepEqual(titles(src), ['Aspirin', 'Ibuprofen', 'Naproxen'], 'discard restores deleted rows');
  assert.equal(src.isDirty.value, false);
  assert.equal(src.rows.byKey(draft.id), undefined, 'the draft is gone');
  assert.equal(src.currentRow.value.id, 'i3', 'the row now at its place is current');

  src.currentRow.value = src.rows.byKey('i2');
  src.session.discard();
  assert.equal(src.currentRow.value.id, 'i2', 'a discard with nothing pending changes nothing');
  src.dispose();
});

source('currentRow.<col>: a two-way step over whatever row is current, writes through the edit state', async () => {
  backends.domain = backend();
  const src = await issues();
  const current = src.bindStep('currentRow');
  const title = current.bindStep('title');
  assert.equal(title.value, undefined, 'no current row yet');
  assert.deepEqual(current.bindProps().map((p) => p.name).slice(5, 7), ['project_id', 'number']);
  assert.equal(current.bindProps().find((p) => p.name === 'id').writable, false);
  assert.equal(current.bindProps().find((p) => p.name === 'title').writable, true);

  src.currentRow.value = src.rows.byKey('i3');
  assert.equal(title.value, 'Naproxen');
  title.value = 'Naproxen 250';
  assert.equal(src.rows.byKey('i3').title, 'Naproxen 250');
  assert.equal(src.df.value.get('title', 2), 'Naproxen 250', 'the writer wrote the frame');
  assert.equal(src.edit.value.isChanged('i3', 'title'), true, 'and tracked it');
  assert.equal(src.isDirty.value, true);
  assert.equal(src.bindStep('isDirty').value, true);

  src.currentRow.value = src.rows.byKey('i1');
  assert.equal(title.value, 'Aspirin', 'follows the current row');
  src.rows.byKey('i1').title = 'Aspirin 100';
  assert.equal(title.value, 'Aspirin 100', 'and its edits');
  assert.equal(src.bindStep('currentRow').bindStep(''), null);
  assert.equal(src.bindStep('nope'), null);
  assert.equal(src.bindStep('access').value, src.access.value);
  assert.equal(src.bindStep('state').value, 'ready');
  src.dispose();
});

source('the frame\'s current row and currentRow are one thing, written in either direction', async () => {
  backends.domain = backend();
  const src = await issues();
  const df = src.df.value;
  assert.equal(df.currentRowIdx, -1);
  assert.equal(src.currentRow.value, null);
  df.currentRowIdx = 1;
  assert.equal(src.currentRow.value.id, 'i2');
  src.currentRow.value = src.rows.byKey('i1');
  assert.equal(df.currentRowIdx, 0, 'and back');
  src.currentRow.value = null;
  assert.equal(df.currentRowIdx, 0, 'no row is not a frame index');
  src.dispose();
});

source('access: the backend\'s answer gates the fields; withAccess: false takes everything for granted', async () => {
  backends.domain = backend({access: {can: {view: true, insert: true, edit: true, delete: false},
    fields: {project_id: 'editable', title: 'editable', number: 'readonly'}}});
  const src = await issues();
  const access = src.access.value;
  assert.equal(access.can('delete'), false);
  assert.equal(access.field('title'), 'editable');
  assert.equal(access.field('number'), 'readonly');
  assert.equal(access.field('weight'), 'hidden', 'a restricted column is absent from fields');
  assert.equal(src.rows.items.value[0]['~can_delete'], false);
  const pristine = src.newRow({project_id: 'p1', title: 'seed'}, {pristine: true});
  assert.equal(src.isDirty.value, false, 'a pristine draft never arms the gate');
  assert.equal(src.changeCount.value, 0);
  pristine.title = 'touched';
  assert.equal(src.changeCount.value, 1, 'until its first write');
  src.discard();
  const draft = src.newRow({project_id: 'p1', title: 'x', weight: 5});
  assert.equal(await src.save(), true);
  const stored = backends.domain.tableSync('grit.issue').rows.at(-1);
  assert.equal(stored.title, 'x');
  assert.equal(stored.weight, undefined, 'only the editable columns are sent');
  assert.equal(src.currentRow.value.weight, 5, 'the row keeps what it was given');
  assert.equal(draft.id, Rows.draftKey(3), 'a held draft proxy is not re-keyed; the current row is');

  const trusting = await issues({withAccess: false});
  assert.equal(trusting.access.value.can('delete'), false, 'the table-level access is always fetched');
  assert.equal(trusting.access.value.field('weight'), 'hidden');
  assert.equal(trusting.rows.items.value[0]['~can_edit'], undefined, 'only the row columns are off');
  src.dispose();
  trusting.dispose();
});

source('~can_share: null off row mode is not carried; a frame without the column neither; a row-mode row is', async () => {
  const rows = ROWS.issue.map((r) => ({...r}));
  rows[2]['~can_share'] = false;
  backends.domain = backend({rows: {issue: rows, project: ROWS.project.map((r) => ({...r}))},
    access: {can: {view: true, insert: true, edit: true, delete: true, share: true},
      fields: {title: 'editable', project_id: 'editable'}}});
  const src = await issues();
  const [first, , third] = src.rows.items.value;
  assert.equal(first['~can_share'], null, 'the JSON shape off row mode');
  assert.equal(first['~can_edit'], true);
  assert.equal(src.access.value.row(first).can('share'), true, 'not carried: the table\'s answer');
  assert.equal(third['~can_share'], false, 'a row-mode table carries booleans');
  assert.equal(src.access.value.row(third).can('share'), false, 'and the row is the truth');
  src.dispose();

  const bare = await issues({withAccess: false});
  const row = bare.rows.items.value[0];
  assert.equal(bare.df.value.columns.byName('~can_share'), null, 'as a d42 frame off row mode: no column');
  assert.equal(row['~can_share'], undefined);
  assert.equal(bare.access.value.row(row).can('share'), true);
  bare.dispose();
});

source('spec: the tray tag round-trips through dump, a bound query narrows, a bound input edits the current row',
  async () => {
    backends.domain = backend();
    const reg = new Registry();
    registerAll(reg);
    const meta = reg.get('u2-domain-source');
    assert.equal(meta.visual, false);
    assert.equal(meta.category, 'Data');
    assert.equal(typeof meta.usage, 'string');
    assert.deepEqual(meta.props.map((p) => p.name), ['table', 'query', 'pageSize', 'withAccess', 'defaults', 'draft']);
    assert.equal(meta.props.find((p) => p.name === 'query').bindable, true);

    const spec = {
      $schema: 'dg-ui/1',
      components: [{tag: 'u2-domain-source', name: 'issues', props: {table: 'grit.issue', pageSize: 10},
        bind: {query: '$.search'}}],
      root: {tag: 'u2-form', name: 'form', children: [
        {tag: 'u2-text-input', name: 'title', props: {label: 'Title'}, bind: {value: '$.issues.currentRow.title'}},
      ]},
    };
    const ctx = new SpecContext({data: {search: 'done = false'}});
    const instance = renderSpec(spec, ctx, reg);
    document.body.append(instance.root);
    await flush();
    assert.deepEqual(instance.dump(), spec);
    const src = instance.node('issues');
    assert.ok(src instanceof DomainSource);
    assert.equal(src.bindStep('query'), src.query, 'the declared prop stays a step');
    assert.deepEqual(titles(src), ['Ibuprofen', 'Naproxen']);
    assert.deepEqual(src.getFunctions().map((f) => f.name), ['refresh', 'loadMore', 'save', 'discard', 'newRow']);

    const input = instance.root.querySelector('input');
    assert.equal(input.value, '', 'no current row: the field is empty, not "undefined"');
    await fn(src, 'newRow').apply({values: {project_id: 'p1', title: 'From the form'}});
    assert.equal(input.value, 'From the form');
    assert.equal(src.isDirty.value, false, 'a spec-driven New is a pristine draft');
    input.value = 'Edited in the form';
    input.dispatchEvent(new Event('input', {bubbles: true}));
    assert.equal(src.currentRow.value.title, 'Edited in the form');
    assert.equal(src.isDirty.value, true);
    assert.equal(await fn(src, 'save').apply(), true, 'cmd:issues.save goes through the session');
    assert.equal(backends.domain.tableSync('grit.issue').rows.at(-1).title, 'Edited in the form');

    assert.equal(titles(src).length, 3, 'the saved draft stays a row');
    src.currentRow.value.title = 'pending';
    ctx.data.search.value = 'title = "Aspirin"';
    await flush();
    assert.equal(titles(src).length, 3, 'H6: a re-query never drops pending edits silently');
    assert.equal(src.isDirty.value, true);
    src.discard();
    ctx.data.search.value = 'title = "Aspirin" or done = true';
    await flush();
    assert.deepEqual(titles(src), ['Aspirin'], 'a new query is a new collection once nothing is pending');
    assert.equal(src.currentRow.value, null);
    instance.dispose();
  });

source('session: save and discard go through it, it announces and emits; a given session takes over', async () => {
  backends.domain = backend();
  const src = await issues();
  const heard = [];
  const subs = [src.session.onSaved.subscribe(() => heard.push('saved')),
    src.session.onDiscarded.subscribe(() => heard.push('discarded'))];
  assert.equal(src.session.isDirty, src.isDirty);
  assert.equal(src.session.isSaving.value, false);
  src.rows.byKey('i1').title = 'Aspirin 100';
  src.rows.byKey('i2').title = 'Ibuprofen 400';
  const saving = src.save();
  assert.equal(src.session.isSaving.value, true, 'while the transaction is in flight');
  assert.equal(await saving, true);
  assert.equal(src.session.isSaving.value, false);
  assert.deepEqual(heard, ['saved']);
  assert.match(document.body.querySelector('.u2-notify-info')?.textContent ?? '', /2 changes saved/, 'the one balloon');
  notify.closeAll();
  src.rows.byKey('i1').title = 'x';
  assert.equal(await src.save(), true);
  assert.match(document.body.querySelector('.u2-notify-info')?.textContent ?? '', /Issue saved/, 'one change: the row');
  src.rows.byKey('i1').title = 'x';
  src.discard();
  assert.deepEqual(heard, ['saved', 'saved', 'discarded']);
  assert.equal(src.isDirty.value, false);
  assert.equal(await src.save(), false, 'nothing pending: nothing landed');
  assert.deepEqual(heard, ['saved', 'saved', 'discarded'], 'and nothing announced');
  for (const s of subs)
    s.unsubscribe();
  src.dispose();

  const calls = [];
  const session = {isDirty: signal(true), changeCount: signal(1), validity: signal(null), isSaving: signal(false),
    onSaved: {subscribe: () => ({unsubscribe() {}})}, onDiscarded: {subscribe: () => ({unsubscribe() {}})},
    save: async () => calls.push('save') && true, discard: () => calls.push('discard')};
  const shared = await issues({session});
  assert.equal(shared.session, session);
  assert.equal(await shared.save(), true);
  await fn(shared, 'save').apply();
  shared.discard();
  assert.deepEqual(calls, ['save', 'save', 'discard'], 'every path goes through the session given');
  shared.dispose();
});

source('commit: a guard refuses first and names why; the writer\'s validity next; the error is the summary', async () => {
  backends.domain = backend();
  const src = await issues();
  src.rows.byKey('i1').title = 'Aspirin 100';
  let problem = 'Title is required';
  const unguard = src.guard(() => problem);
  assert.equal(await src.save(), false);
  assert.equal(src.error.value.code, 'validation');
  assert.equal(src.error.value.message, 'Cannot save: Title is required');
  assert.equal(src.summary.value, 'Cannot save: Title is required');
  assert.equal(src.isDirty.value, true, 'nothing was written');
  problem = null;
  assert.equal(await src.save(), true, 'the guard passes');
  assert.equal(src.error.value, undefined);
  unguard();
  src.newRow({project_id: 'p1'});
  assert.equal(await src.save(), false, 'the writer refuses a required column');
  assert.equal(src.error.value.message, 'Cannot save: Value can\'t be empty');
  src.dispose();
});

source('H8: a re-query replaces the frame and disposes the old one with its writer; so does disposal', async () => {
  const memory = backend();
  const frames = [];
  backends.domain = counting(memory, frames);
  const src = await issues({pageSize: 2});
  const first = frames.at(-1);
  assert.equal(src.df.value, first.df);
  assert.ok(first.df.subscriberCount > 0, 'the source listens to the frame');
  src.rows.byKey('i1').title = 'x';
  src.discard();
  src.query.value = 'done = false';
  await flush();
  assert.equal(first.disposed, 1, 'the replaced frame was disposed');
  assert.equal(first.df.subscriberCount, 0, 'and let go of');
  const second = frames.at(-1);
  assert.notEqual(second, first);
  assert.equal(src.df.value, second.df);
  assert.deepEqual(titles(src), ['Ibuprofen', 'Naproxen']);
  src.dispose();
  assert.equal(second.disposed, 1, 'disposal disposes the current frame too');
  assert.equal(second.df.subscriberCount, 0);
  assert.equal(src.edit.value, undefined);
});

source('stale loads: a frame a later refresh outran never becomes the collection', async () => {
  const memory = backend();
  const frames = [];
  const gate = [];
  const inner = counting(memory, frames);
  backends.domain = {table: async (address) => {
    const table = await inner.table(address);
    const frame = table.frame;
    table.frame = (spec) => new Promise((resolve) => gate.push(() => resolve(frame(spec))));
    return table;
  }};
  const src = new DomainSource({table: 'grit.issue', pageSize: 10}, env());
  src.start();
  await flush();
  assert.equal(gate.length, 1, 'the first page is waiting');
  src.query.value = 'done = true';
  await flush();
  assert.equal(gate.length, 2);
  gate[1]();
  await flush();
  assert.deepEqual(titles(src), ['Aspirin']);
  assert.equal(src.total.value, 1);
  const kept = src.df.value;
  const edit = src.edit.value;
  gate[0]();
  await flush();
  assert.equal(src.df.value, kept, 'the late frame never became the collection');
  assert.equal(src.edit.value, edit, 'the writer was not replaced');
  assert.equal(frames.at(-1).disposed, 1, 'and was disposed on arrival');
  assert.equal(frames[0].disposed, 0, 'the kept one lives');
  assert.deepEqual(titles(src), ['Aspirin']);
  assert.equal(src.total.value, 1);
  src.dispose();
});

source('containment: an unknown table is an error state; no backend throws at construction', async () => {
  backends.domain = backend();
  const src = await issues({table: 'grit.nope'});
  assert.equal(src.state.value, 'error');
  assert.equal(src.error.value.code, 'not-found');
  assert.equal(src.bindStep('error').value, 'Unknown table "grit.nope"');
  assert.equal(src.summary.value, 'Unknown table "grit.nope"');
  assert.deepEqual(src.rows.items.value, []);
  assert.throws(() => src.newRow(), /not loaded/);
  assert.equal(await src.save(), false, 'nothing to save through');
  src.dispose();

  let attempts = 0;
  const flaky = backend();
  backends.domain = {table: (address) => ++attempts === 1 ? Promise.reject(new Error('offline')) : flaky.table(address)};
  const retried = await issues();
  assert.equal(retried.state.value, 'error');
  await retried.refresh();
  assert.equal(retried.state.value, 'ready', 'a failed acquisition is not cached');
  assert.equal(attempts, 2);
  retried.dispose();

  delete backends.domain;
  assert.throws(() => new DomainSource({table: 'grit.issue'}, env()), /no platform backend/);
});

source('schema round trip: the registry example builds, and the sub-bind and query prop agree', async () => {
  backends.domain = backend();
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'div'}, components: [reg.get('u2-domain-source').example]},
    new SpecContext(), reg);
  await flush();
  const src = instance.node('issues');
  assert.equal(src.table, 'grit.issue');
  assert.equal(src.pageSize, 50);
  assert.equal(src.state.value, 'ready');
  assert.equal(SCHEMA.tables.issue.columns.title.isName, true);
  instance.dispose();
});

source('the `source` step hands the source itself to a bound control', async () => {
  backends.domain = backend();
  const src = await issues();
  const step = src.bindStep('source');
  assert.ok(step instanceof Signal, 'a signal, so the resolver stops there');
  assert.equal(step.peek(), src);
  assert.equal(src.bindProps().find((p) => p.name === 'source').type, 'object');
  const reg = new Registry();
  registerAll(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', components: [reg.get('u2-domain-source').example],
    root: {tag: 'div'}}, new SpecContext(), reg);
  await flush();
  assert.equal(instance.resolveBinding('$.issues.source').signal.peek(), instance.node('issues'));
  instance.dispose();
  src.dispose();
});

source('a draft is written under insert, an existing row under edit as its own columns narrow it', async () => {
  const rows = Object.fromEntries(Object.entries(ROWS).map(([t, list]) => [t, list.map((r) => ({...r}))]));
  rows.issue[0]['~can_edit'] = false;
  const fields = {id: 'readonly', title: 'editable', project_id: 'editable', number: 'editable'};
  backends.domain = new MemoryDomainBackend(SCHEMA, {rows,
    access: {can: {view: true, insert: true, edit: true, delete: false, share: false}, fields}});
  const src = await issues();
  const columns = (s) => s.bindStep('currentRow').bindProps().filter((p) => p.writable).map((p) => p.name);
  src.currentRow.value = src.rows.byKey('i1');
  assert.deepEqual(columns(src), [], 'the row says the caller may not edit it');
  src.currentRow.value = src.rows.byKey('i2');
  assert.deepEqual(columns(src), ['project_id', 'number', 'title']);
  const draft = src.newRow({project_id: 'p1', title: 'Draft'});
  assert.deepEqual(columns(src), ['project_id', 'number', 'title'], 'a draft: insert is what counts');
  draft.number = 9;
  src.rows.byKey('i1').title = 'Aspirin 100';
  src.rows.byKey('i2').title = 'Ibu';
  assert.equal(await src.save(), true);
  const store = backends.domain.tableSync('grit.issue').rows;
  assert.equal(store.length, 4);
  assert.deepEqual([store[3].title, store[3].number, store[3].project_id], ['Draft', 9, 'p1'], 'inserted whole');
  assert.equal(store[0].title, 'Aspirin', 'the caller may not edit i1: the change was dropped');
  assert.equal(store[1].title, 'Ibu');
  src.dispose();

  // an insert-only caller: nothing to write on any existing row, everything on a draft
  backends.domain = new MemoryDomainBackend(SCHEMA, {rows: {issue: ROWS.issue.map((r) => ({...r}))},
    access: {can: {view: true, insert: true, edit: false, delete: false, share: false}, fields}});
  const inserter = await issues();
  inserter.currentRow.value = inserter.rows.byKey('i2');
  assert.deepEqual(columns(inserter), []);
  const fresh = inserter.newRow({project_id: 'p1', title: 'Fresh'});
  assert.deepEqual(columns(inserter), ['project_id', 'number', 'title']);
  fresh.number = 1;
  assert.equal(await inserter.save(), true);
  assert.equal(backends.domain.tableSync('grit.issue').rows.at(-1).title, 'Fresh');
  inserter.dispose();
});

source('a saved draft stays current under the id the server gave it; a none cell reads null', async () => {
  backends.domain = backend();
  const src = await issues({pageSize: 5});
  const draft = src.newRow({title: 'Draft'}, {pristine: true});
  assert.equal(draft.id, Rows.draftKey(3));
  assert.equal(src.currentRow.value, draft);
  assert.equal(src.summary.value, '3 issues', 'a draft is not a row of the table yet');
  draft.project_id = 'p1';
  assert.equal(await src.save(), true);
  assert.match(src.currentRow.value?.id, UUID, 'the same row, re-keyed by its id, is current');
  assert.equal(src.currentRow.value?.title, 'Draft');
  assert.equal(src.summary.value, '4 issues');

  const df = src.df.value;
  df.rows[0].number = null;
  assert.equal(src.rows.byKey('i1').number, null);
  assert.equal(src.rows.byKey('i1').priority, null, 'a column the record lacks reads null, not undefined');
  assert.equal(src.rows.byKey('i2').number, 2);
  src.dispose();
});

source('the summary: singular for one, "New <row>" for a draft source, deletions pending', async () => {
  backends.domain = backend();
  const one = await issues({query: 'title = "Aspirin"'});
  assert.equal(one.summary.value, '1 issue');
  one.edit.value.markDeleted('i1');
  assert.equal(one.summary.value, '1 deletion pending');
  one.edit.value.unmarkDeleted('i1');
  assert.equal(one.summary.value, '1 issue');
  one.dispose();
  const some = await issues();
  some.rows.byKey('i2').title = 'Ibu';
  some.edit.value.markDeleted('i1');
  assert.equal(some.summary.value, '2 unsaved changes', 'a deletion among other changes is one of them');
  some.dispose();

  const draft = new DomainSource({table: 'grit.issue', draft: true, defaults: {title: 'Draft'}}, env());
  draft.start();
  assert.equal(draft.summary.value, 'Loading…');
  await flush();
  assert.equal(draft.isDraft, true);
  assert.equal(draft.rows.items.value.length, 1, 'nothing loaded: the pristine draft alone');
  assert.equal(draft.total.value, 0);
  assert.equal(draft.summary.value, 'New issue', 'never "N of M" for a draft source');
  const row = draft.currentRow.value;
  assert.equal(row.title, 'Draft', 'over the defaults');
  assert.equal(draft.isDirty.value, false, 'a pristine draft is not a change');
  row.number = 5;
  assert.equal(draft.summary.value, '1 unsaved change');
  await draft.loadMore();
  assert.equal(draft.rows.items.value.length, 1, 'a draft source never pages');
  draft.dispose();
});

source('a saved delete leaves the rows and the count; the current row moves on', async () => {
  backends.domain = backend();
  const src = await issues({pageSize: 5});
  src.currentRow.value = src.rows.byKey('i2');
  src.edit.value.markDeleted('i2');
  assert.equal(src.rows.byKey('i2')['~state'], 'deleted', 'kept until the save');
  assert.equal(src.currentRow.value.id, 'i2');
  assert.equal(await src.save(), true);
  assert.deepEqual(titles(src), ['Aspirin', 'Naproxen'], 'the writer removed the row without a frame event');
  assert.equal(src.rows.byKey('i2'), undefined, 'no proxy survives its row');
  assert.equal(src.currentRow.value.id, 'i3', 'the row at its place is current');
  assert.equal(src.summary.value, '2 issues');
  assert.equal(backends.domain.tableSync('grit.issue').rows.length, 2);

  src.currentRow.value = src.rows.byKey('i3');
  src.edit.value.markDeleted('i3');
  assert.equal(await src.save(), true);
  assert.deepEqual(titles(src), ['Aspirin']);
  assert.equal(src.currentRow.value.id, 'i1', 'the last row: the one before it');
  assert.equal(src.summary.value, '1 issue');
  src.dispose();
});

source('a draft source: "New <row>" until the draft is saved, its row afterwards, "New" again on a fresh draft', async () => {
  backends.domain = backend();
  const draft = new DomainSource({table: 'grit.issue', draft: true, defaults: {project_id: 'p1'}}, env());
  draft.start();
  await flush();
  const row = draft.currentRow.value;
  assert.equal(draft.summary.value, 'New issue');
  row.title = 'Draft';
  row.number = 7;
  assert.equal(await draft.save(), true);
  assert.equal(draft.currentRow.value.title, 'Draft', 'the saved row stays current');
  assert.equal(Rows.isDraft(draft.currentRow.value), false);
  assert.equal(draft.summary.value, '1 issue');
  draft.newRow({title: 'Another'}, {pristine: true});
  assert.equal(draft.summary.value, 'New issue');
  draft.dispose();
});
