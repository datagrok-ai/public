import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Rows, notify} from '@datagrok-libraries/u2';
import {DomainSource} from '../src/core/domain-source.js';
import {backend} from './domain-fixtures.mjs';

/** Every test runs over its own backend and must leave the live-scope count where it was. */
function source(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const env = () => ({designTime: false, subBinds: {}, resolve: () => null});

async function rows(table, options = {}) {
  const src = testOwner().build(() => new DomainSource({table, ...options}, env()));
  src.start();
  await flush();
  return src;
}

const issues = (options) => rows('grit.issue', options);
const titles = (src) => src.rows.items.value.map((r) => r.title);

/** Marks the rows deleted through the writer and lands the batch. */
async function drop(src, ...ids) {
  for (const id of ids)
    src.edit.value.markDeleted(id);
  assert.equal(await src.save(), true, 'the delete landed');
  await flush();
}

source('a landed delete is soft: the trash shows the row, a staged restore brings it back on Save', async () => {
  emsBackends.domain = backend();
  const live = await issues();
  const trash = await issues({deleted: 'only'});
  assert.deepEqual(titles(trash), [], 'nothing deleted yet');
  await drop(live, 'i2');
  assert.deepEqual(titles(live), ['Aspirin', 'Naproxen']);

  await trash.refresh();
  assert.deepEqual(titles(trash), ['Ibuprofen']);
  assert.equal(Rows.isDeleted(trash.rows.byKey('i2')), true, 'the row carries ~is_deleted');
  assert.equal(trash.summary.value, '1 deleted issue');

  trash.stageRestore(['i2']);
  await flush();
  assert.equal(trash.isDirty.value, true, 'a staged restore is a pending change');
  assert.equal(trash.changeCount.value, 1);
  assert.equal(trash.summary.value, '1 restore pending');
  assert.deepEqual(titles(trash), ['Ibuprofen'], 'nothing has happened on the server yet');
  // pending() is what the session reads to word its balloon and what check() filters: a restored
  // row that is not in it makes both of them answer about an empty batch (B1)
  assert.deepEqual(trash.pending().map((r) => r.id), ['i2']);
  assert.equal(trash.check(), null, 'a trash source with nothing but restores staged may save');

  notify.closeAll();
  assert.equal(await trash.save(), true);
  assert.equal(document.body.querySelector('.u2-notify-info')?.textContent, 'Issue restored',
    'a landed restore is restored, not deleted');
  await flush();
  assert.deepEqual(titles(trash), [], 'the restored row left the trash');
  assert.equal(trash.summary.value, '0 deleted issues');

  await live.refresh();
  assert.deepEqual(titles(live), ['Aspirin', 'Ibuprofen', 'Naproxen']);
  assert.equal(live.rows.byKey('i2').version, 3, 'the delete and the restore are a version each');
  live.dispose();
  trash.dispose();
});

source('Discard drops a staged restore and the row is still deleted', async () => {
  emsBackends.domain = backend();
  const live = await issues();
  await drop(live, 'i2');
  live.dispose();
  const trash = await issues({deleted: 'only'});
  trash.stageRestore(['i2']);
  await flush();
  assert.equal(trash.changeCount.value, 1);
  trash.discard();
  await flush();
  assert.equal(trash.isDirty.value, false);
  assert.equal(trash.rows.byKey('i2')[Rows.STATE], '');
  assert.equal(Rows.isDeleted(trash.rows.byKey('i2')), true, 'still in the trash');
  trash.dispose();
});

source('two staged restores are ONE transaction; a staged delete beside them is still refused', async () => {
  const be = backend();
  emsBackends.domain = be;
  const live = await issues();
  await drop(live, 'i1', 'i2');
  live.dispose();
  const both = await issues({deleted: 'include'});
  both.stageRestore(['i1', 'i2']);
  both.edit.value.markDeleted('i3');
  await flush();
  assert.equal(both.changeCount.value, 3);
  // a source over deleted rows is read-only for everything but the restores it stages
  assert.equal(both.check(), 'deleted rows are read-only until they are restored');
  assert.equal(await both.save(), false);

  both.edit.value.unmarkDeleted('i3');
  await flush();
  assert.equal(await both.save(), true);
  await flush();
  const table = be.tableSync('grit.issue');
  const undeletes = table.history.filter((line) => line.op === 'undelete');
  assert.equal(undeletes.length, 2);
  assert.equal(undeletes[0].tx_id, undeletes[1].tx_id, 'one transaction');
  assert.deepEqual(both.rows.items.value.map((r) => Rows.isDeleted(r)), [false, false, false]);
  both.dispose();
});

source('a trash source with an EDIT pending is still refused with the unchanged sentence', async () => {
  emsBackends.domain = backend();
  const live = await issues();
  await drop(live, 'i1');
  live.dispose();
  const trash = await issues({deleted: 'only'});
  // the narrowed access keeps a control from doing this; the writer is the last line of defence
  trash.edit.value.setValue('i1', 'title', 'Edited');
  await flush();
  assert.equal(trash.check(), 'deleted rows are read-only until they are restored');
  assert.equal(await trash.save(), false);
  assert.deepEqual(trash.pending().map((r) => r[Rows.STATE]), ['modified'],
    'the edit is what the guard sees, and it is not a restore');
  trash.dispose();
});

source('stageRestore over a backend without restore refuses by name before anything is staged', async () => {
  const memory = backend();
  emsBackends.domain = memory;
  const live = await issues();
  await drop(live, 'i2');
  live.dispose();
  emsBackends.domain = {
    table: async (address) => {
      const t = await memory.table(address);
      return {address: t.address, properties: t.properties, info: t.info, access: () => t.access(),
        query: (spec) => t.query(spec), count: (scope) => t.count(scope),
        transaction: (ops) => t.transaction(ops), frame: (spec) => t.frame(spec)};
    },
    saveAll: (edits) => memory.saveAll(edits),
  };
  const src = await issues({deleted: 'include'});
  assert.throws(() => src.stageRestore(['i2']),
    /grit\.issue: the backend does not support restoring deleted rows/);
  assert.equal(src.isDirty.value, false, 'nothing was staged');
  src.dispose();
});

source('the deleted mode is a signal: one source flips to the trash, to both, and back', async () => {
  emsBackends.domain = backend();
  const src = await issues();
  await drop(src, 'i2');

  src.deleted.value = 'only';
  await flush();
  assert.deepEqual(titles(src), ['Ibuprofen']);
  assert.equal(src.readOnly.value, true);

  src.deleted.value = 'include';
  await flush();
  assert.deepEqual(titles(src), ['Aspirin', 'Ibuprofen', 'Naproxen']);
  assert.deepEqual(src.rows.items.value.map((r) => Rows.isDeleted(r)), [false, true, false]);

  src.deleted.value = 'exclude';
  await flush();
  assert.deepEqual(titles(src), ['Aspirin', 'Naproxen']);
  assert.equal(src.readOnly.value, false);
  assert.equal(src.rows.byKey('i1')[Rows.DELETED], undefined, 'the column is not projected at all');
  src.dispose();
});

source('a trash source is read-only: the narrowed access is an upper bound, and save says why', async () => {
  emsBackends.domain = backend();
  const live = await issues();
  await drop(live, 'i1');
  live.dispose();
  const trash = await issues({deleted: 'only'});
  assert.equal(trash.readOnly.value, true);
  const access = trash.access.value;
  assert.equal(access.can('edit'), false);
  assert.equal(access.can('insert'), false);
  assert.equal(access.can('delete'), true, 'Restore is the delete grant');
  assert.equal(access.field('title'), 'readonly', 'every field degrades to text');
  const row = trash.rows.byKey('i1');
  assert.equal(row['~can_edit'], true, 'the row still carries the table\'s own answer');
  assert.equal(access.row(row).can('edit'), false, 'which cannot lift the bound');
  assert.equal(access.row(row).can('delete'), true);
  assert.equal(trash.check(), null, 'a clean trash source is not in a batch and refuses nothing');
  assert.equal(await trash.save(), false, 'there is nothing to save');
  trash.dispose();
});

source('stageRestoreSelection: every selected row staged, then one Save', async () => {
  emsBackends.domain = backend();
  const live = await issues();
  await drop(live, 'i1', 'i3');
  live.dispose();
  const trash = await issues({deleted: 'only'});
  assert.deepEqual(titles(trash), ['Aspirin', 'Naproxen']);
  trash.df.value.selection.set(0, true);
  trash.df.value.selection.set(1, true);
  await flush();
  assert.equal(trash.selection.value.length, 2);
  trash.stageRestoreSelection();
  await flush();
  assert.equal(trash.changeCount.value, 2);
  assert.equal(trash.summary.value, '2 restores pending');
  assert.equal(await trash.save(), true);
  await flush();
  assert.deepEqual(titles(trash), []);
  trash.dispose();
});

source('a restore the backend refuses fails the batch and the rows stay in the trash', async () => {
  emsBackends.domain = backend();
  const live = await issues();
  await drop(live, 'i3');
  live.dispose();
  const projects = await rows('grit.project');
  await drop(projects, 'p2');
  projects.dispose();

  const trash = await issues({deleted: 'only'});
  trash.stageRestore(['i3']);
  await flush();
  assert.equal(await trash.save(), false);
  await flush();
  assert.match(String(trash.error.value.message), /Column "project_id" references a deleted row in "project"/);
  assert.deepEqual(titles(trash), ['Naproxen'], 'the row stays in the trash');
  trash.dispose();
});

source('a deleted source over a backend that cannot restore is refused by name', async () => {
  const memory = backend();
  // the seam without `restore` — a platform an older server answers, or a backend of one's own
  emsBackends.domain = {
    table: async (address) => {
      const t = await memory.table(address);
      return {address: t.address, properties: t.properties, info: t.info, access: () => t.access(),
        support: {...t.support, restore: false},
        query: (spec) => t.query(spec), count: (f, s, d) => t.count(f, s, d),
        transaction: (ops) => t.transaction(ops), frame: (spec) => t.frame(spec)};
    },
    saveAll: (edits) => memory.saveAll(edits),
  };
  const src = await issues({deleted: 'only'});
  assert.equal(src.state.value, 'error');
  assert.equal(src.error.value.code, 'unsupported');
  assert.match(String(src.error.value.message), /grit\.issue: the backend does not support restoring deleted rows/);
  src.dispose();
  assert.throws(() => DomainSource.requireRestore({address: 'grit.issue'}),
    /grit\.issue: the backend does not support restoring deleted rows/, 'and the guard refuses at construction');
});
