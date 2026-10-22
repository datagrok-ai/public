/* `DomainAppOptions.inline` over the memory backend: a new batch's create page carries a form for a new compound
   it may point at — the reference follows that form's first edit, both rows land in one save, and picking an
   existing compound drops the new one. Also the entity page's ⋯ menu and a ribbon without a search box. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');

function scoped(name, body) {
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
    assert.equal(DomainApp.live.size, 0, 'no app left behind');
  });
}

const LAB = {
  name: 'lab',
  tables: {
    compound: {
      singularName: 'Compound', pluralName: 'Compounds',
      columns: {name: {type: 'string', required: true, isName: true, searchable: true}},
    },
    batch: {
      singularName: 'Batch', pluralName: 'Batches',
      columns: {
        compound_id: {type: 'ref', ref: 'compound', required: true},
        amount: {type: 'float'},
      },
    },
  },
};

async function batches() {
  const backend = new MemoryDomainBackend(LAB, {rows: {compound: [{id: 'c1', name: 'Aspirin'}], batch: []}});
  emsBackends.domain = backend;
  const compounds = await domains.table('lab.compound');
  const table = await domains.table('lab.batch');
  const app = domains.app({table, base: '/apps/T/Batches', inline: {compound_id: {table: compounds}},
    immutable: ['compound_id']});
  document.body.append(app.root);
  await flush();
  return {backend, app};
}

const rowsOf = async (backend, table) => (await backend.table(`lab.${table}`)).query({});

scoped('a new batch and its new compound land in one save', async () => {
  const {backend, app} = await batches();
  await app.create();
  await flush();
  const compound = app.inline.value.get('compound_id');
  assert.ok(compound !== undefined, 'the create page has a new compound');
  const section = app.root.querySelector('.u2-domain-app-inline');
  assert.equal(section.hidden, false);
  assert.equal(app.session.isDirty.value, false, 'an untouched create page asks nothing on leaving');
  assert.equal(section.classList.contains('u2-section-invalid'), false, 'nor is it marked wrong before an edit');
  const draft = compound.currentRow.value;
  draft.name = 'Ibuprofen';
  await flush();
  const batch = app.entitySource.value.currentRow.value;
  assert.equal(batch.compound_id, draft.id, 'the batch points at the new compound');
  batch.amount = 2.5;
  assert.equal(await app.session.save(), true);
  await flush();
  const saved = (await rowsOf(backend, 'compound')).find((r) => r.name === 'Ibuprofen');
  assert.ok(saved !== undefined);
  const [row] = (await rowsOf(backend, 'batch'));
  assert.equal(row.compound_id, saved.id);
  assert.equal(row.amount, 2.5);
  assert.equal(section.hidden, true, 'a saved batch has no new compound to fill in');
  await flush();
  await flush();
  assert.equal(app.root.querySelector('[data-u2-name="compound_id"] [data-u2-part="readonly-value"]')?.textContent,
    'Ibuprofen', 'once saved, its compound reads as the row it named, not a picker');
  assert.equal(document.body.querySelector('.u2-notify-warning'), null, 'the save discards nothing while it runs');
  app.dispose();
});

scoped('picking an existing compound drops the new one', async () => {
  const {backend, app} = await batches();
  await app.create();
  await flush();
  const compound = app.inline.value.get('compound_id');
  compound.currentRow.value.name = 'Draft';
  await flush();
  const batch = app.entitySource.value.currentRow.value;
  batch.compound_id = 'c1';
  await flush();
  assert.equal(compound.isDirty.value, false);
  assert.equal(app.root.querySelector('.u2-domain-app-inline').hidden, true);
  assert.equal(await app.session.save(), true);
  assert.deepEqual((await rowsOf(backend, 'compound')).map((r) => r.name), ['Aspirin']);
  assert.equal((await rowsOf(backend, 'batch'))[0].compound_id, 'c1');
  app.dispose();
});

scoped('a reaction to the new compound\'s own write that picks an existing one drops it after the write', async () => {
  const {backend, app} = await batches();
  await app.create();
  await flush();
  const compound = app.inline.value.get('compound_id');
  const edit = compound.edit.peek();
  let writing = false;
  let discardedMidWrite = false;
  const setValue = edit.setValue.bind(edit);
  edit.setValue = (...args) => {
    writing = true;
    try {
      setValue(...args);
    } finally {
      writing = false;
    }
  };
  const discard = edit.discard.bind(edit);
  edit.discard = () => {
    discardedMidWrite ||= writing;
    discard();
  };
  const batch = app.entitySource.value.currentRow.value;
  const sub = edit.onChanged.subscribe(() => {
    if (compound.currentRow.peek()?.name === 'Aspirin')
      app.entitySource.value.edit.peek().setValue(batch.id, 'compound_id', 'c1');
  });
  compound.currentRow.value.name = 'Aspirin';
  await flush();
  await flush();
  assert.equal(discardedMidWrite, false, 'the platform editor would finish that write on a removed row');
  assert.equal(compound.isDirty.value, false);
  assert.notEqual(compound.currentRow.value, null, 'a pristine new compound is ready again');
  assert.equal(app.session.changeCount.value, 1, 'only the batch is a change');
  assert.equal(await app.session.save(), true);
  assert.deepEqual((await rowsOf(backend, 'compound')).map((r) => r.name), ['Aspirin']);
  sub.unsubscribe();
  app.dispose();
});

scoped('a create page opened again with other values takes them: a second child tab\'s New', async () => {
  const {app} = await batches();
  await app.create({amount: 1});
  await flush();
  assert.equal(app.entitySource.value.currentRow.value.amount, 1);
  await app.create({amount: 2});
  await flush();
  assert.equal(app.entitySource.value.currentRow.value.amount, 2, 'not the first parent\'s draft');
  app.dispose();
});

scoped('the entity page menu holds the row\'s own actions; a table nobody searches has no search box', async () => {
  const {app} = await batches();
  assert.equal(app.ribbon().tools.length, 0, 'batches have no searchable column');
  await app.create();
  await flush();
  assert.deepEqual(app.menuActions(), [], 'a draft has nothing to delete');
  app.dispose();
  const compounds = domains.app({table: await domains.table('lab.compound'), base: '/apps/T/Compounds'});
  document.body.append(compounds.root);
  await flush();
  assert.equal(compounds.ribbon().tools.length, 1);
  await compounds.goTo('entity', 'c1');
  await flush();
  assert.ok(compounds.menuActions().some((a) => a.name === 'Delete'));
  assert.ok(!compounds.menuActions().some((a) => a.name === 'Import…'), 'table-wide actions stay on the list');
  compounds.dispose();
});
