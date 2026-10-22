/* The app's entity page over the memory backend: a create page opened again with other values, the
   page's ⋯ menu and a ribbon without a search box, `immutable` columns through `domains.app`, and a
   reference the client set reading its target's name once saved. */

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
  const backend = new MemoryDomainBackend(LAB,
    {rows: {compound: [{id: 'c1', name: 'Aspirin'}], batch: [{id: 'b1', compound_id: 'c1', amount: 1}]}});
  emsBackends.domain = backend;
  const app = domains.app({table: await domains.table('lab.batch'), base: '/apps/T/Batches',
    immutable: ['compound_id']});
  document.body.append(app.root);
  await flush();
  return {backend, app};
}

const readonlyValue = (app, column) =>
  app.root.querySelector(`[data-u2-name="${column}"] [data-u2-part="readonly-value"]`);

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

scoped('an immutable column of the app is picked on a new row and read-only on a saved one', async () => {
  const {app} = await batches();
  await app.goTo('entity', 'b1');
  await flush();
  await flush();
  assert.equal(app.form.value.input('compound_id') === undefined, true,
    'a saved batch cannot be moved to another compound');
  assert.equal(readonlyValue(app, 'compound_id')?.textContent, 'Aspirin', 'it reads as the compound it names');
  assert.equal(app.form.value.input('amount') !== undefined, true, 'the other columns stay editable');
  await app.create();
  await flush();
  assert.equal(app.form.value.input('compound_id') !== undefined, true, 'a new batch picks its compound');
  app.dispose();
});

scoped('a reference the client set reads as its target\'s name once the row is saved', async () => {
  const {backend, app} = await batches();
  await app.create();
  await flush();
  const batch = app.entitySource.value.currentRow.value;
  batch.compound_id = 'c1';
  batch.amount = 2.5;
  assert.equal(await app.session.save(), true);
  await flush();
  await flush();
  const saved = (await (await backend.table('lab.batch')).query({})).find((r) => r.amount === 2.5);
  assert.equal(saved.compound_id, 'c1');
  assert.equal(readonlyValue(app, 'compound_id')?.textContent, 'Aspirin', 'not blank, not the id');
  app.dispose();
});
