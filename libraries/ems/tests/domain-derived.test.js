/* `table.derived` over the memory backend: columns an app computes in the browser from other
   columns of the same row — derived once per change of the inputs, waited for by every Save,
   refusing the Save while the derivation fails, read-only in a form. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, notify} from '@datagrok-libraries/u2';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');

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
  });
}

const stored = (title) => emsBackends.domain.tableSync('grit.issue').rows.find((r) => r.title === title);

async function issues(derive) {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const calls = [];
  table.derived.add({from: ['title'], to: ['description'], derive: (row) => {
    calls.push(row.title);
    return derive(row);
  }});
  const source = table.source({});
  await flush();
  return {table, source, calls};
}

const failing = (row) => {
  if (row.title.includes('?'))
    throw new Error('Cannot read the title');
  return {description: row.title.toUpperCase()};
};

scoped('an edit to an input derives the outputs once, and the Save sends them — a draft\'s before its insert', async () => {
  const {source, calls} = await issues((row) => ({description: `About ${row.title}`}));
  const row = source.rows.byKey('i2');
  row.title = 'Ibuprofen 200';
  await flush();
  assert.equal(row.description, 'About Ibuprofen 200');
  row.weight = 3;
  source.newRow({project_id: 'p1', title: 'Paracetamol'});
  await flush();
  assert.deepEqual(calls, ['Ibuprofen 200', 'Paracetamol'], 'an edit to another column derives nothing');
  assert.equal(await source.save(), true);
  assert.equal(stored('Ibuprofen 200').description, 'About Ibuprofen 200');
  assert.equal(stored('Paracetamol').description, 'About Paracetamol');
  source.dispose();
});

scoped('a Save waits for the derivation in flight', async () => {
  const {source} = await issues(async (row) => {
    await new Promise((resolve) => setTimeout(resolve, 20));
    return {description: `Slow ${row.title}`};
  });
  source.rows.byKey('i1').title = 'Aspirin 81';
  assert.equal(await source.save(), true);
  assert.equal(stored('Aspirin 81').description, 'Slow Aspirin 81', 'the derived value rode the same Save');
  source.dispose();
});

scoped('a failure shows on its input field as it ends and refuses the Save; so does a refusal of the output', async () => {
  const {table, source} = await issues(failing);
  source.currentRow.value = source.rows.byKey('i2');
  const form = domains.form(source);
  document.body.append(form.root);
  await flush();
  assert.equal(form.input('description'), undefined, 'a derived column is text, not an input');
  const title = form.input('title');
  title.value.value = 'Ibu?';
  await source.settled();
  await flush();
  assert.equal(title.root.querySelector('[data-u2-part="error"]').textContent, 'Cannot read the title',
    'no Save was pressed');
  assert.equal(await source.save(), false);
  assert.equal(source.error.value.message, 'Cannot save: Ibu?: Title: Cannot read the title');
  source.rows.byKey('i2').title = 'Naproxen?';
  await source.settled();
  await flush();
  assert.equal(title.validity.value, 'Cannot read the title', 'an input written past the field fails on it too');
  table.validators.add('description', (value) => String(value).includes('TAKEN') ? 'This title is taken' : null);
  title.value.value = 'Taken';
  await source.settled();
  assert.equal(await source.save(), false);
  assert.equal(source.error.value.message, 'Cannot save: Taken: Title: This title is taken',
    'the output\'s refusal names the input the user can change');
  title.value.value = 'Ibuprofen 400';
  assert.equal(await source.save(), true, 'the next good input takes the refusal back');
  assert.equal(stored('Ibuprofen 400').description, 'IBUPROFEN 400');
  form.dispose();
  source.dispose();
});

scoped('every input cleared clears the outputs and the failure', async () => {
  const {source} = await issues(failing);
  const row = source.rows.byKey('i2');
  row.title = 'Ibuprofen 600';
  await source.settled();
  row.title = '';
  await source.settled();
  assert.equal(row.description, null, 'derived from nothing is nothing');
  row.title = 'Ibu?';
  await source.settled();
  row.title = '';
  await source.settled();
  assert.equal(await source.save(), false, 'the title is required');
  assert.doesNotMatch(source.error.value.message, /Cannot read the title/, 'the cleared input took the failure back');
  source.dispose();
});

scoped('a reverted row typed back to the same input is derived again', async () => {
  const {source, calls} = await issues((row) => ({description: `About ${row.title}`}));
  const row = source.rows.byKey('i2');
  const original = row.description ?? null;
  row.title = 'Naproxen';
  await source.settled();
  source.revert();
  await flush();
  assert.equal(row.description ?? null, original, 'the revert took the derived value back too');
  row.title = 'Naproxen';
  await source.settled();
  assert.equal(row.description, 'About Naproxen');
  assert.deepEqual(calls, ['Naproxen', 'Naproxen']);
  source.dispose();
});
