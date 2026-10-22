/* `table.derived` over the memory backend: columns an app computes in the browser from other
   columns of the same row — derived once per change of the inputs, waited for by every Save,
   refusing the Save while the derivation fails, dropped when the inputs moved on, read-only in a
   form. */

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

const stored = (id) => emsBackends.domain.tableSync('grit.issue').rows.find((r) => r.id === id);
const later = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

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

scoped('an edit to an input derives the outputs once, and the Save sends them', async () => {
  const {source, calls} = await issues((row) => ({description: `About ${row.title}`}));
  const row = source.rows.byKey('i2');
  row.title = 'Ibuprofen 200';
  await flush();
  assert.equal(row.description, 'About Ibuprofen 200');
  row.weight = 3;
  await flush();
  assert.deepEqual(calls, ['Ibuprofen 200'], 'an edit to another column derives nothing');
  assert.equal(await source.save(), true);
  assert.equal(stored('i2').description, 'About Ibuprofen 200');
  source.dispose();
});

scoped('a Save waits for the derivation in flight', async () => {
  const {source} = await issues(async (row) => {
    await later(20);
    return {description: `Slow ${row.title}`};
  });
  source.rows.byKey('i1').title = 'Aspirin 81';
  assert.equal(await source.save(), true);
  assert.equal(stored('i1').description, 'Slow Aspirin 81', 'the derived value rode the same Save');
  source.dispose();
});

scoped('a failing derivation refuses the Save, naming the row and the input', async () => {
  const {source} = await issues((row) => {
    if (row.title.includes('?'))
      throw new Error('Cannot read the title');
    return {description: row.title.toUpperCase()};
  });
  const row = source.rows.byKey('i2');
  row.title = 'Ibu?';
  await flush();
  assert.equal(await source.save(), false);
  assert.equal(source.error.value.message, 'Cannot save: Ibu?: Title: Cannot read the title');
  assert.equal(source.problemRow.value, 'i2');
  row.title = 'Ibuprofen 400';
  assert.equal(await source.save(), true, 'the next good input takes the refusal back');
  assert.equal(stored('i2').description, 'IBUPROFEN 400');
  source.dispose();
});

scoped('a result whose inputs moved on is dropped', async () => {
  const {source} = await issues(async (row) => {
    const title = row.title;
    await later(title === 'First' ? 40 : 5);
    return {description: title};
  });
  const row = source.rows.byKey('i3');
  row.title = 'First';
  await flush();
  row.title = 'Second';
  await source.settled();
  assert.equal(row.description, 'Second', 'the late first result did not overwrite the second');
  source.dispose();
});

scoped('a draft is derived before its insert', async () => {
  const {source} = await issues((row) => ({description: `New: ${row.title}`}));
  const draft = source.newRow({project_id: 'p1', title: 'Paracetamol'});
  await flush();
  assert.equal(await source.save(), true);
  const inserted = emsBackends.domain.tableSync('grit.issue').rows.find((r) => r.title === 'Paracetamol');
  assert.equal(inserted.description, 'New: Paracetamol');
  assert.notEqual(draft.id, undefined);
  source.dispose();
});

scoped('a form shows the derived columns as read-only text', async () => {
  const {source} = await issues((row) => ({description: row.title}));
  source.currentRow.value = source.rows.byKey('i1');
  const form = domains.form(source);
  document.body.append(form.root);
  await flush();
  assert.notEqual(form.input('title'), undefined, 'an input stays an input');
  assert.equal(form.input('description'), undefined, 'a derived column is text, not an input');
  form.dispose();
  source.dispose();
});

scoped('every input cleared clears the outputs and the failure', async () => {
  const {source} = await issues((row) => {
    if (row.title.includes('?'))
      throw new Error('Cannot read the title');
    return {description: row.title.toUpperCase()};
  });
  const row = source.rows.byKey('i2');
  row.title = 'Ibuprofen 600';
  await source.settled();
  assert.equal(row.description, 'IBUPROFEN 600');
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

scoped('its own writes and an unrelated edit in flight derive nothing more', async () => {
  const {source, calls} = await issues(async (row) => {
    await later(20);
    return {description: `Slow ${row.title}`};
  });
  const row = source.rows.byKey('i3');
  row.title = 'Once';
  await flush();
  row.weight = 9;
  await source.settled();
  await flush();
  assert.equal(row.description, 'Slow Once');
  assert.deepEqual(calls, ['Once']);
  source.dispose();
});

scoped('a failing derivation shows on its input field as soon as it ends, not first at Save', async () => {
  const {source} = await issues((row) => {
    if (row.title.includes('?'))
      throw new Error('Cannot read the title');
    return {description: row.title.toUpperCase()};
  });
  source.currentRow.value = source.rows.byKey('i2');
  const form = domains.form(source);
  document.body.append(form.root);
  await flush();
  const title = form.input('title');
  const shown = () => title.root.querySelector('[data-u2-part="error"]').textContent;
  title.value.value = 'Ibu?';
  await source.settled();
  await flush();
  assert.equal(title.validity.value, 'Cannot read the title', 'no Save was pressed');
  assert.equal(shown(), 'Cannot read the title');
  assert.equal(form.problem.value, 'Title: Cannot read the title');
  title.value.value = 'Ibuprofen';
  await source.settled();
  await flush();
  assert.equal(title.validity.value, null, 'the next good run takes it back');
  source.rows.byKey('i2').title = 'Naproxen?';
  await source.settled();
  await flush();
  assert.equal(title.validity.value, 'Cannot read the title', 'an input written past the field fails on it too');
  form.dispose();
  source.dispose();
});

scoped('a create form: the failure is on the field once the derivation ends; the draft keeps it until fixed', async () => {
  const {table, source} = await issues((row) => {
    if (row.title.includes('?'))
      throw new Error('This title is taken by GRIT-1');
    return {description: row.title};
  });
  const draft = table.draft({project_id: 'p1'});
  await flush();
  const form = domains.form(draft);
  document.body.append(form.root);
  await flush();
  const title = form.input('title');
  title.value.value = 'Aspirin?';
  await draft.settled();
  await flush();
  assert.equal(title.validity.value, 'This title is taken by GRIT-1');
  assert.equal(title.root.classList.contains('u2-input-untouched'), false, 'an edited field shows its verdict');
  assert.equal(await draft.save(), false, 'and Save still refuses it');
  form.dispose();
  draft.dispose();
  source.dispose();
});
