import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Registry, registerAll} from '@datagrok-libraries/u2';
import {backend, until} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {registerDomainComponents} = await import('../src/dg/registrations.js');

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
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);

async function issues() {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  return src;
}

function type(input, value) {
  input.value = value;
  fire(input, 'input');
}

scoped('typing searches after the pause; Enter at once; Escape clears; the total follows', async () => {
  const src = await issues();
  const search = testOwner().build(() => domains.search(src, {debounceMs: 20}));
  const input = search.root.querySelector('input');
  assert.equal(search.root.dataset.u2, 'domain-search');
  assert.equal(input.placeholder, 'Search title, description…', 'names what it searches');
  type(input, 'asp');
  assert.equal(src.search.value, '', 'not yet');
  await until(() => src.search.value === 'asp');
  assert.equal(src.search.value, 'asp');
  await flush();
  assert.equal(src.total.value, 1, 'the memory backend searched the title');
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Aspirin']);
  type(input, 'ibu');
  fire(input, 'keydown', {key: 'Enter'});
  // the gate answers at once over a clean session, a microtask later
  await flush();
  assert.equal(src.search.value, 'ibu', 'Enter does not wait for the pause');
  fire(input, 'keydown', {key: 'Escape'});
  assert.equal(input.value, '');
  await flush();
  assert.equal(src.search.value, '');
  await flush();
  assert.equal(src.total.value, 3);
  search.dispose();
  src.dispose();
});

scoped('two-way: a search set in code shows in the box and cancels a pending keystroke', async () => {
  const src = await issues();
  const search = testOwner().build(() => domains.search(src, {debounceMs: 20}));
  const input = search.root.querySelector('input');
  type(input, 'pending');
  src.search.value = 'nap';
  assert.equal(input.value, 'nap');
  await wait(30);
  assert.equal(src.search.value, 'nap', 'the keystroke that was pending never lands');
  search.dispose();
  src.dispose();
});

scoped('spec: u2-domain-search is registered with usage', () => {
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const meta = reg.get('u2-domain-search');
  assert.equal(meta.usage.length > 0, true);
  assert.deepEqual(meta.props.slice(0, 3).map((p) => p.name), ['source', 'placeholder', 'debounceMs']);
});

scoped('a search while the session is dirty asks: cancel puts the text back, discard searches', async () => {
  const src = await issues();
  const search = testOwner().build(() => domains.search(src, {debounceMs: 20}));
  const input = search.root.querySelector('input');
  src.rows.byKey('i1').title = 'Edited';
  assert.equal(src.isDirty.value, true);
  type(input, 'nap');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.notEqual(document.body.querySelector('.u2-dialog'), null, 'the gate asks');
  assert.equal(src.search.value, '', 'nothing written yet');
  fire(buttonNamed('CANCEL'), 'click');
  await flush();
  assert.equal(input.value, '', 'the box shows the search in force again');
  assert.equal(src.isDirty.value, true, 'the edit is kept');
  type(input, 'nap');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(src.search.value, 'nap');
  assert.equal(src.isDirty.value, false);
  await flush();
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Naproxen'], 're-queried');
  search.dispose();
  src.dispose();
});

scoped('the box carries the class its sheet selects, and a search during a write-back is refused', async () => {
  const src = await issues();
  const search = testOwner().build(() => domains.search(src, {debounceMs: 20}));
  assert.equal(search.root.classList.contains('u2-domain-search'), true);
  const input = search.root.querySelector('input');
  let release;
  const real = emsBackends.domain.saveAll.bind(emsBackends.domain);
  emsBackends.domain.saveAll = async (edits) => {
    const landed = await real(edits);
    await new Promise((resolve) => release = resolve);
    return landed;
  };
  src.rows.byKey('i1').title = 'Edited';
  const saving = src.save();
  await flush();
  assert.equal(src.isDirty.value, false, 'the batch landed: the rows read clean already');
  type(input, 'nap');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(src.search.value, '', 'the gate refuses through the write-back, clean or not');
  assert.match(document.body.querySelector('.u2-notify-warning')?.textContent ?? '', /Wait for the batch/);
  release();
  assert.equal(await saving, true);
  search.dispose();
  src.dispose();
});

// the source over the issue table, its frame reads counted, the current row i3, one edit pending, a search typed
async function searchedOverEdits() {
  emsBackends.domain = backend();
  const inner = emsBackends.domain.table.bind(emsBackends.domain);
  const reads = {frames: 0};
  emsBackends.domain.table = async (address) => {
    const t = await inner(address);
    return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {frame: (spec) => (reads.frames++, t.frame(spec))});
  };
  const table = await domains.table('grit.issue');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  await wait(10);
  src.currentRow.value = src.rows.byKey('i3');
  src.rows.byKey('i1').title = 'Edited';
  src.search.value = 'Ibuprofen';
  await flush();
  await wait(20);
  reads.frames = 0;
  return {src, reads};
}

scoped('a search changed over pending edits applies when Discard makes the source clean (EQ-008)', async () => {
  const {src, reads} = await searchedOverEdits();
  assert.equal(src.rows.items.value.some((r) => r.title === 'Edited'), true, 'held while dirty');
  src.discard();
  await flush();
  await wait(20);
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Ibuprofen'],
    `after discard the deferred search ran: rows ${JSON.stringify(src.rows.items.value.map((r) => r.title))}`);
  assert.equal(reads.frames, 1, 'one read');
  const current = src.currentRow.value?.id;
  assert.equal(current === undefined || src.rows.indexOf(current) >= 0, true, 'the current row is in the rows shown or none');
  src.dispose();
});

scoped('a search changed over pending edits is applied by the Save\'s own re-base: one read, the current row kept (EQ-008)', async () => {
  const {src, reads} = await searchedOverEdits();
  src.currentRow.value = src.rows.byKey('i2');
  assert.equal(await src.save(), true);
  await flush();
  await wait(30);
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Ibuprofen'], 'the pending search applies');
  assert.equal(reads.frames, 1, 'exactly one read of the rows, no second refresh');
  assert.equal(src.currentRow.value?.id, 'i2', 'the current row is kept');
  // the re-base recorded the scope it read: a later Discard has no search left to apply
  src.rows.byKey('i2').title = 'Edited again';
  src.discard();
  await flush();
  await wait(30);
  assert.equal(reads.frames, 1, 'Discard after the save reads nothing: the search was applied once');
  assert.equal(src.currentRow.value?.id, 'i2');
  src.dispose();
});

for (const [what, change, expected] of [
  ['search', (src) => src.search.value = 'Ibuprofen', ['Ibuprofen']],
  ['query', (src) => src.query.value = 'title = "Naproxen"', ['Naproxen']],
]) {
  scoped(`a ${what} changed while a Save's re-base read is in flight is applied when it lands`, async () => {
    emsBackends.domain = backend();
    const inner = emsBackends.domain.table.bind(emsBackends.domain);
    const reads = {frames: 0, hold: false, waiting: []};
    emsBackends.domain.table = async (address) => {
      const t = await inner(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {frame: async (spec) => {
        reads.frames++;
        const f = await t.frame(spec);
        if (reads.hold)
          await new Promise((go) => reads.waiting.push(go));
        return f;
      }});
    };
    const table = await domains.table('grit.issue');
    const src = testOwner().build(() => table.source({pageSize: 10}));
    await flush();
    await wait(10);
    src.rows.byKey('i1').title = 'Edited';
    reads.frames = 0;
    reads.hold = true;
    const saving = src.save();
    await until(() => reads.waiting.length > 0);
    change(src);
    await flush();
    reads.hold = false;
    reads.waiting.splice(0).forEach((go) => go());
    assert.equal(await saving, true);
    await flush();
    await wait(30);
    assert.deepEqual(src.rows.items.value.map((r) => r.title), expected, `the rows of the new ${what}`);
    assert.equal(reads.frames, 2, 'the re-base read, then one read for the new scope');
    src.dispose();
  });
}

scoped('a search on a table with no searchable column is dropped and the list read once', async () => {
  emsBackends.domain = backend();
  const inner = emsBackends.domain.table.bind(emsBackends.domain);
  const specs = [];
  emsBackends.domain.table = async (address) => {
    const t = await inner(address);
    return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {info: {...t.info, searchableColumns: []},
      frame: (spec) => (specs.push(spec), t.frame(spec))});
  };
  const table = await domains.table('grit.issue');
  const src = testOwner().build(() => table.source({pageSize: 10, search: 'abc'}));
  await until(() => src.state.value === 'ready');
  await wait(10);
  assert.equal(src.search.value, '', 'the search is dropped');
  assert.equal(specs.length, 1, 'the rows are read once');
  src.dispose();
});

scoped('a search typed and cleared while its read is in flight leaves the unsearched rows shown', async () => {
  emsBackends.domain = backend();
  const inner = emsBackends.domain.table.bind(emsBackends.domain);
  const gate = {hold: false, waiting: []};
  emsBackends.domain.table = async (address) => {
    const t = await inner(address);
    return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {frame: async (spec) => {
      const f = await t.frame(spec);
      if (gate.hold)
        await new Promise((go) => gate.waiting.push(go));
      return f;
    }});
  };
  const table = await domains.table('grit.issue');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  await wait(10);
  gate.hold = true;
  src.search.value = 'Ibuprofen';
  await until(() => gate.waiting.length === 1);
  src.search.value = '';
  await flush();
  gate.hold = false;
  gate.waiting.splice(0).forEach((go) => go());
  await wait(20);
  assert.deepEqual(src.rows.items.value.map((r) => r.id), ['i1', 'i2', 'i3'], 'the rows match the empty search');
  src.dispose();
});
