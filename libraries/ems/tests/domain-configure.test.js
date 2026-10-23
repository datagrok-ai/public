/* `domains.configure` on the memory backend: what an app declares on a table — validators,
   renderer, actions — dresses every handle `domains.table` makes for the address, the ones the
   library opens itself (a child tab) included; an async `fn` is awaited, a throw rejects the
   handle, and the returned function unregisters. */

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

/** What Grit declares on issues: a rule that reads another column, a caption, an action. */
async function dress(issues) {
  await Promise.resolve();
  issues.validators.add('done', (value, row) => value === true && !row.priority ? 'Prioritise before closing' : null);
  issues.renderer = {...issues.renderer, caption: (row) => `#${row.number} ${row.title}`};
  issues.actions.add({name: 'Close', run: (row) => { row.done = true; }});
}

scoped('every handle over the address is configured, the child tab\'s too; unregistering stops it', async () => {
  const memory = backend();
  emsBackends.domain = memory;
  const unregister = domains.configure('grit.issue', dress);
  try {
    const first = await domains.table('grit.issue');
    const second = await domains.table('grit.issue');
    assert.notEqual(first, second, 'a handle per call');
    for (const handle of [first, second]) {
      assert.deepEqual(handle.validators.columns, ['done'], 'once each, never doubled');
      assert.equal(handle.renderer.caption({id: 'i1', number: 1, title: 'Aspirin'}), '#1 Aspirin');
    }
    assert.deepEqual((await domains.table('grit.project')).validators.columns, [], 'another table is left alone');

    const projects = await domains.table('grit.project');
    const parent = projects.source({});
    await flush();
    parent.currentRow.value = parent.rows.byKey('p1');
    const children = domains.children(parent, {mode: 'list'});
    await flush();
    await flush();
    const tab = children.entries.value.find((e) => e.table.address === 'grit.issue');
    assert.deepEqual(tab.table.validators.columns, ['done'], 'the tab the library opened carries the rule');
    assert.deepEqual(tab.table.actions.for({id: 'i2', done: false}).map((a) => a.name), ['Close']);

    const issues = children.child(tab.id);
    await flush();
    issues.edit.value.setValue('i2', 'done', true);
    assert.equal(await issues.session.save(), false, 'a write from the child tab is refused by the app\'s rule');
    assert.match(String(issues.error.value), /#2 Ibuprofen: Done: Prioritise before closing/);
    assert.deepEqual(memory.tableSync('grit.issue').rows.map((r) => r.done), [true, false, false]);
    issues.session.discard();
    children.dispose();
    parent.dispose();
  } finally {
    unregister();
  }
  const bare = await domains.table('grit.issue');
  assert.deepEqual(bare.validators.columns, [], 'unregistered: the next handle is bare');
  assert.equal(bare.renderer, bare.defaultRenderer);
});

scoped('one registered again under the same id replaces the earlier: applied once', async () => {
  emsBackends.domain = backend();
  const first = domains.configure('grit.issue', dress, {id: 'grit'});
  // the same bundle evaluated again: an equal function, not the same one
  const second = domains.configure('grit.issue', (issues) => dress(issues), {id: 'grit'});
  try {
    const handle = await domains.table('grit.issue');
    assert.deepEqual(handle.validators.columns, ['done']);
    assert.deepEqual(handle.actions.for({id: 'i2', done: false}).map((a) => a.name), ['Close'], 'one action, not two');
    first();
    assert.deepEqual((await domains.table('grit.issue')).validators.columns, ['done'],
      'the replaced entry\'s unregister leaves the one that replaced it');
  } finally {
    second();
  }
  assert.deepEqual((await domains.table('grit.issue')).validators.columns, []);
});

scoped('a configure that throws rejects the handle', async () => {
  emsBackends.domain = backend();
  const unregister = domains.configure('grit.issue', async () => {
    throw new Error('lookups did not load');
  });
  try {
    await assert.rejects(domains.table('grit.issue'), /lookups did not load/);
  } finally {
    unregister();
  }
  assert.equal((await domains.table('grit.issue')).address, 'grit.issue');
});
