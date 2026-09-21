/* `live: true` (WO 3-11) on fake timers: one probe per interval — the collection counted and
   dated in ONE call under the source's own filter, search and `deleted` mode — refreshing the
   rows while there is nothing to lose, marking the source stale while the session has changes of
   its own, and asking nothing at all while the tab is hidden, while a batch is being written
   back, or after the source is disposed. */

import {test, mock} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {DomainSource} from '../src/sources/domain-source.js';
import {backend} from './domain-fixtures.mjs';

const env = {designTime: false, subBinds: {}, resolve: () => null};

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    mock.timers.enable({apis: ['setInterval']});
    try {
      await body();
    } finally {
      mock.timers.reset();
      document.hidden = false;
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** The memory backend with its probe counted — what a live source is supposed to ask, and how often. */
function counted(options = {}) {
  const inner = backend(options);
  const probes = [];
  return {probes, backend: {
    table: async (address) => {
      const t = await inner.table(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {
        probe: (spec) => {
          probes.push(spec);
          return t.probe(spec);
        },
      });
    },
    saveAll: (edits) => inner.saveAll(edits),
  }};
}

/** Starts the source over `be` and waits for its first load; the backend stays in place, since
 * `save()` goes looking for it through the session. */
async function started(be, options = {}) {
  backends.domain = be.backend;
  const src = new DomainSource({table: 'grit.issue', pageSize: 10, live: true, liveMs: 1000, ...options}, env);
  src.start();
  await flush();
  return src;
}

/** One interval, plus the microtasks the probe's promise chain needs. */
async function tick(times = 1, ms = 1000) {
  for (let i = 0; i < times; i++) {
    mock.timers.tick(ms);
    await flush();
    await flush();
  }
}

scoped('the first probe only baselines; a moved pair reloads a clean source', async () => {
  const be = counted();
  const src = await started(be);
  assert.equal(src.rows.items.value.length, 3);
  assert.deepEqual(be.probes, [], 'nothing is asked before the first interval');

  await tick();
  assert.equal(be.probes.length, 1, 'ONE call per interval, not a page of rows');
  assert.deepEqual(Object.keys(be.probes[0]).sort(), ['filter'], 'the one read scope, the default trash mode left off');
  assert.equal(src.stale.value, false);

  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 9, title: 'New one'}}]);
  await tick();
  assert.equal(be.probes.length, 2);
  await flush();
  assert.equal(src.rows.items.value.length, 4, 'the collection moved, so the rows were read again');
  assert.equal(src.stale.value, false, 'a source that reloaded is not stale');
  src.dispose();
});

scoped('a row edited in place moves the pair too — the count alone would miss it', async () => {
  const be = counted();
  const src = await started(be);
  await tick();
  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'update', table: 'issue', id: 'i1', values: {priority: 'high'}}]);
  await tick();
  await flush();
  assert.equal(src.rows.byKey('i1').priority, 'high');
  src.dispose();
});

scoped('a session with unsaved changes is marked stale instead, and a reload clears it', async () => {
  const be = counted();
  const src = await started(be);
  await tick();
  src.rows.items.value[0].title = 'Mine';
  await flush();
  assert.equal(src.isDirty.value, true);

  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 8, title: 'Theirs'}}]);
  await tick();
  assert.equal(src.stale.value, true, 'the rows are behind the server');
  assert.equal(src.rows.items.value.length, 3, 'and they were NOT dropped from under the edit');
  assert.equal(src.rows.byKey('i1').title, 'Mine');

  // a reload that the pending changes refuse leaves both the rows and the mark alone
  await src.refresh();
  await flush();
  assert.equal(src.stale.value, true);
  assert.equal(src.rows.items.value.length, 3);

  src.discard();
  await src.refresh();
  await flush();
  assert.equal(src.stale.value, false);
  assert.equal(src.rows.items.value.length, 4);
  src.dispose();
});

scoped('a save clears the stale mark: the re-base has just read the collection again', async () => {
  const be = counted();
  const src = await started(be);
  await tick();
  src.rows.items.value[0].title = 'Mine';
  await flush();
  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 6, title: 'Theirs'}}]);
  await tick();
  assert.equal(src.stale.value, true);

  assert.equal(await src.save(), true);
  await flush();
  assert.equal(src.stale.value, false, 'a save re-reads the window, so the mark cannot outlive it');
  assert.equal(src.rows.byKey('i1').title, 'Mine');
  src.dispose();
});

scoped('a probe that fails keeps the baseline, and a run of failures stops the timer', async () => {
  const be = counted();
  let fail = 0;
  const inner = be.backend.table;
  be.backend.table = async (address) => {
    const t = await inner(address);
    return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {
      probe: (spec) => fail > 0 ? (fail--, Promise.reject(new Error('offline'))) : t.probe(spec),
    });
  };
  const src = await started(be);
  await tick();
  const issue = await be.backend.table('grit.issue');

  // one bad tick over a change: the baseline stands, so the next good probe still sees it
  fail = 1;
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 5, title: 'Missed'}}]);
  await tick();
  assert.equal(src.rows.items.value.length, 3, 'the failed tick changed nothing');
  await tick();
  await flush();
  assert.equal(src.rows.items.value.length, 4, 'and the change was not baselined away by the outage');

  fail = DomainSource.liveFailures + 2;
  await tick(DomainSource.liveFailures);
  const asked = be.probes.length;
  await tick(3);
  assert.equal(be.probes.length, asked, 'after three failures in a row the timer gives up');
  src.live.value = false;
  src.live.value = true;
  await flush();
  fail = 0;
  await tick();
  assert.equal(be.probes.length > asked, true, 'and re-setting `live` starts a new one');
  src.dispose();
});

scoped('nothing is asked while the tab is hidden, and nothing after dispose', async () => {
  const be = counted();
  const src = await started(be);
  document.hidden = true;
  await tick(3);
  assert.deepEqual(be.probes, [], 'a background tab is throttled and has no one reading it');
  document.hidden = false;
  await tick();
  assert.equal(be.probes.length, 1);

  src.dispose();
  await tick(3);
  assert.equal(be.probes.length, 1, 'the timer went with the source');
});

scoped('live is a signal: turning it off stops the timer, turning it on starts one', async () => {
  const be = counted();
  const src = await started(be, {live: false});
  await tick(2);
  assert.deepEqual(be.probes, [], 'a source that is not live never probes');
  src.live.value = true;
  await flush();
  await tick();
  assert.equal(be.probes.length, 1);
  src.live.value = false;
  await flush();
  await tick(2);
  assert.equal(be.probes.length, 1);
  src.dispose();
});

scoped('the probe carries the source\'s own filter, search and deleted mode', async () => {
  const be = counted();
  const src = await started(be, {query: 'done = false', search: 'pro', deleted: 'include'});
  await tick();
  assert.deepEqual(be.probes[0], {filter: 'done = false', search: 'pro', deleted: 'include'});
  src.dispose();
});

scoped('a backend that cannot probe is simply not polled', async () => {
  const inner = backend();
  backends.domain = {
    table: async (address) => {
      const t = await inner.table(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {probe: undefined});
    },
    saveAll: (edits) => inner.saveAll(edits),
  };
  const src = new DomainSource({table: 'grit.issue', pageSize: 10, live: true, liveMs: 1000}, env);
  src.start();
  await flush();
  await tick(3);
  assert.equal(src.stale.value, false);
  assert.equal(src.rows.items.value.length, 3, 'the rows are the ones that were loaded');
  src.dispose();
});

/* The app's half of WO 3-11: the status bar carries the stale mark, and `refresh()` is what the
   word "Refresh" in it stands for — through the unsaved gate, since a reload drops the changes. */
register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {DomainApp} = await import('../src/dg/domain/app.js');
const {DgDomainBackend} = await import('../src/dg/domain/backend.js');
const grok = await import('datagrok-api/grok');

/** A `DgDomainTable` over a stubbed js-api client, with what the poll asked for counted. */
async function platformTable(support = {}, can = {}) {
  const calls = {version: 0, aggregate: []};
  const client = {
    access: async () => ({can: {view: true, ...can}, fields: {}, support: {
      systemColumns: ['id', 'version', 'created_on', 'updated_on', 'author_id'], writes: true,
      deleted: true, restore: true, audit: true, ancestors: false, probe: true, version: true, watch: true,
      updateWhere: true, captions: true, transaction: true, concurrency: 'version', filters: 'full',
      batch: {upsert: true, partial: true, validate: true, skipDuplicates: true},
      ...support}}),
    version: async () => ({seq: ++calls.version, at: '2026-09-16T10:05:13Z'}),
    aggregate: async (spec) => (calls.aggregate.push(spec), [{count: 3, last: '2026-09-16T10:00:00Z'}]),
  };
  grok.dapi.domains.table = () => client;
  grok.dapi.domains.registry = {rowProperties: async () => [], tableInfo: async () => ({
    nameColumn: 'title', singularName: 'issue', pluralName: 'issues', businessKey: [],
    searchableColumns: [], constraints: [], refFilters: {}, permissions: [], childTables: []})};
  return {calls, table: await new DgDomainBackend().table('grit.issue')};
}

/** The platform's probe adapter mirrored over the memory backend, so the source can be driven
 * against the token shape headlessly; the adapter itself is pinned by the test above. */
function tokened(options = {}) {
  const inner = backend(options);
  const calls = {token: 0, aggregate: 0};
  const state = {seq: 1};
  return {calls, state, backend: {
    table: async (address) => {
      const t = await inner.table(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {
        probe: async (scope) => {
          if ((scope.filter ?? '') === '' && (scope.search ?? '') === '' && (scope.deleted ?? 'exclude') === 'exclude') {
            calls.token++;
            return {count: -1, last: String(state.seq)};
          }
          calls.aggregate++;
          return t.probe(scope);
        },
      });
    },
    saveAll: (edits) => inner.saveAll(edits),
  }};
}

scoped('the platform poll asks the change token over the whole table and aggregates over anything narrower',
  async () => {
    const {calls, table} = await platformTable();
    try {
      assert.deepEqual(await table.probe(), {count: -1, last: '1'}, 'not counted: the token answers alone');
      assert.deepEqual(await table.probe({filter: undefined, search: undefined, deleted: undefined}),
        {count: -1, last: '2'});
      assert.deepEqual(await table.probe({filter: [], deleted: 'exclude'}), {count: -1, last: '3'},
        'an empty condition tree is what an empty filter builder compiles to — still the whole table');
      assert.deepEqual(calls.aggregate, [], 'and none of it costs a scan');

      assert.deepEqual(await table.probe({filter: 'done = false'}), {count: 3, last: '2026-09-16T10:00:00Z'});
      await table.probe({search: 'pro'});
      await table.probe({deleted: 'only'});
      assert.equal(calls.version, 3, 'a narrowed read never asks the token: it says nothing about a subset');
      assert.equal(calls.aggregate.length, 3);
      assert.deepEqual(calls.aggregate[0].measures,
        [{fn: 'count'}, {fn: 'max', column: 'updated_on', as: 'last'}]);
    } finally {
      delete grok.dapi.domains.table;
      delete grok.dapi.domains.registry;
    }
  });

scoped('the token is skipped where it cannot answer: no table-level View, or a token the engine does not move',
  async () => {
    for (const [support, can, why] of [
      [{}, {view: false}, 'a reader who reaches rows through row grants only'],
      [{version: false}, {}, 'a registration the platform writes, whose token never moves'],
    ]) {
      const {calls, table} = await platformTable(support, can);
      try {
        assert.deepEqual(await table.probe(), {count: 3, last: '2026-09-16T10:00:00Z'}, why);
        assert.equal(calls.version, 0, `${why}: the token is never asked`);
        assert.equal(calls.aggregate.length, 1, `${why}: the aggregate over what the caller sees answers`);
      } finally {
        delete grok.dapi.domains.table;
        delete grok.dapi.domains.registry;
      }
    }
  });

scoped('a token refused mid-session hands the poll to the aggregate for good', async () => {
  const {calls, table} = await platformTable();
  try {
    assert.deepEqual(await table.probe(), {count: -1, last: '1'});
    grok.dapi.domains.table().version = async () => { throw new Error('403 forbidden'); };
    assert.deepEqual(await table.probe(), {count: 3, last: '2026-09-16T10:00:00Z'}, 'the same poll answers');
    await table.probe();
    assert.equal(calls.aggregate.length, 2, 'and the token is not asked again');
  } finally {
    delete grok.dapi.domains.table;
    delete grok.dapi.domains.registry;
  }
});

scoped('a change token that moves refreshes a clean source and marks a dirty one stale', async () => {
  const be = tokened();
  const src = await started(be);
  await tick();
  assert.equal(be.calls.token, 1, 'one token read per interval');
  await tick();
  assert.equal(src.stale.value, false, 'a token that did not move is not a change');

  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 9, title: 'Theirs'}}]);
  be.state.seq++;
  await tick();
  await flush();
  assert.equal(src.rows.items.value.length, 4, 'the token moved, so the rows were read again');

  // the reload re-baselines, so the next poll is the one that can see a move again
  await tick();
  src.rows.items.value[0].title = 'Mine';
  await flush();
  be.state.seq++;
  await tick();
  assert.equal(src.stale.value, true, 'a session with changes of its own is marked, not reloaded');
  assert.equal(be.calls.aggregate, 0, 'nothing was counted along the way');
  src.dispose();
});

scoped('narrowing the query re-baselines: the shapes never meet in one comparison', async () => {
  const be = tokened();
  const src = await started(be);
  await tick();
  assert.equal(be.calls.token, 1);

  src.query.value = 'done = false';
  await flush();
  await tick();
  assert.equal(be.calls.aggregate, 1, 'a narrowed source counts');
  assert.equal(src.stale.value, false, 'and the first poll after the change only baselines');
  const rows = src.rows.items.value.length;
  await tick();
  assert.equal(src.rows.items.value.length, rows, 'no false change out of the shape swap');

  src.query.value = '';
  await flush();
  await tick();
  assert.equal(be.calls.token, 2, 'back to the whole table, back to the token');
  assert.equal(src.stale.value, false);
  src.dispose();
});

const dialogButton = (text) =>
  [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);

scoped('the app says "Data changed — Refresh" and its Refresh goes through the unsaved gate', async () => {
  const be = counted();
  backends.domain = be.backend;
  const table = await domains.table('grit.issue');
  const app = domains.app({table, base: '/apps/T/Issues', pageSize: 10});
  document.body.append(app.root);
  await flush();
  app.listSource.live.value = true;
  await flush();
  // the app's own source keeps the 30 s default
  await tick(1, DomainSource.liveMs);

  app.listSource.rows.items.value[0].title = 'Mine';
  await flush();
  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 7, title: 'Theirs'}}]);
  await tick(1, DomainSource.liveMs);
  assert.equal(app.stale.value, true);
  assert.equal(app.summary.value, '3 issues — 1 unsaved change — Data changed — Refresh');

  const done = app.refresh();
  await flush();
  fire(dialogButton('DISCARD'), 'click');
  assert.equal(await done, true);
  await flush();
  assert.equal(app.stale.value, false);
  assert.equal(app.listSource.rows.items.value.length, 4);
  assert.equal(app.summary.value, '4 issues');
  app.dispose();
  assert.equal(DomainApp.live.size, 0, 'no app left behind');
});

scoped('DomainAppOptions.live reaches the list source — the option the apps actually pass', async () => {
  const be = counted();
  backends.domain = be.backend;
  const table = await domains.table('grit.issue');
  const app = domains.app({table, base: '/apps/T/Issues', pageSize: 10, live: true, liveMs: 1000});
  document.body.append(app.root);
  await flush();
  assert.equal(app.listSource.live.value, true, 'the app asked for it, so the list polls');
  await tick();
  assert.equal(be.probes.length, 1);

  // the simulated user's scenario: a draft pending here, an insert from another session there
  app.listSource.newRow({title: 'Mine', project_id: 'p1'});
  await flush();
  const issue = await be.backend.table('grit.issue');
  await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', number: 4, title: 'Theirs'}}]);
  await tick();
  assert.equal(app.stale.value, true);
  assert.match(app.summary.value, /Data changed — Refresh/);
  app.dispose();
  assert.equal(DomainApp.live.size, 0, 'no app left behind');
});

scoped('the caption follows the values a save writes back into the row', async () => {
  const be = counted();
  backends.domain = be.backend;
  const table = await domains.table('grit.issue');
  const app = domains.app({table, base: '/apps/T/Issues', pageSize: 10});
  document.body.append(app.root);
  await flush();
  assert.equal(await app.goTo('entity', 'i1'), true);
  await flush();
  const crumb = () => [...app.breadcrumbs.root.querySelectorAll('.u2-breadcrumbs-item, .u2-breadcrumbs-current')]
    .map((el) => el.textContent).pop();
  assert.equal(crumb(), 'Aspirin');

  // the row is a keyed proxy: a write-back moves its cells without `currentRow` ever changing,
  // and the caption used to keep whatever it read first (a uuid, for a row saved with no name yet)
  const row = app.entitySource.value.currentRow.value;
  row.title = 'Aspirin 100';
  await flush();
  assert.equal(crumb(), 'Aspirin 100');
  assert.equal(await app.session.save(), true);
  await flush();
  assert.equal(crumb(), 'Aspirin 100', 'and it still says so after the batch and its re-base');
  assert.equal(app.summary.value, 'Aspirin 100', 'the status line reads the same row');
  app.dispose();
});

scoped('the trash reads newest-deleted first, and leaving it gives the order back', async () => {
  const be = counted();
  backends.domain = be.backend;
  const table = await domains.table('grit.issue');
  const app = domains.app({table, base: '/apps/T/Issues', pageSize: 10});
  document.body.append(app.root);
  await flush();
  assert.equal(app.listSource.sort.value, '', 'a live list is in the table order');

  assert.equal(await app.setTrash(true), true);
  await flush();
  assert.equal(app.listSource.sort.value, '!updated_on', 'a soft delete stamps updated_on');
  assert.equal(app.listSource.deleted.value, 'only');

  assert.equal(await app.setTrash(false), true);
  await flush();
  assert.equal(app.listSource.sort.value, '');
  app.dispose();
  assert.equal(DomainApp.live.size, 0, 'no app left behind');
});
