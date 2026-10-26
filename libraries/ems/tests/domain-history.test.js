/* `domains.history` (2-8) over the memory backend's audit: the current row's entries newest first,
   an update line naming the changed captions (never the version stamp), the refresh after the
   session saves, a draft's "Not saved yet" until it is, the hint without a row, and the
   `u2-domain-history` tag. `DG` and `grok` come from tests/dg-stub.mjs. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, signal, Registry, SpecContext, renderSpec, registerAll, notify} from '@datagrok-libraries/u2';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainHistory} = await import('../src/dg/history.js');
const {registerDomainComponents} = await import('../src/dg/registrations.js');
const grok = await import('datagrok-api/grok');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    const find = grok.dapi.users.find;
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      grok.dapi.users.find = find;
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** The fetch is a timer away, the render another: two flushes settle the view. */
async function settle() {
  await flush();
  await flush();
}

function lines(history) {
  return [...history.view.root.querySelectorAll('.u2-domain-history-line')];
}

/** An entry's change lines as "<caption>: <before> → <after>". */
function changes(line) {
  return [...line.querySelectorAll('.u2-domain-history-change')].map((c) => {
    const [field, before, , after] = [...c.children].map((e) => e.textContent);
    return `${field}: ${before} → ${after}`;
  });
}

scoped('entries from the audit, newest first; an update names the changed captions; refreshed on save', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  await flush();
  src.currentRow.value = src.rows.byKey('i1');
  const history = domains.history(src);
  assert.equal(history.root.dataset.u2, 'domain-history');
  assert.equal(history.row, src.currentRow, 'the source\'s current row by default');
  await settle();
  assert.equal(history.root.querySelector('[data-u2-part="hint"]').hidden, true);
  assert.equal(history.view.root.querySelector('.u2-async-empty')?.textContent, 'No history yet');

  src.rows.byKey('i1').title = 'Aspirin 2';
  src.rows.byKey('i1').weight = 3;
  assert.equal(await src.save(), true);
  await settle();
  let rows = lines(history);
  assert.equal(rows.length, 1, 'one update');
  assert.equal(rows[0].querySelector('.u2-badge').textContent, 'updated');
  assert.equal(rows[0].querySelector('.u2-domain-history-actor').textContent, 'system');
  assert.equal(rows[0].querySelector('.u2-timestamp') !== null, true);
  assert.deepEqual(changes(rows[0]), ['Title: Aspirin → Aspirin 2', 'Weight: 1.5 → 3'],
    'one line per changed column; the version stamp is not a change');

  src.rows.byKey('i1').done = false;
  assert.equal(await src.save(), true);
  await settle();
  rows = lines(history);
  assert.equal(rows.length, 2);
  assert.deepEqual(changes(rows[0]), ['Done: true → false'], 'newest first');
  assert.equal(history.changesOf({op: 'insert', before: null, after: {title: 'x'}}), '');
  assert.equal(history.changesOf({op: 'update', before: {title: 'A', version: 1, priority: null},
    after: {title: 'B', version: 2, priority: 'high', nosuch: 1}}), 'title: A → B, Priority: (empty) → high');
  history.dispose();
  src.dispose();
});

scoped('the actor is the user\'s name, looked up once per id; a draft is not saved yet; no row asks for one', async () => {
  emsBackends.domain = backend();
  const memory = emsBackends.domain;
  const looked = [];
  grok.dapi.users.find = async (id) => {
    looked.push(id);
    return id === 'u1' ? {friendlyName: 'Ann'} : null;
  };
  const table = await domains.table('grit.issue');
  const draft = table.draft({project_id: 'p1', title: 'Draft'});
  await flush();
  const history = domains.history(draft);
  await settle();
  const hint = history.root.querySelector('[data-u2-part="hint"]');
  assert.equal(hint.textContent, 'Not saved yet');
  assert.equal(history.view.root.hidden, true);

  draft.currentRow.value.number = 7;
  assert.equal(await draft.save(), true);
  const id = draft.currentRow.value.id;
  memory.tableSync('grit.issue').history.push(
    {id, tx_id: 't2', op: 'update', actor_id: 'u1', ts: '2026-09-13T10:00:00Z', before: {title: 'Draft'}, after: {title: 'D2'}},
    {id, tx_id: 't3', op: 'update', actor_id: 'u2', ts: '2026-09-13T11:00:00Z', before: {title: 'D2'}, after: {title: 'D3'}},
    {id, tx_id: 't4', op: 'promote', actor_id: 'u1', ts: '2026-09-13T12:00:00Z', before: null, after: null});
  history.view.refresh();
  await settle();
  assert.equal(hint.hidden, true, 'the saved row has a history');
  const rows = lines(history);
  assert.deepEqual(rows.map((r) => r.querySelector('.u2-badge').textContent), ['shared', 'updated', 'updated', 'created']);
  assert.deepEqual(rows.map((r) => r.querySelector('.u2-domain-history-actor').textContent), ['Ann', 'u2', 'Ann', 'system'],
    'a known user by name, an unknown one by id');
  assert.deepEqual(looked.sort(), ['u1', 'u2'], 'one lookup per actor');

  draft.currentRow.value = null;
  await flush();
  assert.equal(hint.textContent, 'Select a issue to see its history.');
  assert.equal(history.view.root.hidden, true);
  history.dispose();
  draft.dispose();
});

scoped('the trail loads 50 entries at a time; "Show older" loads the older ones down to the first', async () => {
  emsBackends.domain = backend();
  const issue = emsBackends.domain.tableSync('grit.issue');
  for (let i = 0; i < 60; i++) {
    issue.history.push({id: 'i1', tx_id: `t${i}`, op: 'update', actor_id: null, ts: '2026-09-13T10:00:00Z',
      before: {title: `T${i}`}, after: {title: `T${i + 1}`}});
  }
  const pages = [];
  const audit = issue.audit.bind(issue);
  Object.defineProperty(issue, 'audit', {configurable: true, value: (id, page) => {
    pages.push(page);
    return audit(id, page);
  }});
  const src = (await domains.table('grit.issue')).source({pageSize: 10});
  await flush();
  src.currentRow.value = src.rows.byKey('i1');
  const history = domains.history(src);
  await settle();
  assert.deepEqual(pages, [{limit: 50}]);
  const more = history.view.root.querySelector('[data-u2-part="more"]');
  const shown = () => [lines(history).length, more.hidden];
  assert.deepEqual(shown(), [50, false]);

  fire(more, 'click');
  await settle();
  assert.deepEqual([pages[1], shown()], [{limit: 50, before: 10}, [60, true]]);
  history.dispose();
  src.dispose();
});

scoped('a table that keeps no history says so instead of showing an empty list', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support, audit: false};
  Object.defineProperty(issue, 'audit', {value: undefined, configurable: true});
  emsBackends.domain = memory;
  const table = await domains.table('grit.issue');
  const source = table.source();
  await flush();
  source.currentRow.value = source.rows.items.value[0];
  const history = domains.history(source);
  await settle();
  const hint = history.root.querySelector('[data-u2-part="hint"]');
  assert.equal(hint.textContent, 'This table keeps no history.');
  assert.equal(history.view.root.hidden, true, 'and no empty list nothing will ever fill');
  history.dispose();
  source.dispose();
});

scoped('spec: u2-domain-history over a bound source', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  assert.equal(reg.get('u2-domain-history').props[0].name, 'source', 'the source, then the shared appearance props');
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-domain-history', name: 'h',
    bind: {source: '$.issues'}}}, new SpecContext({data: {issues: signal(src)}}), reg);
  await flush();
  assert.equal(instance.node('h') instanceof DomainHistory, true);
  assert.equal(instance.node('h').source, src);
  instance.dispose();
  src.dispose();
});

scoped('a reference cell shows the name it points at, not the uuid; the pane carries a History header', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  await flush();
  src.currentRow.value = src.rows.byKey('i1');
  const history = domains.history(src);
  await settle();
  assert.equal(history.root.querySelector('.u2-section-title').textContent, 'History');

  src.rows.byKey('i1').project_id = 'p2';
  assert.equal(await src.save(), true);
  await settle();
  assert.deepEqual(changes(lines(history)[0]), ['Project id: Grit → Datagrok'],
    'a ref cell shows the name of the row it points at, not the uuid');
  history.dispose();
  src.dispose();
});

scoped('HI-1: "(empty)" for a missing value, dates as dates, derived columns by name, long entries fold, folded summary', async () => {
  emsBackends.domain = backend();
  const issue = emsBackends.domain.tableSync('grit.issue');
  issue.history.push({id: 'i1', tx_id: 't1', op: 'update', actor_id: null, ts: '2026-09-13T10:00:00Z',
    before: {title: 'A', description: null, due: null, weight: 1, done: false, number: 1, tags: null, priority: 'low'},
    after: {title: 'B', description: 'd', due: '2026-10-02T00:00:00Z', weight: 2, done: true, number: 2, tags: ['x'],
      priority: 'high'}});
  const src = (await domains.table('grit.issue')).source({pageSize: 10});
  await flush();
  src.currentRow.value = src.rows.byKey('i1');
  const history = domains.history(src, undefined, {include: ['title', 'description', 'due', 'weight', 'done', 'number']});
  await settle();
  assert.equal(history.section.expanded.value, false, 'folded by default');
  assert.equal(history.root.querySelector('.u2-section-hint').textContent.startsWith('updated '), true,
    'the folded header sums up the newest entry');
  const [line] = lines(history);
  const all = changes(line);
  assert.equal(all.length, 6);
  assert.equal(all.includes('Description: (empty) → d'), true, all.join('; '));
  const due = all.find((c) => c.startsWith('Due: '));
  assert.equal(/^Due: \(empty\) → [^:]*2026$/.test(due), true, `a date at midnight UTC is a date, no time: ${due}`);
  assert.equal(line.querySelector('.u2-domain-history-derived').textContent, 'Also changedTags, Priority',
    'the columns outside the form by name only');
  const folded = [...line.querySelectorAll('.u2-domain-history-change, .u2-domain-history-derived')].filter((l) => l.hidden);
  assert.equal(folded.length, 2, 'past five lines the rest fold');
  const expand = line.querySelector('.u2-domain-history-expand');
  assert.equal(expand.textContent, '2 more');
  fire(expand, 'click');
  assert.equal(folded.every((l) => !l.hidden) && !expand.isConnected, true);
  history.dispose();
  src.dispose();
});

scoped('VR-39: a failed load says "History is unavailable" with the server\'s words as its tooltip', async () => {
  emsBackends.domain = backend();
  const issue = emsBackends.domain.tableSync('grit.issue');
  Object.defineProperty(issue, 'audit', {configurable: true, value: async () => {
    throw new Error('NoSuchMethodError: method not found: \'call\'');
  }});
  const src = (await domains.table('grit.issue')).source({pageSize: 10});
  await flush();
  src.currentRow.value = src.rows.byKey('i1');
  const history = domains.history(src);
  await settle();
  assert.equal(history.view.root.querySelector('.u2-async-error-message').textContent, 'History is unavailable');
  assert.equal(history.view.root.title, 'NoSuchMethodError: method not found: \'call\'');
  history.dispose();
  src.dispose();
});
