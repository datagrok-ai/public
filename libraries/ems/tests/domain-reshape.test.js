import {emsBackends} from '../src/core/domain-backend.js';
import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {SCHEMA, backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');
const {DomainViews} = await import('../src/dg/views.js');
const grok = await import('datagrok-api/grok');

const BASE = '/apps/T/Issues';
const NOTICE = '1 unsaved change — Properties changed: save or discard to reload';
const DROPPED = 'Priority was removed in Edit properties…: undo your change to it, or Discard all changes';
const buttonNamed = (text) => [...document.body.querySelectorAll('button')].find((b) => b.textContent === text);

/** The schema as Edit properties… left it: the issue table's columns with `change` applied. */
function reshaped(change) {
  const issue = SCHEMA.tables.issue;
  const schema = {...SCHEMA, tables: {...SCHEMA.tables, issue: {...issue, columns: change({...issue.columns})}}};
  const rows = backend().tableSync('grit.issue').rows;
  return new MemoryDomainBackend(schema, {rows: {project: backend().tableSync('grit.project').rows, issue: rows}});
}

/** A docked app over grit.issue, the views the reload docks, and the cleanup. */
async function docked(address = 'grit.issue', options = {path: BASE}) {
  emsBackends.domain = backend();
  const table = await domains.table(address);
  const view = table.app(options);
  document.body.append(view.root);
  await flush();
  const added = [];
  const saved = grok.shell.addView;
  grok.shell.addView = (v) => (added.push(v), v);
  const done = () => {
    grok.shell.addView = saved;
    for (const v of [view, ...added])
      DomainApp.of(v)?.dispose();
  };
  return {app: DomainApp.of(view), view, added, done};
}

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

scoped('reshaped: a pending change in a children tab to a column its table lost is never sent; Discard rebuilds', async () => {
  const {app, view, added, done} = await docked('grit.project', {path: '/apps/T/Projects', children: {mode: 'list'}});
  try {
    assert.equal(await app.goTo('entity', 'p1'), true);
    await flush();
    app.session.sources.value.find((s) => s.table === 'grit.issue').rows.byKey('i1').priority = 'low';
    await flush();
    const next = reshaped(({priority, ...columns}) => columns);
    let sent = 0;
    const real = next.saveAll.bind(next);
    next.saveAll = (edits) => (sent++, real(edits));
    emsBackends.domain = next;
    await DomainViews.reload('grit');
    await flush();
    assert.equal(app.status.value.notice, NOTICE, 'the page says its columns changed');
    assert.equal(app.session.validity.value, DROPPED, 'Save is off and says why');
    assert.equal(await app.session.save(), false);
    assert.equal(sent, 0, 'the change to the dropped child column is not sent');
    app.session.discard();
    await flush();
    assert.equal(added.length, 1, 'rebuilt once the change is discarded');
    assert.equal(view.dart.closed, 1);
  } finally {
    done();
  }
});

scoped('reshaped: Refresh asks to discard the pending changes, then rebuilds the page with the new columns', async () => {
  const {app, added, done} = await docked();
  try {
    app.listSource.rows.byKey('i1').title = 'Edited';
    await flush();
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    await DomainViews.reload('grit');
    await flush();
    const reload = app.ribbon().main.find((c) => c.root.dataset.u2 === 'refresh-button');
    assert.equal(reload.root.hidden, false, 'Refresh is offered: the page is stale');
    fire(reload.root.querySelector('button') ?? reload.root, 'click');
    await flush();
    assert.equal(added.length, 0, 'nothing rebuilt while the gate asks');
    fire(buttonNamed('DISCARD'), 'click');
    await flush();
    assert.equal(added.length, 1, 'rebuilt once, after the discard');
    assert.ok(DomainApp.of(added[0]).table.properties.some((p) => p.name === 'severity'));
  } finally {
    done();
  }
});

scoped('reshaped: an edit undone by hand leaves the session clean, and the page is rebuilt then', async () => {
  const {app, added, done} = await docked();
  try {
    app.listSource.rows.byKey('i1').title = 'Edited';
    await flush();
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    await DomainViews.reload('grit');
    await flush();
    assert.equal(app.status.value.notice, NOTICE);
    app.listSource.rows.byKey('i1').title = 'Aspirin';
    await flush();
    assert.equal(app.session.isDirty.value, false, 'nothing left to save or discard');
    assert.equal(added.length, 1, 'rebuilt once the session is clean');
  } finally {
    done();
  }
});

scoped('reshaped: clearing a draft value in a removed column turns Save back on', async () => {
  const {app, done} = await docked();
  try {
    const draft = app.listSource.newRow({project_id: 'p1', number: 9, title: 'New', priority: 'low'});
    await flush();
    emsBackends.domain = reshaped(({priority, ...columns}) => columns);
    await DomainViews.reload('grit');
    await flush();
    assert.equal(app.session.validity.value, DROPPED);
    draft.priority = null;
    await flush();
    assert.equal(app.session.changeCount.value, 1, 'still one change: the draft');
    assert.equal(app.session.validity.value, null, 'nothing removed is pending any more');
  } finally {
    done();
  }
});

scoped('reshaped: a rebuild that fails says why, and Refresh tries again', async () => {
  const {app, done} = await docked();
  try {
    app.listSource.rows.byKey('i1').title = 'Edited';
    await flush();
    let tries = 0;
    app.reshape(new Map(), () => (tries++, Promise.reject(new Error('The table did not load'))));
    app.session.discard();
    await flush();
    assert.match(document.body.querySelector('[role="alert"]')?.textContent ?? '', /The table did not load/);
    assert.equal(app.status.value.notice, 'Data changed — Refresh', 'nothing left to save or discard');
    assert.equal(await app.refresh(), true);
    assert.equal(tries, 2, 'the page is still marked: Refresh rebuilds again');
  } finally {
    done();
  }
});

scoped('Refresh rebuilds the page when the table\'s columns changed meanwhile: a newly shared column shows', async () => {
  const {added, app, view, done} = await docked();
  try {
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    assert.equal(await app.refresh(), true);
    await flush();
    assert.equal(view.dart.closed, 1);
    assert.ok(DomainApp.of(added[0]).table.properties.some((p) => p.name === 'severity'));
  } finally {
    done();
  }
});
