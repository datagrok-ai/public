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

scoped('a pending change to a column Edit properties… removed is never sent; Discard rebuilds the page', async () => {
  const {app, added, done} = await docked('grit.project', {path: '/apps/T/Projects', children: {mode: 'list'}});
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
    assert.match(app.session.validity.value, /^Priority was removed/);
    assert.equal(await app.session.save(), false);
    assert.equal(sent, 0);
    app.session.discard();
    await flush();
    assert.equal(added.length, 1);
  } finally {
    done();
  }
});

scoped('Refresh over changed properties asks to discard the pending changes before it rebuilds the page', async () => {
  const {app, added, done} = await docked();
  try {
    app.listSource.rows.byKey('i1').title = 'Edited';
    await flush();
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    await DomainViews.reload('grit');
    await flush();
    const reload = app.ribbon().main.find((c) => c.root.dataset.u2 === 'refresh-button');
    fire(reload.root.querySelector('button') ?? reload.root, 'click');
    await flush();
    assert.equal(added.length, 0, 'nothing rebuilt while the gate asks');
    fire(buttonNamed('DISCARD'), 'click');
    await flush();
    assert.equal(added.length, 1);
    assert.ok(DomainApp.of(added[0]).table.properties.some((p) => p.name === 'severity'));
  } finally {
    done();
  }
});
