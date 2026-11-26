import {emsBackends} from '../src/core/domain-backend.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {SCHEMA, backend, scoped, until} from './domain-fixtures.mjs';

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
    assert.equal(DomainApp.live.size, 0, 'no app left behind');
  };
  return {app: DomainApp.of(view), view, added, done};
}

scoped('Edit properties… rebuilds a page with its own options; one with pending changes once, when they are saved', async () => {
  const {app, view, added, done} = await docked('grit.issue', {path: BASE, name: 'Mine', query: 'done = false'});
  try {
    assert.equal(await app.goTo('entity', 'i2'), true);
    await flush();
    app.entitySource.value.currentRow.value.title = 'Edited';
    await flush();
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    await DomainViews.reload('grit');
    await flush();
    assert.deepEqual([added.length, app.stale.value, app.session.changeCount.value], [0, true, 1],
      'the edit stays, and the page says it is behind');
    const refreshing = app.refresh();
    await flush();
    fire(buttonNamed('SAVE'), 'click');
    assert.equal(await refreshing, true);
    await until(() => added.length === 1);
    await flush();
    assert.equal(added.length, 1, 'Refresh and the save rebuild it once');
    assert.equal(view.dart.closed, 1);
    const rebuilt = DomainApp.of(added[0]);
    assert.ok(rebuilt.table.properties.some((p) => p.name === 'severity'));
    assert.deepEqual([added[0].name, rebuilt.listSource.query.value, rebuilt.entity.value], ['Mine', 'done = false', 'i2'],
      'the page keeps its own options and its row');
  } finally {
    done();
  }
});

scoped('a second Edit properties landing while the rebuild reads the table rebuilds again', async () => {
  const {added, done} = await docked();
  try {
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    const first = DomainViews.reload('grit');
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}, impact: {type: 'string'}}));
    await Promise.all([first, DomainViews.reload('grit')]);
    assert.equal(added.length, 1);
    assert.equal(DomainApp.of(added[0]).table.properties.some((p) => p.name === 'impact'), true);
  } finally {
    done();
  }
});
