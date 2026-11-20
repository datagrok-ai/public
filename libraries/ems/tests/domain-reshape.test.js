import {emsBackends} from '../src/core/domain-backend.js';
import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {SCHEMA, backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');
const {DomainViews} = await import('../src/dg/views.js');
const grok = await import('datagrok-api/grok');

const BASE = '/apps/T/Issues';

/** The schema as Edit properties… left it: the issue table's columns with `change` applied. */
function reshaped(change) {
  const issue = SCHEMA.tables.issue;
  const schema = {...SCHEMA, tables: {...SCHEMA.tables, issue: {...issue, columns: change({...issue.columns})}}};
  const rows = backend().tableSync('grit.issue').rows;
  return new MemoryDomainBackend(schema, {rows: {project: backend().tableSync('grit.project').rows, issue: rows}});
}

/** A docked app over grit.issue, the views the reload docks, and the cleanup. */
async function docked() {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const view = table.app({path: BASE});
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

scoped('reshaped: a page with unsaved changes is marked when Edit properties… adds a column, and rebuilt with it after its Save', async () => {
  const {app, view, added, done} = await docked();
  try {
    app.listSource.rows.byKey('i1').title = 'Edited';
    await flush();
    emsBackends.domain = reshaped((columns) => ({...columns, severity: {type: 'string'}}));
    await DomainViews.reload('grit');
    await flush();
    assert.equal(app.status.value.notice, '1 unsaved change — Properties changed — save or discard, then the page reloads',
      'a dirty page says its columns changed');
    assert.equal(app.stale.value, true, 'and is stale: Refresh is offered');
    assert.equal(added.length, 0, 'not rebuilt under the pending change');

    assert.equal(await app.session.save(), true);
    await flush();
    assert.equal(added.length, 1, 'rebuilt once the change is saved');
    assert.equal(view.dart.closed, 1, 'in place of the old view');
    const fresh = DomainApp.of(added[0]);
    assert.ok(fresh.table.properties.some((p) => p.name === 'severity'), 'the new view has the new column');
    assert.equal(emsBackends.domain.tableSync('grit.issue').rows.find((r) => r.id === 'i1').title, 'Edited');
  } finally {
    done();
  }
});

scoped('reshaped: a pending change to a column Edit properties… dropped is never sent; Discard rebuilds the page', async () => {
  const {app, view, added, done} = await docked();
  try {
    app.listSource.rows.byKey('i1').priority = 'low';
    await flush();
    const next = reshaped(({priority, ...columns}) => columns);
    let sent = 0;
    const real = next.saveAll.bind(next);
    next.saveAll = (edits) => (sent++, real(edits));
    emsBackends.domain = next;
    await DomainViews.reload('grit');
    await flush();
    assert.equal(app.status.value.notice, '1 unsaved change — Properties changed — save or discard, then the page reloads',
      'a dirty page says its columns changed');
    assert.equal(app.session.validity.value, 'Priority was removed in Edit properties… — discard its change',
      'Save is off and says why');
    assert.equal(await app.session.save(), false);
    assert.equal(sent, 0, 'the change to the dropped column is not sent');
    assert.equal(added.length, 0);

    app.session.discard();
    await flush();
    assert.equal(added.length, 1, 'rebuilt once the change is discarded');
    assert.equal(view.dart.closed, 1);
    assert.equal(DomainApp.of(added[0]).table.properties.some((p) => p.name === 'priority'), false,
      'the new view has no dropped column');
  } finally {
    done();
  }
});
