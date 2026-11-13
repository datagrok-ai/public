import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Rows, notify} from '@datagrok-libraries/u2';
import {SharedSession} from '../src/core/session.js';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');

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

const BASE = '/apps/T/Issues';

/** Records every row of the frame at every event it fires, so the order a draft becomes
 * observable in can be asserted: stamped, then ONE notification. */
function watched(domain, seen) {
  const table = domain.table.bind(domain);
  domain.table = async (address) => {
    const handle = await table(address);
    const frame = handle.frame.bind(handle);
    handle.frame = async (spec) => {
      const made = await frame(spec);
      const observe = () => {
        for (let i = 0; i < made.df.rowCount; i++)
          seen.push(`${made.df.get('id', i) ?? ''}:${made.df.get(Rows.STATE, i) ?? ''}`);
      };
      const subs = [made.df.onRowsAdded.subscribe(observe), made.df.onValuesChanged.subscribe(observe),
        made.edit.onChanged.subscribe(observe)];
      return {...made, dispose: () => {
        for (const s of subs)
          s.unsubscribe();
        made.dispose();
      }};
    };
    return handle;
  };
  return domain;
}

async function app(options = {}) {
  const table = await domains.table('grit.issue');
  const a = testOwner().build(() => domains.app({table, base: BASE, pageSize: 10, ...options}));
  document.body.append(a.root);
  await flush();
  return a;
}

scoped('a draft is observable only stamped: its id and its state are there at the first event', async () => {
  const seen = [];
  emsBackends.domain = watched(backend(), seen);
  const a = await app();
  assert.equal(await a.create(), true);
  await flush();
  assert.equal(seen.some((row) => row.startsWith(Rows.DRAFT_PREFIX)), true, 'the draft was observed');
  assert.deepEqual(seen.filter((row) => row.startsWith(Rows.DRAFT_PREFIX) && !row.endsWith(':new')), [],
    'never half-stamped: no event carries a draft that is not new');
  assert.deepEqual(seen.filter((row) => row.startsWith(':')), [], 'and none carries an unkeyed row');
  assert.equal(a.entity.value, DomainApp.NEW, 'the transient unkeyed row never became the entity');
  assert.equal(a.path.value, `${BASE}?entity=new`);
  const source = a.entitySource.value;
  const row = source.currentRow.value;
  assert.equal(Rows.isDraft(row), true);
  assert.equal(row[Rows.STATE], 'new');
  const view = source.access.value.row(row);
  assert.equal(view.isDraft, true);
  assert.equal(view.field('title'), 'editable');
  assert.notEqual(a.form.value.input('title'), undefined, 'the create form is editable, not text');
  a.dispose();
});

scoped('a cold ?entity=new opens the create page', async () => {
  emsBackends.domain = backend();
  const a = await app();
  location.search = '?entity=new';
  assert.equal(await a.open(), true);
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, DomainApp.NEW);
  assert.equal(a.path.value, `${BASE}?entity=new`);
  a.dispose();
});

scoped('the table validators refuse a save no form ever saw', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const off = table.validators.add('priority',
    (value, row) => value === 'high' && !row.reporter ? 'Assign before escalating' : null);
  const source = testOwner().build(() => table.source({}));
  await flush();
  const row = source.rows.byKey('i2');
  row.priority = 'high';
  await flush();
  assert.equal(source.isDirty.value, true);
  assert.equal(await source.save(), false);
  assert.equal(source.error.value.message, 'Cannot save: Ibuprofen: Priority: Assign before escalating');
  assert.equal(source.summary.value, 'Cannot save: Ibuprofen: Priority: Assign before escalating',
    'the refusal is what the status bar says, not the change count');
  assert.equal(source.problemRow.value, 'i2', 'and it names the row');
  const balloons = document.body.querySelectorAll('.u2-notify-error');
  assert.equal(balloons.length, 1, 'one balloon for a refusal no form ever saw');
  assert.match(balloons[0].textContent, /Assign before escalating/);
  assert.equal(await source.save(), false, 'refused again');
  assert.equal(document.body.querySelectorAll('.u2-notify-error').length, 1, 'the new refusal replaces the balloon');
  row.reporter = 'u1';
  await flush();
  assert.equal(source.error.value, undefined, 'the next edit takes the refusal back');
  assert.equal(document.body.querySelector('.u2-notify-error'), null, 'and its balloon');
  assert.match(source.summary.value, /unsaved change/);
  assert.equal(await source.save(), true);
  assert.equal(emsBackends.domain.tableSync('grit.issue').rows.find((r) => r.id === 'i2').priority, 'high');
  off();
  source.dispose();
});

scoped('the writer\'s own cell refusal names the row and the field, and marks the row', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const source = testOwner().build(() => table.source({}));
  await flush();
  const draft = source.newRow({project_id: 'p1', number: 7});
  await flush();
  assert.equal(source.validity.value, 'Value can\'t be empty', 'the writer says only what is wrong');
  assert.equal(await source.save(), false);
  assert.equal(source.error.value.message, 'Cannot save: New issue: Title: Value can\'t be empty',
    'the refusal names the row and the field the writer refused over');
  assert.equal(source.summary.value, 'Cannot save: New issue: Title: Value can\'t be empty');
  assert.equal(source.problemRow.value, draft.id, 'and marks the row');
  draft.title = 'Titled';
  await flush();
  assert.equal(await source.save(), true);
  source.dispose();
});

scoped('the validators run over a row the frame filter hides — the writer sends it all the same', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const off = table.validators.add('priority', (value) => value === 'high' ? 'Escalation needs a reviewer' : null);
  const source = testOwner().build(() => table.source({}));
  await flush();
  source.rows.byKey('i2').priority = 'high';
  const at = source.rows.items.value.findIndex((r) => r.id === 'i2');
  const df = source.df.value;
  df.filter = {get: (i) => i !== at};
  df.onFilterChanged.fire(undefined);
  await flush();
  assert.equal(source.rows.byKey('i2'), undefined, 'the edited row is out of sight');
  assert.equal(source.isDirty.value, true, 'and still pending');
  assert.equal(await source.save(), false, 'a hidden edit is checked like any other');
  assert.match(source.error.value.message, /Escalation needs a reviewer/);
  assert.equal(source.problemRow.value, 'i2');
  off();
  source.dispose();
});

scoped('a pristine parent the batch inserts is checked too: its validator refuses the save', async () => {
  emsBackends.domain = backend();
  const session = new SharedSession();
  const projects = await domains.table('grit.project');
  const issues = await domains.table('grit.issue');
  const off = projects.validators.add('name', (value) => value === 'Pristine' ? 'Pick another name' : null);
  const parentSource = testOwner().build(() => projects.source({session}));
  const childSource = testOwner().build(() => issues.source({session}));
  await flush();
  const parent = parentSource.newRow({key: 'NEW', name: 'Pristine'}, {pristine: true});
  childSource.newRow({project_id: parent.id, title: 'Child'});
  assert.equal(parentSource.isDirty.value, false, 'the parent is in the batch without being dirty');
  assert.equal(await session.save(), false, 'and its guards decide the batch all the same');
  assert.match(parentSource.error.value.message, /Pick another name/);
  assert.equal(emsBackends.domain.tableSync('grit.project').rows.length, 2, 'nothing landed');
  off();
  parentSource.dispose();
  childSource.dispose();
});
