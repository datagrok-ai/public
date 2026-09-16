/* `domains.import` (WO 3-9) on the memory backend: the wizard's four steps, the mapping
   auto-matched by name against the columns the caller may write, every blocking problem the Dart
   wizard defines, a skipped column absent from what is posted, the preview's own verdicts, upsert
   merging by the business key, and the report reaching the caller. The memory backend has no bulk
   endpoint, so the posting path here is the transaction fallback — the same rows either way. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {backends} from '../src/sources/backends.js';
import {notify} from '../src/components/display/notify.js';
import {SCHEMA, backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const DG = await import('datagrok-api/dg');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')]
  .find((b) => b.textContent === text);
const named = (name) => Control.forElement(document.querySelector(`[data-u2-name="${name}"]`));
const reason = () => document.body.querySelector('.u2-wizard-reason')?.textContent ?? '';
const mapping = () => document.body.querySelector('.u2-domain-import-mapping');
const issues = (memory) => memory.tableSync('grit.issue').rows;

/** A source frame of the given columns: every value a string, as a CSV import delivers them. */
function frame(columns, rows, name = 'import.csv') {
  return new DG.DataFrame(columns.map((c) => ({name: c, type: 'string'})), rows, name);
}

/** Opens the wizard over `grit.issue` and returns the promise plus the mapping table, ready to
 * drive: the first step is shown, the frame is picked, and the mapping step is on screen. */
async function opened(memory, source) {
  backends.domain = memory;
  const handle = await domains.table('grit.issue');
  const running = domains.import(handle);
  await flush();
  named('source').value.value = source;
  await flush();
  // the mapping table has no layout in the shim: give it a viewport so its rows render
  return {running, handle};
}

async function toMapping() {
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  mapping().clientHeight = 400;
  fire(mapping(), 'scroll');
  await flush();
}

scoped('the mapping is auto-matched by name and caption; a system column is never a target', async () => {
  const memory = backend();
  const {running} = await opened(memory, frame(['title', 'Priority', 'nonesuch'], [{title: 'A', Priority: 'low'}]));
  await toMapping();
  assert.equal(named('title').value.value, 'title');
  assert.equal(named('Priority').value.value, 'priority', 'matched on the caption, case-insensitively');
  assert.equal(named('nonesuch').value.value, '(skip)');
  const labels = named('title').items.map((x) => typeof x === 'string' ? x : x.label);
  assert.equal(labels.some((l) => l === 'Id' || l === 'Version'), false, 'the system columns are not offered');
  assert.equal(labels.includes('Title *'), true, 'a target the table will not take empty carries the asterisk');
  assert.equal(labels.includes('Description'), true, 'and one it will takes none');
  assert.equal(named('title').root.classList.contains('u2-input-required'), true,
    'the choice standing on a required target wears the form convention');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('a frame handed in as an option is mapped while the wizard is still being built', async () => {
  const memory = backend();
  backends.domain = memory;
  const handle = await domains.table('grit.issue');
  const running = domains.import(handle,
    {source: frame(['project_id', 'title'], [{project_id: 'p1', title: 'Handed in'}])});
  await flush();
  await toMapping();
  assert.equal(named('title').value.value, 'title', 'the mapping is there without a pick');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('every blocking problem the wizard gates on', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'name2'], [{project_id: 'p1', title: 'A', name2: 'B'}]));
  await toMapping();
  const proceed = async () => {
    fire(buttonNamed('NEXT'), 'click');
    await flush();
    return reason();
  };
  named('name2').value.value = 'title';
  await flush();
  assert.match(await proceed(), /Two source columns are mapped to "title"/);

  for (const column of ['project_id', 'title', 'name2'])
    named(column).value.value = '(skip)';
  await flush();
  assert.match(await proceed(), /Map at least one column/);

  named('title').value.value = 'description';
  await flush();
  assert.match(await proceed(), /Required column "Project id" is not mapped/,
    'insert mode names every required column left unmapped');

  named('project_id').value.value = 'project_id';
  await flush();
  assert.match(await proceed(), /Required column "Title" is not mapped/);

  named('title').value.value = 'title';
  named('mode').value.value = 'upsert';
  await flush();
  assert.match(await proceed(), /Upsert merges by the business key — map a column to "Number"/);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('the preview reports the schema\'s own verdicts over the mapped cells', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'number', 'priority'],
      [{project_id: 'p1', title: 'A', number: '3', priority: 'low'},
        {project_id: 'p1', title: 'B', number: '3.0', priority: 'urgent'}]));
  await toMapping();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  const preview = document.body.querySelector('.u2-domain-import-preview');
  assert.match(preview.textContent, /1 of 2 checked rows have problems/);
  assert.match(preview.textContent, /Integer value expected, passed: "3\.0"/);
  assert.match(preview.textContent, /Must be one of: low, high/);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('only the mapped columns are posted, under their target names', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'note'], [{project_id: 'p1', title: 'Imported', note: 'skipped'}]));
  await toMapping();
  assert.equal(named('note').value.value, '(skip)');
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  const landed = issues(memory).find((r) => r.title === 'Imported');
  assert.notEqual(landed, undefined, 'the row was inserted');
  assert.equal('note' in landed, false, 'a skipped source column never reached the backend');
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent,
    /1 inserted, 0 updated, 0 skipped, 0 failed/);
  fire(buttonNamed('CLOSE'), 'click');
  const result = await running;
  assert.equal(result.inserted, 1);
  assert.deepEqual(result.rows.map((r) => r.status), ['inserted']);
});

scoped('upsert merges by the business key instead of inserting again', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'number', 'title'], [{project_id: 'p1', number: 2, title: 'Renamed'}]));
  await toMapping();
  named('mode').value.value = 'upsert';
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  assert.equal(issues(memory).length, 3, 'nothing was inserted');
  assert.equal(issues(memory).find((r) => r.id === 'i2').title, 'Renamed');
  fire(buttonNamed('CLOSE'), 'click');
  const result = await running;
  assert.equal(result.updated, 1);
  assert.equal(result.inserted, 0);
});

scoped('a refused row leaves nothing behind and the report says so', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Fine'}, {project_id: 'p1', title: ''}]));
  await toMapping();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  assert.equal(issues(memory).length, 3, 'the transaction rolled back');
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent, /column "title"/);
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal(await running, null, 'a failed import reports nothing');
});

scoped('a table-level insert denial still opens the wizard; no visible column refuses it', async () => {
  // the row-mode deviation: `can.insert` is a false negative there, so column security decides
  backends.domain = backend({access: {can: {view: true, insert: false, edit: false, delete: false, share: false},
    fields: {project_id: 'editable', title: 'editable'}}});
  const open = await domains.table('grit.issue');
  const running = domains.import(open,
    {source: frame(['project_id', 'title'], [{project_id: 'p1', title: 'A'}])});
  await flush();
  await toMapping();
  assert.equal(named('title').value.value, 'title', 'the writable columns are offered');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);

  backends.domain = backend({access: {can: {view: true, insert: true, edit: true, delete: true, share: true},
    fields: {}}});
  const hidden = await domains.table('grit.issue');
  assert.equal(await domains.import(hidden), null);
  assert.match(document.body.querySelector('.u2-notify-warning')?.textContent ?? '',
    /no column of grit.issue you may write/);
});

scoped('the preview draws the rows as they would land, refused cells marked', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'number'], [{project_id: 'p1', title: 'A', number: '3'},
      {project_id: 'p1', title: 'B', number: '3.0'}]));
  await toMapping();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  const rows = document.body.querySelector('.u2-domain-import-rows');
  rows.clientHeight = 200;
  fire(rows, 'scroll');
  await flush();
  assert.deepEqual([...rows.querySelectorAll('.u2-data-table-head')].map((c) => c.textContent),
    ['Project id', 'Title', 'Number'], 'the target captions, not the DB names');
  const cells = [...rows.querySelectorAll('.u2-data-table-row')]
    .map((r) => [...r.children].map((c) => c.textContent).join('|'));
  assert.deepEqual(cells, ['p1|A|3', 'p1|B|3.0'], 'the mapped columns as they would be sent');
  const bad = rows.querySelector('.u2-cell-error');
  assert.notEqual(bad, null, 'the cell the schema refuses is marked');
  assert.match(bad.title, /Integer value expected/);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('the report step is terminal: no BACK, no CANCEL, and the button reads CLOSE', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Terminal'}]));
  await toMapping();
  assert.notEqual(buttonNamed('BACK'), undefined, 'BACK is there while the work can still be undone');
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  assert.equal(buttonNamed('FINISH'), undefined, 'the last step does not offer to finish again');
  assert.notEqual(buttonNamed('CLOSE'), undefined);
  assert.equal(buttonNamed('BACK').style.display, 'none', 'nothing to go back to: the rows are written');
  assert.equal(buttonNamed('CANCEL').style.display, 'none', 'and nothing left to cancel');
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).inserted, 1);
});

scoped('a source column is matched the way a person reads it, not byte for byte', async () => {
  const memory = backend();
  // 'Chemical title' → title, 'Project ID' → project_id, 'Weight (kg)' → weight, and a column
  // naming nothing is skipped — case, spaces and punctuation are not identity
  const {running} = await opened(memory, frame(['Project ID', 'Chemical title', 'Weight (kg)', 'Batch'],
    [{'Project ID': 'p1', 'Chemical title': 'A', 'Weight (kg)': '2', 'Batch': 'x'}]));
  await toMapping();
  assert.equal(named('Project ID').value.value, 'project_id');
  assert.equal(named('Chemical title').value.value, 'title');
  assert.equal(named('Weight (kg)').value.value, 'weight');
  assert.equal(named('Batch').value.value, '(skip)');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});
