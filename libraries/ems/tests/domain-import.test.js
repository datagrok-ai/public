/* `domains.import` (WO 3-9) on the memory backend: the wizard's four steps, the mapping
   auto-matched by name against the columns the caller may write, every blocking problem the Dart
   wizard defines, a skipped column absent from what is posted, the preview's own verdicts, upsert
   merging by the business key, and the report reaching the caller. The memory backend has no bulk
   endpoint, so the posting path here is the transaction fallback — the same rows either way. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Control, notify} from '@datagrok-libraries/u2';
import {SCHEMA, backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const DG = await import('datagrok-api/dg');

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

const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')]
  .find((b) => b.textContent === text);
/** The preview's committing button: it reads "Import", with the count where the dry run saw every row. */
const importButton = () => [...document.body.querySelectorAll('.u2-dialog button')]
  .find((b) => /^Import\b/.test(b.textContent));
const named = (name) => Control.forElement(document.querySelector(`[data-u2-name="${name}"]`));
const reason = () => document.body.querySelector('.u2-wizard-reason')?.textContent ?? '';
const mapping = () => document.body.querySelector('.u2-domain-import-mapping');
const issues = (memory) => memory.tableSync('grit.issue').rows;

/** A source frame of the given columns: every value a string, as a CSV import delivers them. */
function frame(columns, rows, name = 'import.csv') {
  return new DG.DataFrame(columns.map((c) => ({name: c, type: 'string'})), rows, name);
}

/** The storage answers no dry run: the preview shows the mapped rows and the commit is the only check. */
function noDryRun(memory) {
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support, batch: {...issue.support.batch, validate: false}};
  issue.validate = undefined;
  return issue;
}

/** Opens the wizard over `grit.issue` and returns the promise plus the mapping table, ready to
 * drive: the first step is shown, the frame is picked, and the mapping step is on screen. */
async function opened(memory, source) {
  emsBackends.domain = memory;
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

/** Into the preview step — the dry run is a round trip, so the render is a flush behind it. */
async function toPreview() {
  fire(buttonNamed('NEXT'), 'click');
  await flush();
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

scoped('an immutable key column is a target, and auto-matched by name: an import inserts', async () => {
  // mapping only: the memory backend declares no `immutable` field, so the insert itself is not proven here
  const memory = backend({access: {can: {view: true, insert: true, edit: true, delete: true, share: true},
    fields: {project_id: 'editable', number: 'immutable', title: 'editable', description: 'editable'}}});
  const {running} = await opened(memory, frame(['number', 'title'], [{number: '7', title: 'Keyed'}]));
  await toMapping();
  assert.equal(named('number').value.value, 'number', 'the key typed once on insert is matched by name');
  const labels = named('title').items.map((x) => typeof x === 'string' ? x : x.label);
  assert.equal(labels.includes('Number'), true, 'the key is offered as a target');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('a frame handed in as an option is mapped while the wizard is still being built', async () => {
  const memory = backend();
  emsBackends.domain = memory;
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

scoped('the preview is the backend\'s dry run: its verdicts, its counts, its issues', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'number', 'priority'],
      [{project_id: 'p1', title: 'A', number: '3', priority: 'low'},
        {project_id: 'p1', title: 'B', number: '3.0', priority: 'urgent'}]));
  await toMapping();
  await toPreview();
  const preview = document.body.querySelector('.u2-domain-import-preview');
  assert.match(preview.textContent, /1 will be added, 0 updated, 0 skipped, 1 row has errors\./);
  assert.match(reason(), /1 row has errors\. Uncheck "All or nothing" on Source to import the 1 valid row/,
    'the counts are per row; all-or-nothing makes one bad row the verdict on all of them');
  assert.match(preview.textContent, /Must be one of: low, high/,
    'the issue lines are the report\'s, not a second implementation of the rules');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('the dry run is posted once per entry into the step, not per keystroke', async () => {
  const memory = backend();
  const posts = [];
  const table = memory.tableSync('grit.issue');
  const real = table.validate.bind(table);
  table.validate = (rows, options) => {
    posts.push({rows, options});
    return real(rows, options);
  };
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'note'], [{project_id: 'p1', title: 'A', note: 'x'}]));
  await toMapping();
  await toPreview();
  assert.equal(posts.length, 1, 'one post on the way in');
  assert.deepEqual(posts[0].rows, [{project_id: 'p1', title: 'A'}], 'the mapped columns only');
  assert.equal(posts[0].options.mode, 'insert');

  fire(buttonNamed('BACK'), 'click');
  await flush();
  named('note').value.value = 'description';
  await flush();
  named('note').value.value = '(skip)';
  await flush();
  assert.equal(posts.length, 1, 'typing on another step posts nothing');

  await toPreview();
  assert.equal(posts.length, 1, 'the same mapping is the same dry run');
  fire(buttonNamed('BACK'), 'click');
  await flush();
  named('note').value.value = 'description';
  await flush();
  await toPreview();
  assert.equal(posts.length, 2, 'a mapping the run did not cover is posted again');
  assert.deepEqual(posts[1].rows, [{project_id: 'p1', title: 'A', description: 'x'}]);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('a dry run that cannot reach the backend shows why and never blocks the import', async () => {
  const memory = backend();
  memory.tableSync('grit.issue').validate = () => Promise.reject(new Error('the stand is down'));
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Regardless'}]));
  await toMapping();
  await toPreview();
  const preview = document.body.querySelector('.u2-domain-import-preview');
  assert.match(preview.textContent, /the stand is down/);
  assert.equal(reason(), '', 'the commit is the authority: NEXT is not gated on the preview');
  fire(importButton(), 'click');
  await flush();
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent,
    /1 inserted, 0 updated, 0 skipped, 0 failed/);
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).inserted, 1);
});

scoped('a storage with no dry run previews the mapped rows without verdicts, and the import still posts', async () => {
  const memory = backend();
  noDryRun(memory);
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'note'], [{project_id: 'p1', title: 'Unchecked', note: 'x'}]));
  await toMapping();
  await toPreview();
  const preview = document.body.querySelector('.u2-domain-import-preview');
  assert.match(preview.textContent, /Rows are checked when imported\./);
  const rows = preview.querySelector('.u2-domain-import-rows');
  rows.clientHeight = 200;
  fire(rows, 'scroll');
  await flush();
  assert.deepEqual([...rows.querySelectorAll('.u2-data-table-head')].map((c) => c.textContent),
    ['Project id', 'Title'], 'the mapped columns under their captions, no verdict column');
  assert.deepEqual([...rows.querySelectorAll('.u2-data-table-row')]
    .map((r) => [...r.children].map((c) => c.textContent).join('|')), ['p1|Unchecked']);
  assert.equal(reason(), '', 'nothing gates the import');
  assert.equal(importButton().textContent, 'Import', 'nothing counted the rows that would land');
  fire(importButton(), 'click');
  await flush();
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent,
    /1 inserted, 0 updated, 0 skipped, 0 failed/);
  fire(buttonNamed('CLOSE'), 'click');
  const report = await running;
  assert.equal(report.inserted, 1);
  assert.notEqual(issues(memory).find((r) => r.title === 'Unchecked'), undefined, 'the row landed through batch');
});

scoped('only the mapped columns are posted, under their target names', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'note'], [{project_id: 'p1', title: 'Imported', note: 'skipped'}]));
  await toMapping();
  assert.equal(named('note').value.value, '(skip)');
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(importButton(), 'click');
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
  fire(importButton(), 'click');
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
  noDryRun(memory);
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Fine'}, {project_id: 'p1', title: ''}]));
  await toMapping();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(importButton(), 'click');
  await flush();
  assert.equal(issues(memory).length, 3, 'nothing was committed');
  const shown = document.body.querySelector('.u2-domain-import-report').textContent;
  assert.match(shown, /Import aborted — 1 row has errors/);
  assert.match(shown, /Uncheck "All or nothing" to import the 1 valid row and skip the rest\./,
    'an aborted all-or-nothing run says what to change');
  assert.match(shown, /Title: Value can't be empty/, 'the per-row report names the field the row was refused on');
  assert.match(shown, /Source row/, 'the issue lines are headed by the source row they are about');
  fire(buttonNamed('CLOSE'), 'click');
  const report = await running;
  assert.equal(report.error, 'validation', 'an allOrNothing abort answers the report, not a throw');
  assert.equal(report.errorCount, 1);
  assert.deepEqual(report.rows.map((r) => [r.index, r.status]), [[1, 'error']]);
});

scoped('a table-level insert denial still opens the wizard; no visible column refuses it', async () => {
  // the row-mode deviation: `can.insert` is a false negative there, so column security decides
  emsBackends.domain = backend({access: {can: {view: true, insert: false, edit: false, delete: false, share: false},
    fields: {project_id: 'editable', title: 'editable'}}});
  const open = await domains.table('grit.issue');
  const running = domains.import(open,
    {source: frame(['project_id', 'title'], [{project_id: 'p1', title: 'A'}])});
  await flush();
  await toMapping();
  assert.equal(named('title').value.value, 'title', 'the writable columns are offered');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);

  emsBackends.domain = backend({access: {can: {view: true, insert: true, edit: true, delete: true, share: true},
    fields: {}}});
  const hidden = await domains.table('grit.issue');
  assert.equal(await domains.import(hidden), null);
  assert.match(document.body.querySelector('.u2-notify-warning')?.textContent ?? '',
    /no column of grit.issue you may write/);
});

scoped('the preview draws the rows as they would land, each under its verdict, refused cells marked', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'priority'], [{project_id: 'p1', title: 'A', priority: 'low'},
      {project_id: 'p1', title: 'B', priority: 'urgent'}]));
  await toMapping();
  await toPreview();
  const rows = document.body.querySelector('.u2-domain-import-rows');
  rows.clientHeight = 200;
  fire(rows, 'scroll');
  await flush();
  assert.deepEqual([...rows.querySelectorAll('.u2-data-table-head')].map((c) => c.textContent),
    ['#', 'Result', 'Project id', 'Title', 'Priority', 'Note'], 'the verdict leads, then the target captions and why');
  const cells = [...rows.querySelectorAll('.u2-data-table-row')]
    .map((r) => [...r.children].map((c) => c.textContent).join('|'));
  assert.deepEqual(cells, ['1|Add|p1|A|low|', '2|Error|p1|B|urgent|Priority: Must be one of: low, high'],
    'what the backend said it would do with the row, beside the row as it would be sent');
  const bad = rows.querySelector('.u2-cell-error');
  assert.notEqual(bad, null, 'the cell the backend refuses is marked, by the column it named');
  assert.match(bad.title, /Must be one of: low, high/);
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
  fire(importButton(), 'click');
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

scoped('an aborted all-or-nothing run is not a dead end: BACK, the flag off, and the valid rows land', async () => {
  const memory = backend();
  noDryRun(memory);
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Fine'}, {project_id: 'p1', title: ''}]));
  await toMapping();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(importButton(), 'click');
  await flush();
  const back = [...document.body.querySelectorAll('.u2-domain-import-report button')]
    .find((b) => b.textContent === 'BACK');
  assert.notEqual(back, undefined, 'the report step offers the way back');
  fire(back, 'click');
  await flush();
  assert.notEqual(named('allOrNothing'), null, 'the source step is back');
  named('allOrNothing').value.value = false;
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  assert.equal(named('title').value.value, 'title', 'and the mapping it was left with stands');
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(importButton(), 'click');
  await flush();
  assert.equal(issues(memory).length, 4, 'the valid row lands once all-or-nothing is off');
  fire(buttonNamed('CLOSE'), 'click');
  const report = await running;
  assert.equal(report.inserted, 1);
  assert.equal(report.errorCount, 1);
});

scoped('the options offered are the ones the storage declares, and an option not offered is not sent', async () => {
  // an external binding: all-or-nothing is the only way it imports, a duplicate is always an
  // error, and this one has no upsert either
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support,
    batch: {...issue.support.batch, upsert: false, partial: false, skipDuplicates: false}};
  const posts = [];
  const real = issue.batch.bind(issue);
  issue.batch = (rows, options) => {
    posts.push(options);
    return real(rows, options);
  };
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Declared'}]));
  assert.equal(document.body.querySelector('[data-u2-name="allOrNothing"]'), null, 'no "All or nothing"');
  assert.equal(document.body.querySelector('[data-u2-name="errorOnDuplicate"]'), null, 'no "Error on duplicate"');
  assert.equal(document.body.querySelector('[data-u2-name="mode"]'), null, 'no upsert: no mode to choose');
  await toMapping();
  await toPreview();
  fire(importButton(), 'click');
  await flush();
  assert.deepEqual(posts, [{mode: 'insert'}], 'nothing the storage refuses is on the wire');
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).inserted, 1);
});

scoped('an aborted run on a storage without partial imports never recommends the flag it has not got', async () => {
  const memory = backend();
  const issue = noDryRun(memory);
  issue.support = {...issue.support, batch: {...issue.support.batch, partial: false}};
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Fine'}, {project_id: 'p1', title: ''}]));
  await toMapping();
  await toPreview();
  fire(importButton(), 'click');
  await flush();
  const shown = document.body.querySelector('.u2-domain-import-report').textContent;
  assert.match(shown, /Import aborted — 1 row has errors/);
  assert.doesNotMatch(shown, /Uncheck "All or nothing"/);
  assert.notEqual([...document.body.querySelectorAll('.u2-domain-import-report button')]
    .find((b) => b.textContent === 'BACK'), undefined, 'the way back is still the mapping');
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).error, 'validation');
});

scoped('an upsert on a storage that cannot tell insert from update reports merged rows', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.batch = async (rows) => ({inserted: 0, updated: 0, merged: rows.length, skipped: 0, errorCount: 0,
    rows: rows.map((_, index) => ({index, id: `k${index}`, status: 'merged'}))});
  const {running} = await opened(memory,
    frame(['project_id', 'number', 'title'], [{project_id: 'p1', number: 1, title: 'One'},
      {project_id: 'p1', number: 2, title: 'Two'}]));
  await toMapping();
  named('mode').value.value = 'upsert';
  await flush();
  await toPreview();
  fire(importButton(), 'click');
  await flush();
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent,
    /0 inserted, 2 merged, 0 skipped, 0 failed/);
  assert.match(document.body.querySelector('.u2-notify-info').textContent, /0 inserted, 2 merged/);
  fire(buttonNamed('CLOSE'), 'click');
  const report = await running;
  assert.equal(report.merged, 2);
  assert.deepEqual(report.rows.map((r) => r.status), ['merged', 'merged']);
});

scoped('a warehouse refusal — the failing rows and no totals — still renders the report and the way back', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.batch = async () => ({error: 'validation', rows: [{index: 1, id: null, status: 'error',
    errors: [{column: 'title', code: 'unique', message: 'duplicate key value violates unique constraint'}]}]});
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'One'}, {project_id: 'p1', title: 'One'}]));
  await toMapping();
  await toPreview();
  fire(importButton(), 'click');
  await flush();
  const shown = document.body.querySelector('.u2-domain-import-report').textContent;
  assert.match(shown, /Import aborted — 1 row has errors; nothing was committed\./);
  assert.match(shown, /Title is already used/, 'the row error is listed, a unique refusal in words');
  assert.notEqual([...document.body.querySelectorAll('.u2-domain-import-report button')]
    .find((b) => b.textContent === 'BACK'), undefined, 'the way back stands');
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).error, 'validation');
});

scoped('an importer maps, derives and keeps rows out, and makes the preview its own; the report says why', async () => {
  const memory = backend();
  emsBackends.domain = memory;
  const projects = await domains.table('grit.project');
  projects.importer = {derives: ['description'],
    match: (column) => column.name === 'code' ? 'key' : undefined,
    prepare: async (rows) => rows.map((r, i) =>
      i === 1 ? {existing: 'p1', label: 'Exists', note: 'Already there'} : i === 2 ? {sameAs: 0, note: 'Same as row 1'} :
        i === 3 ? {error: 'Cannot read it'} : {values: {description: `About ${r.name}`}}),
    columns: (defaults) => [...defaults, {name: 'loud', header: 'Loud', render: (row) => String(row.values.name).toUpperCase()}]};
  const running = domains.import(projects, {source: frame(['code', 'name'],
    [{code: 'NEW', name: 'New'}, {code: 'GRIT', name: 'Grit'}, {code: 'NEW', name: 'New'}, {code: 'BAD', name: 'Bad'}])});
  await flush();
  named('allOrNothing').value.value = false;
  named('writeBack').value.value = false;
  await toMapping();
  assert.equal(named('name').items.some((x) => x?.label === 'Description'), false, 'what the importer derives is not offered');
  assert.equal(named('code').value.value, 'key', 'the importer\'s match comes first');
  await toPreview();
  const rows = document.body.querySelector('.u2-domain-import-rows');
  rows.clientHeight = 400;
  fire(rows, 'scroll');
  await flush();
  const preview = document.body.querySelector('.u2-domain-import-preview').textContent;
  assert.match(preview, /1 will be added, 0 updated, 2 skipped, 1 row has errors\./);
  assert.match(preview, /Exists.*Already there.*Same as row 1.*Cannot read it/, 'the verdict\'s label on the badge');
  assert.match(preview, /Loud.*NEW.*GRIT/, 'and the importer\'s own column');
  fire(importButton(), 'click');
  await flush();
  await flush();
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent, /1 inserted, 0 updated, 2 skipped, 1 failed\./);
  fire(buttonNamed('CLOSE'), 'click');
  const report = await running;
  assert.deepEqual(report.rows.map((r) => [r.index, r.status]), [[0, 'inserted'], [1, 'skipped'], [2, 'skipped'], [3, 'error']]);
  assert.equal(memory.tableSync('grit.project').rows.find((r) => r.key === 'NEW')?.description, 'About New');
});

scoped('the preview\'s button names what lands, and is held with the reason while all-or-nothing would abort', async () => {
  const memory = backend();
  const {running} = await opened(memory,
    frame(['project_id', 'number', 'title', 'priority'], [{project_id: 'p1', number: '5', title: 'A', priority: 'low'},
      {project_id: 'p1', number: '6', title: 'B', priority: 'urgent'},
      {project_id: 'p1', number: '7', title: 'C', priority: 'high'}]));
  await toMapping();
  await toPreview();
  assert.equal(buttonNamed('NEXT'), undefined, 'the committing button is not one more NEXT');
  assert.equal(importButton().textContent, 'Import 2 issues', 'the rows the dry run says would land');
  assert.equal(importButton().disabled, true, 'all-or-nothing would abort: nothing to press');
  assert.equal(reason(), '1 row has errors. Uncheck "All or nothing" on Source to import the 2 valid rows ' +
    'and skip the rest.', 'the reason is where the mapping step gives its own');
  assert.doesNotMatch(document.body.querySelector('.u2-domain-import-preview').textContent, /Nothing will be imported/,
    'the body keeps the counts alone');

  fire(buttonNamed('BACK'), 'click');
  await flush();
  fire(buttonNamed('BACK'), 'click');
  await flush();
  named('allOrNothing').value.value = false;
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  await toPreview();
  assert.equal(importButton().disabled, false);
  assert.equal(reason(), '');
  fire(importButton(), 'click');
  await flush();
  assert.match(document.body.querySelector('.u2-domain-import-report').textContent, /2 inserted, 0 updated, 0 skipped, 1 failed/);
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).inserted, 2);
});

scoped('a storage without partial imports holds the import without naming a flag it has not got', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support, batch: {...issue.support.batch, partial: false}};
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Fine'}, {project_id: 'p1', title: ''}]));
  await toMapping();
  await toPreview();
  assert.equal(importButton().disabled, true);
  assert.equal(reason(), '1 row has errors, and this import is all or nothing.');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('the import waits for the dry run, and reads plain "Import" past the rows it checks', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  let answer;
  issue.validate = (rows) => new Promise((resolve) => answer = () => resolve({validateOnly: true,
    rowCount: rows.length, willInsert: rows.length, willUpdate: 0, willSkip: 0, errorCount: 0,
    rows: rows.map((_, index) => ({index, predicted: 'insert'}))}));
  const rows = Array.from({length: 1001}, (_, i) => ({project_id: 'p1', title: `T${i}`}));
  const {running} = await opened(memory, frame(['project_id', 'title'], rows));
  await toMapping();
  await toPreview();
  assert.equal(importButton().disabled, true, 'held, silently, while the rows are checked');
  assert.equal(reason(), '');
  answer();
  await flush();
  assert.equal(importButton().disabled, false);
  assert.equal(importButton().textContent, 'Import', 'the dry run saw the first 1,000 rows: no count');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('a row the dry run skips says why in the preview', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.validate = async (rows) => ({validateOnly: true, rowCount: rows.length, willInsert: 1, willUpdate: 0,
    willSkip: 1, errorCount: 0, rows: [{index: 0, predicted: 'insert'}, {index: 1, predicted: 'skip', existingId: 'i2'}]});
  const {running} = await opened(memory,
    frame(['project_id', 'number', 'title'], [{project_id: 'p1', number: 9, title: 'New'},
      {project_id: 'p1', number: 2, title: 'Ibuprofen'}]));
  await toMapping();
  await toPreview();
  const table = document.body.querySelector('.u2-domain-import-rows');
  table.clientHeight = 200;
  fire(table, 'scroll');
  await flush();
  const [first, second] = [...table.querySelectorAll('.u2-data-table-row')];
  assert.equal([...second.children].at(-1).textContent, 'Duplicate business key');
  assert.equal(first.children[4].title, 'New', 'a cell a narrow column cuts reads in full on hover');
  assert.equal(importButton().textContent, 'Import 1 issue');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});
