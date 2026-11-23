// the memory backend has no bulk endpoint: the posting path here is the transaction fallback, the same rows either way
import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {layout} from '@datagrok-libraries/u2/tests/layout.mjs';
import {backends, Scope, Control, notify} from '@datagrok-libraries/u2';
import {SCHEMA, ROWS, backend} from './domain-fixtures.mjs';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';

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
const reportText = () => document.body.querySelector('.u2-domain-import-report').textContent;
const issues = (memory) => memory.tableSync('grit.issue').rows;
/** The outcome a step leads with, one entry per count. */
const counts = (step) => [...document.body.querySelectorAll(`.u2-domain-import-${step} .u2-domain-import-count`)]
  .map((c) => c.textContent);
/** The report's BACK, in the wizard's footer — shown only after an aborted run. */
const reportBack = () => [...document.body.querySelectorAll('.u2-dialog .u2-wizard-buttons button')]
  .find((b) => b.textContent === 'BACK' && b.style.display !== 'none');

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
  const running = testOwner().build(() => domains.import(handle));
  await flush();
  named('source').value = source;
  await flush();
  // the mapping table has no layout in the shim: give it a viewport so its rows render
  return {running, handle};
}

async function toMapping() {
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  layout(mapping(), {clientHeight: 400});
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
  assert.equal(named('title').value, 'title');
  assert.equal(named('Priority').value, 'priority', 'matched on the caption, case-insensitively');
  assert.equal(named('nonesuch').value, '(skip)');
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
  assert.equal(named('number').value, 'number', 'the key typed once on insert is matched by name');
  const labels = named('title').items.map((x) => typeof x === 'string' ? x : x.label);
  assert.equal(labels.includes('Number'), true, 'the key is offered as a target');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('a frame handed in as an option is mapped while the wizard is still being built', async () => {
  const memory = backend();
  emsBackends.domain = memory;
  const handle = await domains.table('grit.issue');
  const running = testOwner().build(() => domains.import(handle,
    {source: frame(['project_id', 'title'], [{project_id: 'p1', title: 'Handed in'}])}));
  await flush();
  await toMapping();
  assert.equal(named('title').value, 'title', 'the mapping is there without a pick');
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
  named('name2').value = 'title';
  await flush();
  assert.match(await proceed(), /Columns "title" and "name2" are both mapped to "Title"/,
    'the target by its caption, both columns by name');

  for (const column of ['project_id', 'title', 'name2'])
    named(column).value = '(skip)';
  await flush();
  assert.match(await proceed(), /Map at least one column/);

  named('title').value = 'description';
  await flush();
  assert.match(await proceed(), /Required column "Project id" is not mapped/,
    'insert mode names every required column left unmapped');

  named('project_id').value = 'project_id';
  await flush();
  assert.match(await proceed(), /Required column "Title" is not mapped/);

  named('title').value = 'title';
  named('mode').value = 'upsert';
  await flush();
  assert.match(await proceed(), /Add or update matches rows by Number — map a column to "Number"/);
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
  assert.deepEqual(counts('preview'), ['1 to add', '1 error'], 'one vocabulary, zeros left out');
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
  named('note').value = 'description';
  await flush();
  named('note').value = '(skip)';
  await flush();
  assert.equal(posts.length, 1, 'typing on another step posts nothing');

  await toPreview();
  assert.equal(posts.length, 1, 'the same mapping is the same dry run');
  fire(buttonNamed('BACK'), 'click');
  await flush();
  named('note').value = 'description';
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
  assert.deepEqual(counts('report'), ['1 added']);
  fire(buttonNamed('CLOSE'), 'click');
  assert.equal((await running).inserted, 1);
});

scoped('a storage with no dry run previews the mapped rows without verdicts, and the import still posts', async () => {
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support, batch: {...issue.support.batch, validate: false}};
  issue.validate = undefined;
  const {running} = await opened(memory,
    frame(['project_id', 'title', 'note'], [{project_id: 'p1', title: 'Unchecked', note: 'x'}]));
  await toMapping();
  await toPreview();
  const preview = document.body.querySelector('.u2-domain-import-preview');
  assert.match(preview.textContent, /Rows are checked when imported\./);
  const rows = preview.querySelector('.u2-domain-import-rows');
  layout(rows, {clientHeight: 200});
  fire(rows, 'scroll');
  await flush();
  assert.deepEqual([...rows.querySelectorAll('.u2-data-table-head')].map((c) => c.textContent),
    ['Project id', 'Title'], 'the mapped columns under their captions, no verdict column');
  assert.deepEqual([...rows.querySelectorAll('.u2-data-table-row')]
    .map((r) => [...r.children].map((c) => c.textContent).join('|')), ['p1|Unchecked']);
  assert.equal(reason(), '', 'nothing gates NEXT');
  fire(importButton(), 'click');
  await flush();
  assert.deepEqual(counts('report'), ['1 added']);
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
  assert.equal(named('note').value, '(skip)');
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  fire(importButton(), 'click');
  await flush();
  const landed = issues(memory).find((r) => r.title === 'Imported');
  assert.notEqual(landed, undefined, 'the row was inserted');
  assert.equal('note' in landed, false, 'a skipped source column never reached the backend');
  assert.deepEqual(counts('report'), ['1 added']);
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
  named('mode').value = 'upsert';
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
  const shown = reportText();
  assert.match(shown, /Import aborted — 1 row has errors/);
  assert.match(shown, /Uncheck "All or nothing" to import the 1 valid row and skip the rest\./,
    'an aborted all-or-nothing run says what to change');
  assert.match(shown, /Title: Value can't be empty/, 'the per-row report names the field the row was refused on');
  const table = document.body.querySelector('.u2-domain-import-report .u2-domain-import-rows');
  assert.deepEqual([...table.querySelectorAll('.u2-data-table-head')].map((c) => c.textContent), ['#', 'Result', 'Note'],
    'the preview\'s headers, and no key column where no row has a key');
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
  const running = testOwner().build(() => domains.import(open,
    {source: frame(['project_id', 'title'], [{project_id: 'p1', title: 'A'}])}));
  await flush();
  await toMapping();
  assert.equal(named('title').value, 'title', 'the writable columns are offered');
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
  layout(rows, {clientHeight: 200});
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
  assert.equal(named('Project ID').value, 'project_id');
  assert.equal(named('Chemical title').value, 'title');
  assert.equal(named('Weight (kg)').value, 'weight');
  assert.equal(named('Batch').value, '(skip)');
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
  const back = reportBack();
  assert.notEqual(back, undefined, 'the report step offers the way back, in its footer');
  fire(back, 'click');
  await flush();
  assert.notEqual(named('allOrNothing'), null, 'the source step is back');
  named('allOrNothing').value = false;
  await flush();
  fire(buttonNamed('NEXT'), 'click');
  await flush();
  assert.equal(named('title').value, 'title', 'and the mapping it was left with stands');
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
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support, batch: {...issue.support.batch, partial: false}};
  const {running} = await opened(memory,
    frame(['project_id', 'title'], [{project_id: 'p1', title: 'Fine'}, {project_id: 'p1', title: ''}]));
  await toMapping();
  await toPreview();
  assert.equal(reason(), '1 row has errors, and this import is all or nothing.', 'the fixed flag reads as on, unnamed');
  noDryRun(memory);
  fire(buttonNamed('BACK'), 'click');
  await flush();
  await toPreview();
  fire(importButton(), 'click');
  await flush();
  const shown = reportText();
  assert.match(shown, /Import aborted — 1 row has errors/);
  assert.doesNotMatch(shown, /Uncheck "All or nothing"/);
  assert.notEqual(reportBack(), undefined, 'the way back is still the mapping');
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
  named('mode').value = 'upsert';
  await flush();
  await toPreview();
  fire(importButton(), 'click');
  await flush();
  assert.deepEqual(counts('report'), ['2 merged']);
  assert.equal(document.body.querySelector('.u2-notify-info'), null, 'the open report is the answer, not a balloon');
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
  const shown = reportText();
  assert.match(shown, /Import aborted — 1 row has errors; nothing was committed\./);
  assert.match(shown, /Title is already used/, 'the row error is listed, a unique refusal in words');
  assert.notEqual(reportBack(), undefined, 'the way back stands');
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
  const running = testOwner().build(() => domains.import(projects, {source: frame(['code', 'name'],
    [{code: 'NEW', name: 'New'}, {code: 'GRIT', name: 'Grit'}, {code: 'NEW', name: 'New'}, {code: 'BAD', name: 'Bad'}])}));
  await flush();
  named('allOrNothing').value = false;
  named('writeBack').value = false;
  await toMapping();
  assert.equal(named('name').items.some((x) => x?.label === 'Description'), false, 'what the importer derives is not offered');
  assert.equal(named('code').value, 'key', 'the importer\'s match comes first');
  await toPreview();
  const rows = document.body.querySelector('.u2-domain-import-rows');
  layout(rows, {clientHeight: 400});
  fire(rows, 'scroll');
  await flush();
  const preview = document.body.querySelector('.u2-domain-import-preview').textContent;
  assert.deepEqual(counts('preview'), ['1 to add', '2 to skip', '1 error']);
  assert.match(preview, /Skip.*Already there.*Skip.*Same as row 1.*Cannot read it/,
    'the badge says the verb the counts use, the verdict says why');
  assert.match(preview, /Loud.*NEW.*GRIT/, 'and the importer\'s own column');
  fire(importButton(), 'click');
  await flush();
  await flush();
  assert.deepEqual(counts('report'), ['1 added', '2 skipped', '1 error']);
  fire(buttonNamed('CLOSE'), 'click');
  const report = await running;
  assert.deepEqual(report.rows.map((r) => [r.index, r.status]), [[0, 'inserted'], [1, 'skipped'], [2, 'skipped'], [3, 'error']]);
  assert.equal(memory.tableSync('grit.project').rows.find((r) => r.key === 'NEW')?.description, 'About New');
});

/** The server itemizes at most 1,000 rows of a batch report (`maxBatchReportRows` in `batch_loader.dart`) —
 * errors, then duplicates, then updates, then inserts — and counts them all. */
function capped(table) {
  const batch = table.batch.bind(table);
  const rank = {error: 0, duplicate: 1, updated: 2, merged: 2, inserted: 3};
  table.batch = async (rows, options) => {
    const report = await batch(rows, options);
    return {...report, rows: [...report.rows].sort((a, b) => rank[a.status] - rank[b.status] || a.index - b.index)
      .slice(0, 1000)};
  };
}

/** Imports `source` into `grit.project` and answers the report's text and the column written back. */
async function importProjects(memory, source, allOrNothing, mode = 'insert') {
  // the frame double has no string column to add: the write-back's column is kept here
  const written = [];
  source.columns.getUnusedName = (name) => source.columns.names().some((n) => n.toLowerCase() === name.toLowerCase()) ?
    `${name} (2)` : name;
  source.columns.addNewString = (name) => ({name, init(value) {
    for (let i = 0; i < source.rowCount; i++)
      written.push(value(i));
    return this;
  }});
  emsBackends.domain = memory;
  const projects = await domains.table('grit.project');
  const running = testOwner().build(() => domains.import(projects, {source}));
  await flush();
  named('allOrNothing').value = allOrNothing;
  await toMapping();
  named('mode').value = mode;
  await flush();
  await toPreview();
  fire(importButton(), 'click');
  for (let i = 0; i < 5; i++)
    await flush();
  const shown = reportText();
  fire(buttonNamed('CLOSE'), 'click');
  await running;
  return {shown, written};
}

const projectRows = (count, prefix) => Array.from({length: count}, (_, i) =>
  ({key: `${prefix}-${String(i + 1).padStart(4, '0')}`, name: `Project ${i + 1}`}));

scoped('a commit past the 1,000 rows the server itemizes still writes back every row\'s key', async () => {
  const memory = backend();
  capped(memory.tableSync('grit.project'));
  // the existing GRIT is itemized first (a duplicate), so the cap leaves out more of the inserts
  const rows = [...projectRows(1200, 'BIG'), {key: 'GRIT', name: 'Grit again'}];
  const {shown, written} = await importProjects(memory, frame(['key', 'name'], rows), false);
  assert.match(shown, /1,200 added/);
  assert.equal(written.filter((key, i) => key !== rows[i].key).length, 0,
    'every imported row has its key written back, past the 1,000 rows the server itemizes');
});

const UID = '5d1f2a4e-0b6c-4c1e-9a53-7a3f0e6b9c10';

// the server's guid check: dashed, either case, untrimmed; an empty or missing cell is no reference
scoped('a user column takes ids as the server reads them: a login or a padded id is refused by name', async () => {
  const row = (number, reporter) => ({project_id: 'p1', number: String(number), title: `T${number}`, reporter});
  const {running} = await opened(backend(), frame(['project_id', 'number', 'title', 'reporter'],
    [row(11, 'admin'), row(12, ` ${UID} `), row(13, UID.toUpperCase()), row(14, ''), row(15, null)]));
  await toMapping();
  await toPreview();
  assert.deepEqual(counts('preview'), ['3 to add', '2 errors'], 'a login or a padded id is refused in the preview');
  assert.match(document.body.textContent, /Reporter: expected a user id \(uuid\), got "admin"/);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

scoped('a group column takes group ids: a group name is refused by name, a group id passes', async () => {
  const issue = SCHEMA.tables.issue;
  const memory = new MemoryDomainBackend({...SCHEMA, tables: {...SCHEMA.tables, issue: {...issue,
    columns: {...issue.columns, team: {type: 'group'}}}}},
  {rows: Object.fromEntries(Object.entries(ROWS).map(([t, list]) => [t, list.map((r) => ({...r}))]))});
  const {running} = await opened(memory, frame(['project_id', 'title', 'team'],
    [{project_id: 'p1', title: 'A', team: 'Developers'}, {project_id: 'p1', title: 'B', team: UID}]));
  await toMapping();
  await toPreview();
  assert.deepEqual(counts('preview'), ['1 to add', '1 error'], 'a group name is refused, a group id passes');
  assert.match(document.body.textContent, /Team: expected a group id \(uuid\), got "Developers"/);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
});

/** The server's answer to a batch a table CHECK fails (`repository.dart` `_mapPgError`): the statement aborts, so
 * the dry run and the commit both get HTTP 409 `{error: 'validation', rows: [...]}` pinned on row 0, which the
 * js-api client resolves rather than rejects — no `validateOnly`, no counts, no `status` or `predicted`. */
function checkRefusing(memory, {constraint = 'weight_positive', message = 'Weight is positive',
  when = (rows) => rows.some((r) => Number(r.weight) < 0)} = {}) {
  const issue = memory.tableSync('grit.issue');
  const refuse = (rows) => when(rows) ? {error: 'validation', rows: [{index: 0, id: null,
    errors: [{column: null, code: 'check', message, constraint}]}]} : null;
  const validate = issue.validate?.bind(issue);
  const batch = issue.batch.bind(issue);
  if (validate !== undefined)
    issue.validate = async (rows, options) => refuse(rows) ?? validate(rows, options);
  issue.batch = async (rows, options) => refuse(rows) ?? batch(rows, options);
}

const weights = (...list) => frame(['project_id', 'title', 'weight'],
  list.map((weight, i) => ({project_id: 'p1', title: `QA-${i + 1}`, weight})));

async function rowCells(step) {
  const rows = document.body.querySelector(`.u2-domain-import-${step} .u2-domain-import-rows`);
  layout(rows, {clientHeight: 400});
  fire(rows, 'scroll');
  await flush();
  return [...rows.querySelectorAll('.u2-data-table-row')].map((r) => [...r.children].map((c) => c.textContent).join('|'));
}

/** The preview of `source` over a backend whose dry run is refused as `refusal` says; answers what it shows. */
async function refusedPreview(source, refusal, memory = backend()) {
  checkRefusing(memory, refusal);
  const {running} = await opened(memory, source);
  await toMapping();
  await toPreview();
  const shown = {cells: await rowCells('preview'), counts: counts('preview'), text: document.body.textContent,
    disabled: importButton().disabled, button: importButton().textContent};
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  return shown;
}

scoped('a dry run the server refuses whole over a CHECK puts it on the rows that break it, as SQL reads NULL',
  async () => {
    const shown = await refusedPreview(weights('1', '-1', '', '-3'));
    assert.doesNotMatch(shown.button, /NaN/, 'the import button never counts NaN rows');
    assert.deepEqual(shown.cells, ['1||p1|QA-1|1|', '2|Error|p1|QA-2|-1|Weight is positive', '3||p1|QA-3||',
      '4|Error|p1|QA-4|-3|Weight is positive'], 'the message is under the rows that break it, not on row 1');
    assert.equal(shown.disabled, true, 'nothing would land: the button is off');
  });

// a CHECK the browser cannot place stays on the batch: never on the server's arbitrary row 0, never a button on
scoped('a refused CHECK no row can be named for stays on the batch and keeps the button off', async () => {
  const issue = SCHEMA.tables.issue;
  const floored = new MemoryDomainBackend({...SCHEMA, tables: {...SCHEMA.tables, issue: {...issue,
    constraints: {...issue.constraints, floor: {expr: 'weight >= $min', message: 'Under the floor'}}}}},
  {rows: Object.fromEntries(Object.entries(ROWS).map(([t, list]) => [t, list.map((r) => ({...r}))]))});
  // undeclared, unevaluable ($min is unbound), and one no row breaks
  for (const [weight, refusal, memory] of [['-1', {constraint: 'nonesuch'}], ['-1', {constraint: 'floor'}, floored],
    ['2', {when: () => true}]]) {
    const shown = await refusedPreview(weights('1', weight), refusal, memory);
    assert.equal(shown.cells.some((c) => /Error/.test(c)), false, `${refusal.constraint ?? 'no row breaks it'}: no row is guessed`);
    assert.match(shown.text, /it does not say which row, and no row here breaks it/);
    assert.equal(shown.disabled, true);
  }
});

scoped('a commit the server refuses whole reports each row that breaks the CHECK as an error with its message',
  async () => {
    const memory = backend();
    noDryRun(memory);
    checkRefusing(memory);
    const {running} = await opened(memory, weights('1', '-1', '2', '-3'));
    named('allOrNothing').value = false;
    await toMapping();
    await toPreview();
    fire(importButton(), 'click');
    for (let i = 0; i < 3; i++)
      await flush();
    assert.deepEqual((await rowCells('report')).map((r) => r.split('|').filter((c) => c !== '').join('|')),
      ['1|Not imported|Another row has errors', '2|Error|Weight is positive', '3|Not imported|Another row has errors',
        '4|Error|Weight is positive'], 'a row carrying errors is an error, whatever its status');
    fire(buttonNamed('CLOSE'), 'click');
    await running;
  });

scoped('a commit refused whole names the offender past the 1,000 rows previewed, not "uncheck All or nothing"',
  async () => {
    const memory = backend();
    checkRefusing(memory);
    const {running} = await opened(memory, weights(...Array.from({length: 1005}, (_, i) => i === 1003 ? '-2' : '1')));
    await toMapping();
    await toPreview();
    fire(importButton(), 'click');
    for (let i = 0; i < 3; i++)
      await flush();
    assert.doesNotMatch(reportText(), /Uncheck "All or nothing"/, 'unchecking it cannot help a refused statement');
    fire(buttonNamed('CLOSE'), 'click');
    const report = await running;
    assert.deepEqual(report.rows.filter((r) => r.errors?.length).map((r) => r.index), [1003]);
  });

/** Opens `grit.project` with `source`, goes to the preview and answers the button and the reason. */
async function projectPreview(source) {
  emsBackends.domain = backend();
  const running = testOwner().build(async () => domains.import(await domains.table('grit.project'), {source}));
  await flush();
  await toMapping();
  await toPreview();
  const shown = {disabled: importButton().disabled, counts: counts('preview')};
  fire(buttonNamed('CANCEL'), 'click');
  await running;
  return shown;
}

const existing = (count) => Array.from({length: count}, () => ({key: 'GRIT', name: 'Grit again'}));

scoped('every previewed row skipped keeps the button off, up to the preview cap and not past it', async () => {
  const atCap = await projectPreview(frame(['key', 'name'], existing(1000)));
  assert.deepEqual(atCap.counts, ['1,000 to skip']);
  assert.equal(atCap.disabled, true, '1,000 rows are all judged and none would land: the button is off');
  const past = await projectPreview(frame(['key', 'name'], existing(1001)));
  assert.equal(past.disabled, false, 'the row past the 1,000 previewed is unjudged: the commit decides');
});

scoped('a commit that lands nothing adds no write-back column', async () => {
  const {shown, written} = await importProjects(backend(), frame(['key', 'name'], existing(1001)), true);
  assert.match(shown, /1,001 skipped/);
  assert.deepEqual(written, [], 'no column was added to the user\'s table');
});
