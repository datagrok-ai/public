/* `domains.bulkEdit` (WO 3-8) on the memory backend: only targets that can be written are offered
   (never "all N matching" while a search the server cannot take is set), one aligned form whose
   field picker decides what is sent, a muted line under each chosen field (the rows' value, how
   many differ, or "will be cleared"), the button caption following the target, whatever blocks it
   its tooltip, the table's validators over the rows as they would land, refusals in the dialog
   rather than a balloon, the unsaved gate before it opens, and the report in the table's words. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Control, Tooltip, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {SCHEMA, backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');

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
const control = (name) => Control.forElement(document.querySelector(`.u2-dialog [data-u2-name="${name}"]`));
const field = (name) => control(name);
const shown = (name) => document.querySelector(`.u2-dialog [data-u2-name="${name}"]`)?.closest('[hidden]') === null;
const hint = (name) => field(name).root.querySelector('.u2-domain-bulk-current').textContent;
const balloon = (kind) => document.body.querySelector(`.u2-notify-${kind}`)?.textContent ?? '';
const issues = (memory) => memory.tableSync('grit.issue').rows;

/** The field picker: what the edit changes, by column name. */
const picker = () => control('fields');
const offered = () => picker().getItems({all: true}).list.map((i) => i.key);
function choose(...names) {
  picker().value.value = names;
}

/** "Apply to": a choice where there are two targets, plain text where there is one. */
const target = () => control('target');
const targets = () => target().items;
const targetText = () => document.querySelector('.u2-dialog [data-u2-name="target"]').textContent;
/** The primary button — its caption follows the target. */
const ok = () => [...document.querySelectorAll('.u2-dialog [data-u2-part="footer"] button')].at(-1);
const okDisabled = () => ok().disabled === true;
/** A refused write, said inside the dialog. */
const refusal = () => document.querySelector('.u2-dialog [data-u2-part="error-row"]').textContent;
/** What blocks the button: its tooltip, never a red line in the body. */
async function okReason() {
  fire(ok(), 'pointerenter');
  await new Promise((r) => setTimeout(r, 700));
  const text = Tooltip.root.textContent;
  fire(ok(), 'pointerleave');
  return text;
}

/** The dialog up and rendered; the promise is handed back WRAPPED — an async function
 * returning it would chain on it and never come back until the dialog is answered. */
async function opened(source) {
  const running = domains.bulkEdit(source);
  await flush();
  return {running};
}

async function table(memory, options = {}) {
  emsBackends.domain = memory;
  const handle = await domains.table('grit.issue');
  const source = handle.source(options);
  await flush();
  return {handle, source};
}

/** `n` issues on one project, all open — enough of them to reach the server's 1000-row cap. */
function many(n) {
  return new MemoryDomainBackend(SCHEMA, {rows: {
    project: [{id: 'p1', key: 'GRIT', name: 'Grit'}],
    issue: Array.from({length: n}, (_, i) =>
      ({id: `x${i}`, project_id: 'p1', number: i + 1, title: `T${i}`, done: false})),
  }});
}

function select(source, keep) {
  const df = source.df.value;
  df.selection = {get: keep};
  df.onSelectionChanged.fire(undefined);
}

scoped('the chosen columns are written to the selection; one not chosen is left alone', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'number > 0'});
  select(source, (i) => i !== 1);
  assert.deepEqual(source.selection.value.map((r) => r.id), ['i1', 'i3']);

  const {running} = await opened(source);
  assert.deepEqual(targets(), ['2 selected issues', 'All 3 issues matching the filter']);
  assert.equal(target().value.value, '2 selected issues', 'a selection is the default target');
  assert.equal(ok().textContent, 'Update 2 issues', 'the button says what it does');
  assert.equal(offered().includes('id'), false, 'a system column is not offered');
  assert.equal(offered().includes('number'), false, 'a business key is not offered for two rows');
  assert.equal(offered().includes('project_id'), false);
  assert.equal(shown('done'), false, 'a field not chosen is not in the form');
  assert.equal(okDisabled(), true);
  assert.equal(await okReason(), 'Choose the fields to change');
  assert.equal(refusal(), '', 'nothing red before the user did anything');
  assert.equal(document.querySelector('.u2-domain-bulk-include, .u2-domain-bulk-reason'), null);

  choose('done');
  await flush();
  assert.equal(shown('done'), true, 'a chosen field shows its editor');
  assert.equal(shown('title'), false);
  assert.equal(okDisabled(), false);
  assert.equal(hint('done'), '2 different values', 'what the write overwrites');
  field('done').value.value = true;
  field('title').value.value = 'NOT WRITTEN';
  fire(ok(), 'click');
  assert.equal(await running, 1, 'Aspirin is already done: only the row the values change is written');
  assert.match(balloon('info'), /Updated 1 issue\./, 'the table\'s words, never "rows"');

  const rows = issues(memory);
  assert.deepEqual(rows.map((r) => r.done), [true, false, true], 'the two selected rows only');
  assert.deepEqual(rows.map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen'],
    'the column not chosen was not in the request');
  source.dispose();
});

scoped('under each chosen field: the shared value, how many differ, or that it will be cleared', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'number > 0'});
  select(source, (i) => i !== 0);
  const {running} = await opened(source);
  choose('title', 'done', 'weight');
  await flush();
  assert.equal(hint('weight'), 'Will be cleared', 'a chosen field left empty clears the column');
  assert.equal(hint('title'), 'Required: it cannot be cleared', 'unless the column is required');
  assert.equal(field('title').root.querySelector('.u2-input-error')?.textContent ?? '', '',
    'nothing red in the body: the button\'s tooltip says it');
  field('title').value.value = 'Renamed';
  await flush();
  assert.equal(hint('title'), '2 different values');
  assert.equal(hint('done'), 'Currently false', 'both selected issues are open');

  target().value.value = 'All 3 issues matching the filter';
  await flush();
  assert.equal(ok().textContent, 'Update 3 issues', 'the caption follows the target');
  assert.equal(hint('title'), '', 'the matching rows are not loaded: no current value to show');
  assert.equal(hint('weight'), 'Will be cleared');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  assert.deepEqual(issues(memory).map((r) => r.version), [1, 1, 1], 'nothing was written');
  source.dispose();
});

scoped('one target is plain text, and one row is offered its business key', async () => {
  const memory = backend();
  const {source} = await table(memory);
  select(source, (i) => i === 1);
  const {running} = await opened(source);
  assert.equal(targetText(), 'Apply to1 selected issue', 'a caption and its text, no choice of one');
  assert.equal(document.querySelector('.u2-dialog [data-u2="radio-input"]'), null);
  assert.equal(offered().includes('number'), true, 'a business key can change on one row');
  assert.equal(ok().textContent, 'Update 1 issue');
  choose('title');
  await flush();
  field('title').value.value = 'X';
  await flush();
  assert.equal(hint('title'), 'Currently Ibuprofen');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('only a chosen field holds the button: one taken back with a bad value no longer blocks it', async () => {
  const memory = backend();
  const {source} = await table(memory);
  select(source, (i) => i === 1);
  const {running} = await opened(source);
  choose('number', 'done');
  await flush();
  field('number').value.value = 0;
  await flush();
  assert.equal(okDisabled(), true, 'a chosen field with a bad value holds the button');
  choose('done');
  await flush();
  assert.equal(shown('number'), false);
  assert.equal(okDisabled(), false, 'the field taken back is not in the request, nor in the verdict');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('a selection past the cap disables the button and says how many are selected', async () => {
  const memory = many(1001);
  const {source} = await table(memory, {pageSize: 1100});
  select(source, () => true);
  assert.equal(source.selection.value.length, 1001);
  const {running} = await opened(source);
  choose('done');
  await flush();
  field('done').value.value = true;
  assert.equal(okDisabled(), true);
  assert.match(await okReason(), /at most 1000 issues at a time; 1001 are selected/);
  assert.equal(issues(memory).filter((r) => r.done === true).length, 0, 'nothing was written');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('"all matching" writes the list\'s filter, capped, and reports that more match', async () => {
  const memory = many(1001);
  const {source} = await table(memory, {query: 'done = false'});
  const {running} = await opened(source);
  assert.equal(targetText(), 'Apply toAll 1,001 issues matching the filter', 'the filter is the only target');
  assert.equal(ok().textContent, 'Update 1000 of 1,001 issues', 'one edit writes at most 1000');
  choose('done');
  await flush();
  field('done').value.value = true;
  fire(ok(), 'click');
  assert.equal(await running, 1000);
  assert.match(balloon('info'), /Updated 1000 issues — more match/);
  assert.equal(issues(memory).filter((r) => r.done === true).length, 1000);
  source.dispose();
});

scoped('the table\'s validators run over the selection as it would land; a refusal stays in the dialog', async () => {
  const memory = backend();
  const {handle, source} = await table(memory, {query: 'number > 0'});
  handle.validators.add('done', (value, row) => value === true && !row.priority ? 'Prioritise before closing' : null);
  select(source, (i) => i !== 0);
  const {running} = await opened(source);
  choose('done');
  await flush();
  field('done').value.value = true;
  fire(ok(), 'click');
  await flush();
  assert.match(refusal(), /Ibuprofen: Done: Prioritise before closing/, 'the row and the field, named');
  assert.equal(balloon('error'), '', 'never a balloon');
  assert.deepEqual(issues(memory).map((r) => r.version), [1, 1, 1], 'nothing was written');

  choose('done', 'priority');
  await flush();
  field('priority').value.value = 'high';
  fire(ok(), 'click');
  assert.equal(await running, 2, 'the written priority is what the rule reads');
  assert.deepEqual(issues(memory).map((r) => r.done), [true, true, true]);
  source.dispose();
});

scoped('"all matching" reads the rows to check them, and refuses past what one edit writes', async () => {
  const memory = backend();
  const {handle, source} = await table(memory, {query: 'done = false'});
  handle.validators.add('done', (value, row) => value === true && !row.priority ? 'Prioritise before closing' : null);
  const {running} = await opened(source);
  choose('done');
  await flush();
  field('done').value.value = true;
  fire(ok(), 'click');
  await flush();
  assert.match(refusal(), /Ibuprofen: Done: Prioritise before closing/);
  assert.deepEqual(issues(memory).map((r) => r.version), [1, 1, 1]);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();

  const big = await table(many(1001), {query: 'done = false'});
  big.handle.validators.add('done', () => null);
  const second = await opened(big.source);
  choose('done');
  await flush();
  assert.equal(okDisabled(), true);
  assert.match(await okReason(), /matches more than 1000 issues/);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await second.running, null);
  big.source.dispose();
});

scoped('a server refusal is the dialog\'s error line, and keeps what the user chose', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'number > 0'});
  select(source, () => true);
  const {running} = await opened(source);
  choose('title');
  await flush();
  fire(ok(), 'click');
  await flush();
  await flush();
  assert.notEqual(refusal(), '', 'the server refused clearing a required column');
  assert.equal(balloon('error'), '');
  assert.equal(shown('title'), true, 'the dialog is still up, the field still chosen');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('nothing chosen keeps the button disabled; CANCEL resolves null', async () => {
  const memory = backend();
  const {source} = await table(memory);
  select(source, () => true);
  const {running} = await opened(source);
  assert.equal(okDisabled(), true);
  assert.equal(balloon('warning'), '', 'a refusal is not a balloon on the far edge of the screen');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  assert.deepEqual(issues(memory).map((r) => r.version), [1, 1, 1], 'nothing was written');
  source.dispose();
});

scoped('a search with a selection still offers the selection alone', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'done = false'});
  source.search.value = 'ibu';
  await flush();
  select(source, () => true);
  const {running} = await opened(source);
  assert.equal(targetText(), 'Apply to1 selected issue', 'the filter target is not offered under a search');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('no target the server would take: the dialog never opens, and says the way out', async () => {
  const memory = backend();
  const {source} = await table(memory);
  assert.equal(await domains.bulkEdit(source), null, 'nothing selected and nothing filtered');
  assert.match(balloon('warning'), /Select the issues to write, or filter the list/);
  assert.equal(picker(), undefined, 'no dialog promising a target it would then refuse');

  notify.closeAll();
  source.query.value = 'done = false';
  source.search.value = 'asp';
  await flush();
  assert.equal(await domains.bulkEdit(source), null);
  assert.match(balloon('warning'), /clear the search box/,
    'a search the update endpoint cannot take is said BEFORE the dialog, not after OK');
  source.dispose();
});

scoped('the unsaved gate runs before the dialog opens', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'done = false'});
  source.edit.value.setValue('i2', 'title', 'Edited');
  assert.equal(source.isDirty.value, true);
  const running = domains.bulkEdit(source);
  await flush();
  assert.equal(picker(), undefined, 'the bulk dialog is not up yet');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null, 'a cancelled gate cancels the whole thing');

  const second = domains.bulkEdit(source);
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(source.isDirty.value, false);
  assert.notEqual(picker(), undefined, 'and a discard lets it through');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await second, null);
  source.dispose();
});

scoped('the list\'s own selection is the collection\'s, and only a user makes one', async () => {
  const memory = backend();
  emsBackends.domain = memory;
  const handle = await domains.table('grit.issue');
  const source = handle.source({});
  const list = domains.list(source);
  document.body.append(list.root);
  list.root.clientHeight = 400;
  await flush();
  assert.deepEqual(source.selection.value, [], 'a cold list selects nothing: the load put the lead there');
  assert.notEqual(source.currentRow.value, null, 'though it has a current row');

  const row = (id) => list.root.querySelector(`[data-u2-row="${id}"]`);
  fire(row('i1'), 'click');
  await flush();
  assert.deepEqual(source.selection.value.map((r) => r.id), ['i1'], 'a click is a selection');
  fire(row('i3'), 'click', {ctrlKey: true});
  await flush();
  assert.deepEqual(source.selection.value.map((r) => r.id), ['i1', 'i3']);

  const {running} = await opened(source);
  assert.equal(targetText(), 'Apply to2 selected issues');
  choose('done');
  await flush();
  field('done').value.value = true;
  fire(ok(), 'click');
  assert.equal(await running, 1, 'Aspirin is already done');
  assert.deepEqual(issues(memory).map((r) => r.done), [true, false, true]);

  assert.deepEqual(source.selection.value.map((r) => r.id), ['i1', 'i3'],
    'the rows a bulk edit just wrote are still the rows the user picked: the refresh keeps them');

  fire(list.list.root, 'keydown', {key: 'Escape'});
  await flush();
  assert.deepEqual(source.selection.value, [], 'Escape clears it');

  // a query that answers none of the picked rows is a new collection, with nothing selected
  fire(row('i1'), 'click');
  await flush();
  assert.deepEqual(source.selection.value.map((r) => r.id), ['i1']);
  source.query.value = 'title = "Naproxen"';
  await flush();
  assert.deepEqual(source.selection.value, [], 'and a collection holding none of them selects none');
  list.dispose();
  source.dispose();
});

scoped('a trash source is refused: deleted rows are read-only until they are restored', async () => {
  const memory = backend();
  const {source: live} = await table(memory);
  live.edit.value.markDeleted('i2');
  assert.equal(await live.session.save(), true);
  await flush();
  live.dispose();

  const {source} = await table(memory, {deleted: 'only'});
  assert.deepEqual(source.rows.items.value.map((r) => r.title), ['Ibuprofen']);
  assert.equal(await domains.bulkEdit(source), null);
  assert.match(balloon('warning'), /deleted rows are read-only/);
  assert.equal(picker(), undefined, 'no dialog, and nothing posted');
  assert.deepEqual(issues(memory).map((r) => r.version), [1, 2, 1]);
  source.dispose();
});

scoped('a backend that declares no writes offers no bulk edit, and says so by name', async () => {
  const memory = backend();
  // the method lives on the prototype: shadowed here, as a table whose declared support has no
  // writes — the member and the flag always answer together
  const issue = memory.tableSync('grit.issue');
  Object.defineProperty(issue, 'updateWhere', {value: undefined, configurable: true});
  issue.support = {...issue.support, writes: false};
  const {source} = await table(memory);
  assert.equal(await domains.bulkEdit(source), null);
  assert.match(balloon('error'), /does not support bulk edits/);
  source.dispose();
});

scoped('on a filtered list, the one row the user clicked is not the target: the filter is', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'number > 0'});
  select(source, (i) => i === 1);
  source.currentRow.value = source.selection.value[0];
  const {running} = await opened(source);
  assert.deepEqual(targets(), ['1 selected issue', 'All 3 issues matching the filter'], 'both are still offered');
  assert.equal(target().value.value, 'All 3 issues matching the filter');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);

  select(source, (i) => i !== 1);
  const again = await opened(source);
  assert.equal(target().value.value, '2 selected issues', 'a selection the user made is the target');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await again.running, null);
  source.dispose();
});

scoped('a required field wears the asterisk, and choosing it empty holds the button: it cannot be cleared', async () => {
  const memory = backend();
  const {source} = await table(memory);
  select(source, (i) => i === 1);
  const {running} = await opened(source);
  choose('title');
  await flush();
  assert.equal(field('title').root.classList.contains('u2-input-required'), true);
  field('title').value.value = '';
  await flush();
  assert.equal(hint('title'), 'Required: it cannot be cleared');
  assert.equal(okDisabled(), true);
  assert.match(await okReason(), /Title is required/);
  field('title').value.value = 'Kept';
  await flush();
  assert.equal(okDisabled(), false);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});
