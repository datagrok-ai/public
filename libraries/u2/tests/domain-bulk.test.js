/* `domains.bulkEdit` (WO 3-8) on the memory backend: only targets that can be written are offered
   (never "all N matching" while a search the server cannot take is set), whatever blocks OK is a
   line above the buttons with OK disabled, the include checkbox — outside the input it enables —
   decides what is sent, nothing in the dialog is "required" (a checked field left empty CLEARS
   the column), the unsaved gate runs before it opens, and the report speaks the table's words. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {backends} from '../src/sources/backends.js';
import {notify} from '../src/components/display/notify.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {SCHEMA, backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');

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
const includeBox = (name) => document.querySelector(`[data-u2-include="${name}"]`);
const field = (name) => Control.forElement(includeBox(name).parentElement.querySelector('.u2-input-root'));
const balloon = (kind) => document.body.querySelector(`.u2-notify-${kind}`)?.textContent ?? '';
const issues = (memory) => memory.tableSync('grit.issue').rows;

/** A real click, as a user's is: the box activates itself and the input follows. */
function check(name) {
  includeBox(name).click();
}

/** The target picker: its value is `selected` or `matching`. */
const target = () => Control.forElement(document.querySelector('[data-u2-name="target"]'));
const targets = () => [...target().root.querySelectorAll('option')].map((o) => o.textContent);
/** What blocks OK, said in the dialog rather than across the screen. */
const reason = () => document.body.querySelector('.u2-domain-bulk-reason')?.textContent ?? '';
const okDisabled = () => buttonNamed('OK').disabled === true;

/** The dialog up and rendered; the promise is handed back WRAPPED — an async function
 * returning it would chain on it and never come back until the dialog is answered. */
async function opened(source) {
  const running = domains.bulkEdit(source);
  await flush();
  return {running};
}

async function table(memory, options = {}) {
  backends.domain = memory;
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

scoped('the checked columns are written to the selection; an unchecked one is left alone', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'number > 0'});
  select(source, (i) => i !== 1);
  assert.deepEqual(source.selection.value.map((r) => r.id), ['i1', 'i3']);

  const {running} = await opened(source);
  assert.equal(target().value.value, 'selected', 'a selection is the default target');
  assert.deepEqual(targets(), ['2 selected', 'all 3 matching'], 'no empty option in front of the targets');
  assert.equal(includeBox('id'), null, 'a system column is not offered');
  // the box is OUTSIDE the input's root: inside a disabled one it would be dead (BULK-CHECKBOX)
  assert.equal(includeBox('done').parentElement.classList.contains('u2-domain-bulk-row'), true);
  assert.equal(field('done').enabled, false, 'an unchecked column has no editor');
  assert.equal(okDisabled(), true);
  assert.match(reason(), /Check the fields to write/);
  assert.equal(field('title').root.querySelector('.u2-input-error')?.textContent ?? '', '',
    'nothing in a bulk dialog is required: an empty checked field CLEARS the column');
  check('done');
  await flush();
  assert.equal(field('done').enabled, true, 'a real click on the box enables the input');
  assert.equal(okDisabled(), false);
  assert.equal(reason(), '');
  field('done').value.value = true;
  field('title').value.value = 'NOT WRITTEN';
  fire(buttonNamed('OK'), 'click');
  assert.equal(await running, 2);
  assert.match(balloon('info'), /Updated 2 issues/, 'the table\'s words, never "rows"');

  const rows = issues(memory);
  assert.deepEqual(rows.map((r) => r.done), [true, false, true], 'the two selected rows only');
  assert.deepEqual(rows.map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen'],
    'the unchecked column was not in the request');
  source.dispose();
});

scoped('nothing checked keeps OK disabled and says so inline; CANCEL resolves null', async () => {
  const memory = backend();
  const {source} = await table(memory);
  select(source, () => true);
  const {running} = await opened(source);
  assert.equal(okDisabled(), true);
  assert.match(reason(), /Check the fields to write/);
  assert.equal(balloon('warning'), '', 'a refusal is not a balloon on the far edge of the screen');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  assert.deepEqual(issues(memory).map((r) => r.version), [1, 1, 1], 'nothing was written');
  source.dispose();
});

scoped('a selection past the cap disables OK and says how many are selected', async () => {
  const memory = many(1001);
  const {source} = await table(memory, {pageSize: 1100});
  select(source, () => true);
  assert.equal(source.selection.value.length, 1001);
  const {running} = await opened(source);
  check('done');
  await flush();
  field('done').value.value = true;
  assert.equal(okDisabled(), true);
  assert.match(reason(), /at most 1000 issues at a time; 1001 are selected/);
  assert.equal(issues(memory).filter((r) => r.done === true).length, 0, 'nothing was written');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('"all matching" writes the list\'s filter, capped, and reports that more match', async () => {
  const memory = many(1001);
  const {source} = await table(memory, {query: 'done = false'});
  const {running} = await opened(source);
  assert.deepEqual(targets(), ['all 1001 matching'], 'nothing selected: the filter is the only target');
  check('done');
  await flush();
  field('done').value.value = true;
  fire(buttonNamed('OK'), 'click');
  assert.equal(await running, 1000);
  assert.match(balloon('info'), /Updated 1000 issues — more match/);
  assert.equal(issues(memory).filter((r) => r.done === true).length, 1000);
  source.dispose();
});

scoped('no target the server would take: the dialog never opens, and says the way out', async () => {
  const memory = backend();
  const {source} = await table(memory);
  assert.equal(await domains.bulkEdit(source), null, 'nothing selected and nothing filtered');
  assert.match(balloon('warning'), /Select the issues to write, or filter the list/);
  assert.equal(includeBox('done'), null, 'no dialog promising a target it would then refuse');

  notify.closeAll();
  source.query.value = 'done = false';
  source.search.value = 'asp';
  await flush();
  assert.equal(await domains.bulkEdit(source), null);
  assert.match(balloon('warning'), /clear the search box/,
    'a search the update endpoint cannot take is said BEFORE the dialog, not after OK');
  source.dispose();
});

scoped('a search with a selection still offers the selection alone', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'done = false'});
  source.search.value = 'ibu';
  await flush();
  select(source, () => true);
  const {running} = await opened(source);
  assert.deepEqual(targets(), ['1 selected'], 'the filter target is not offered under a search');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null);
  source.dispose();
});

scoped('the unsaved gate runs before the dialog opens', async () => {
  const memory = backend();
  const {source} = await table(memory, {query: 'done = false'});
  source.edit.value.setValue('i2', 'title', 'Edited');
  assert.equal(source.isDirty.value, true);
  const running = domains.bulkEdit(source);
  await flush();
  assert.equal(includeBox('done'), null, 'the bulk dialog is not up yet');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await running, null, 'a cancelled gate cancels the whole thing');

  const second = domains.bulkEdit(source);
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(source.isDirty.value, false);
  assert.notEqual(includeBox('done'), null, 'and a discard lets it through');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await second, null);
  source.dispose();
});

scoped('the list\'s own selection is the collection\'s, and only a user makes one', async () => {
  const memory = backend();
  backends.domain = memory;
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
  assert.deepEqual(targets(), ['2 selected']);
  check('done');
  await flush();
  field('done').value.value = true;
  fire(buttonNamed('OK'), 'click');
  assert.equal(await running, 2);
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
  assert.equal(includeBox('done'), null, 'no dialog, and nothing posted');
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
