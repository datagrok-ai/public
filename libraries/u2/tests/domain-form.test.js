/* `domains.form` (WO-9) and the Save/Discard buttons over the memory backend: fields follow the
   access, an edit writes through the source, the form follows the current row, a create form's
   draft, the table's validators, the server's column errors on the fields, the version-conflict
   dialog both ways, the reference editors, and the `u2-domain-form` tag. `DG` and `grok` come
   from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {backends} from '../src/sources/backends.js';
import {DomainSource} from '../src/sources/domain-source.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {notify} from '../src/components/display/notify.js';
import {backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {DomainForm} = await import('../src/dg/domain/form.js');
const {Rows} = await import('../src/sources/rows-like.js');
const {registerDomainComponents} = await import('../src/dg/domain/registrations.js');
const grok = await import('datagrok-api/grok');

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

const CAN = {view: true, insert: true, edit: true, delete: true, share: true};

/** A started source over the issues, through the handle, with its first page loaded — over the
 * memory backend the test installed, or a fresh one (importing dg installs the platform one). */
async function issues(options = {}) {
  if (!(backends.domain instanceof MemoryDomainBackend))
    backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10, ...options});
  await flush();
  return {table, src};
}

scoped('fields follow the access: hidden absent, readonly text, editable an input', async () => {
  backends.domain = backend({access: {can: CAN,
    fields: {id: 'readonly', title: 'editable', number: 'readonly', project_id: 'editable'}}});
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  assert.equal(form.root.dataset.u2, 'domain-form');
  assert.deepEqual(form.form.getWidgetStatus().inputs.map((f) => [f.name, f.access]),
    [['project_id', 'editable'], ['number', 'readonly'], ['title', 'editable']],
    'the system columns are not fields');
  assert.equal(form.input('number'), undefined, 'readonly is text, not an input');
  assert.equal(form.root.querySelector('[data-u2-name="number"] [data-u2-part="readonly-value"]').textContent, '1');
  assert.equal(form.input('title').value.value, 'Aspirin');
  assert.equal(form.form.root.dataset.u2Row, 'i1');
  form.dispose();
  src.dispose();
});

scoped('a row the caller may not edit renders as text; a draft is written under insert', async () => {
  const rows = backend().tableSync('grit.issue').rows.map((r) => ({...r}));
  rows[1]['~can_edit'] = false;
  backends.domain = backend({rows: {issue: rows, project: backend().tableSync('grit.project').rows},
    access: {can: {...CAN, edit: false}, fields: {title: 'editable', number: 'editable'}}});
  const {table, src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  assert.deepEqual(form.form.getWidgetStatus().inputs.map((f) => f.access), ['readonly', 'readonly'],
    'no table-wide edit and no row saying otherwise: text');
  src.currentRow.value = src.rows.byKey('i2');
  await flush();
  assert.deepEqual(form.form.getWidgetStatus().inputs.map((f) => f.access), ['readonly', 'readonly'],
    'a row narrowing edit to false');
  form.dispose();
  src.dispose();
  const draft = table.draft({title: 'New'});
  await flush();
  const create = domains.form(draft);
  await flush();
  assert.deepEqual(create.form.getWidgetStatus().inputs.map((f) => f.access), ['editable', 'editable'],
    'insert is granted: the draft is editable');
  create.dispose();
  draft.dispose();
});

scoped('an edit writes through the source; discard restores the field', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  form.input('title').value.value = 'Aspirin 100';
  await flush();
  assert.equal(src.isDirty.value, true);
  assert.equal(src.changeCount.value, 1);
  assert.equal(src.rows.byKey('i1').title, 'Aspirin 100');
  assert.equal(src.edit.value.isChanged('i1', 'title'), true);
  src.discard();
  await flush();
  assert.equal(src.isDirty.value, false);
  assert.equal(form.input('title').value.value, 'Aspirin', 'the form re-read the row');
  form.dispose();
  src.dispose();
});

scoped('the form follows the current row; without one it says so', async () => {
  const {src} = await issues();
  const form = domains.form(src.currentRow);
  await flush();
  assert.equal(form.form, null);
  assert.equal(form.root.querySelector('[data-u2-part="empty"]').textContent, 'Select an issue to edit.');
  src.currentRow.value = src.rows.byKey('i1');
  await flush();
  const first = form.form;
  assert.equal(form.input('title').value.value, 'Aspirin');
  src.currentRow.value = src.rows.byKey('i2');
  await flush();
  assert.notEqual(form.form, first, 'a new row is a new form');
  assert.equal(form.input('title').value.value, 'Ibuprofen');
  assert.equal(form.form.root.dataset.u2Row, 'i2');
  src.currentRow.value = null;
  await flush();
  assert.equal(form.form, null);
  assert.equal(form.validity.value, null);
  form.dispose();
  src.dispose();
});

scoped('a bare row or a foreign signal needs its source', async () => {
  backends.domain = backend();
  const src = new DomainSource({table: 'grit.issue'});
  src.start();
  await flush();
  assert.throws(() => domains.form(src.currentRow), /pass `source`/);
  const row = src.rows.byKey('i3');
  assert.throws(() => domains.form(row), /pass `source`/);
  const form = domains.form(row, {source: src});
  await flush();
  assert.equal(form.input('title').value.value, 'Naproxen');
  const bound = domains.form(src.currentRow, {source: src});
  assert.equal(bound.source, src);
  bound.dispose();
  form.dispose();
  src.dispose();
});

scoped('a create form: the draft is pristine until touched, save inserts it', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const draft = table.draft({project_id: 'p1', title: 'Draft'});
  await flush();
  assert.equal(draft.state.value, 'ready');
  assert.equal(draft.isDraft, true);
  assert.equal(draft.rows.items.value.length, 1, 'the draft alone: nothing is loaded');
  assert.notEqual(draft.currentRow.value, null);
  assert.equal(draft.isDirty.value, false, 'pristine');
  const form = domains.form(draft);
  await flush();
  assert.equal(form.input('title').value.value, 'Draft');
  form.input('number').value.value = 4;
  await flush();
  assert.equal(draft.isDirty.value, true);
  assert.equal(await draft.save(), true);
  const store = backends.domain.tableSync('grit.issue').rows;
  assert.equal(store.length, 4);
  assert.equal(store[3].title, 'Draft');
  assert.equal(store[3].number, 4);
  assert.equal(store[3].project_id, 'p1');
  form.dispose();
  draft.dispose();
});

scoped('the table\'s validators reach the field, after the form\'s own rules', async () => {
  const {table, src} = await issues();
  table.validators.add('title', (v) => String(v ?? '').length < 5 ? 'At least 5 characters' : null);
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  assert.equal(form.validity.value, null);
  form.input('title').value.value = 'Abc';
  await flush();
  assert.equal(form.input('title').validity.value, 'At least 5 characters');
  assert.equal(form.validity.value, 'At least 5 characters');
  assert.equal(form.validate(), false);
  form.input('title').value.value = '';
  await flush();
  assert.equal(form.input('title').validity.value, 'Value can\'t be empty', 'required first');
  form.input('title').value.value = 'Abcdef';
  await flush();
  assert.equal(form.validity.value, null);
  form.dispose();
  src.dispose();
});

scoped('a refused save is one balloon; a failed load is the hint\'s to show, not a balloon', async () => {
  const {src} = await issues();
  backends.domain.transaction = async () => {
    throw Object.assign(new Error('Row rejected'), {code: 'validation'});
  };
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  form.input('title').value.value = 'Ab';
  await flush();
  assert.equal(await src.save(), false);
  await flush();
  assert.match(document.body.querySelector('.u2-notify-error').textContent, /Row rejected/);
  assert.equal(form.input('title').value.value, 'Ab', 'the edit stays pending');
  assert.equal(src.isDirty.value, true);
  form.dispose();
  src.dispose();
  notify.closeAll();

  backends.domain = backend();
  const broken = new DomainSource({table: 'grit.nope'});
  const hint = domains.form(broken);
  broken.start();
  await flush();
  assert.match(hint.root.querySelector('[data-u2-part="empty"]').textContent, /Unknown table/);
  assert.equal(document.body.querySelector('.u2-notify-error'), null, 'no balloon for a load failure');
  hint.dispose();
  broken.dispose();
});

scoped('a 403 drops the access caches on top of the balloon', async () => {
  const {src} = await issues();
  backends.domain.transaction = async () => {
    throw Object.assign(new Error('Forbidden'), {code: 'forbidden'});
  };
  const form = domains.form(src);
  const before = grok.dapi.domains.invalidated;
  src.rows.byKey('i1').title = 'Ab';
  assert.equal(await src.save(), false);
  await flush();
  assert.equal(grok.dapi.domains.invalidated, before + 1, 'the form reports the refused save');
  assert.match(document.body.querySelector('.u2-notify-error').textContent, /Forbidden/);
  form.dispose();
  src.dispose();
});

scoped('a reference column is picked from its table, a user column through the user picker', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  const pick = form.input('project_id');
  assert.equal(pick.root.dataset.u2, 'domain-pick');
  assert.equal(pick.value.value, 'p1');
  await flush();
  assert.equal(pick.root.querySelector('input').value, 'Grit', 'the id resolved to the project\'s name');
  const reporter = form.input('reporter');
  assert.equal(reporter.root.dataset.u2, 'pick-input');
  assert.equal(reporter.value.value ?? '', '', 'the form seeds an empty string cell with ""');
  assert.equal(reporter.typeAhead.selected.value, null);
  assert.equal(src.isDirty.value, false, 'rendering the form writes nothing');
  form.dispose();
  src.dispose();
});

scoped('New: hidden without insert, disabled until loaded, the draft made current counts as a change', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source();
  const add = domains.newButton(src, {project_id: 'p1', title: 'Fresh'});
  assert.equal(add.root.dataset.u2, 'new-button');
  assert.equal(add.root.disabled, true, 'nothing to add to until the table is loaded');
  await flush();
  assert.equal(add.root.disabled, false);
  assert.equal(add.root.hidden, false);
  fire(add.root, 'click');
  await flush();
  const draft = src.currentRow.value;
  assert.equal(draft.title, 'Fresh');
  assert.equal(src.isDirty.value, true, 'the press is a change of its own: Save and the leave gate see it');
  assert.equal(src.summary.value, '1 unsaved change');
  assert.equal(src.rows.items.value.length, 4);
  add.dispose();
  src.dispose();

  backends.domain = backend({access: {can: {view: true, insert: false, edit: true, delete: true, share: true},
    fields: {title: 'editable'}}});
  const reader = (await domains.table('grit.issue')).source();
  const denied = domains.newButton(reader);
  await flush();
  assert.equal(denied.root.hidden, true, 'permission ⇒ hidden');
  denied.dispose();
  reader.dispose();
});

scoped('Save and Discard follow the state: disabled while clean, while saving, and after', async () => {
  const {src} = await issues();
  const save = domains.saveButton(src);
  const discard = domains.discardButton(src);
  await flush();
  assert.equal(save.root.dataset.u2, 'save-button');
  assert.equal(save.button.textContent, 'Save');
  assert.equal(save.button.disabled, true);
  assert.equal(discard.button.disabled, true);
  src.rows.byKey('i1').title = 'X';
  await flush();
  assert.equal(save.button.disabled, false);
  assert.equal(discard.button.disabled, false);
  const store = backends.domain.tableSync('grit.issue');
  const real = backends.domain.transaction.bind(backends.domain);
  let release;
  backends.domain.transaction = (ops) => new Promise((resolve) => release = () => resolve(real(ops)));
  fire(save.button, 'click');
  await flush();
  assert.equal(save.button.disabled, true, 'a save in flight');
  assert.equal(discard.button.disabled, true, 'locked together');
  release();
  await flush();
  await flush();
  assert.equal(src.isDirty.value, false);
  assert.equal(save.button.disabled, true);
  assert.equal(store.rows[0].title, 'X');
  src.rows.byKey('i1').title = 'Y';
  await flush();
  fire(discard.button, 'click');
  await flush();
  assert.equal(src.isDirty.value, false);
  assert.equal(src.rows.byKey('i1').title, 'X');
  save.dispose();
  discard.dispose();
  src.dispose();
});

scoped('spec: u2-domain-form over a bound source edits its current row; unbound it is a placeholder', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  assert.equal(reg.get('u2-domain-form').usage.length > 0, true);
  const ctx = new SpecContext({data: {issues: signal(src)}});
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-splitter', children: [
    {tag: 'u2-domain-form', name: 'form', bind: {source: '$.issues'}, props: {include: ['title', 'number']}},
    {tag: 'u2-domain-form', name: 'orphan'},
  ]}}, ctx, reg);
  await flush();
  const form = instance.node('form');
  assert.equal(form instanceof DomainForm, true);
  assert.deepEqual(form.form.properties.map((p) => p.name), ['title', 'number']);
  form.input('title').value.value = 'Spec';
  await flush();
  assert.equal(src.rows.byKey('i1').title, 'Spec');
  assert.match(instance.root.querySelector('.u2-spec-error').textContent, /bind "source"/);
  instance.dispose();
  src.dispose();
});

scoped('spec round trip: a form and a list bound through `$.issues.source` over a `u2-domain-source`', async () => {
  backends.domain = backend();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const spec = {
    $schema: 'dg-ui/1',
    components: [{tag: 'u2-domain-source', name: 'issues', props: {table: 'grit.issue', pageSize: 10}}],
    root: {tag: 'u2-splitter', name: 'page', children: [
      {tag: 'u2-domain-list', name: 'list', bind: {source: '$.issues.source'}, props: {mode: 'brief'}},
      {tag: 'u2-domain-form', name: 'form', bind: {source: '$.issues.source'}, props: {include: ['title', 'number']}},
    ]},
  };
  const instance = renderSpec(spec, new SpecContext(), reg);
  document.body.append(instance.root);
  await flush();
  assert.deepEqual(instance.dump(), spec);
  const src = instance.node('issues');
  const list = instance.node('list');
  const form = instance.node('form');
  assert.equal(list.source, src);
  assert.equal(form.source, src);
  assert.equal(instance.root.querySelector('.u2-spec-error'), null, 'both bound, nothing placeholdered');
  list.list.selectedIndex.value = 2;
  await flush();
  assert.equal(form.input('title').value.value, 'Naproxen', 'the list\'s selection is the form\'s row');
  form.input('title').value.value = 'Naproxen 500';
  await flush();
  assert.equal(src.rows.byKey('i3').title, 'Naproxen 500');
  assert.equal(src.isDirty.value, true);
  instance.dispose();
});

scoped('a readonly reference shows the target\'s name, not its uuid — the id until it resolves, and after a refresh', async () => {
  backends.domain = backend({access: {can: CAN, fields: {title: 'editable', project_id: 'readonly', reporter: 'readonly'}}});
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  const value = () => form.root.querySelector('[data-u2-name="project_id"] [data-u2-part="readonly-value"]').textContent;
  assert.equal(value(), 'p1', 'the id stands in until the name arrives');
  await flush();
  assert.equal(value(), 'Grit', 'the project\'s name column');
  form.input('title').value.value = 'Aspirin 100';
  await flush();
  assert.equal(value(), 'Grit', 'a refresh re-reads the row and the caption is put back');
  assert.equal(form.root.querySelector('[data-u2-name="reporter"] [data-u2-part="readonly-value"]').textContent, '',
    'an empty reference stays empty');
  src.currentRow.value = src.rows.byKey('i3');
  await flush();
  assert.equal(value(), 'Datagrok', 'a new row resolves its own reference');
  form.dispose();
  src.dispose();
});

scoped('system columns: out by default, a muted footer on request — captions, local times, never a draft\'s key', async () => {
  const {table, src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const plain = domains.form(src);
  await flush();
  assert.equal(plain.root.querySelector('[data-u2-part="system"]'), null);
  assert.equal(plain.form.properties.some((p) => ['id', 'version', 'created_on'].includes(p.name)), false);
  plain.dispose();
  const form = domains.form(src, {system: 'footer'});
  await flush();
  const footer = form.root.querySelector('[data-u2-part="system"]');
  const row = (name) => footer.querySelector(`[data-u2-name="${name}"]`);
  assert.deepEqual([...footer.querySelectorAll('.u2-input-label')].map((el) => el.textContent),
    ['Id', 'Version', 'Created', 'Updated', 'Author']);
  assert.equal(row('id').querySelector('[data-u2-part="readonly-value"]').textContent, 'i1');
  assert.equal(row('version').querySelector('[data-u2-part="readonly-value"]').textContent, '1');
  assert.equal(row('created_on').querySelector('.u2-timestamp') !== null, true, 'a local short date, the full one on hover');
  src.currentRow.value = null;
  form.dispose();
  const draft = table.draft({title: 'Draft'});
  await flush();
  const create = domains.form(draft, {system: 'footer'});
  await flush();
  assert.equal(create.root.querySelector('[data-u2-part="system"] [data-u2-name="id"] [data-u2-part="readonly-value"]').textContent,
    'assigned on save');
  assert.deepEqual([...create.root.querySelectorAll('[data-u2-part="system"] .u2-input-label')]
    .map((el) => el.textContent), ['Id'], 'a draft has no version, times or author yet');
  assert.equal(create.input('title').root.classList.contains('u2-input-required'), true,
    'a column that must hold a value is marked on its caption');
  assert.equal(create.input('description').root.classList.contains('u2-input-required'), false);
  create.dispose();
  draft.dispose();
  src.dispose();
});

scoped('a pristine draft is not red: verdicts show once a field is edited or left, or on validate()', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const draft = table.draft({project_id: 'p1'});
  await flush();
  const form = domains.form(draft);
  await flush();
  const title = form.input('title');
  assert.equal(title.validity.value, 'Value can\'t be empty', 'the verdict exists');
  assert.equal(title.root.classList.contains('u2-input-untouched'), true, 'but is not shown yet');
  // DOM doubles are cyclic: compare identity as a boolean, never hand them to assert's inspector
  assert.equal(document.activeElement === title.root.querySelector('input'), true,
    'the draft starts with its first empty required field focused');
  fire(title.root.querySelector('input'), 'focusout');
  assert.equal(title.root.classList.contains('u2-input-untouched'), false, 'leaving the field shows it');
  const number = form.input('number');
  assert.equal(number.root.classList.contains('u2-input-untouched'), true);
  number.value.value = 3;
  await flush();
  assert.equal(number.root.classList.contains('u2-input-untouched'), false, 'editing the field shows it');
  const due = form.input('due');
  assert.equal(due.root.classList.contains('u2-input-untouched'), true);
  assert.equal(form.validate(), false);
  assert.equal(due.root.classList.contains('u2-input-untouched'), false, 'validate() shows every verdict');
  form.dispose();
  draft.dispose();
});

scoped('an edited field carries the amber edge until the change is saved or discarded', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  const title = form.input('title');
  assert.equal(title.root.classList.contains('u2-input-changed'), false);
  title.value.value = 'Aspirin 100';
  await flush();
  assert.equal(title.root.classList.contains('u2-input-changed'), true);
  assert.equal(form.input('number').root.classList.contains('u2-input-changed'), false);
  src.discard();
  await flush();
  assert.equal(title.root.classList.contains('u2-input-changed'), false);
  form.dispose();
  src.dispose();
});

scoped('Save names what it saved and returns the focus to the form; Discard returns it too', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  const save = domains.saveButton(src);
  const discard = domains.discardButton(src);
  assert.equal(save.session, src.session, 'paired through the source: nothing is passed twice');
  await flush();
  form.input('title').value.value = 'Aspirin 100';
  await flush();
  document.body.append(save.root, discard.root, form.root);
  save.button.focus();
  fire(save.button, 'click');
  await flush();
  await flush();
  assert.equal(src.isDirty.value, false);
  assert.match(document.body.querySelector('.u2-notify-info')?.textContent ?? '', /Issue saved/);
  assert.equal(document.activeElement === form.input('project_id').root.querySelector('input'), true,
    'the focus is back on the first field');
  form.input('title').value.value = 'Zzz';
  await flush();
  discard.button.focus();
  fire(discard.button, 'click');
  await flush();
  assert.equal(src.isDirty.value, false);
  assert.equal(document.activeElement === form.input('project_id').root.querySelector('input'), true);
  save.button.focus();
  src.activate.value++;
  await flush();
  assert.equal(document.activeElement === form.input('project_id').root.querySelector('input'), true,
    'a list\'s Enter, through the source, focuses the form');
  save.dispose();
  discard.dispose();
  form.dispose();
  src.dispose();
});

scoped('the focus goes to the first empty required field after New and after Save, not to a pre-filled picker', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const draft = table.draft({project_id: 'p1'});
  await flush();
  const form = domains.form(draft);
  document.body.append(form.root);
  await flush();
  await flush();
  const focused = () => [...form.form.inputs].find((i) => i.root.contains(document.activeElement))?.name;
  assert.equal(focused(), 'title', 'the project is filled by the defaults: the title is what needs typing');
  form.input('title').value.value = 'Typed';
  await flush();
  assert.equal(await draft.save(), true);
  await flush();
  await flush();
  assert.equal(focused(), 'project_id', 'every required field filled: the first one');
  draft.newRow({project_id: 'p2'}, {pristine: true});
  await flush();
  await flush();
  assert.equal(focused(), 'title', 'New: past the pre-filled picker again');
  form.dispose();
  draft.dispose();
});

scoped('Save carries its shortcut as the tooltip and is the session\'s primary button; Tab from the last field lands on it',
  async () => {
    const {src} = await issues();
    src.currentRow.value = src.rows.byKey('i1');
    const form = domains.form(src);
    const save = domains.saveButton(src);
    const discard = domains.discardButton(src);
    document.body.append(form.root, save.root, discard.root);
    await flush();
    assert.equal(save.button.title, 'Save (Ctrl+S)');
    assert.equal(discard.button.title, 'Discard changes');
    assert.equal(src.session.primaryButton, save.button);
    const last = form.form.inputs.at(-1);
    const editor = last.root.querySelector('.u2-input-editor');
    const focusable = editor.matches('input, select, textarea, button') ? editor :
      editor.querySelector('input, select, textarea, button') ?? editor;
    focusable.focus();
    let prevented = !fire(focusable, 'keydown', {key: 'Tab', cancelable: true});
    assert.equal(prevented, false, 'nothing to save: Tab goes where the browser takes it');
    form.input('title').value.value = 'Aspirin 100';
    await flush();
    prevented = !fire(focusable, 'keydown', {key: 'Tab', cancelable: true});
    assert.equal(prevented, true);
    assert.equal(document.activeElement, save.button, 'Tab from the last field reaches Save');
    focusable.focus();
    prevented = !fire(focusable, 'keydown', {key: 'Tab', shiftKey: true, cancelable: true});
    assert.equal(prevented, false, 'Shift+Tab walks back through the fields');
    const first = form.form.inputs[0].root.querySelector('.u2-input-editor');
    prevented = !fire(first, 'keydown', {key: 'Tab', cancelable: true});
    assert.equal(prevented, false, 'only the last field hands over');
    save.dispose();
    assert.equal(src.session.primaryButton, undefined, 'released with the button');
    discard.dispose();
    form.dispose();
    src.dispose();
  });

scoped('a form built after the list bumped activate does not steal the focus; a later bump moves it', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  src.activate.value++;
  const other = document.createElement('button');
  document.body.append(other);
  other.focus();
  const form = domains.form(src);
  document.body.append(form.root);
  await flush();
  await flush();
  assert.equal(document.activeElement === other, true, 'the earlier bump is not a request to this form');
  src.activate.value++;
  await flush();
  assert.equal(document.activeElement === form.input('project_id').root.querySelector('input'), true);
  form.dispose();
  src.dispose();
});

scoped('a refused save names the field by its caption; the verdicts are shown', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const draft = table.draft({project_id: 'p1'});
  await flush();
  const form = domains.form(draft);
  const save = domains.saveButton(draft);
  await flush();
  draft.currentRow.value.number = 2;
  await flush();
  assert.equal(form.problem.value, 'Title is required');
  fire(save.button, 'click');
  await flush();
  await flush();
  assert.match(document.body.querySelector('.u2-notify-error')?.textContent ?? '', /Cannot save: Title is required/);
  assert.equal(draft.summary.value, 'Cannot save: Title is required', 'the guard\'s refusal is the source\'s error');
  assert.equal(form.input('title').root.classList.contains('u2-input-untouched'), false, 'the verdict is shown');
  assert.equal(backends.domain.tableSync('grit.issue').rows.length, 3, 'nothing was written');
  form.input('title').value.value = 'Titled';
  await flush();
  assert.equal(form.problem.value, null);
  form.dispose();
  assert.equal(await draft.save(), true, 'the guard left with the form');
  save.dispose();
  draft.dispose();
});

scoped('the keyboard path: Ctrl+S and Ctrl+Enter save a dirty, valid form; Esc does nothing destructive', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  assert.deepEqual(form.getWidgetStatus().shortcuts, {'Ctrl+S': 'Save', 'Ctrl+Enter': 'Save'});
  const store = backends.domain.tableSync('grit.issue');
  form.input('title').value.value = 'Aspirin 100';
  await flush();
  fire(form.root, 'keydown', {key: 'Escape'});
  await flush();
  assert.equal(src.isDirty.value, true, 'Esc discards nothing');
  const notPrevented = fire(form.root, 'keydown', {key: 's', ctrlKey: true, cancelable: true});
  assert.equal(notPrevented, false, 'the browser\'s Save page is not what was meant');
  await flush();
  await flush();
  assert.equal(src.isDirty.value, false);
  assert.equal(store.rows[0].title, 'Aspirin 100');
  form.input('title').value.value = 'Aspirin 200';
  await flush();
  fire(form.root, 'keydown', {key: 'Enter', ctrlKey: true, cancelable: true});
  await flush();
  await flush();
  assert.equal(store.rows[0].title, 'Aspirin 200');
  assert.match(document.body.querySelector('.u2-notify-info')?.textContent ?? '', /Issue saved/);
  form.dispose();
  src.dispose();
});

scoped('New with a function carries fields over from the row that was current', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i3');
  const add = domains.newButton(src, (last) => ({project_id: last?.project_id ?? 'p1', title: ''}));
  await flush();
  fire(add.root, 'click');
  await flush();
  assert.equal(src.currentRow.value.project_id, 'p2', 'Naproxen\'s project');
  add.dispose();
  src.dispose();
});

scoped('a schema constraint refuses Save with its message; fixing the value clears it; an empty cell passes', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  assert.equal(form.validity.value, null);
  form.input('weight').value.value = -1;
  await flush();
  assert.equal(form.input('weight').validity.value, 'Weight is positive', 'on the first column the expr names');
  assert.equal(form.problem.value, 'Weight: Weight is positive');
  assert.equal(await src.save(), false);
  assert.match(src.error.value.message, /Weight is positive/);
  assert.equal(src.isDirty.value, true);
  form.input('weight').value.value = 2;
  await flush();
  assert.equal(form.input('weight').validity.value, null);
  assert.equal(form.problem.value, null);
  form.input('weight').value.value = null;
  await flush();
  assert.equal(form.problem.value, null, 'NULL passes a CHECK');
  form.dispose();
  src.dispose();
});

scoped('a constraint naming a column the caller may not see is left to the server', async () => {
  backends.domain = backend({access: {can: CAN, fields: {title: 'editable', number: 'editable', project_id: 'editable'}}});
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  src.rows.byKey('i1').weight = -5;
  await flush();
  assert.equal(form.problem.value, null);
  assert.equal(form.validate(), true);
  form.dispose();
  src.dispose();
});

scoped('a ref column with a schema filter picks among the rows its siblings allow', async () => {
  backends.domain = new MemoryDomainBackend({name: 'grit', tables: {
    project: {friendlyName: 'Projects', columns: {name: {type: 'string', isName: true}}},
    issue: {friendlyName: 'Issues', columns: {
      project_id: {type: 'ref', ref: 'project', friendlyName: 'Project'},
      title: {type: 'string', isName: true},
      parent_id: {type: 'ref', ref: 'issue', filter: 'project_id = $project_id', friendlyName: 'Parent'},
    }},
  }}, {rows: {
    project: [{id: 'p1', name: 'Grit'}, {id: 'p2', name: 'Datagrok'}],
    issue: [{id: 'i1', project_id: 'p1', title: 'Aspirin'}, {id: 'i2', project_id: 'p1', title: 'Ibuprofen'},
      {id: 'i3', project_id: 'p2', title: 'Naproxen'}],
  }});
  const table = await domains.table('grit.issue');
  assert.deepEqual(table.info.refFilters, {parent_id: 'project_id = $project_id'});
  const src = table.source();
  await flush();
  src.currentRow.value = src.rows.byKey('i3');
  const form = domains.form(src, {overrides: {parent_id: {debounceMs: 0}}});
  document.body.append(form.root);
  await flush();
  const input = form.input('parent_id').root.querySelector('input');
  input.focus();
  await flush();
  await flush();
  const options = () => [...document.body.querySelectorAll('.u2-typeahead-option')].map((el) => el.textContent);
  assert.deepEqual(options(), ['Naproxen'], 'the issues of the row\'s project');
  input.blur();
  form.input('project_id').value.value = null;
  await flush();
  input.focus();
  await flush();
  await flush();
  assert.deepEqual(options(), []);
  assert.equal(document.body.querySelector('.u2-typeahead-empty')?.textContent, 'Pick a project first');
  form.dispose();
  src.dispose();
});

scoped('the system footer follows the row: a saved draft shows the number, version, times and author', async () => {
  backends.domain = backend();
  const find = grok.dapi.users.find;
  grok.dapi.users.find = async (id) => ({id, friendlyName: 'Ada Lovelace'});
  try {
    const table = await domains.table('grit.issue');
    const draft = table.draft({project_id: 'p1', title: 'Fresh'});
    await flush();
    const form = domains.form(draft, {system: 'footer'});
    await flush();
    const labels = () => [...form.root.querySelectorAll('[data-u2-part="system"] .u2-input-label')]
      .map((el) => el.textContent);
    assert.deepEqual(labels(), ['Id'], 'a draft has none of them yet');
    form.input('title').value.value = 'Fresh issue';
    await flush();
    assert.equal(await draft.save(), true);
    await flush();
    await flush();
    assert.equal(Rows.isDraft(draft.currentRow.value), false, 'the draft is a saved row now');
    assert.deepEqual(labels(), ['Id', 'Version', 'Created', 'Updated', 'Author'],
      'and the footer says so without the form being rebuilt');
    const value = (name) => form.root
      .querySelector(`[data-u2-part="system"] [data-u2-name="${name}"] [data-u2-part="readonly-value"]`).textContent;
    assert.notEqual(value('id'), 'assigned on save');
    // the insert answers `{id, version}` alone: the rest is the writer's post-save re-read, and it
    // has to be in the frame before the footer is built
    const row = draft.currentRow.value;
    assert.equal(row.version, 1);
    assert.equal(typeof row.created_on, 'string');
    assert.equal(row.author_id, 'me');
    assert.equal(value('version'), '1');
    const full = {month: 'short', day: 'numeric', year: 'numeric', hour: '2-digit', minute: '2-digit'};
    assert.equal(value('created_on'), new Date(row.created_on).toLocaleString(undefined, full),
      'the full date, the way a grid cell spells it — never a bare "Sep 13"');
    await flush();
    assert.equal(value('author_id'), 'Ada Lovelace', 'the author is a name, not the uuid the insert stamped');
    form.dispose();
    draft.dispose();
  } finally {
    grok.dapi.users.find = find;
  }
});

scoped('a create form arrives with the first field focused and no dropdown over it', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const draft = table.draft();
  await flush();
  const form = domains.form(draft, {system: 'footer'});
  await flush();
  await flush();
  const pick = form.input('project_id');
  assert.equal(document.activeElement, pick.typeAhead.root.querySelector('input'),
    'the first field still empty is where typing starts');
  assert.equal(pick.typeAhead.isOpen.value, false, 'a dropdown the user did not open is in the way');
  form.dispose();
  draft.dispose();
});

scoped('a save closes the picker its Ctrl+S came from, so nothing covers the form that comes back', async () => {
  const {table, src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  const pick = form.input('project_id');
  pick.typeAhead.root.querySelector('input').focus();
  await flush();
  assert.equal(pick.typeAhead.isOpen.value, true, 'a focused ref picker offers its rows');
  src.rows.byKey('i1').title = 'Aspirin 100';
  await flush();
  assert.equal(await src.session.save(), true);
  await flush();
  await flush();
  assert.equal(pick.typeAhead.isOpen.value, false);
  form.dispose();
  src.dispose();
});

scoped('a constraint over two columns is read with SQL null semantics: an empty column is not a free pass', async () => {
  backends.domain = new MemoryDomainBackend({name: 'grit', tables: {
    issue: {friendlyName: 'Issues',
      constraints: {span: {expr: 'low >= 0 and high >= 0', message: 'Both bounds must be positive'}},
      columns: {title: {type: 'string', isName: true}, low: {type: 'int'}, high: {type: 'int'}}},
  }}, {rows: {issue: [{id: 'i1', title: 'A', low: 1, high: 1}]}});
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i1');
  const form = domains.form(src);
  await flush();
  assert.equal(form.problem.value, null);
  form.input('high').value.value = null;
  await flush();
  assert.equal(form.problem.value, null, 'one bound empty, the other true: unknown, and admitted');
  form.input('low').value.value = -1;
  await flush();
  assert.equal(form.problem.value, 'Low: Both bounds must be positive', 'false and unknown is false, as Postgres reads it');
  form.dispose();
  src.dispose();
});
