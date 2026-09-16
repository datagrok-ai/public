/* `domains.list` (WO-10) over the memory backend: rows by id through the table's renderer (the
   handler where one claims the rows, the schema otherwise), Open and Delete gated by access per
   row, selection and the current row one thing, the next page near the bottom, the loading /
   empty / error area, and the `u2-domain-list` tag. `DG` comes from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {backends} from '../src/sources/backends.js';
import {DomainSource} from '../src/sources/domain-source.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {ROWS, backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {Rows} = await import('../src/sources/rows-like.js');
const {DomainList} = await import('../src/dg/domain/list.js');
const {allowedActions} = await import('../src/components/actions/actions.js');
const {registerDomainComponents} = await import('../src/dg/domain/registrations.js');
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
      DG.ObjectHandler.registered.length = 0;
      DG.DomainObjectHandler.opened.length = 0;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

async function issues(options = {}) {
  if (!(backends.domain instanceof MemoryDomainBackend))
    backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10, ...options});
  await flush();
  return {table, src};
}

/** Mounted with a viewport, so the virtual list renders rows. */
function mount(list, height = 300) {
  document.body.append(list.root);
  list.list.root.clientHeight = height;
  fire(list.list.root, 'scroll');
  return [...list.root.querySelectorAll('.u2-list-row')];
}

const names = (list, options) => list.actionsFor(options).map((a) => a.name);

scoped('rows render by id; the schema names them without a handler, the handler where one claims them', async () => {
  const {src} = await issues();
  const list = domains.list(src);
  assert.equal(list.root.dataset.u2, 'domain-list');
  const rows = mount(list);
  assert.deepEqual(rows.map((r) => r.dataset.u2Row), ['i1', 'i2', 'i3']);
  assert.deepEqual(rows.map((r) => r.querySelector('.u2-domain-list-content').textContent),
    ['Aspirin', 'Ibuprofen', 'Naproxen']);
  assert.equal(rows[0].querySelector('.u2-domain-list-name') !== null, true, 'the schema renderer');
  assert.equal(rows[0].querySelectorAll('.u2-row-actions button').length, 2, 'Open and Delete');
  list.dispose();
  src.dispose();
});

scoped('a handler that claims the rows renders them, resolved once per handle', async () => {
  DG.ObjectHandler.register(new DG.DomainObjectHandler('grit.issue'));
  const {table, src} = await issues();
  const handled = domains.list(src);
  const items = mount(handled);
  assert.deepEqual(items.map((r) => r.querySelector('.test-handler-item')?.textContent),
    ['Aspirin', 'Ibuprofen', 'Naproxen'], 'the handler\'s list item, built off the row\'s values');
  assert.equal(table.renderer.caption(src.rows.byKey('i1')), 'Aspirin', 'the name column first');
  handled.dispose();
  src.dispose();
});

scoped('cards: the schema card carries the name, up to three text fields and the creation time', async () => {
  const {src} = await issues();
  const list = domains.list(src, {mode: 'cards'});
  const rows = mount(list);
  assert.equal(list.root.classList.contains('u2-domain-list-cards'), true);
  assert.equal(rows[0].style.height, '64px', 'room for a title, a description and the time');
  const card = rows[2].querySelector('.u2-domain-card');
  assert.equal(card.querySelector('.u2-domain-card-title').textContent, 'Naproxen');
  assert.equal(card.querySelector('.u2-domain-card-description').textContent, 'high',
    'one muted line: the description column, else the first text field with a value');
  assert.equal(card.querySelector('.u2-domain-card-time') !== null, true);
  const draft = src.newRow({title: ''}, {pristine: true});
  await flush();
  const draftCard = [...list.root.querySelectorAll('.u2-list-row')].find((r) => r.dataset.u2Row === draft.id);
  assert.equal(draftCard.querySelector('.u2-domain-card-title').textContent, 'New issue');
  assert.equal(draftCard.querySelector('.u2-domain-card-title').classList.contains('u2-domain-draft'), true);
  list.dispose();
  const custom = domains.list(src, {render: (row) => Object.assign(document.createElement('b'), {textContent: row.id})});
  assert.deepEqual(mount(custom).map((r) => r.querySelector('b').textContent).slice(0, 3), ['i1', 'i2', 'i3']);
  custom.dispose();
  src.dispose();
});

scoped('actions: Open on saved rows through the handler, Delete under the row\'s own access, the table\'s own', async () => {
  const rows = ROWS.issue.map((r) => ({...r}));
  rows[1]['~can_delete'] = false;
  backends.domain = backend({rows: {issue: rows, project: ROWS.project.map((r) => ({...r}))}});
  const {table, src} = await issues();
  table.actions.add({name: 'Escalate', icon: 'arrow-up', requires: 'edit', when: (r) => r.title !== 'Naproxen',
    run: (r) => r.title = `${r.title}!`});
  const list = domains.list(src, {actions: [{name: 'Copy id', run: () => {}}]});
  assert.deepEqual(names(list, src.rows.byKey('i1')), ['Open', 'Delete', 'Escalate', 'Copy id']);
  assert.deepEqual(names(list, src.rows.byKey('i2')), ['Open', 'Escalate', 'Copy id'], 'no ~can_delete');
  assert.deepEqual(names(list, src.rows.byKey('i3')), ['Open', 'Delete', 'Copy id'], 'when: false');
  const draft = src.newRow({title: 'Draft'});
  assert.deepEqual(names(list, draft), ['Delete', 'Escalate', 'Copy id'], 'a draft has no address to open');
  list.actionsFor(src.rows.byKey('i1')).find((a) => a.name === 'Escalate').run();
  assert.equal(src.rows.byKey('i1').title, 'Aspirin!', 'bound to the row');
  list.actionsFor(src.rows.byKey('i1')).find((a) => a.name === 'Open').run();
  assert.equal(DG.DomainObjectHandler.opened.length, 1);
  assert.equal(DG.DomainObjectHandler.opened[0].values.title, 'Aspirin!');
  list.actionsFor(src.rows.byKey('i1')).find((a) => a.name === 'Delete').run();
  await flush();
  assert.deepEqual(src.rows.items.value.map((r) => r.title), ['Aspirin!', 'Ibuprofen', 'Naproxen', 'Draft'],
    'a row marked deleted stays, struck through, until the save');
  assert.equal(src.isDirty.value, true);
  assert.deepEqual(names(list, src.rows.byKey('i1')), ['Restore'], 'the one action on it');
  const rendered = mount(list);
  assert.equal(rendered[0].classList.contains('u2-domain-list-deleted'), true);
  assert.equal(rendered[1].classList.contains('u2-domain-list-deleted'), false);
  assert.equal(src.summary.value, '2 unsaved changes', 'the deletion and the draft added above');
  list.actionsFor(src.rows.byKey('i1')).find((a) => a.name === 'Restore').run();
  await flush();
  assert.equal(src.summary.value, '2 unsaved changes', 'the Escalate edit is back, and the draft');
  assert.deepEqual(names(list, src.rows.byKey('i1')), ['Open', 'Delete', 'Escalate', 'Copy id']);
  list.dispose();
  src.dispose();

  const bare = new DomainSource({table: 'grit.issue'});
  bare.start();
  await flush();
  const plain = domains.list(bare);
  assert.deepEqual(names(plain, bare.rows.byKey('i1')), ['Delete'], 'no handle: nothing to open with');
  plain.dispose();
  bare.dispose();
});

scoped('permission ⇒ hidden: without the delete capability Delete is gone everywhere', async () => {
  backends.domain = backend({access: {can: {view: true, insert: true, edit: true, delete: false, share: false},
    fields: {title: 'editable'}}});
  const {src} = await issues();
  const list = domains.list(src);
  assert.deepEqual(names(list, src.rows.byKey('i1')), ['Open']);
  const rows = mount(list);
  assert.equal(rows[0].querySelectorAll('.u2-row-actions button').length, 1);
  list.dispose();
  src.dispose();
});

scoped('selection and the current row are one thing, in both directions', async () => {
  const {src} = await issues();
  const list = domains.list(src);
  mount(list);
  assert.equal(list.list.selectedIndex.value, 0, 'a loaded list starts on its first row');
  assert.equal(src.currentRow.value.id, 'i1');
  list.list.selectedIndex.value = 1;
  await flush();
  assert.equal(src.currentRow.value.id, 'i2');
  src.currentRow.value = src.rows.byKey('i3');
  await flush();
  assert.equal(list.list.selectedIndex.value, 2);
  src.currentRow.value = null;
  await flush();
  assert.equal(list.list.selectedIndex.value, -1);
  const rows = [...list.root.querySelectorAll('.u2-list-row')];
  fire(rows[0], 'click');
  await flush();
  assert.equal(src.currentRow.value.id, 'i1');
  list.dispose();
  src.dispose();
});

scoped('near the bottom the next page is loaded into the same collection', async () => {
  const {src} = await issues({pageSize: 2});
  assert.equal(src.rows.items.value.length, 2);
  const list = domains.list(src);
  mount(list, 100);
  const root = list.list.root;
  root.scrollHeight = 56;
  root.scrollTop = 0;
  fire(root, 'scroll');
  await flush();
  assert.equal(src.rows.items.value.length, 3);
  assert.equal(src.summary.value, '3 issues');
  list.dispose();
  src.dispose();
});

scoped('the status area: loading, then empty or the failure with a retry', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({query: 'title = "zzz"'});
  const list = domains.list(src);
  const status = list.root.querySelector('[data-u2-part="status"]');
  assert.equal(status.querySelector('.u2-loader') !== null, true, 'loading');
  await flush();
  assert.equal(status.textContent, 'No issues match the filter.', 'an empty result says why it is empty');
  list.dispose();
  src.dispose();

  const broken = new DomainSource({table: 'grit.nope'});
  broken.start();
  const failing = domains.list(broken, {empty: 'Nothing here'});
  await flush();
  const area = failing.root.querySelector('[data-u2-part="status"]');
  assert.match(area.textContent, /Unknown table/);
  const retry = area.querySelector('button');
  assert.equal(retry.textContent, 'Retry');
  backends.domain = backend();
  fire(retry, 'click');
  await flush();
  assert.match(area.textContent, /Unknown table/, 'the table handle is resolved once per source');
  failing.dispose();
  broken.dispose();
});

scoped('spec: u2-domain-list over a bound source, its rows selectable', async () => {
  const {src} = await issues();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-domain-list', name: 'list',
    bind: {source: '$.issues'}, props: {mode: 'cards', itemHeight: 40}}},
  new SpecContext({data: {issues: signal(src)}}), reg);
  await flush();
  const list = instance.node('list');
  assert.equal(list instanceof DomainList, true);
  assert.equal(list.mode, 'cards');
  document.body.append(instance.root);
  list.list.root.clientHeight = 200;
  fire(list.list.root, 'scroll');
  const rows = [...instance.root.querySelectorAll('.u2-list-row')];
  assert.equal(rows.length, 3);
  assert.equal(rows[0].style.height, '40px');
  fire(rows[1], 'click');
  await flush();
  assert.equal(src.currentRow.value.id, 'i2');
  instance.dispose();
  src.dispose();
});

scoped('a list over a table starts on its first row; a draft source does not', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source();
  const list = domains.list(src);
  await flush();
  assert.equal(src.currentRow.value?.id, 'i1', 'the first row is current once loaded');
  assert.equal(list.list.selectedIndex.value, 0);
  list.dispose();
  src.dispose();
  const draft = table.draft({title: 'D'});
  const over = domains.list(draft);
  await flush();
  assert.equal(draft.currentRow.value.id.startsWith('~'), true, 'the draft stays current');
  over.dispose();
  draft.dispose();
});

scoped('keyboard: Enter hands the row to the form through the source, Delete deletes, actions rove with the selection', async () => {
  DG.ObjectHandler.register(new DG.DomainObjectHandler('grit.issue'));
  const {src} = await issues();
  const list = domains.list(src);
  const rows = mount(list);
  const buttons = (row) => [...row.querySelectorAll('.u2-row-actions button')].map((b) => b.tabIndex);
  assert.deepEqual(buttons(rows[1]), [-1, -1], 'an unselected row: no action in the tab order');
  assert.deepEqual(buttons(rows[0]), [0, 0], 'the selected (first) row\'s actions are tabbable');
  list.list.selectedIndex.value = 1;
  await flush();
  assert.deepEqual(buttons(rows[1]), [0, 0], 'the selection roves');
  assert.deepEqual(buttons(rows[0]), [-1, -1]);
  assert.equal(src.activate.value, 0);
  fire(list.list.root, 'keydown', {key: 'Enter'});
  assert.equal(src.activate.value, 1, 'Enter bumps the source: the paired form takes the focus');
  assert.equal(src.currentRow.value.id, 'i2');
  fire(list.list.root, 'keydown', {key: 'Delete'});
  await flush();
  assert.equal(src.rows.byKey('i2')['~state'], 'deleted');
  list.dispose();
  src.dispose();
});

scoped('a list over a source with a current row keeps it', async () => {
  const {src} = await issues();
  src.currentRow.value = src.rows.byKey('i2');
  const list = domains.list(src);
  await flush();
  assert.equal(src.currentRow.value?.id, 'i2', 'construction did not clear it');
  assert.equal(list.list.selectedIndex.value, 1, 'and the selection shows it');
  list.dispose();
  src.dispose();
});

scoped('allowedActions and Rows.isDraft', () => {
  const access = {can: (c) => c !== 'delete', row() { return this; }};
  assert.deepEqual(allowedActions([{name: 'a', run() {}}, {name: 'b', requires: 'delete', run() {}},
    {name: 'c', requires: 'edit', run() {}}], {access, row: {}}).map((a) => a.name), ['a', 'c']);
  assert.equal(Rows.isDraft({id: '~new:3'}), true);
  assert.equal(Rows.isDraft({id: 'i1'}), false);
});

scoped('cards: the recipe card unless a registered handler subclass paints its own', async () => {
  DG.ObjectHandler.register(new DG.DomainObjectHandler('grit.issue'));
  const {src} = await issues();
  const plain = domains.list(src, {mode: 'cards'});
  const rows = mount(plain);
  assert.equal(rows[0].querySelector('.u2-domain-card-title')?.textContent, 'Aspirin',
    'the reflective default paints the generic table: the recipe card instead');
  plain.dispose();
  src.dispose();
  DG.ObjectHandler.registered.length = 0;

  class IssueHandler extends DG.DomainObjectHandler {
    renderCard(x) {
      const el = document.createElement('div');
      el.className = 'test-issue-card';
      el.textContent = `card: ${x.values.title}`;
      return el;
    }
  }
  DG.ObjectHandler.register(new IssueHandler('grit.issue'));
  const {src: own} = await issues();
  const custom = domains.list(own, {mode: 'cards'});
  const cards = mount(custom);
  assert.equal(cards[0].querySelector('.test-issue-card')?.textContent, 'card: Aspirin', 'a subclass\'s own card');
  custom.dispose();
  own.dispose();
});

scoped('a bare table: the searchable columns beside the name, and the row a refusal names is marked', async () => {
  const {table, src} = await issues();
  src.rows.byKey('i1').description = 'Acetylsalicylic acid';
  await flush();
  const list = domains.list(src);
  const rows = mount(list);
  assert.equal(rows[0].querySelector('.u2-domain-list-name').textContent, 'Aspirin');
  assert.equal(rows[0].querySelector('.u2-domain-list-details').textContent, 'Acetylsalicylic acid',
    'what a search would match the row by');
  assert.equal(rows[1].querySelector('.u2-domain-list-details'), null, 'nothing to add for a row without one');

  const off = table.validators.add('priority',
    (value, row) => value === 'high' && !row.reporter ? 'Assign before escalating' : null);
  src.rows.byKey('i2').priority = 'high';
  assert.equal(await src.save(), false);
  await flush();
  assert.equal(src.problemRow.value, 'i2');
  const marked = list.list.root.querySelector('[data-u2-row="i2"]');
  assert.equal(marked.classList.contains('u2-domain-list-invalid'), true);
  assert.equal(marked.getAttribute('aria-invalid'), 'true');
  assert.match(marked.title, /Assign before escalating/);
  src.discard();
  await flush();
  assert.equal(list.list.root.querySelector('[data-u2-row="i2"]').classList.contains('u2-domain-list-invalid'),
    false, 'a discard takes the mark back with the refusal');
  off();
  list.dispose();
  src.dispose();
});

scoped('a table whose app brought its own renderer says what that renderer says', async () => {
  const {table, src} = await issues();
  src.rows.byKey('i1').description = 'Acetylsalicylic acid';
  await flush();
  table.renderer = {...table.renderer, caption: (row) => `#${row.number}`};
  const list = domains.list(src);
  const rows = mount(list);
  assert.equal(rows[0].querySelector('.u2-domain-list-details'), null);
  list.dispose();
  src.dispose();
});

scoped('a double-click on a row activates it, the way Enter does; a single click only selects', async () => {
  const {src} = await issues();
  const list = domains.list(src);
  const rows = mount(list);
  fire(rows[1].querySelector('.u2-domain-list-content'), 'click');
  await flush();
  assert.equal(src.currentRow.value.id, 'i2', 'a click selects');
  assert.equal(src.activate.value, 0, 'and nothing more');
  fire(rows[1].querySelector('.u2-domain-list-content'), 'dblclick');
  assert.equal(src.activate.value, 1, 'a double-click hands the row over, as Enter does');
  assert.equal(src.currentRow.value.id, 'i2');
  list.dispose();
  src.dispose();
});

scoped('no class inside a row is also on the list root: a row rule reaching the host collapses it', async () => {
  const {src} = await issues();
  src.rows.byKey('i1').description = 'Acetylsalicylic acid';
  await flush();
  const classes = (el, out = []) => {
    out.push(...el.className.split(' ').filter((c) => c !== ''));
    for (const child of el.children)
      classes(child, out);
    return out;
  };
  for (const mode of ['brief', 'cards']) {
    const list = domains.list(src, {mode});
    const rows = mount(list);
    const host = new Set(list.root.className.split(' '));
    assert.equal(host.has(`u2-domain-list-${mode}`), true, 'the mode is a class on the host');
    for (const cls of classes(rows[0]))
      assert.equal(host.has(cls), false, `${cls} is on both a ${mode} row and the list root`);
    list.dispose();
  }
  const brief = domains.list(src);
  const row = mount(brief)[0];
  assert.equal(row.querySelector('.u2-domain-list-line .u2-domain-list-details').textContent,
    'Acetylsalicylic acid');
  brief.dispose();
  src.dispose();
});

scoped('an empty result says what it is empty OF, and a search offers the way out of it', async () => {
  const {src} = await issues({query: 'number > 100'});
  const list = domains.list(src);
  mount(list);
  await flush();
  const status = list.root.querySelector('[data-u2-part="status"]');
  assert.equal(list.root.classList.contains('u2-domain-list-blank'), true,
    'the message takes the rows own area and the empty scroller goes');
  assert.match(status.textContent, /No issues match the filter\./);
  assert.equal(status.querySelector('button'), null, 'nothing to clear: the filter is not the search box');

  src.query.value = '';
  src.search.value = 'nonesuch';
  await flush();
  assert.match(status.textContent, /No issues match "nonesuch"\./);
  fire(status.querySelector('button'), 'click');
  await flush();
  assert.equal(src.search.value, '', 'Clear search brings the rows back');
  assert.equal(list.root.classList.contains('u2-domain-list-blank'), false);
  list.dispose();
  src.dispose();
});

scoped('a trash row says WHEN it was deleted, and Restore says so', async () => {
  const {src: live} = await issues();
  live.edit.value.markDeleted('i2');
  assert.equal(await live.session.save(), true);
  await flush();
  live.dispose();
  const {src} = await issues({deleted: 'only'});
  const list = domains.list(src);
  mount(list);
  await flush();
  const row = list.root.querySelector('[data-u2-row="i2"]');
  const detail = row.querySelector('.u2-domain-list-details');
  assert.notEqual(detail, null, 'a trash row carries a detail line');
  assert.equal(detail.classList.contains('u2-timestamp'), true,
    'when it was deleted, read in the reader own zone — never the wire string');
  assert.equal(detail.textContent.includes('GMT'), false);
  assert.match(detail.title, /2026/, 'the full moment is the title');
  list.actionsFor(src.rows.byKey('i2'))[0].run();
  await flush();
  assert.match([...document.body.querySelectorAll('.u2-notify-info')].map((e) => e.textContent).join('|'),
    /Restored "Ibuprofen"/, 'a restore says what came back');
  list.dispose();
  src.dispose();
});
