/* `domainList` (WO-10) over the memory backend: rows by id through the table's renderer (the
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
const {domainList, DomainList} = await import('../src/dg/domain/list.js');
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
  const list = domainList(src);
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
  const handled = domainList(src);
  const items = mount(handled);
  assert.deepEqual(items.map((r) => r.querySelector('.test-handler-item')?.textContent),
    ['Aspirin', 'Ibuprofen', 'Naproxen'], 'the handler\'s list item, built off the row\'s values');
  assert.equal(table.renderer.caption(src.rows.byKey('i1')), 'Aspirin', 'the name column first');
  handled.dispose();
  src.dispose();
});

scoped('cards: the schema card carries the name, up to three text fields and the creation time', async () => {
  const {src} = await issues();
  const list = domainList(src, {mode: 'cards'});
  const rows = mount(list);
  assert.equal(list.root.classList.contains('u2-domain-list-cards'), true);
  assert.equal(rows[0].style.height, '96px', 'room for a handler\'s card');
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
  const custom = domainList(src, {render: (row) => Object.assign(document.createElement('b'), {textContent: row.id})});
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
  const list = domainList(src, {actions: [{name: 'Copy id', run: () => {}}]});
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
  const plain = domainList(bare);
  assert.deepEqual(names(plain, bare.rows.byKey('i1')), ['Delete'], 'no handle: nothing to open with');
  plain.dispose();
  bare.dispose();
});

scoped('permission ⇒ hidden: without the delete capability Delete is gone everywhere', async () => {
  backends.domain = backend({access: {can: {view: true, insert: true, edit: true, delete: false, share: false},
    fields: {title: 'editable'}}});
  const {src} = await issues();
  const list = domainList(src);
  assert.deepEqual(names(list, src.rows.byKey('i1')), ['Open']);
  const rows = mount(list);
  assert.equal(rows[0].querySelectorAll('.u2-row-actions button').length, 1);
  list.dispose();
  src.dispose();
});

scoped('selection and the current row are one thing, in both directions', async () => {
  const {src} = await issues();
  const list = domainList(src);
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
  const list = domainList(src);
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
  const list = domainList(src);
  const status = list.root.querySelector('[data-u2-part="status"]');
  assert.equal(status.querySelector('.u2-loader') !== null, true, 'loading');
  await flush();
  assert.equal(status.textContent, 'No issues.');
  list.dispose();
  src.dispose();

  const broken = new DomainSource({table: 'grit.nope'});
  broken.start();
  const failing = domainList(broken, {empty: 'Nothing here'});
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
  const list = domainList(src);
  await flush();
  assert.equal(src.currentRow.value?.id, 'i1', 'the first row is current once loaded');
  assert.equal(list.list.selectedIndex.value, 0);
  list.dispose();
  src.dispose();
  const draft = table.draft({title: 'D'});
  const over = domainList(draft);
  await flush();
  assert.equal(draft.currentRow.value.id.startsWith('~'), true, 'the draft stays current');
  over.dispose();
  draft.dispose();
});

scoped('keyboard: Enter hands the row to the form through the source, Delete deletes, actions rove with the selection', async () => {
  DG.ObjectHandler.register(new DG.DomainObjectHandler('grit.issue'));
  const {src} = await issues();
  const list = domainList(src);
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
  const list = domainList(src);
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
  assert.equal(Rows.isDraft({id: '~row:3'}), true);
  assert.equal(Rows.isDraft({id: 'i1'}), false);
});

scoped('cards: the recipe card unless a registered handler subclass paints its own', async () => {
  DG.ObjectHandler.register(new DG.DomainObjectHandler('grit.issue'));
  const {src} = await issues();
  const plain = domainList(src, {mode: 'cards'});
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
  const custom = domainList(own, {mode: 'cards'});
  const cards = mount(custom);
  assert.equal(cards[0].querySelector('.test-issue-card')?.textContent, 'card: Aspirin', 'a subclass\'s own card');
  custom.dispose();
  own.dispose();
});
