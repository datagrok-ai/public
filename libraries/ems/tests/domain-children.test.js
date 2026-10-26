/* `domains.children` (2-8) over the memory backend: a tab per table referring to the parent (two
   FKs to one parent → two labelled tabs), the shown tab's source in the parent's session and
   queried by the FK with New pre-filled, the pane following the current row, create-related-inline
   (a draft parent and a child draft as one transaction, the child's FK the parent's real id, the
   tab re-queried by it), the `u2-domain-children` tag — and `DomainPick.resolve` naming a draft
   from its live source without a query. `DG` comes from tests/dg-stub.mjs. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, signal, Rows, Registry, SpecContext, renderSpec, registerAll, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {Property, WidgetDescriptor, platform} from '@datagrok-libraries/u2/tests/platform-doubles.mjs';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainChildren} = await import('../src/dg/children.js');
const {DomainPick} = await import('../src/dg/pick.js');
const {registerDomainComponents} = await import('../src/dg/registrations.js');

/** Companies with contacts, and deals that point at a company twice (seller and buyer). */
const SCHEMA = {
  name: 'crm',
  tables: {
    company: {friendlyName: 'Companies', singularName: 'Company',
      columns: {name: {type: 'string', required: true, isName: true}}},
    contact: {friendlyName: 'Contacts', columns: {
      company_id: {type: 'ref', ref: 'company', required: true},
      name: {type: 'string', required: true, isName: true},
    }},
    deal: {friendlyName: 'Deals', columns: {
      seller_id: {type: 'ref', ref: 'company', friendlyName: 'Seller'},
      buyer_id: {type: 'ref', ref: 'company', friendlyName: 'Buyer'},
      title: {type: 'string', isName: true},
    }},
  },
};
const ROWS = {
  company: [{id: 'c1', name: 'Acme'}, {id: 'c2', name: 'Globex'}],
  contact: [{id: 'k1', company_id: 'c1', name: 'Ann'}, {id: 'k2', company_id: 'c2', name: 'Bob'}],
  deal: [{id: 'd1', seller_id: 'c1', buyer_id: 'c2', title: 'Big'}],
};
const crm = () => new MemoryDomainBackend(SCHEMA,
  {rows: Object.fromEntries(Object.entries(ROWS).map(([t, list]) => [t, list.map((r) => ({...r}))]))});
const UUID = /^[0-9a-f-]{36}$/;
const CONTACTS = 'crm.contact.company_id';

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    WidgetDescriptor.registry = [new WidgetDescriptor('Grid', [new Property('allowEdit', 'bool', {defaultValue: false})])];
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      WidgetDescriptor.registry = [];
      platform.reset();
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** The parent's load, then the child handles. */
async function settle() {
  await flush();
  await flush();
}

scoped('a tab per child table, two FKs told apart; the shown tab\'s source shares the session, queried by the FK', async () => {
  emsBackends.domain = crm();
  const companies = await domains.table('crm.company');
  const parent = companies.source({pageSize: 10});
  await flush();
  parent.currentRow.value = parent.rows.byKey('c1');
  const children = domains.children(parent, {mode: 'list'});
  assert.equal(children.root.dataset.u2, 'domain-children');
  assert.deepEqual(children.entries.value, [], 'the child tables come with the handles');
  await settle();
  assert.deepEqual(children.entries.value.map((e) => e.label), ['Contacts', 'Deals (Seller)', 'Deals (Buyer)']);
  assert.deepEqual(children.entries.value.map((e) => e.id), [CONTACTS, 'crm.deal.seller_id', 'crm.deal.buyer_id']);
  assert.deepEqual(children.entries.value.map((e) => e.fk), ['company_id', 'seller_id', 'buyer_id']);
  assert.equal(children.tabs.activeTab.value, CONTACTS);
  const contacts = children.child(CONTACTS);
  assert.equal(contacts.session, parent.session, 'one session');
  assert.equal(contacts.query.value, 'company_id = "c1"');
  assert.deepEqual(contacts.defaults, {company_id: 'c1'}, 'New is pre-filled with the parent');
  assert.equal(contacts.isEmpty, false);
  await flush();
  assert.deepEqual(contacts.rows.items.value.map((r) => r.name), ['Ann']);
  assert.equal(children.child('crm.deal.seller_id'), undefined, 'a tab not shown builds nothing');
  assert.equal(children.root.querySelectorAll('.u2-tabs-tab').length, 3);
  const panes = [...children.root.querySelectorAll('[data-u2-part="pane"]')];
  assert.equal(panes.length, 1, 'only the shown tab has a pane');
  assert.equal(panes[0].querySelector('[data-u2="new-button"]') !== null, true, 'a tab-local New');
  assert.equal(panes[0].querySelector('.u2-domain-list') !== null, true);
  assert.equal(panes[0].querySelector('.u2-domain-form') !== null, true, 'list mode: a list beside a form');

  children.tabs.activeTab.value = 'crm.deal.buyer_id';
  await flush();
  const bought = children.child('crm.deal.buyer_id');
  assert.equal(bought.query.value, 'buyer_id = "c1"');
  await flush();
  assert.deepEqual(bought.rows.items.value, []);

  parent.currentRow.value = parent.rows.byKey('c2');
  await settle();
  assert.equal(contacts.scope.isDisposed, true, 'a new parent row is a new child source');
  assert.equal(children.child(CONTACTS).query.value, 'company_id = "c2"');
  assert.deepEqual(children.child(CONTACTS).rows.items.value.map((r) => r.name), ['Bob']);
  assert.deepEqual(children.child('crm.deal.buyer_id').rows.items.value.map((r) => r.title), ['Big']);
  assert.equal(parent.session.sources.value.length, 3, 'the parent and the two shown children');

  parent.currentRow.value = null;
  await flush();
  assert.equal(children.child(CONTACTS), undefined);
  assert.equal(panes[0].textContent, 'Select a company.');
  children.dispose();
  parent.dispose();
  assert.deepEqual(parent.session.sources.value, []);
});

scoped('a child table that does not open fails its own tab; the others work', async () => {
  emsBackends.domain = crm();
  const unregister = domains.configure('crm.deal', () => {
    throw new Error('deal rules did not load');
  }, {id: 'test'});
  try {
    const companies = await domains.table('crm.company');
    const parent = companies.source();
    await flush();
    parent.currentRow.value = parent.rows.byKey('c1');
    const children = domains.children(parent, {mode: 'list'});
    await settle();
    assert.equal(children.root.querySelectorAll('.u2-tabs-tab').length, 3, 'every child table keeps its tab');
    assert.deepEqual(children.entries.value.map((e) => e.id), [CONTACTS], 'only the tables that opened');
    assert.equal(children.root.querySelector('[data-u2-part="status"]').textContent, '', 'no error over all tabs');
    await flush();
    assert.deepEqual(children.child(CONTACTS).rows.items.value.map((r) => r.name), ['Ann']);
    children.tabs.activeTab.value = 'crm.deal.seller_id';
    await flush();
    const pane = [...children.root.querySelectorAll('[data-u2-part="pane"]')].at(-1);
    assert.equal(pane.querySelector('.u2-domain-children-error').textContent, 'deal rules did not load');
    children.dispose();
    parent.dispose();
  } finally {
    unregister();
  }
});

scoped('a failed tab first is named in the plural and gives way to a tab with rows', async () => {
  emsBackends.domain = crm();
  const unregister = domains.configure('crm.deal', () => {
    throw new Error('deal rules did not load');
  }, {id: 'test'});
  try {
    const companies = await domains.table('crm.company');
    const parent = companies.source();
    await flush();
    parent.currentRow.value = parent.rows.byKey('c1');
    const children = domains.children(parent, {tables: ['deal', 'contact']});
    await settle();
    const labels = [...children.root.querySelectorAll('.u2-tabs-tab')].map((t) => t.textContent);
    assert.deepEqual(labels.slice(0, 2), ['Deals (Seller)', 'Deals (Buyer)']);
    assert.equal(children.tabs.activeTab.value, CONTACTS, 'the tab with rows, not the failed first one');
    children.dispose();
    parent.dispose();
  } finally {
    unregister();
  }
});

scoped('create-related-inline: a draft parent and its child draft save as one transaction; the tab re-queries by the real id', async () => {
  const memory = crm();
  emsBackends.domain = memory;
  const companies = await domains.table('crm.company');
  const parent = companies.draft();
  await flush();
  const children = domains.children(parent, {mode: 'list'});
  await settle();
  const contacts = children.child(CONTACTS);
  assert.equal(contacts.isEmpty, true, 'a draft parent: nothing to query');
  assert.equal(contacts.query.value, '');
  assert.equal(Rows.isDraft(contacts.defaults.company_id), true, 'the defaults carry the draft id');
  await flush();
  assert.equal(contacts.state.value, 'ready');
  const kid = contacts.newRow({name: 'Peter'});
  assert.equal(kid.company_id, parent.currentRow.value.id);
  parent.currentRow.value.name = 'Initech';
  assert.equal(parent.session.changeCount.value, 2);
  assert.equal(await parent.session.save(), true);
  await settle();
  const saved = parent.currentRow.value;
  assert.match(saved.id, UUID, 'the parent has its real id');
  const stored = memory.tableSync('crm.contact').rows.find((r) => r.name === 'Peter');
  assert.equal(stored.company_id, saved.id, 'the child\'s FK is the parent\'s real id');
  assert.equal(parent.session.isDirty.value, false);
  const fresh = children.child(CONTACTS);
  assert.notEqual(fresh, contacts, 'the tab re-queried by the real id');
  assert.equal(fresh.query.value, `company_id = "${saved.id}"`);
  assert.equal(fresh.isEmpty, false);
  await flush();
  assert.deepEqual(fresh.rows.items.value.map((r) => [r.name, r.company_id]), [['Peter', saved.id]]);
  assert.equal(fresh.isDirty.value, false, 'reloaded from the backend, clean');
  children.dispose();
  parent.dispose();
});

scoped('grid mode hosts a domains.grid per pane; `tables` narrows the tabs', async () => {
  emsBackends.domain = crm();
  const companies = await domains.table('crm.company');
  const parent = companies.source();
  await flush();
  parent.currentRow.value = parent.rows.byKey('c1');
  const children = domains.children(parent, {tables: ['deal']});
  await settle();
  assert.deepEqual(children.entries.value.map((e) => e.label), ['Deals (Seller)', 'Deals (Buyer)']);
  const pane = children.root.querySelector('[data-u2-part="pane"]');
  assert.equal(pane.querySelector('[data-u2="domain-grid"]') !== null, true);
  // the chain the height rules hang off: tab strip → panel → pane → grid host
  assert.equal(children.tabs.root.classList.contains('u2-domain-children-tabs'), true);
  assert.equal(pane.parentElement.classList.contains('u2-tabs-panel'), true);
  assert.equal(pane.querySelector('[data-u2="new-button"]') !== null, true);
  assert.equal(children.child('crm.deal.seller_id').query.value, 'seller_id = "c1"');
  children.dispose();
  parent.dispose();
});

scoped('spec: u2-domain-children over a bound source', async () => {
  emsBackends.domain = crm();
  const companies = await domains.table('crm.company');
  const src = companies.source();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  assert.deepEqual(reg.get('u2-domain-children').props.slice(0, 3).map((p) => p.name), ['source', 'tables', 'mode']);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-domain-children', name: 'kids',
    bind: {source: '$.companies'}, props: {mode: 'list', tables: ['contact']}}},
  new SpecContext({data: {companies: signal(src)}}), reg);
  await settle();
  const kids = instance.node('kids');
  assert.equal(kids instanceof DomainChildren, true);
  assert.equal(kids.mode, 'list');
  assert.deepEqual(kids.entries.value.map((e) => e.label), ['Contacts']);
  instance.dispose();
  src.dispose();
});

scoped('DomainPick.resolve: a draft id is the draft\'s name from its live source, no query issued', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.project');
  const src = table.source();
  await flush();
  const t = emsBackends.domain.tableSync('grit.project');
  const queries = [];
  const query = t.query.bind(t);
  t.query = (spec) => {
    queries.push(spec);
    return query(spec);
  };
  const named = src.newRow({key: 'NP', name: 'New Project X'});
  assert.deepEqual(await DomainPick.resolve('grit.project', named.id), {id: named.id, name: 'New Project X'});
  const unnamed = src.newRow({key: 'NN'});
  assert.deepEqual(await DomainPick.resolve('grit.project', unnamed.id), {id: unnamed.id, name: 'New project'});
  assert.equal(queries.length, 0, 'resolved from the live source');
  const gone = `${Rows.DRAFT_PREFIX}nobody`;
  assert.deepEqual(await DomainPick.resolve('grit.project', gone), {id: gone, name: gone}, 'a draft nobody holds: itself');
  assert.equal(queries.length, 1);
  src.dispose();
});

scoped('a tab label is the plural as a caption: underscores out, first letter up', async () => {
  emsBackends.domain = new MemoryDomainBackend({name: 'crm', tables: {
    company: {friendlyName: 'Companies', columns: {name: {type: 'string', required: true, isName: true}}},
    order_line: {columns: {company_id: {type: 'ref', ref: 'company', required: true}}},
  }}, {rows: {company: [{id: 'c1', name: 'Acme'}], order_line: []}});
  const companies = await domains.table('crm.company');
  const parent = companies.source({pageSize: 10});
  await flush();
  parent.currentRow.value = parent.rows.byKey('c1');
  const children = domains.children(parent, {mode: 'list'});
  await settle();
  assert.deepEqual(children.entries.value.map((e) => e.label), ['Order lines']);
  children.dispose();
  parent.dispose();
});

scoped('an empty child collection says so, and the caller names the tab order', async () => {
  emsBackends.domain = crm();
  const companies = await domains.table('crm.company');
  const parent = companies.source({pageSize: 10});
  await flush();
  parent.currentRow.value = parent.rows.byKey('c1');
  const children = domains.children(parent, {tables: ['deal', 'contact'], mode: 'list'});
  await settle();
  assert.deepEqual(children.entries.value.map((e) => e.label), ['Deals (Seller)', 'Deals (Buyer)', 'Contacts'],
    'the tabs are in the order the caller named the tables');
  const pane = children.root.querySelector('[data-u2-part="pane"]');
  await flush();
  assert.deepEqual(children.child('crm.deal.seller_id').rows.items.value.map((r) => r.title), ['Big']);
  assert.equal(pane.querySelector('.u2-domain-children-empty').hidden, true, 'the collection is not empty');

  children.tabs.activeTab.value = 'crm.deal.buyer_id';
  await settle();
  const empty = [...children.root.querySelectorAll('.u2-domain-children-empty')];
  assert.equal(empty.some((el) => !el.hidden && el.textContent === 'No deals (buyer).'), true,
    'an empty child collection says so instead of showing an empty collection');
  children.dispose();
  parent.dispose();
});

scoped('a row opens on the first tab that has rows; the first one when none has', async () => {
  emsBackends.domain = crm();
  emsBackends.domain.tableSync('crm.company').rows.push({id: 'c3', name: 'Initech'}, {id: 'c4', name: 'Umbrella'});
  emsBackends.domain.tableSync('crm.deal').rows.push({id: 'd2', seller_id: 'c3', title: 'Small'});
  const companies = await domains.table('crm.company');
  const parent = companies.source({pageSize: 10});
  await flush();
  parent.currentRow.value = parent.rows.byKey('c3');
  const children = domains.children(parent, {mode: 'list'});
  await settle();
  await flush();
  assert.equal(children.tabs.activeTab.value, 'crm.deal.seller_id', 'Contacts has none for this one');
  assert.deepEqual(children.child('crm.deal.seller_id').rows.items.value.map((r) => r.title), ['Small']);
  children.dispose();

  const bare = domains.children(parent, {mode: 'list'});
  parent.currentRow.value = parent.rows.byKey('c4');
  await settle();
  await flush();
  assert.equal(bare.tabs.activeTab.value, CONTACTS, 'nothing anywhere: the first tab stands');
  bare.dispose();

  const filled = domains.children(parent, {mode: 'list'});
  parent.currentRow.value = parent.rows.byKey('c1');
  await settle();
  await flush();
  assert.equal(filled.tabs.activeTab.value, CONTACTS, 'and the first tab when it is the one with rows');
  filled.dispose();
  parent.dispose();
});

scoped('a New from the children toolbar is an unsaved change of its own: the gate and Save see it', async () => {
  emsBackends.domain = crm();
  const companies = await domains.table('crm.company');
  const parent = companies.source({pageSize: 10});
  await flush();
  parent.currentRow.value = parent.rows.byKey('c1');
  const children = domains.children(parent, {mode: 'list'});
  await settle();
  const contacts = children.child(CONTACTS);
  assert.equal(parent.session.isDirty.value, false);
  children.root.querySelector('button[data-u2="new-button"]').click();
  await flush();
  assert.equal(Rows.isDraft(contacts.currentRow.value), true, 'the draft is current');
  assert.equal(parent.session.isDirty.value, true, 'the press is the change — Ctrl+S and the leave gate answer');
  assert.equal(parent.session.summary.value, '1 unsaved change');
  children.dispose();
  parent.dispose();
});

scoped('VR-35/36, VS-10: counts in the labels, a heading for one table, "Add <singular>", dirty and refused tabs marked', async () => {
  emsBackends.domain = crm();
  const companies = await domains.table('crm.company');
  const parent = companies.source({pageSize: 10});
  await flush();
  parent.currentRow.value = parent.rows.byKey('c1');
  const children = domains.children(parent, {mode: 'list'});
  await settle();
  await flush();
  const labels = () => [...children.root.querySelectorAll('.u2-tabs-tab .u2-tabs-label')].map((t) => t.textContent);
  assert.deepEqual(labels(), ['Contacts (1)', 'Deals (Seller) (1)', 'Deals (Buyer) (0)'], 'every tab says how many');
  assert.equal(children.root.classList.contains('u2-domain-children-single'), false);
  const add = children.root.querySelector('button[data-u2="new-button"]');
  assert.equal(add.textContent, 'Add contact', 'the child button says what it adds');
  const tab = children.root.querySelector(`.u2-tabs-tab[data-id="${CONTACTS}"]`);
  add.click();
  await flush();
  assert.equal(tab.querySelector('.u2-tabs-mark').classList.contains('u2-tabs-dirty'), true, 'a pending child marks its tab');
  assert.equal(await parent.session.save(), false, 'the draft contact has no name');
  await flush();
  assert.equal(tab.classList.contains('u2-tabs-invalid'), true, 'the refused child marks its tab');
  parent.session.discard();
  await flush();
  assert.equal(tab.classList.contains('u2-tabs-invalid') || tab.querySelector('.u2-tabs-mark').classList.contains('u2-tabs-dirty'),
    false, 'and the marks go with the changes');

  children.tabs.activeTab.value = 'crm.deal.buyer_id';
  await settle();
  assert.equal(children.root.classList.contains('u2-domain-children-blank'), true, 'an empty tab takes only its own height');
  children.dispose();

  const one = domains.children(parent, {tables: ['contact'], mode: 'grid'});
  await settle();
  await flush();
  assert.equal(one.root.classList.contains('u2-domain-children-single'), true);
  assert.equal(one.root.querySelector('.u2-domain-children-heading .u2-section-title').textContent, 'Contacts (1)',
    'one child table is a heading with its count');
  one.dispose();
  parent.dispose();
});
