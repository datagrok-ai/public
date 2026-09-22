/* ManifestEditor over the Northwind draft: the tree's rows, checkboxes and badges; the three
   panels for the three selection kinds; a checkbox in the tree and an "include" link in the
   panel driving the model; a diagnostic landing on its row and its panel; view mode as text;
   plan() resolved to logical names; the three tags registered. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {ManifestEditor} from '../src/dg/domain/authoring/manifest-editor.js';
import {authoring} from '../src/dg/domain/authoring/index.js';

const DRAFT = JSON.parse(readFileSync(new URL('./fixtures/authoring/northwind-draft.json', import.meta.url), 'utf8'));
const draft = () => JSON.parse(JSON.stringify(DRAFT));

function ui(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

async function editor(options = {}) {
  const e = new ManifestEditor(draft(), {context: {mode: 'create', storage: 'external'},
    groups: ['Sales', 'Developers'], friendlyName: 'Northwind sales', ...options});
  document.body.append(e.root);
  e.tree.tree.root.querySelector('.u2-list').clientHeight = 800;
  await flush();
  return e;
}

const rows = (e) => e.tree.root.querySelectorAll('.u2-tree-row');
const names = (e) => rows(e).map((r) => r.querySelector('.u2-manifest-node-name').textContent);
const rowOf = (e, name) => rows(e).find((r) => r.querySelector('.u2-manifest-node-name').textContent === name);
const badges = (row) => row.querySelectorAll('.u2-badge').map((b) => b.textContent);
const panel = (e) => e.panel.root;
const title = (e) => panel(e).querySelector('.u2-manifest-panel-title').textContent;
const inputNames = (e) => panel(e).querySelectorAll('.u2-input-root').map((i) => i.dataset.u2Name);

async function select(e, name) {
  fire(rowOf(e, name), 'click');
  await flush();
}

ui('the tree: schema › tables › columns, the first included table open; checkboxes, locks and badges', async () => {
  const e = await editor();
  assert.equal(e.root.dataset.u2, 'manifest-editor');
  assert.deepEqual(names(e), ['northwind', 'orders', 'orderid', 'customerid', 'employeeid', 'orderdate', 'requireddate',
    'shippeddate', 'shipvia', 'freight', 'shipname', 'shipaddress', 'shipcity', 'shipcountry', 'tracking_number',
    'order_details', 'order_summary', 'audit_log']);
  assert.equal(e.tree.root.querySelector('.u2-manifest-tree-count').textContent, '2 of 2 bindable tables');

  const box = (name) => rowOf(e, name).querySelector('.u2-tree-check');
  assert.equal(rowOf(e, 'northwind').querySelector('.u2-tree-check'), null, 'the schema has no checkbox');
  assert.deepEqual([box('orders').checked, box('orders').disabled], [true, false]);
  assert.deepEqual([box('orderid').checked, box('orderid').disabled], [true, true], 'a key column is locked');
  assert.deepEqual([box('tracking_number').checked, box('tracking_number').disabled], [false, true], 'unsupported');
  assert.deepEqual([box('order_summary').checked, box('order_summary').disabled], [false, true], 'a view');
  assert.deepEqual([box('audit_log').checked, box('audit_log').disabled], [false, true], 'keyless');
  assert.equal(rowOf(e, 'order_summary').classList.contains('u2-tree-row-disabled'), true);
  assert.match(rowOf(e, 'audit_log').title, /no primary key/);
  assert.match(rowOf(e, 'tracking_number').title, /not supported/);

  assert.deepEqual(badges(rowOf(e, 'orders')), ['orderid']);
  assert.deepEqual(badges(rowOf(e, 'orderid')), ['key']);
  assert.deepEqual(badges(rowOf(e, 'customerid')), ['→ customers']);
  assert.equal(rowOf(e, 'customerid').querySelector('.u2-badge-warning').textContent, '→ customers', 'plain');
  assert.deepEqual(badges(rowOf(e, 'shipname')), ['name']);
  assert.deepEqual(badges(rowOf(e, 'order_summary')), ['view']);
  assert.deepEqual(badges(rowOf(e, 'audit_log')), ['no key']);
  assert.equal(rowOf(e, 'shipvia').querySelector('.u2-manifest-node-hint').textContent, 'int');
  e.dispose();
});

ui('the schema panel: name, friendly name, writable, the access grid with the creator inherited', async () => {
  const e = await editor();
  assert.deepEqual(e.selected.value, {kind: 'schema'});
  assert.equal(title(e), 'Schemanorthwind');
  assert.deepEqual(inputNames(e), ['name', 'friendlyName', 'writable', 'access-schema']);
  const grid = panel(e).querySelector('.u2-access-grid');
  assert.deepEqual(grid.querySelectorAll('tbody tr').map((r) => r.querySelector('.u2-access-grid-principal').textContent),
    ['You (creator)']);
  assert.deepEqual(grid.querySelectorAll('tbody .u2-access-grid-check').map((b) => b.disabled), [true, true, true]);
  const picker = grid.querySelector('.u2-access-grid-add');
  picker.value = 'Sales';
  fire(picker, 'change');
  assert.deepEqual(e.access.grants.value, [{scope: {kind: 'schema'}, group: {id: 'Sales', label: 'Sales'},
    view: true, edit: false, delete: false}]);
  assert.deepEqual(grid.querySelectorAll('tbody tr')[1].querySelectorAll('.u2-access-grid-check').map((b) => b.disabled),
    [false, true, true], 'Edit and Delete need a writable binding');

  const writable = panel(e).querySelector('[data-u2-name="writable"] .u2-input-checkbox');
  writable.click();
  await flush();
  assert.equal(e.model.writable.value, true);
  assert.deepEqual(panel(e).querySelector('.u2-access-grid').querySelectorAll('tbody tr')[1]
    .querySelectorAll('.u2-access-grid-check').map((b) => b.disabled), [false, false, false]);
  e.dispose();
});

ui('the table panel: names, the key as text, the pickers, relationships with their status, the inherited access', async () => {
  const e = await editor();
  await select(e, 'order_details');
  assert.deepEqual(e.selected.value, {kind: 'table', table: 'order_details'});
  assert.equal(title(e), 'Tableorder_details');
  assert.deepEqual(inputNames(e), ['logical', 'friendlyName', 'nameColumn', 'searchable', 'access-order_details']);
  assert.equal(panel(e).querySelector('[data-u2-name="key"] .u2-form-readonly-value').textContent,
    'orderid, productid — the primary key; the row id encodes it');
  const relations = panel(e).querySelectorAll('.u2-manifest-relation');
  assert.deepEqual(relations.map((r) => r.querySelector('.u2-badge').textContent), ['ref', 'plain value']);
  assert.equal(relations[0].querySelector('.u2-link'), null, 'nothing to fix on a ref');
  assert.equal(relations[1].querySelector('.u2-link'), null, 'products is not in the draft — nothing to include');

  fire(rowOf(e, 'orders').querySelector('.u2-tree-check'), 'click');
  await flush();
  assert.equal(e.model.table('orders').included, false);
  assert.deepEqual(e.selected.value, {kind: 'table', table: 'order_details'}, 'the checkbox did not move the selection');
  const plain = panel(e).querySelectorAll('.u2-manifest-relation')[0];
  assert.equal(plain.querySelector('.u2-badge').textContent, 'plain value');
  assert.equal(plain.querySelector('.u2-link').textContent, 'include orders');
  assert.equal(e.model.column('order_details', 'orderid').type, 'int');
  fire(plain.querySelector('.u2-link'), 'click');
  await flush();
  assert.equal(e.model.table('orders').included, true);
  assert.equal(panel(e).querySelectorAll('.u2-manifest-relation')[0].querySelector('.u2-badge').textContent, 'ref');
  e.dispose();
});

ui('the column panel: logical name, the type beside the warehouse type, required, name, searchable, reference, visibility', async () => {
  const e = await editor();
  await select(e, 'shipname');
  assert.equal(title(e), 'Columnorders.shipname');
  assert.deepEqual(inputNames(e), ['logical', 'required', 'isName', 'searchable', 'visibility', 'visibleTo']);
  assert.equal(panel(e).querySelector('[data-u2-name="type"] .u2-form-readonly-value').textContent, 'string');
  assert.equal(panel(e).querySelector('[data-u2-name="isName"] .u2-input-checkbox').checked, true);
  assert.equal(panel(e).querySelector('.u2-manifest-relation'), null, 'no foreign key, no Reference section');
  const chips = panel(e).querySelectorAll('.u2-chip');
  assert.deepEqual(chips.map((c) => c.textContent), ['Sales', 'Developers']);
  assert.equal(chips[0].disabled, true, 'everyone may see the row until "only these groups" is chosen');
  // the shim does not uncheck a radio's siblings on click
  const radios = panel(e).querySelectorAll('[data-u2-name="visibility"] input');
  radios[0].checked = false;
  radios[1].click();
  await flush();
  fire(panel(e).querySelectorAll('.u2-chip')[0], 'click');
  assert.deepEqual(e.access.visibility.value, [{table: 'orders', column: 'shipname',
    groups: [{id: 'Sales', label: 'Sales'}]}]);

  await select(e, 'orderid');
  assert.equal(panel(e).querySelector('[data-u2-name="required"] .u2-input-checkbox').disabled, true, 'a key is always required');
  assert.match(panel(e).textContent, /A key column is visible to everyone/);
  await select(e, 'tracking_number');
  assert.equal(title(e), 'Columnorders.tracking_numbernot bindablebigint is not supported: lossy through the platform');
  assert.deepEqual(inputNames(e), []);
  await select(e, 'shipvia');
  assert.equal(panel(e).querySelector('.u2-manifest-relation').textContent,
    'foreign key to shippersplain valueshippers is not in the draft');
  e.dispose();
});

ui('a rename commits through the panel; a refused name keeps the old one and shows why', async () => {
  const e = await editor();
  await select(e, 'orders');
  const input = panel(e).querySelector('[data-u2-name="logical"] input');
  input.value = 'Orders';
  fire(input, 'change');
  await flush();
  assert.equal(e.model.table('orders').logical, 'orders');
  assert.match(panel(e).querySelector('[data-u2-name="logical"] .u2-input-error').textContent, /lowercase/);
  const again = panel(e).querySelector('[data-u2-name="logical"] input');
  again.value = 'order';
  fire(again, 'change');
  await flush();
  assert.equal(e.model.table('orders').logical, 'order');
  assert.equal(rowOf(e, 'orders').querySelector('.u2-manifest-node-hint').textContent, '· order');
  assert.deepEqual(Object.keys(e.model.toJSON().tables), ['order', 'order_details']);
  e.dispose();
});

ui('a diagnostic addressed by manifest path marks its row and shows in its panel', async () => {
  const e = await editor();
  e.diagnostics.value = [{path: 'tables.orders.columns.shipvia', code: 'external-ref-key', message: 'shippers.shipperid is not a single-column key'},
    {code: 'external-unreachable', message: 'the warehouse did not answer'}];
  await flush();
  assert.equal(rowOf(e, 'shipvia').querySelector('.u2-manifest-node').classList.contains('u2-manifest-node-problem'), true);
  assert.match(rowOf(e, 'shipvia').title, /single-column key/);
  assert.equal(rowOf(e, 'northwind').querySelector('.u2-manifest-node').classList.contains('u2-manifest-node-problem'), true);
  assert.deepEqual(panel(e).querySelectorAll('.u2-manifest-panel-diagnostic .u2-badge').map((b) => b.textContent),
    ['external-unreachable'], 'the schema panel lists the schema-wide finding');
  await e.select({kind: 'column', table: 'orders', column: 'shipvia'});
  await flush();
  assert.deepEqual(e.selected.value, {kind: 'column', table: 'orders', column: 'shipvia'});
  assert.deepEqual(panel(e).querySelectorAll('.u2-manifest-panel-diagnostic').map((d) => d.textContent),
    ['external-ref-keyshippers.shipperid is not a single-column key']);
  e.dispose();
});

ui('view mode: every field is text, the checkboxes are locked, the access grid cannot be edited', async () => {
  const e = await editor({context: {mode: 'view', storage: 'external'}});
  assert.equal(e.offer.editable, false);
  assert.deepEqual(panel(e).querySelectorAll('[data-u2="readonly-field"]').map((f) => f.dataset.u2Name),
    ['name', 'friendlyName', 'writable']);
  assert.equal(panel(e).querySelector('.u2-access-grid-add').disabled, true);
  assert.equal(rowOf(e, 'orders').querySelector('.u2-tree-check').disabled, true);
  assert.equal(e.tree.root.querySelector('.u2-manifest-tree-header .u2-link'), null, 'no Check all / Clear');
  await select(e, 'shipname');
  assert.deepEqual(inputNames(e), []);
  assert.match(panel(e).textContent, /Everyone who may see the row/);
  e.dispose();
});

ui('plan() fans a schema row out over every included table, merges it with the table rows, resolves logical names', async () => {
  const e = await editor({groups: [{id: 'g-sales', label: 'Sales'}, 'Developers']});
  const sales = {id: 'g-sales', label: 'Sales'};
  const dev = {id: 'Developers', label: 'Developers'};
  e.access.addGroup({kind: 'schema'}, sales);
  e.access.addGroup({kind: 'table', table: 'orders'}, dev);
  e.access.addGroup({kind: 'table', table: 'orders'}, sales);
  e.access.setGrant({kind: 'table', table: 'orders'}, 'g-sales', 'edit', true);
  e.access.addGroup({kind: 'table', table: 'order_details'}, dev);
  e.access.setVisibility('orders', 'freight', [sales]);
  e.access.setVisibility('orders', 'shipcity', [sales]);
  e.model.renameTable('orders', 'order');
  e.model.renameColumn('orders', 'freight', 'freight_cost');
  e.model.includeColumn('orders', 'shipcity', false);
  e.model.name.value = 'northwind_sales';
  let plan = e.plan();
  assert.equal(plan.name, 'northwind_sales');
  assert.equal(plan.friendlyName, 'Northwind sales');
  assert.equal(plan.manifest.name, 'northwind_sales');
  assert.deepEqual(plan.grants, [
    {table: 'order', group: sales, view: true, edit: true, delete: false},
    {table: 'order_details', group: sales, view: true, edit: false, delete: false},
    {table: 'order', group: dev, view: true, edit: false, delete: false},
    {table: 'order_details', group: dev, view: true, edit: false, delete: false},
  ], 'one grant per table and group; the schema row and the table row for Sales merged');
  assert.deepEqual(plan.restrictions, [{table: 'order', column: 'freight_cost', groups: [sales]}]);
  e.model.includeTable('order_details', false);
  plan = e.plan();
  assert.deepEqual(Object.keys(plan.manifest.tables), ['order']);
  assert.deepEqual(plan.grants.map((g) => [g.table, g.group.id]), [['order', 'g-sales'], ['order', 'Developers']],
    'an excluded table gets nothing');
  e.dispose();
});

ui('the namespace and the tags', async () => {
  assert.equal(authoring.ManifestEditor, ManifestEditor);
  assert.deepEqual(Object.keys(authoring).sort(), ['AccessModel', 'ManifestContextPanel', 'ManifestEditor',
    'ManifestModel', 'ManifestRules', 'ManifestTree', 'fieldOffer']);
});

ui('the editor tag is registered next to the other domain controls; the panes alone are not', async () => {
  const {register} = await import('node:module');
  register('./dg-stub.mjs', import.meta.url);
  const {registerDomainComponents} = await import('../src/dg/domain/registrations.js');
  const {Registry} = await import('../src/spec/registry.js');
  const reg = new Registry();
  registerDomainComponents(reg);
  const meta = reg.get('u2-manifest-editor');
  assert.equal(meta.category, 'Inputs');
  assert.equal(reg.get('u2-manifest-tree'), undefined, 'a tree and a panel over separate models cannot share state');
  assert.equal(reg.get('u2-manifest-panel'), undefined);
  const built = meta.create({draft: draft(), groups: ['Sales']});
  assert.equal(built.root.dataset.u2, 'manifest-editor');
  built.dispose();
  const {domains} = await import('../src/dg/domain/index.js');
  assert.equal(domains.authoring, authoring);
});
