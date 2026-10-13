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

async function editor(options = {}, env = draft()) {
  const e = new ManifestEditor(env, {context: {mode: 'create', storage: 'external'},
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
  assert.deepEqual([box('orderid').checked, box('orderid').disabled, box('orderid').getAttribute('aria-disabled')],
    [true, false, 'true'], 'a key column is locked, painted checked');
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
  assert.deepEqual(inputNames(e), ['friendlyName', 'name', 'writable', 'access-schema']);
  // every hint is its input's postfix — on the editor's line, never on a line of its own
  assert.equal(panel(e).querySelector('[data-u2-name="name"] .u2-input-postfix').textContent,
    'registered as ext_northwind');
  assert.equal(panel(e).querySelector('[data-u2-name="writable"] .u2-input-postfix').textContent,
    'users with Edit on a table may insert, update and delete rows in the warehouse');
  assert.equal(panel(e).querySelector('.u2-manifest-panel-hint'), null, 'no hint spans while editable');
  const friendly = panel(e).querySelector('[data-u2-name="friendlyName"] input');
  friendly.value = 'Northwind Sales 2';
  fire(friendly, 'change');
  await flush();
  assert.equal(e.model.name.value, 'northwind_sales_2', 'the identifier follows the friendly name');
  assert.equal(panel(e).querySelector('[data-u2-name="name"] input').value, 'northwind_sales_2');
  assert.equal(panel(e).querySelector('[data-u2-name="name"] .u2-input-postfix').textContent,
    'registered as ext_northwind_sales_2');
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

ui('a schema row inherited by a table is labelled, with no group list offered up front', async () => {
  const e = await editor({groups: undefined});
  e.access.addGroup({kind: 'schema'}, {id: 'g-sales', label: 'Sales'});
  await select(e, 'order_details');
  assert.deepEqual(panel(e).querySelectorAll('.u2-access-grid-principal').map((p) => p.textContent),
    ['You (creator)', 'Sales']);
  e.dispose();
});

ui('the locked boxes show what the plan grants: nothing beyond View on a read-only binding, the creator row included', async () => {
  const e = await editor();
  e.access.addGroup({kind: 'schema'}, {id: 'g-sales', label: 'Sales'});
  e.access.setGrant({kind: 'schema'}, 'g-sales', 'edit', true);
  await select(e, 'orders');
  const checked = () => panel(e).querySelectorAll('.u2-access-grid tbody tr')
    .map((r) => r.querySelectorAll('.u2-access-grid-check').map((b) => b.checked));
  assert.deepEqual(checked(), [[true, false, false], [true, false, false]],
    'the creator and the inherited Sales row, as plan() grants');
  assert.deepEqual(e.plan().grants.map((g) => [g.table, g.edit]), [['orders', false], ['order_details', false]]);
  e.model.setWritable(true);
  await flush();
  assert.deepEqual(checked(), [[true, true, true], [true, true, false]]);
  assert.deepEqual(e.plan().grants.map((g) => [g.table, g.edit]), [['orders', true], ['order_details', true]]);
  e.dispose();
});

ui('the Writable toggle rebuilds the panel with focus kept: the box itself, or a grid cell', async () => {
  const e = await editor();
  e.access.addGroup({kind: 'schema'}, {id: 'g-sales', label: 'Sales'});
  await flush();
  const writable = () => panel(e).querySelector('[data-u2-name="writable"] .u2-input-checkbox');
  const box = writable();
  box.focus();
  box.click();
  await flush();
  assert.equal(e.model.writable.value, true);
  assert.equal(writable() !== box, true, 'the panel was rebuilt');
  assert.equal(document.activeElement === writable(), true, 'a second Space toggles it back');
  const cells = () => panel(e).querySelectorAll('.u2-access-grid tbody tr')[1].querySelectorAll('.u2-access-grid-check');
  cells()[0].focus();
  e.model.setWritable(false);
  await flush();
  assert.equal(document.activeElement === cells()[0], true, 'Sales › View, in the rebuilt grid');
  cells()[1].focus();
  e.model.setWritable(true);
  await flush();
  assert.equal(document.activeElement === cells()[1], true, 'Sales › Edit, unlocked again');
  e.dispose();
});

ui('Tab out of a field that commits on change: the rebuild waits for the Tab to land, then keeps its target', async () => {
  const e = await editor();
  await select(e, 'orders');
  const friendly = panel(e).querySelector('[data-u2-name="friendlyName"] input');
  friendly.focus();
  fire(friendly, 'keydown', {key: 'Tab'});
  // what the browser reports inside `change` on a Tab: nothing focused yet
  document.activeElement = document.body;
  friendly.value = 'Sales orders';
  fire(friendly, 'change');
  assert.equal(e.model.table('orders').friendlyName, 'Sales orders', 'committed');
  assert.equal(panel(e).querySelector('[data-u2-name="friendlyName"] input') === friendly, true, 'not rebuilt yet');
  const next = panel(e).querySelector('[data-u2-name="nameColumn"] select');
  next.focus();
  await flush();
  const rebuilt = panel(e).querySelector('[data-u2-name="nameColumn"] select');
  assert.equal(rebuilt !== next, true, 'rebuilt once the Tab landed');
  assert.equal(document.activeElement === rebuilt, true, 'on the field the Tab went to');
  e.dispose();
});

ui('the table panel: names, the key as text, the pickers, relationships with their status, the inherited access', async () => {
  const e = await editor();
  await select(e, 'order_details');
  assert.deepEqual(e.selected.value, {kind: 'table', table: 'order_details'});
  assert.equal(title(e), 'Tableorder_details');
  assert.deepEqual(inputNames(e), ['logical', 'friendlyName', 'nameColumn', 'searchable', 'access-order_details']);
  assert.equal(panel(e).querySelector('[data-u2-name="key"] .u2-form-readonly-value').textContent,
    'orderid, productidthe primary key; the row id encodes it');
  assert.equal(panel(e).querySelector('[data-u2-name="key"] .u2-manifest-panel-hint').textContent,
    'the primary key; the row id encodes it', 'the hint on the key line');
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
  assert.equal(panel(e).querySelector('[data-u2-name="isName"] .u2-input-postfix').textContent,
    'shown wherever a row is referred to');
  assert.equal(panel(e).querySelector('[data-u2-name="searchable"] .u2-input-postfix').textContent,
    'one per table');
  assert.equal(panel(e).querySelectorAll('.u2-manifest-panel-note').length, 1, 'only the visibility note stands alone');
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

ui('a principal picker feeds the access grid and the visibility chips; the value stays {id, label}', async () => {
  const picks = [];
  const e = await editor({groups: [], principalPicker: (onPick) => {
    const b = document.createElement('button');
    b.className = 'test-pick';
    b.addEventListener('click', () => onPick({id: 'g-sales', label: 'Sales'}));
    picks.push(b);
    return b;
  }});
  assert.equal(panel(e).querySelector('.u2-access-grid-add'), null, 'no select when a picker is given');
  panel(e).querySelector('.u2-access-grid tfoot .test-pick').click();
  await flush();
  assert.deepEqual(e.access.grants.value, [{scope: {kind: 'schema'}, group: {id: 'g-sales', label: 'Sales'},
    view: true, edit: false, delete: false}]);
  assert.deepEqual(panel(e).querySelectorAll('.u2-access-grid-principal').map((p) => p.textContent),
    ['You (creator)', 'Sales'], 'a picked principal is shown by its label');
  await select(e, 'shipname');
  const visiblePicker = () => panel(e).querySelector('[data-u2-name="visibleTo"] .u2-input-box .test-pick');
  assert.equal(visiblePicker().style.display, 'none', 'visible to everyone: no picker to restrict with');
  const radios = panel(e).querySelectorAll('[data-u2-name="visibility"] input');
  radios[0].checked = false;
  radios[1].click();
  await flush();
  assert.equal(visiblePicker().style.display, '', 'only these groups: the picker is offered');
  assert.deepEqual(panel(e).querySelectorAll('.u2-chip').map((c) => c.textContent), [], 'nothing offered up front');
  visiblePicker().click();
  await flush();
  assert.deepEqual(panel(e).querySelectorAll('.u2-chip').map((c) => [c.textContent, c.getAttribute('aria-pressed')]),
    [['Sales', 'true']]);
  assert.deepEqual(e.access.visibilityOf('orders', 'shipname'), [{id: 'g-sales', label: 'Sales'}]);
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
    ['friendlyName', 'name', 'writable']);
  assert.equal(panel(e).querySelector('[data-u2-name="name"] .u2-form-readonly-value .u2-manifest-panel-hint').textContent,
    'registered as ext_northwind');
  assert.equal(panel(e).querySelector('.u2-access-grid-add').disabled, true);
  assert.equal(rowOf(e, 'orders').querySelector('.u2-tree-check').getAttribute('aria-disabled'), 'true');
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
  e.model.setWritable(true);
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

ui('a bindable table the draft did not request: counted in the header, badged "not in this draft", locked', async () => {
  const env = draft();
  env.inventory.tables.push({remote: 'customers', logical: 'customers', bindable: true, key: ['customerid']});
  const e = await editor({}, env);
  assert.equal(e.tree.root.querySelector('.u2-manifest-tree-count').textContent, '2 of 3 bindable tables');
  const row = rowOf(e, 'customers');
  assert.deepEqual(badges(row), ['not in this draft']);
  assert.deepEqual([row.querySelector('.u2-tree-check').checked, row.querySelector('.u2-tree-check').disabled], [false, true]);
  assert.match(row.title, /not in this draft/);
  await select(e, 'customers');
  assert.match(title(e), /Tablecustomersnot in this draft/);
  assert.deepEqual(inputNames(e), []);
  e.dispose();
});

ui('a column with two foreign keys lists both; when the foreign keys could not be read the table panel says so', async () => {
  const env = draft();
  env.inventory.relations.push({table: 'orders', column: 'shipvia', targetTable: 'carriers', targetColumn: 'carrierid',
    status: 'plain', code: 'external-ref-ambiguous', message: 'two foreign keys'});
  const e = await editor({}, env);
  assert.deepEqual(badges(rowOf(e, 'shipvia')), ['→ shippers', '→ carriers']);
  await select(e, 'shipvia');
  assert.deepEqual(panel(e).querySelectorAll('.u2-manifest-relation').map((r) => r.textContent),
    ['foreign key to shippersplain valueshippers is not in the draft', 'foreign key to carriersplain valuetwo foreign keys']);
  e.dispose();

  const unread = draft();
  unread.inventory.relations = [];
  unread.diagnostics = [{code: 'external-relations-unavailable', message: 'The foreign keys of "public" could not be read'}];
  const e2 = await editor({}, unread);
  await select(e2, 'orders');
  assert.match(panel(e2).querySelector('.u2-manifest-panel-note').textContent, /foreign keys of this schema could not be read/);
  assert.equal(panel(e2).querySelector('.u2-manifest-relation'), null);
  e2.dispose();
});

ui('WO-A5.1 #9: plan() grants no Edit or Delete on a read-only binding or a read-only table — the lock the panel shows', async () => {
  const e = await editor();
  const sales = {id: 'g-sales', label: 'Sales'};
  e.access.addGroup({kind: 'schema'}, sales);
  e.access.setGrant({kind: 'schema'}, 'g-sales', 'edit', true);
  e.access.setGrant({kind: 'schema'}, 'g-sales', 'delete', true);
  const rights = () => e.plan().grants.map((g) => [g.table, g.view, g.edit, g.delete]);
  assert.deepEqual(rights(), [['orders', true, false, false], ['order_details', true, false, false]],
    'a read-only binding: View only');
  e.model.setWritable(true);
  e.model.setReadOnly('order_details', true);
  assert.deepEqual(rights(), [['orders', true, true, true], ['order_details', true, false, false]],
    'a writable binding: Edit and Delete where the table is not read-only');
  e.dispose();
});
