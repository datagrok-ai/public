/* The manifest editor over a REGISTERED binding (edit mode): the baseline round-trips as a no-op
   apply, registered names are locked by membership, the draft supplies candidates and the
   catalog facts that mark drift, declared refs stay authoritative and a new foreign key is a
   suggestion, dormant opt-outs and unexposed keys survive, the access snapshot round-trips and
   edits become exact permission-triple deltas, and the three-way rebase keeps what still
   applies, reports conflicts and drops what vanished; the panels in edit mode. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {ManifestEditor} from '../src/dg/domain/authoring/manifest-editor.js';
import {ManifestModel, AccessModel} from '../src/dg/domain/authoring/manifest-model.js';

const FIXTURE = JSON.parse(readFileSync(new URL('./fixtures/authoring/northwind-registered.json', import.meta.url), 'utf8'));
const copy = (x) => JSON.parse(JSON.stringify(x));
const baseline = () => copy(FIXTURE.manifest);
const draft = () => copy(FIXTURE.draft);
const snapshot = () => copy(FIXTURE.snapshot);
const TOKENS = {ifVersion: '3', ifIncarnation: '2026-09-20T10:00:00.000Z'};
const ORDERS = {kind: 'table', table: 'orders'};
const SALES = {id: 'g-sales', label: 'Sales'};
const ME = {id: 'g-me', label: 'askalkin'};
const DEV = {id: 'g-dev', label: 'Developers'};

const model = (d = draft()) => new ManifestModel(d, {baseline: baseline(), friendlyName: 'Northwind sales', description: 'Sales data'});
const editor = (options = {}, d = draft()) => new ManifestEditor(d, {context: {mode: 'edit', storage: 'external'},
  baseline: baseline(), snapshot: snapshot(), friendlyName: 'Northwind sales', description: 'Sales data', ...options});
const column = (m, table, remote) => m.columns(table).value.find((c) => c.remote === remote);
const ids = (changes) => changes.map((c) => c.id);
const bare = (ops) => ops.map(({change: _change, ...op}) => op);

test('open without an edit is a no-op apply; edited and reverted, it is a no-op again', () => {
  const e = editor();
  assert.equal(e.model.editing, true);
  assert.deepEqual([e.model.version, e.model.incarnation], ['3', TOKENS.ifIncarnation]);
  assert.deepEqual(e.model.tables.value.map((t) => [t.remote, t.logical, t.included, t.registered]), [
    ['orders', 'orders', true, true], ['order_details', 'order_lines', true, true], ['shippers', 'shippers', true, true],
    ['employees', 'employees', false, false], ['customers', 'customers', false, false], ['audit_log', 'audit_log', false, false]]);
  let plan = e.editPlan();
  assert.deepEqual(plan.payload, TOKENS, 'a redundant remote name in the registered descriptor is no change');
  assert.deepEqual(plan.changes, []);
  assert.equal(e.model.toJSON().tables.order_lines.table, 'order_details');
  assert.deepEqual(e.model.toJSON().tables.order_lines.columns.order_id, {type: 'ref', ref: 'orders', required: true, column: 'orderid'});

  e.model.setRequired('orders', 'freight', true);
  e.model.includeTable('orders', false);
  e.model.includeColumn('order_details', 'discount', true);
  e.model.renameTable('employees', 'staff');
  e.model.setSchemaFriendlyName('Other');
  e.model.setDescription('');
  e.model.setWritable(true);
  e.access.setGrant(ORDERS, 'g-sales', 'edit', false);
  e.access.setVisibility('orders', 'shipname', [SALES]);
  e.access.setVisibility('orders', 'freight', null);
  assert.notDeepEqual(e.editPlan().payload, TOKENS);
  e.model.setRequired('orders', 'freight', false);
  e.model.includeTable('orders', true);
  e.model.includeColumn('order_details', 'discount', false);
  e.model.renameTable('employees', 'employees');
  e.model.setSchemaFriendlyName('Northwind sales');
  e.model.setDescription('Sales data');
  e.model.setWritable(false);
  e.access.setGrant(ORDERS, 'g-sales', 'edit', true);
  e.access.setVisibility('orders', 'shipname', null);
  e.access.setVisibility('orders', 'freight', [SALES, ME]);
  plan = e.editPlan();
  assert.deepEqual(plan.payload, TOKENS);
  assert.deepEqual(plan.changes, []);
  assert.throws(() => e.plan(), /editPlan/);
  e.dispose();
});

test('registered names are locked by baseline membership, after uncheck and recheck too; candidates rename; the identifier is locked', () => {
  const m = model();
  assert.match(m.renameTable('orders', 'sales_orders'), /registered/);
  m.includeTable('orders', false);
  m.includeTable('orders', true);
  assert.match(m.renameTable('orders', 'sales_orders'), /registered/);
  assert.match(m.checkColumnName('orders', 'freight', 'cost'), /registered/);
  assert.equal(m.renameTable('orders', 'orders'), null, 'the name it has passes');
  assert.equal(m.table('orders').logical, 'orders');
  assert.equal(m.renameTable('employees', 'staff'), null);
  assert.equal(m.renameColumn('orders', 'tracking_number', 'tracking'), null, 'a candidate column on a registered table');
  assert.match(m.renameTable('customers', 'staff'), /already named/);
  m.setSchemaName('other');
  assert.equal(m.name.value, 'northwind_sales');
  m.setSchemaFriendlyName('Sales');
  assert.equal(m.name.value, 'northwind_sales', 'the friendly name does not drive a registered identifier');
  assert.equal(m.friendlyName.value, 'Sales');
});

test('the draft supplies candidates, unchecked: tables and supported columns the baseline does not carry', () => {
  const d = draft();
  d.manifest.tables.orders.columns.tracking_number.isName = true;
  d.manifest.tables.orders.columns.tracking_number.searchable = true;
  const m = model(d);
  const tracking = column(m, 'orders', 'tracking_number');
  assert.deepEqual([tracking.included, tracking.registered, tracking.supported, tracking.type], [false, false, true, 'string']);
  assert.deepEqual([tracking.isName, tracking.searchable], [false, false], 'the registered table has its name column');
  assert.equal(column(m, 'order_details', 'discount').included, false);
  assert.equal(m.toJSON().tables.orders.columns.tracking_number, undefined);
  m.includeColumn('orders', 'tracking_number', true);
  m.includeTable('employees', true);
  const json = m.toJSON();
  assert.deepEqual(json.tables.orders.columns.tracking_number, {type: 'string'});
  assert.deepEqual(json.tables.employees, {businessKey: ['employeeid'],
    columns: {employeeid: {type: 'int', required: true}, lastname: {type: 'string'}}});
  assert.deepEqual(ids(m.changes()), ['column:orders.tracking_number', 'table:employees'], 'in manifest order');
  assert.deepEqual(m.changes().map((c) => c.text), ['Column orders.tracking_number added', 'Table employees added']);
  const patch = m.patch();
  assert.deepEqual(Object.keys(patch.tables), ['orders', 'employees'], 'changed tables travel whole; order_lines does not');
  assert.equal(patch.dropTables, undefined);
  assert.equal(m.table('audit_log').bindable, false, 'an unbindable one stays shown');
});

test('drift on registered items: missing remotely, a changed type, a changed key — each with its reason; blockers gate Validate', () => {
  const m = model();
  assert.equal(m.catalog, 'known');
  const shippers = m.table('shippers');
  assert.deepEqual([shippers.drift.kind, shippers.drift.blocks, shippers.included], ['missing', true, true]);
  assert.match(shippers.reason, /missing remotely/);
  const address = column(m, 'orders', 'shipaddress');
  assert.deepEqual([address.drift.kind, address.drift.blocks], ['missing', true]);
  const city = column(m, 'orders', 'shipcity');
  assert.deepEqual([city.drift.kind, city.drift.blocks], ['type', true]);
  assert.match(city.drift.reason, /declared string; the warehouse column is int4 \(int\)/);
  const freight = column(m, 'orders', 'freight');
  assert.deepEqual([freight.drift.kind, freight.drift.blocks], ['type', false], 'float over int reads under a read-only binding');
  assert.match(freight.drift.reason, /read as float/);
  m.setWritable(true);
  assert.equal(column(m, 'orders', 'freight').drift.blocks, true, 'a float written over an int is lossy');
  assert.equal(m.blockers.value.some((b) => b.path === 'tables.orders.columns.freight'), true);
  m.setWritable(false);
  assert.equal(column(m, 'orders', 'freight').drift.blocks, false);
  assert.equal(column(m, 'orders', 'orderid').drift, undefined);
  assert.equal(column(m, 'order_details', 'orderid').drift, undefined, 'a ref is judged by the key type it points at');
  const lines = m.table('order_details');
  assert.deepEqual([lines.drift.kind, lines.drift.blocks], ['key', false]);
  assert.match(lines.drift.reason, /the warehouse key is orderid; registered with orderid, productid.*external-key-mismatch/);
  assert.equal(m.table('employees').drift, undefined, 'a candidate has nothing to drift from');

  const coded = draft();
  coded.inventory.tables.push({remote: 'shippers', bindable: false, code: 'external-table-missing', message: 'gone'});
  coded.inventory.columns.find((c) => c.remote === 'shipcity').code = 'external-param-collision';
  const byCode = model(coded);
  assert.equal(byCode.table('shippers').drift.kind, 'missing', 'the draft\'s own "missing" code');
  assert.deepEqual([column(byCode, 'orders', 'shipcity').drift.kind, column(byCode, 'orders', 'shipcity').drift.blocks], ['type', true],
    'another code is no type verdict: the declared type is judged against the mapped one');

  assert.deepEqual(m.blockers.value.map((b) => [b.path, b.code]), [
    ['tables.orders.columns.shipaddress', 'external-column-missing'],
    ['tables.orders.columns.shipcity', 'external-column-type'],
    ['tables.shippers', 'external-table-missing']]);
  m.includeTable('shippers', false);
  m.includeColumn('orders', 'shipaddress', false);
  assert.deepEqual(m.blockers.value.map((b) => b.path), ['tables.orders.columns.shipcity'], 'removed explicitly: no longer kept');
  const patch = m.patch();
  assert.deepEqual(patch.dropTables, ['shippers']);
  assert.equal(patch.tables.orders.columns.shipaddress, undefined);
  assert.deepEqual(m.changes().map((c) => [c.id, c.removes]), [['column:orders.shipaddress', true], ['table:shippers', true]]);
});

test('a failed or partial catalog read is "unknown": nothing missing, no blockers, no candidates', () => {
  const none = new ManifestModel(null, {baseline: baseline()});
  assert.equal(none.catalog, 'unknown');
  assert.equal(none.tables.value.length, 3);
  assert.equal(none.table('shippers').drift.kind, 'unknown');
  assert.equal(column(none, 'orders', 'shipaddress').drift.kind, 'unknown');
  assert.deepEqual(none.blockers.value, []);
  assert.deepEqual(none.patch(), {});
  const unreachable = draft();
  unreachable.diagnostics = [{code: 'external-unreachable', message: 'Connection could not be introspected'}];
  const partial = model(unreachable);
  assert.equal(partial.catalog, 'unknown');
  assert.equal(partial.table('shippers').drift.kind, 'unknown', 'never "missing" off a read that did not answer');
  const foreign = draft();
  foreign.inventory.tables = foreign.inventory.tables.filter((t) => !['orders', 'order_details'].includes(t.remote));
  const listsNone = model(foreign);
  assert.equal(listsNone.catalog, 'unknown', 'a read that lists none of the registered tables is a grants problem, not a warehouse that lost them all');
  assert.equal(listsNone.table('orders').drift.kind, 'unknown');
  const relationsOnly = draft();
  relationsOnly.inventory.relations = [];
  relationsOnly.diagnostics = [{code: 'external-relations-unavailable', message: 'no foreign keys'}];
  const noRefs = model(relationsOnly);
  assert.equal(noRefs.catalog, 'known', 'the tables and columns were read');
  assert.equal(noRefs.table('shippers').drift.kind, 'missing');
  assert.equal(column(noRefs, 'order_details', 'orderid').type, 'ref', 'a declared ref needs no warehouse foreign key');
});

test('declared refs are authoritative: excluding the target demotes and says so, including it restores; a new foreign key is a suggestion', () => {
  const m = model();
  m.includeTable('orders', false);
  const line = column(m, 'order_details', 'orderid');
  assert.deepEqual([line.type, line.relations[0].canFix, line.relations[0].ref], ['int', true, false]);
  assert.deepEqual(m.changes().map((c) => c.text),
    ['Table orders removed', 'Column order_lines.order_id: ref to orders demoted to a plain int']);
  m.includeTable('orders', true);
  assert.equal(column(m, 'order_details', 'orderid').type, 'ref');
  assert.deepEqual(m.changes(), []);

  let fk = m.relationOf('orders', 'employeeid');
  assert.deepEqual([fk.suggested, fk.canFix, fk.ref, column(m, 'orders', 'employeeid').type], [false, false, false, 'int'],
    'the target is not included: nothing to suggest yet');
  m.includeTable('employees', true);
  fk = m.relationOf('orders', 'employeeid');
  assert.deepEqual([fk.suggested, fk.ref], [true, false], 'including the target promotes nothing by itself');
  assert.match(fk.reason, /the registered manifest does not declare/);
  m.setRef('orders', 'employeeid', true);
  fk = m.relationOf('orders', 'employeeid');
  assert.deepEqual([fk.ref, fk.promoted, fk.suggested], [true, true, false]);
  assert.deepEqual(m.toJSON().tables.orders.columns.employeeid, {type: 'ref', ref: 'employees'});
  assert.deepEqual(ids(m.changes()), ['column:orders.employeeid:ref', 'table:employees']);
  assert.equal(m.changes()[0].text, 'Column orders.employeeid: now a ref to employees');
  m.includeTable('employees', false);
  assert.deepEqual([column(m, 'orders', 'employeeid').type, m.relationOf('orders', 'employeeid').canFix], ['int', true]);
  m.includeTable('employees', true);
  assert.equal(column(m, 'orders', 'employeeid').type, 'ref');
  m.setRef('orders', 'employeeid', false);
  assert.deepEqual([column(m, 'orders', 'employeeid').type, m.relationOf('orders', 'employeeid').suggested], ['int', true]);
  m.setRef('order_details', 'orderid', false);
  assert.equal(column(m, 'order_details', 'orderid').type, 'ref', 'a declared ref cannot be demoted by hand');
});

test('a read-only opt-out stays dormant across writable off and on; unexposed keys round-trip in the whole descriptor', () => {
  const e = editor();
  const lines = () => e.model.toJSON().tables.order_lines;
  assert.deepEqual([e.model.writable.value, lines().writable, lines().singularName], [false, false, 'order line']);
  e.model.setWritable(true);
  assert.equal(lines().writable, false);
  assert.equal(e.model.writes(e.model.table('order_details')), false);
  let plan = e.editPlan();
  assert.deepEqual(plan.payload, {...TOKENS, storage: {writable: true}});
  assert.deepEqual(plan.changes, [{id: 'schema:writable', text: 'Writable: on', removes: false}]);
  e.model.setWritable(false);
  assert.equal(lines().writable, false);
  assert.deepEqual(e.editPlan().payload, TOKENS);
  e.model.setFriendlyName('order_details', 'Lines');
  e.model.setSchemaFriendlyName('Northwind');
  e.model.setDescription('EU sales');
  plan = e.editPlan();
  assert.deepEqual(plan.payload, {...TOKENS, friendlyName: 'Northwind', description: 'EU sales', tables: {order_lines: {
    friendlyName: 'Lines', singularName: 'order line', writable: false, table: 'order_details',
    businessKey: ['order_id', 'productid'],
    columns: {order_id: {type: 'ref', ref: 'orders', required: true, column: 'orderid'}, productid: {type: 'int', required: true},
      unitprice: {type: 'float'}, quantity: {type: 'int'}},
  }}});
  assert.deepEqual(ids(plan.changes), ['schema:friendlyName', 'schema:description', 'table:order_lines:friendlyName']);
  e.dispose();
});

test('access from the snapshot: per-table direct grants with complete sets, column states, read-only targets; the delta is exact', () => {
  const e = editor();
  const access = e.access;
  assert.equal(access.editing, true);
  assert.deepEqual(access.grantsOf(ORDERS).value, [
    {scope: ORDERS, group: ME, view: true, edit: true, delete: true, other: ['Share']},
    {scope: ORDERS, group: SALES, view: true, edit: true, delete: false}]);
  assert.deepEqual(access.grantsOf({kind: 'schema'}).value, [], 'no schema-scope rows are reconstructed');
  assert.deepEqual(access.grantsOf({kind: 'table', table: 'order_details'}).value, [], 'unreadable: none, not empty');
  assert.deepEqual([access.canEdit(ORDERS), access.canEdit({kind: 'table', table: 'order_details'}),
    access.canEdit({kind: 'table', table: 'employees'})], [true, false, true], 'a table this edit adds is editable');
  assert.deepEqual([access.canEditColumn('orders', 'freight'), access.canEditColumn('order_details', 'unitprice'),
    access.canEditColumn('order_details', 'discount')], [true, false, false]);
  assert.deepEqual(access.entityOf('orders'), 't-orders');
  assert.deepEqual(access.columnOf('orders', 'freight'), {table: 'orders', column: 'freight', groups: [SALES, ME], edit: [ME],
    other: [{group: ME, permission: 'Share'}]});
  assert.deepEqual(access.columnOf('order_details', 'unitprice'), {table: 'order_details', column: 'unitprice', groups: [], edit: [], unknown: true});
  assert.equal(access.columnOf('orders', 'shipname'), undefined, 'unrestricted');
  assert.deepEqual(access.delta(), {grant: [], revoke: [], restrict: [], unrestrict: []});
  assert.deepEqual(new AccessModel(access.toJSON()).toJSON(), access.toJSON());

  access.setGrant(ORDERS, 'g-sales', 'edit', false);
  assert.deepEqual(bare(access.delta().revoke), [{table: 'orders', group: SALES, permission: 'Edit'}], 'one permission; View stays');
  let plan = e.editPlan();
  assert.deepEqual(plan.payload.access, {grant: [], revoke: [{table: 'orders', group: 'g-sales', permission: 'Edit'}],
    restrict: [], unrestrict: []});
  assert.deepEqual(plan.changes, [{id: 'access:revoke:orders:g-sales:Edit', table: 'orders', removes: false,
    text: 'orders: Edit revoked from Sales'}]);
  access.setGrant(ORDERS, 'g-me', 'delete', false);
  assert.deepEqual(access.grantsOf(ORDERS).value[0].other, ['Share'], 'Share is never touched');
  access.removeGroup(ORDERS, 'g-sales');
  assert.deepEqual(access.delta().revoke.map((r) => `${r.group.id}:${r.permission}`), ['g-me:Delete', 'g-sales:View', 'g-sales:Edit'],
    'a removed row revokes the permissions it held, one triple each');
  e.dispose();
});

test('the "every table" bulk action writes per-table rows; a read-only table yields none; a dropped table\'s ops are omitted', () => {
  const e = editor();
  const access = e.access;
  const tables = ['orders', 'order_details', 'shippers'];
  assert.deepEqual(access.everyTable(tables), [{group: SALES, view: true, edit: false, delete: false}],
    'held alike on the tables that can be edited: View on both, Edit on orders alone');
  access.setEveryTable(tables, [...access.everyTable(tables), {group: DEV, view: true, edit: false, delete: false}]);
  assert.deepEqual(access.grantsOf(ORDERS).value.map((g) => [g.group.id, g.view, g.edit]), [['g-me', true, true], ['g-sales', true, true], ['g-dev', true, false]],
    'the bulk row does not narrow a table\'s own wider grant');
  assert.deepEqual(access.grantsOf({kind: 'table', table: 'shippers'}).value.map((g) => g.group.id), ['g-sales', 'g-dev']);
  assert.deepEqual(access.grantsOf({kind: 'table', table: 'order_details'}).value, []);
  assert.deepEqual(bare(access.delta().grant), [{table: 'orders', group: DEV, permission: 'View'}, {table: 'shippers', group: DEV, permission: 'View'}]);
  assert.deepEqual(access.everyTable(tables).map((g) => g.group.id), ['g-sales', 'g-dev']);
  e.model.includeTable('shippers', false);
  const plan = e.editPlan();
  assert.deepEqual(plan.payload.dropTables, ['shippers']);
  assert.deepEqual(plan.payload.access.grant, [{table: 'orders', group: 'g-dev', permission: 'View'}], 'the purge owns the dropped table');
  e.model.includeTable('shippers', true);
  access.setEveryTable(tables, [{group: DEV, view: true, edit: false, delete: false}]);
  assert.deepEqual(access.grantsOf(ORDERS).value.map((g) => [g.group.id, g.view, g.edit]), [['g-me', true, true], ['g-sales', false, true], ['g-dev', true, false]],
    'Sales gone from the bulk view loses on orders only the View the view showed; its Edit there stays');
  assert.deepEqual(access.grantsOf({kind: 'table', table: 'shippers'}).value.map((g) => g.group.id), ['g-dev'], 'nothing left: the row goes');
  assert.deepEqual(e.editPlan().payload.access.revoke, [{table: 'orders', group: 'g-sales', permission: 'View'},
    {table: 'shippers', group: 'g-sales', permission: 'View'}]);
  e.dispose();
});

test('column access: restrict with the groups let in (Edit where they may edit the table), narrow by exact triples, unrestrict; keys and read-only columns yield nothing', () => {
  const e = editor();
  const access = e.access;
  access.setVisibility('orders', 'shipname', [SALES]);
  assert.deepEqual(access.columnOf('orders', 'shipname'), {table: 'orders', column: 'shipname', groups: [SALES], edit: [SALES]});
  assert.deepEqual(bare(access.delta().restrict), [{table: 'orders', column: 'shipname',
    grant: [{group: SALES, permission: 'View'}, {group: SALES, permission: 'Edit'}], revoke: []}]);
  access.setVisibility('orders', 'freight', [ME]);
  assert.deepEqual(access.columnOf('orders', 'freight').edit, [ME]);
  access.setVisibility('orders', 'shipcountry', []);
  access.setVisibility('orders', 'orderid', [SALES]);
  access.setVisibility('order_details', 'unitprice', [SALES]);
  let plan = e.editPlan();
  assert.deepEqual(plan.payload.access, {grant: [], revoke: [], unrestrict: [], restrict: [
    {table: 'orders', column: 'shipname', grant: [{group: 'g-sales', permission: 'View'}, {group: 'g-sales', permission: 'Edit'}], revoke: []},
    {table: 'orders', column: 'freight', grant: [], revoke: [{group: 'g-sales', permission: 'View'}]},
    {table: 'orders', column: 'shipcountry', grant: [], revoke: []},
  ]}, 'a key column and a column the caller cannot share are left out');
  assert.deepEqual(plan.changes.map((c) => c.text), ['orders.shipname restricted — visible to Sales; Edit for Sales',
    'orders.freight: View revoked from Sales', 'orders.shipcountry restricted — visible to nobody else']);
  access.setVisibility('orders', 'freight', null);
  plan = e.editPlan();
  assert.deepEqual(plan.payload.access.unrestrict, [{table: 'orders', column: 'freight'}]);
  assert.equal(plan.changes.find((c) => c.id === 'access:unrestrict:orders.freight').text, 'orders.freight visible to everyone again');
  e.model.includeColumn('orders', 'shipname', false);
  assert.equal(e.editPlan().payload.access.restrict.some((r) => r.column === 'shipname'), false, 'a column that is out gets no op');
  e.dispose();

  const mine = editor({author: ME});
  mine.access.setVisibility('orders', 'shipname', [SALES]);
  mine.access.setVisibility('orders', 'shipcountry', [ME]);
  mine.access.setVisibility('orders', 'freight', [SALES]);
  const kept = mine.editPlan();
  assert.deepEqual(kept.payload.access.restrict, [
    {table: 'orders', column: 'shipname', grant: [{group: 'g-sales', permission: 'View'}, {group: 'g-sales', permission: 'Edit'},
      {group: 'g-me', permission: 'View'}, {group: 'g-me', permission: 'Edit'}], revoke: []},
    {table: 'orders', column: 'shipcountry', grant: [{group: 'g-me', permission: 'View'}, {group: 'g-me', permission: 'Edit'}], revoke: []},
    {table: 'orders', column: 'freight', grant: [], revoke: [{group: 'g-me', permission: 'View'}, {group: 'g-me', permission: 'Edit'}]},
  ], 'the author keeps a first-restricted column explicitly unless already on it; narrowing an existing restriction adds nobody');
  assert.deepEqual(kept.changes.map((c) => c.text), ['orders.shipname restricted — visible to Sales and you; Edit for Sales and you',
    'orders.shipcountry restricted — visible to askalkin; Edit for askalkin', 'orders.freight: View revoked from askalkin, Edit revoked from askalkin']);
  mine.dispose();

  const blind = editor({snapshot: null});
  assert.deepEqual([blind.access.editing, blind.access.canEdit(ORDERS), blind.access.canEditColumn('orders', 'shipname')], [true, false, false]);
  blind.access.addGroup(ORDERS, DEV);
  blind.access.setVisibility('orders', 'shipname', [DEV]);
  assert.deepEqual(blind.editPlan().payload, TOKENS, 'a snapshot that could not be read changes nothing');
  blind.dispose();
});

test('rebase: disjoint edits kept on both sides, a same-field change is a conflict, a vanished item is dropped, absorbed edits count as applied', async () => {
  const e = editor();
  document.body.append(e.root);
  e.tree.tree.root.querySelector('.u2-list').clientHeight = 800;
  await flush();
  e.model.setRequired('orders', 'freight', true);
  e.model.setFriendlyName('orders', 'Sales orders');
  e.model.setRequired('shippers', 'companyname', true);
  e.model.includeTable('employees', true);
  e.model.setSchemaFriendlyName('NW');
  e.access.setGrant(ORDERS, 'g-sales', 'edit', false);
  e.access.addGroup(ORDERS, DEV);
  e.access.addGroup({kind: 'table', table: 'shippers'}, DEV);
  e.access.setVisibility('orders', 'shipname', [SALES]);
  assert.equal(e.stale, false);

  const nb = baseline();
  nb.version = '4';
  nb.tables.orders.friendlyName = 'Orders (EU)';
  delete nb.tables.shippers;
  const ns = snapshot();
  ns.version = '4';
  delete ns.tables.shippers;
  for (const key of Object.keys(ns.columns).filter((k) => k.startsWith('shippers.')))
    delete ns.columns[key];
  ns.tables.orders.grants[1].permissions = ['View', 'Delete'];
  ns.tables.orders.coreSchema.canShare = false;
  for (const key of Object.keys(ns.columns).filter((k) => k.startsWith('orders.')))
    ns.columns[key].canShare = false;
  Object.assign(ns.columns['orders.freight'], {view: null, edit: null, other: null});
  // shippers is unregistered on the server but still in the warehouse: a candidate again
  const nd = draft();
  nd.manifest.tables.shippers = baseline().tables.shippers;
  nd.inventory.tables.push({remote: 'shippers', logical: 'shippers', bindable: true, key: ['shipperid']});
  nd.inventory.columns.push({table: 'shippers', remote: 'shipperid', dbType: 'int4', type: 'int'},
    {table: 'shippers', remote: 'companyname', dbType: 'text', type: 'string'});
  const report = e.rebase(nb, ns, nd, {friendlyName: 'NW', description: 'Sales data'});
  assert.deepEqual(ids(report.applied), ['table:employees', 'schema:friendlyName', 'column:orders.freight:required',
    'access:revoke:orders:g-sales:Edit', 'access:grant:orders:g-dev:View']);
  assert.deepEqual(report.conflicts, [{id: 'table:orders:friendlyName', table: 'orders', removes: false,
    text: 'Table orders: friendly name "Sales orders"', from: 'Orders', to: 'Sales orders', server: 'Orders (EU)'}]);
  assert.deepEqual(report.dropped.map((d) => [d.id, d.to]), [['column:shippers.companyname:required', true],
    ['access:grant:shippers:g-dev:View', true], ['access:restrict:orders.shipname', ['g-sales']]],
    'an edit on a table that is a candidate again, and one on columns the caller can no longer share');
  assert.equal(e.model.table('shippers').registered, false);
  assert.equal(e.model.table('orders').friendlyName, 'Orders (EU)', 'the server\'s value stands on a conflict');
  assert.deepEqual([e.model.version, e.model.friendlyName.value, e.stale], ['4', 'NW', false]);
  assert.equal(e.access.columnOf('orders', 'shipname'), undefined);
  const plan = e.editPlan();
  assert.deepEqual(plan.payload, {ifVersion: '4', ifIncarnation: TOKENS.ifIncarnation, tables: {
    orders: plan.payload.tables.orders, employees: plan.payload.tables.employees},
  access: {grant: [{table: 'orders', group: 'g-dev', permission: 'View'}], revoke: [], restrict: [], unrestrict: []}});
  assert.equal(plan.payload.tables.orders.columns.freight.required, true);
  assert.deepEqual(e.access.grantsOf(ORDERS).value.map((g) => [g.group.id, g.view, g.edit, g.delete]),
    [['g-me', true, true, true], ['g-sales', true, false, true], ['g-dev', true, false, false]]);
  await flush();
  const names = e.tree.root.querySelectorAll('.u2-tree-row').map((r) => r.querySelector('.u2-manifest-node-name').textContent);
  assert.equal(names.includes('employees'), true, 'the tree follows the new baseline');
  assert.deepEqual(badges(rowOf(e, 'shippers')), ['new']);
  e.dispose();
  resetDom();
  await flush();

  const later = snapshot();
  later.version = '4';
  const behind = editor({snapshot: later});
  assert.equal(behind.stale, true, 'the snapshot was read at another version than the manifest');
  behind.dispose();
});

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

async function shown(options = {}) {
  const e = editor(options);
  document.body.append(e.root);
  e.tree.tree.root.querySelector('.u2-list').clientHeight = 800;
  await flush();
  return e;
}

const rows = (e) => e.tree.root.querySelectorAll('.u2-tree-row');
const rowOf = (e, name) => rows(e).find((r) => r.querySelector('.u2-manifest-node-name').textContent === name);
const badges = (row) => row.querySelectorAll('.u2-badge').map((b) => b.textContent);
const panel = (e) => e.panel.root;
const inputNames = (e) => panel(e).querySelectorAll('.u2-input-root').map((i) => i.dataset.u2Name);
const readonlyNames = (e) => panel(e).querySelectorAll('[data-u2="readonly-field"]').map((f) => f.dataset.u2Name);
const notes = (e) => panel(e).querySelectorAll('.u2-manifest-panel-note').map((n) => n.textContent);

async function select(e, name) {
  fire(rowOf(e, name), 'click');
  await flush();
}

ui('edit mode panels: the identifier as text, a description, no creator row, the bulk grid over every table', async () => {
  const e = await shown();
  assert.deepEqual(inputNames(e), ['friendlyName', 'description', 'writable', 'access-schema']);
  assert.deepEqual(readonlyNames(e), ['name']);
  const grid = panel(e).querySelector('.u2-access-grid');
  assert.deepEqual(grid.querySelectorAll('tbody .u2-access-grid-principal').map((p) => p.textContent), ['Sales'],
    'the groups every shareable table grants alike; no "You (creator)" row');
  assert.match(notes(e).join(' '), /written on each of them/);
  grid.querySelectorAll('tbody .u2-access-grid-check')[0].click();
  await flush();
  assert.deepEqual(bare(e.access.delta().revoke), [{table: 'orders', group: SALES, permission: 'View'},
    {table: 'shippers', group: SALES, permission: 'View'}], 'a bulk change lands on each table');
  e.dispose();
});

ui('edit mode rows and panels: drift badges with their reasons, "new" candidates, locked names, read-only access, the ref suggestion', async () => {
  const e = await shown();
  assert.deepEqual(badges(rowOf(e, 'shippers')), ['shipperid', 'missing remotely']);
  assert.deepEqual(badges(rowOf(e, 'shipaddress')), ['missing remotely']);
  assert.deepEqual(badges(rowOf(e, 'shipcity')), ['type changed']);
  assert.deepEqual(badges(rowOf(e, 'employees')), ['new']);
  assert.deepEqual(badges(rowOf(e, 'tracking_number')), ['new']);
  assert.deepEqual(badges(rowOf(e, 'order_details')), ['orderid, productid', 'key changed']);
  assert.match(rowOf(e, 'shippers').title, /missing remotely/);
  assert.deepEqual([rowOf(e, 'shippers').querySelector('.u2-tree-check').checked, rowOf(e, 'shippers').querySelector('.u2-tree-check').disabled],
    [true, false], 'stays checked; removable explicitly');

  await select(e, 'orders');
  assert.deepEqual(readonlyNames(e), ['logical', 'key']);
  assert.match(panel(e).querySelector('[data-u2-name="logical"]').textContent, /for life/);
  assert.deepEqual(panel(e).querySelectorAll('.u2-access-grid tbody .u2-access-grid-principal').map((p) => p.textContent), ['askalkin', 'Sales']);
  assert.match(notes(e).join(' '), /askalkin also holds Share/);
  await select(e, 'shippers');
  assert.deepEqual(panel(e).querySelectorAll('.u2-manifest-panel-title .u2-badge').map((b) => b.textContent), ['missing remotely']);
  await select(e, 'employees');
  assert.equal(inputNames(e).includes('logical'), true, 'a candidate renames');
  await select(e, 'order_details');
  assert.equal(panel(e).querySelector('.u2-access-grid-add').disabled, true);
  assert.match(notes(e).join(' '), /cannot share this table/);
  await select(e, 'shipcity');
  assert.match(panel(e).querySelector('[data-u2-name="type"]').textContent, /declared string/);
  await e.select({kind: 'column', table: 'order_details', column: 'unitprice'});
  await flush();
  assert.match(notes(e).join(' '), /cannot be read here/);
  assert.equal(inputNames(e).includes('visibility'), false);

  await select(e, 'employeeid');
  assert.equal(panel(e).querySelector('.u2-manifest-relation .u2-link'), null, 'employees is not included: nothing to suggest');
  fire(rowOf(e, 'employees').querySelector('.u2-tree-check'), 'click');
  await flush();
  const suggestion = panel(e).querySelector('.u2-manifest-relation .u2-link');
  assert.equal(suggestion.textContent, 'make it a ref');
  fire(suggestion, 'click');
  await flush();
  assert.equal(e.model.column('orders', 'employeeid').type, 'ref');
  assert.equal(panel(e).querySelector('.u2-manifest-relation .u2-link').textContent, 'keep it a plain value');
  e.dispose();
});
