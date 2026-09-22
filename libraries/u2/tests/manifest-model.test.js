/* ManifestModel over the Northwind draft fixture: the untouched draft round-trips, including
   and excluding a table flips the dependent relation between ref and plain, renames are checked
   with the server's rules and refused names keep the old one, one name column per table, an
   excluded column is dropped, key columns cannot go, writable and the read-only opt-out, the
   inventory's unbindable tables and unsupported columns are shown and never submitted; the
   access model round-trips; the field offer's arms. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {ManifestModel, AccessModel} from '../src/dg/domain/authoring/manifest-model.js';
import {ManifestRules} from '../src/dg/domain/authoring/manifest-rules.js';
import {fieldOffer} from '../src/dg/domain/authoring/editor-context.js';

const DRAFT = JSON.parse(readFileSync(new URL('./fixtures/authoring/northwind-draft.json', import.meta.url), 'utf8'));

/** A fresh copy per test: the model keeps the declarations it was handed. */
const draft = () => JSON.parse(JSON.stringify(DRAFT));
const column = (model, table, remote) => model.columns(table).value.find((c) => c.remote === remote);

test('toJSON of the untouched draft is the draft manifest', () => {
  const model = new ManifestModel(draft());
  assert.deepEqual(model.toJSON(), DRAFT.manifest);
  assert.deepEqual(model.tables.value.map((t) => [t.remote, t.included, t.bindable]),
    [['orders', true, true], ['order_details', true, true], ['order_summary', false, false], ['audit_log', false, false]]);
  assert.equal(model.tables.value[2].code, 'external-table-view');
  assert.match(model.tables.value[3].reason, /no primary key/);
});

test('excluding a table flips the dependent relation to plain; including it back makes the ref again', () => {
  const model = new ManifestModel(draft());
  const before = column(model, 'order_details', 'orderid');
  assert.equal(before.type, 'ref');
  assert.equal(before.relation.ref, true);
  assert.equal(before.relation.reason, 'ref: orders');

  model.includeTable('orders', false);
  const plain = column(model, 'order_details', 'orderid');
  assert.equal(plain.type, 'int', 'the scalar type of the key it pointed at');
  assert.equal(plain.relation.ref, false);
  assert.equal(plain.relation.canFix, true);
  assert.equal(plain.relation.reason, 'target not included — stays a plain value');
  const json = model.toJSON();
  assert.equal(json.tables.orders, undefined);
  assert.deepEqual(json.tables.order_details.columns.orderid, {type: 'int', required: true});

  model.includeTable('orders', true);
  assert.deepEqual(model.toJSON(), DRAFT.manifest);
  // a relation the warehouse reported but the draft did not qualify stays plain whatever is included
  const customer = column(model, 'orders', 'customerid');
  assert.equal(customer.type, 'string');
  assert.equal(customer.relation.canFix, false);
  assert.equal(customer.relation.reason, 'customers is not in the draft');
});

test('a rename is checked with the mirrored rules; a refused name keeps the old one', () => {
  const model = new ManifestModel(draft());
  assert.match(model.renameTable('orders', 'Orders'), /lowercase/);
  assert.match(model.renameTable('orders', 'transaction'), /reserved/);
  assert.match(model.renameTable('orders', 'order_details'), /already named/);
  assert.equal(model.table('orders').logical, 'orders');
  assert.match(model.renameColumn('orders', 'shipvia', 'id'), /system column/);
  assert.match(model.renameColumn('orders', 'shipvia', 'x_ship'), /reserved prefix/);
  assert.match(model.renameColumn('orders', 'shipvia', 'freight'), /already named/);
  assert.match(model.renameColumn('orders', 'shipvia', ''), /required/);
  assert.equal(column(model, 'orders', 'shipvia').logical, 'shipvia');
  assert.match(model.checkSchemaName('schemas'), /reserved/);
  assert.equal(model.checkSchemaName('northwind_sales'), null);
  assert.equal(ManifestRules.IDENTIFIER.source, '^[a-z][a-z0-9_]*$');

  assert.equal(model.renameTable('orders', 'order'), null);
  assert.equal(model.renameColumn('orders', 'shipvia', 'ship_via'), null);
  assert.equal(model.renameColumn('order_details', 'orderid', 'order_id'), null);
  const json = model.toJSON();
  assert.deepEqual(Object.keys(json.tables), ['order', 'order_details']);
  assert.equal(json.tables.order.table, 'orders', 'the remote name follows the rename');
  assert.deepEqual(json.tables.order.columns.ship_via, {type: 'int', column: 'shipvia'});
  assert.deepEqual(json.tables.order_details.columns.order_id, {type: 'ref', ref: 'order', required: true, column: 'orderid'});
  assert.deepEqual(json.tables.order_details.businessKey, ['order_id', 'productid'], 'the key follows the rename');
  assert.equal(model.renameColumn('orders', 'shipcountry', 'ship_country'), null);
  assert.deepEqual(model.toJSON().tables.order.filters, [{column: 'ship_country'}], 'a preserved filter follows the rename');
  model.includeColumn('orders', 'shipcountry', false);
  assert.equal(model.toJSON().tables.order.filters, undefined, 'and goes with the excluded column');
});

test('one name column and at most one searchable per table, strings only', () => {
  const model = new ManifestModel(draft());
  model.setNameColumn('orders', 'shipcountry');
  assert.deepEqual(model.columns('orders').value.filter((c) => c.isName).map((c) => c.remote), ['shipcountry']);
  model.setNameColumn('orders', 'freight');
  assert.deepEqual(model.columns('orders').value.filter((c) => c.isName).map((c) => c.remote), ['shipcountry'],
    'a float cannot name a row');
  model.setSearchable('orders', null);
  assert.equal(model.columns('orders').value.some((c) => c.searchable), false);
  const json = model.toJSON();
  assert.equal(json.tables.orders.columns.shipname.isName, undefined);
  assert.deepEqual(json.tables.orders.columns.shipcountry, {type: 'string', isName: true});
});

test('an excluded column is dropped; a key column stays; required, writable and the read-only opt-out', () => {
  const model = new ManifestModel(draft());
  model.includeColumn('orders', 'shipaddress', false);
  model.includeColumn('orders', 'orderid', false);
  model.setRequired('orders', 'shipname', true);
  model.setRequired('orders', 'orderid', false);
  let json = model.toJSON();
  assert.equal(json.tables.orders.columns.shipaddress, undefined);
  assert.deepEqual(json.tables.orders.columns.orderid, {type: 'int', required: true});
  assert.equal(json.tables.orders.columns.shipname.required, true);
  assert.equal(column(model, 'orders', 'shipaddress').included, false, 'still shown, unchecked');

  model.setReadOnly('order_details', true);
  assert.equal(model.toJSON().tables.order_details.writable, undefined, 'no opt-out under a read-only storage');
  model.setWritable(true);
  json = model.toJSON();
  assert.equal(json.storage.writable, true);
  assert.equal(json.tables.order_details.writable, false);
  model.setWritable(false);
  assert.deepEqual(model.toJSON().storage, DRAFT.manifest.storage);
});

test('the inventory is shown, never submitted: unsupported columns, unbindable tables, warehouse types', () => {
  const model = new ManifestModel(draft());
  const tracking = column(model, 'orders', 'tracking_number');
  assert.equal(tracking.supported, false);
  assert.equal(tracking.dbType, 'bigint');
  assert.match(tracking.reason, /not supported/);
  model.includeColumn('orders', 'tracking_number', true);
  model.includeTable('audit_log', true);
  assert.equal(model.toJSON().tables.orders.columns.tracking_number, undefined);
  assert.equal(model.toJSON().tables.audit_log, undefined);
  assert.equal(model.table('audit_log').included, false);
  model.includeTables(false);
  assert.deepEqual(model.toJSON().tables, {});
  model.includeTables(true);
  assert.deepEqual(model.tables.value.map((t) => t.included), [true, true, false, false]);
});

test('resolvePath addresses a diagnostic to its node by the current logical names', () => {
  const model = new ManifestModel(draft());
  model.renameTable('orders', 'order');
  assert.deepEqual(model.resolvePath('tables.order.columns.shipvia'), {kind: 'column', table: 'orders', column: 'shipvia'});
  assert.deepEqual(model.resolvePath('tables.order.businessKey'), {kind: 'table', table: 'orders'});
  assert.deepEqual(model.resolvePath('tables.orders'), {kind: 'schema'}, 'the old name owns nothing');
  assert.deepEqual(model.resolvePath('storage.connection'), {kind: 'schema'});
  assert.deepEqual(model.resolvePath(undefined), {kind: 'schema'});
});

test('the access model round-trips and edits by scope', () => {
  const access = new AccessModel();
  access.addGroup('schema', 'Sales');
  access.addGroup('schema', 'Sales');
  access.addGroup('orders', 'Chemists');
  access.setGrant('orders', 'Chemists', 'edit', true);
  access.setVisibility('orders', 'freight', ['Sales']);
  access.setVisibility('orders', 'shipname', null);
  const json = access.toJSON();
  assert.deepEqual(json, {
    grants: [{scope: 'schema', group: 'Sales', view: true, edit: false, delete: false},
      {scope: 'orders', group: 'Chemists', view: true, edit: true, delete: false}],
    visibility: [{table: 'orders', column: 'freight', groups: ['Sales']}],
  });
  assert.deepEqual(new AccessModel(json).toJSON(), json);
  assert.deepEqual(access.grantsOf('orders').value.map((g) => g.group), ['Chemists']);
  access.removeGroup('schema', 'Sales');
  access.setGrants('orders', [{group: 'Developers', view: true, edit: false, delete: false}]);
  assert.deepEqual(access.grants.value.map((g) => [g.scope, g.group]), [['orders', 'Developers']]);
  assert.deepEqual(access.visibilityOf('orders', 'freight'), ['Sales']);
  assert.equal(access.visibilityOf('orders', 'shipname'), null);
});

test('the field offer: the external vocabulary under create and view; the other arms refuse by name', () => {
  const create = fieldOffer({mode: 'create', storage: 'external'});
  assert.equal(create.editable, true);
  assert.equal(create.type, 'readonly');
  assert.equal(create.refs, 'relations');
  assert.deepEqual([create.defaultValue, create.autoNumber, create.immutable, create.unique, create.choices],
    [false, false, false, false, false]);
  assert.equal(fieldOffer({mode: 'view', storage: 'external'}).editable, false);
  assert.throws(() => fieldOffer({mode: 'create', storage: 'domain'}), /not built yet/);
  assert.throws(() => fieldOffer({mode: 'edit', storage: 'external'}), /"edit" mode is not built yet/);
});
