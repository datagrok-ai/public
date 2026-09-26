/* The "Create domain schema" dialog over stubbed dapi calls: the draft is read once the
   connection and the schema are picked (both preset: the dialog opens on Design, the Connection
   step a BACK away), the Design step is built over it (one table included when the caller named
   one; the identifier proposed from the connection and the schema, clear of the registered
   names; the platform group picker in the access grids), the review shows the manifest, VALIDATE is the dry run bound to the exact
   payload it checked (a refusal lands on the rows and blocks CREATE; a payload changed since is
   not validated), CREATE is the real one followed by the column restrictions and then the table
   grants — a schema row is the same grant on every included table, never a schema grant — and
   the dialog ends on Created with what went through and what did not, RETRY ACCESS over the
   failed steps alone; the promise resolves to the name with the access outcome. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {notify} from '../src/components/display/notify.js';

register('./dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {domains} = await import('../src/dg/domain/index.js');

const DRAFT = JSON.parse(readFileSync(new URL('./fixtures/authoring/northwind-draft.json', import.meta.url), 'utf8'));
const CONN = {id: 'c1', nqName: 'NorthwindBinding:PostgresNorthwind', friendlyName: 'PostgresNorthwind',
  dataSource: 'Postgres', isDatabase: true};
const NOT_BINDABLE = [{id: 'c-s3', friendlyName: 'Bucket', dataSource: 'S3', isDatabase: false},
  {id: 'c-dom', friendlyName: 'Domains', dataSource: 'Domain', isDatabase: true}];
const SALES = {id: 'g-sales', label: 'Sales'};
const DEV = {id: 'g-dev', label: 'Developers'};
/** The groups the platform knows; a test takes one out to model a group deleted since the pick. */
const GROUPS = ['g-sales', 'g-dev', 'u-sam'];
const SCHEMA = {kind: 'schema'};
const ORDERS = {kind: 'table', table: 'orders'};

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** What the registry answers a probe by name: the manifest of a registered schema, a typed
 * not-found otherwise. */
const REGISTERED = new Map([['postgresnorthwind_public', {storage: {kind: 'external'}, tables: {}}]]);
const notFound = (name) => Object.assign(new Error(`Domain schema "${name}" not found`), {status: 404, code: 'not-found'});

/** Every dapi call the dialog makes, recorded; the dry run answers what `dryRun` says; a step
 * named in `failing` refuses with its message. */
function stub(calls, dryRun, failing = {}) {
  grok.shell.settings = {enableDomainDatabases: true};
  globalThis.grok_Dapi_Domains_SchemaCreated = (_dart, name) => calls.push(['announce', name]);
  globalThis.grok_Dapi_Domains_SchemaAltered = (_dart, name) => calls.push(['altered', name]);
  grok.dapi.connections = {list: async () => [CONN, ...NOT_BINDABLE], getSchemas: async () => ['public', 'audit']};
  grok.dapi.permissions = {check: async (_c, right) => right !== 'DataConnection.RemoveRows',
    checkGlobal: async () => true};
  grok.dapi.groups = {getGroupsLookup: async (query) => [{id: 'g-sales', friendlyName: 'Sales', personal: false},
    {id: 'u-sam', friendlyName: 'Sam Sales', personal: true}].filter((g) => g.friendlyName.toLowerCase().includes(query)),
  list: async ({filter}) => {
    calls.push(['groups', filter]);
    return GROUPS.filter((id) => filter.includes(`"${id}"`)).map((id) => ({id}));
  }};
  const refuse = (step) => {
    if (failing[step])
      throw new Error(failing[step]);
  };
  Object.assign(grok.dapi.domains, {
    schemas: {list: async () => {
      throw new Error('the registry is never listed');
    }},
    draft: async (body) => {
      calls.push(['draft', body]);
      return JSON.parse(JSON.stringify(DRAFT));
    },
    createSchema: async (name, options) => {
      calls.push(['create', name, options]);
      return options.dryRun ? dryRun() : {id: 'id-1', name, pgSchema: `ext_${name}`, version: '1'};
    },
    schema: (name) => ({
      grant: async (group, permission) => calls.push(['schema.grant', name, group, permission]),
      manifest: async () => {
        calls.push(['registry', name]);
        if (!REGISTERED.has(name))
          throw notFound(name);
        return REGISTERED.get(name);
      },
    }),
    table: (address) => ({
      grant: async (group, permission) => {
        calls.push(['table.grant', address, group, permission]);
        refuse('grant');
      },
      shareColumn: async (column, group, permission) => {
        calls.push(['shareColumn', address, column, group, permission]);
        refuse('shareColumn');
      },
      restrictColumn: async (column) => calls.push(['restrictColumn', address, column]),
    }),
  });
}

const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')]
  .find((b) => b.textContent === text);
const reason = () => document.querySelector('.u2-wizard-reason').textContent;
const status = () => document.querySelector('.u2-wizard-status');
const accessCalls = (calls) => calls.filter((c) => c[0].endsWith('grant') || c[0].endsWith('Column')).map((c) => c.join(' '));

/** Design (opened on, the name set) › Review, ready to validate. */
async function toReview(dialog, name) {
  await flush();
  dialog.editor.model.setSchemaName(name);
  dialog.wizard.next();
  await flush();
}

/** A pick in the platform group picker: the look-up is typed into (past its debounce), the first
 * row — highlighted as the candidates land — taken with Enter. */
async function pick(input, query) {
  input.focus();
  input.value = query;
  fire(input, 'input');
  await new Promise((r) => setTimeout(r, 200));
  await flush();
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
}

scoped('connection › design › review › created: draft, dry run, create, restrictions then fanned-out grants', async () => {
  const calls = [];
  let refuse = true;
  stub(calls, () => {
    if (refuse)
      throw Object.assign(new Error('Manifest refused'), {code: 'manifest-validation',
        errors: [{path: 'tables.orders.columns.orderid', code: 'external-column-type', message: 'not this'}]});
    return {status: 'ok', issues: []};
  });
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public', table: 'orders'});
  const done = dialog.open();
  await flush();
  assert.deepEqual(calls.filter((c) => c[0] === 'draft').map((c) => c[1]),
    [{connection: CONN.nqName, schema: 'public', catalog: undefined}], 'the draft is read once both are picked');
  assert.match(document.querySelector('.u2-binding-facts').textContent, /Postgres · public: \d+ of \d+ tables bindable/);
  assert.equal(document.querySelector('.u2-binding-access').textContent,
    'You may introspect and query this connection; writes need RemoveRows on the connection, which you lack');
  assert.equal(dialog.wizard.currentStep.value, 'design', 'connection and schema preset: opened on Design');
  const editor = dialog.editor;
  assert.ok(editor, 'the editor is built over the draft');
  assert.deepEqual(editor.model.tables.value.filter((t) => t.included).map((t) => t.remote), ['orders'],
    'only the named table starts included');
  assert.equal(editor.model.friendlyName.value, 'PostgresNorthwind public 2', 'the connection and the remote schema, numbered like the identifier');
  assert.equal(editor.model.name.value, 'postgresnorthwind_public_2', 'harmonized, past the registered one');
  const panel = editor.panel.root;
  assert.equal(panel.querySelector('[data-u2-name="name"] input').value, 'postgresnorthwind_public_2');
  assert.equal(panel.querySelector('[data-u2-name="name"] .u2-input-postfix').textContent,
    'registered as ext_postgresnorthwind_public_2');
  const writable = panel.querySelector('[data-u2-name="writable"] .u2-input-checkbox');
  assert.equal(writable.disabled, true, 'the Writable switch is off with the reason');
  assert.equal(panel.querySelector('[data-u2-name="writable"] .u2-input-postfix').textContent,
    'writes need RemoveRows on the connection, which you lack');
  assert.equal(buttonNamed('NEXT').disabled, false, 'the proposed name passes');
  editor.model.setSchemaName('Bad Name');
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.match(reason(), /^Name: /);
  editor.model.setSchemaName('northwind_sales');

  dialog.wizard.back();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'connection', 'the Connection step is a BACK away, pre-filled');
  assert.equal(dialog.connection.value.value, CONN.id);
  assert.equal(dialog.schema.value.value, 'public');
  dialog.wizard.next();
  await flush();
  assert.equal(dialog.editor, editor, 'the same draft: the editor is kept');

  await editor.select({kind: 'schema'});
  await flush();
  assert.equal(editor.panel.root.querySelector('.u2-access-grid-add'), null, 'no fixed list');
  await pick(editor.panel.root.querySelector('.u2-access-grid .u2-typeahead input'), 'sales');
  assert.deepEqual(editor.access.grants.value, [{scope: SCHEMA, group: SALES, view: true, edit: false, delete: false}],
    'the platform look-up: the pick is stored as {id, label}');
  assert.equal(editor.panel.root.querySelector('.u2-access-grid .u2-typeahead input').value, '', 'ready for the next');
  editor.access.addGroup(ORDERS, DEV);
  editor.access.setGrant(ORDERS, 'g-dev', 'edit', true);
  editor.access.setVisibility('orders', 'freight', []);
  await editor.select({kind: 'column', table: 'orders', column: 'freight'});
  await flush();
  await pick(editor.panel.root.querySelector('[data-u2-name="visibleTo"] .u2-typeahead input'), 'sales');
  assert.deepEqual(editor.access.visibilityOf('orders', 'freight'), [{id: 'g-sales', label: 'Sales'}]);
  await pick(editor.panel.root.querySelector('[data-u2-name="visibleTo"] .u2-typeahead input'), 'sam');
  assert.deepEqual(editor.access.visibilityOf('orders', 'freight'), [SALES, {id: 'u-sam', label: 'Sam Sales'}],
    'a user through the personal group');
  assert.deepEqual(editor.panel.root.querySelectorAll('.u2-chip').map((c) => c.textContent), ['Sales', 'Sam Sales']);
  editor.access.setVisibility('orders', 'freight', [SALES]);
  assert.equal(buttonNamed('NEXT').disabled, false);

  dialog.wizard.next();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.match(document.querySelector('.u2-binding-json').textContent, /"name": "northwind_sales"/);
  assert.equal(buttonNamed('CREATE').disabled, true, 'CREATE waits for a green validate');

  buttonNamed('VALIDATE').click();
  await flush();
  const dryRuns = calls.filter((c) => c[0] === 'create' && c[2].dryRun === true);
  assert.equal(dryRuns.length, 1);
  assert.equal(dryRuns[0][1], 'northwind_sales');
  assert.deepEqual(Object.keys(dryRuns[0][2].manifest.tables), ['orders']);
  assert.equal(editor.diagnostics.value.length, 1, 'the refusal lands on the editor');
  assert.equal(document.querySelectorAll('.u2-binding-issue').length, 1);
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.equal(buttonNamed('CREATE').disabled, true);

  refuse = false;
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(editor.diagnostics.value.length, 0);
  assert.equal(status().textContent, 'Validated');
  assert.equal(buttonNamed('CREATE').disabled, false);

  buttonNamed('CREATE').click();
  await flush();
  const create = calls.find((c) => c[0] === 'create' && c[2].dryRun === undefined);
  assert.equal(create[1], 'northwind_sales');
  assert.deepEqual(Object.keys(create[2].manifest.tables), ['orders']);
  // the restriction first, then the grants: the schema row is a grant on every included table
  // (orders alone here), the table row another — by group id, never through schema.grant; the
  // binding is read-only, so the Edit the Developers row asked for is not granted
  assert.deepEqual(accessCalls(calls), ['shareColumn northwind_sales.orders freight g-sales View',
    'table.grant northwind_sales.orders g-sales View', 'table.grant northwind_sales.orders g-dev View']);
  assert.equal(dialog.wizard.currentStep.value, 'created', 'access rows: the dialog stays on the report');
  assert.match(document.querySelector('.u2-binding-created').textContent, /northwind_sales is registered/);
  assert.equal(document.querySelectorAll('.u2-binding-created-applied span').length, 4, 'a title and three lines');
  assert.equal(document.querySelector('.u2-binding-created-failed'), null);
  assert.equal(buttonNamed('RETRY ACCESS').disabled, true, 'nothing to retry');
  assert.equal(buttonNamed('CANCEL').style.display, 'none');
  assert.equal(grok.dapi.domains.invalidated > 0, true);
  assert.deepEqual(calls.filter((c) => c[0] === 'announce'), [], 'an answered create needs no announcement');
  buttonNamed('CLOSE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, 'northwind_sales');
  assert.deepEqual(result.access, {applied: ['orders.freight visible to Sales', 'orders: View for Sales',
    'orders: View for Developers'], failed: []});
  assert.equal(document.querySelector('.u2-dialog'), null, 'the dialog is gone');
});

scoped('a failed restriction withholds that table\'s grants; RETRY ACCESS re-runs the failed steps alone', async () => {
  const calls = [];
  const failing = {shareColumn: 'column schemas are locked'};
  stub(calls, () => ({status: 'ok', issues: []}), failing);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  dialog.editor.access.setVisibility('orders', 'freight', [SALES]);
  buttonNamed('VALIDATE').click();
  await flush();
  buttonNamed('CREATE').click();
  await flush();
  assert.deepEqual(accessCalls(calls), ['shareColumn nw.orders freight g-sales View', 'table.grant nw.order_details g-sales View'],
    'orders got no grant while its column stayed unrestricted; order_details went through');
  assert.equal(dialog.wizard.currentStep.value, 'created');
  assert.deepEqual([...document.querySelectorAll('.u2-binding-created-failed span')].map((s) => s.textContent), [
    '2 access steps failed — retry, or redo them from the schema\'s page',
    'orders.freight visible to Sales: column schemas are locked',
    'orders: View for Sales: withheld — a column restriction of orders failed',
  ]);
  assert.equal(status().textContent, '2 access steps failed');
  assert.equal(buttonNamed('RETRY ACCESS').disabled, false);

  delete failing.shareColumn;
  calls.length = 0;
  buttonNamed('RETRY ACCESS').click();
  await flush();
  assert.deepEqual(accessCalls(calls), ['shareColumn nw.orders freight g-sales View', 'table.grant nw.orders g-sales View'],
    'only the failed steps ran, the restriction first');
  assert.equal(document.querySelector('.u2-binding-created-failed'), null);
  assert.equal(status().textContent, 'Access applied');
  assert.equal(buttonNamed('RETRY ACCESS').disabled, true);
  document.querySelector('.u2-dialog-close').click();
  await flush();
  const result = await done;
  assert.deepEqual(result.access.failed, [], 'a created schema closed with ✕ is reported created, not cancelled');
  assert.deepEqual(result.access.applied, ['order_details: View for Sales', 'orders.freight visible to Sales',
    'orders: View for Sales'], 'in the order they went through');
});

scoped('a validation is bound to the exact payload: a change since is not validated, a late answer is dropped', async () => {
  const calls = [];
  let release;
  stub(calls, () => ({status: 'ok', issues: []}));
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  const createSchema = grok.dapi.domains.createSchema;
  grok.dapi.domains.createSchema = (name, options) => new Promise((r) => release = () => r(createSchema(name, options)));
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(buttonNamed('CREATE').disabled, true, 'the footer waits');
  assert.equal(buttonNamed('CANCEL').disabled, true);
  dialog.editor.model.friendlyName.value = 'Changed meanwhile';
  release();
  await flush();
  assert.equal(status().textContent, 'Changed while validating — validate again');
  assert.equal(buttonNamed('CREATE').disabled, true);
  assert.equal(reason(), 'Validate before creating');

  grok.dapi.domains.createSchema = createSchema;
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(status().textContent, 'Validated');
  assert.equal(buttonNamed('CREATE').disabled, false);
  dialog.wizard.back();
  await flush();
  dialog.editor.model.includeColumn('orders', 'shipcity', false);
  dialog.wizard.next();
  await flush();
  assert.equal(status().textContent, '', 'Review after an edit: not validated');
  assert.equal(buttonNamed('CREATE').disabled, true);
  dialog.editor.model.includeColumn('orders', 'shipcity', true);
  dialog.wizard.back();
  await flush();
  dialog.wizard.next();
  await flush();
  assert.equal(status().textContent, 'Validated', 'the same payload again is what was validated');
  assert.equal(buttonNamed('CREATE').disabled, false);

  // no access rows: the short way — a toast, the app, the promise, no Created step
  const invalidated = grok.dapi.domains.invalidated;
  buttonNamed('CREATE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, 'nw');
  assert.deepEqual(result.access, {applied: [], failed: []});
  assert.equal(document.querySelector('.u2-dialog'), null);
  assert.equal(grok.dapi.domains.invalidated, invalidated + 1, 'the completion boundary, once');
});

scoped('a refused create keeps the dialog open with the refusal named; cancel resolves null', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.domains.createSchema = async (name, options) => {
    if (options.dryRun)
      return {status: 'ok', issues: []};
    throw Object.assign(new Error('Domain schema "x" is already registered'), {code: 'schema-name-taken'});
  };
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'x');
  buttonNamed('VALIDATE').click();
  await flush();
  buttonNamed('CREATE').click();
  await flush();
  assert.notEqual(document.querySelector('.u2-dialog'), null, 'still open');
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.equal(status().textContent, 'Domain schema "x" is already registered');
  assert.equal(dialog.editor.diagnostics.value[0].code, 'schema-name-taken');
  assert.equal(buttonNamed('CREATE').disabled, true, 'validate again after a change');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('a draft the server refuses is named on the connection step and blocks NEXT', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.domains.draft = async () => {
    throw Object.assign(new Error('You don\'t have DataConnection.Query on connection "x"'), {code: 'forbidden'});
  };
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await flush();
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.equal(reason(), 'You don\'t have DataConnection.Query on connection "x"');
  fire(document.querySelector('.u2-dialog'), 'keydown', {key: 'Escape'});
  await flush();
  assert.equal(await done, null);
});

scoped('a connection alone opens on Connection; Design follows once the schema is picked and the draft is in', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const dialog = new domains.authoring.BindingDialog({connection: CONN});
  const done = dialog.open();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'connection');
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.equal(reason(), 'Pick a connection and a schema');
  dialog.schema.value.value = 'audit';
  await flush();
  assert.equal(buttonNamed('NEXT').disabled, false);
  dialog.wizard.next();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'design');
  assert.equal(dialog.editor.model.friendlyName.value, 'PostgresNorthwind audit');
  assert.equal(dialog.editor.model.name.value, 'postgresnorthwind_audit');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('createBinding refuses while domain databases are off', async () => {
  grok.shell.settings = {enableDomainDatabases: false};
  await assert.rejects(domains.authoring.createBinding({connection: CONN}), /Beta feature/);
  assert.equal(document.querySelector('.u2-dialog'), null);
});

scoped('createBinding refuses by name without the CreateDomainSchema privilege, before anything is read', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.permissions.checkGlobal = async (permission) => {
    calls.push(['global', permission]);
    return false;
  };
  await assert.rejects(domains.authoring.createBinding({connection: CONN, schema: 'public'}),
    /Requires the CreateDomainSchema privilege/);
  assert.equal(document.querySelector('.u2-dialog'), null);
  assert.deepEqual(calls, [['global', 'CreateDomainSchema']], 'no draft, no connections');
});

scoped('the table preset matches the remote name whatever its case; an unmatched one is said on the status line', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const upper = new domains.authoring.BindingDialog({connection: CONN, schema: 'public', table: 'ORDERS'});
  const first = upper.open();
  await flush();
  assert.deepEqual(upper.editor.model.tables.value.filter((t) => t.included).map((t) => t.remote), ['orders'],
    'the catalog spells it in lowercase; the caller in uppercase');
  assert.equal(status().textContent, '');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await first, null);

  const unknown = new domains.authoring.BindingDialog({connection: CONN, schema: 'public', table: 'nope'});
  const second = unknown.open();
  await flush();
  assert.deepEqual(unknown.editor.model.tables.value.filter((t) => t.included).map((t) => t.remote),
    ['orders', 'order_details']);
  assert.equal(status().textContent, 'Table nope is not in public — every table starts included');
  assert.equal(status().classList.contains('u2-wizard-status-error'), false);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await second, null);
});

scoped('the catalog preset stays with the preset connection; the schemas are offered sorted, case aside', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const OTHER = {id: 'c2', nqName: 'Other:Scratch', friendlyName: 'Scratch', dataSource: 'Postgres', isDatabase: true};
  grok.dapi.connections.list = async () => [CONN, OTHER];
  grok.dapi.connections.getSchemas = async (conn, catalog) => {
    calls.push(['schemas', conn.id, catalog]);
    return ['public', 'Audit', 'archive'];
  };
  const dialog = new domains.authoring.BindingDialog({connection: CONN, catalog: 'Northwind'});
  const done = dialog.open();
  await flush();
  assert.deepEqual(dialog.schema.items.map((i) => i.value), ['archive', 'Audit', 'public']);
  dialog.schema.value.value = 'public';
  await flush();
  dialog.connection.value.value = OTHER.id;
  await flush();
  dialog.schema.value.value = 'public';
  await flush();
  assert.deepEqual(calls.filter((c) => c[0] === 'schemas'), [['schemas', 'c1', 'Northwind'], ['schemas', 'c2', null]],
    'the other connection is read whole');
  assert.deepEqual(calls.filter((c) => c[0] === 'draft').map((c) => [c[1].connection, c[1].catalog]),
    [[CONN.nqName, 'Northwind'], [OTHER.nqName, undefined]]);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('every connection past the first page is offered; the identifier is probed by name past every registered one', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const many = Array.from({length: 510}, (_, i) => ({id: `c${i}`, nqName: `P:c${i}`, isDatabase: true,
    friendlyName: `conn ${String(i).padStart(3, '0')}`, dataSource: 'Postgres'}));
  grok.dapi.connections.list = async ({pageSize, pageNumber, order}) => {
    calls.push(['connections', pageSize, pageNumber, order]);
    return many.slice((pageNumber - 1) * pageSize, pageNumber * pageSize);
  };
  REGISTERED.set('conn_000_public', {});
  REGISTERED.set('conn_000_public_2', {});
  try {
    const dialog = new domains.authoring.BindingDialog({connection: many[0], schema: 'public'});
    const done = dialog.open();
    await flush();
    assert.deepEqual(calls.filter((c) => c[0] === 'connections'), [['connections', 500, 1, 'id'], ['connections', 500, 2, 'id']],
      'pages in a stable order');
    assert.equal(dialog.connection.items.length, 510);
    assert.deepEqual(calls.filter((c) => c[0] === 'registry').map((c) => c[1]),
      ['conn_000_public', 'conn_000_public_2', 'conn_000_public_3'], 'probed until a name is free');
    assert.equal(dialog.editor.model.name.value, 'conn_000_public_3');
    assert.equal(dialog.editor.model.friendlyName.value, 'conn 000 public 3');
    assert.equal(status().textContent, '');
    buttonNamed('CANCEL').click();
    await flush();
    assert.equal(await done, null);
  } finally {
    REGISTERED.delete('conn_000_public');
    REGISTERED.delete('conn_000_public_2');
  }
});

scoped('a warehouse read that never answers fails by name after the timeout; the dialog stays usable', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const never = () => new Promise(() => {});
  grok.dapi.connections.getSchemas = never;
  const {BindingDialog} = domains.authoring;
  domains.authoring.BindingWizard.readTimeout = 20;
  try {
    const dialog = new BindingDialog({connection: CONN, catalog: 'bogus'});
    const done = dialog.open();
    await new Promise((r) => setTimeout(r, 40));
    await flush();
    assert.equal(document.querySelector('.u2-binding-facts').textContent,
      'Postgres · The schemas could not be read: The connection did not answer within 0.02 s');
    assert.equal(dialog.connection.enabled, true);
    grok.dapi.connections.getSchemas = async () => ['public'];
    grok.dapi.domains.draft = never;
    dialog.connection.value.value = null;
    await flush();
    dialog.connection.value.value = CONN.id;
    await flush();
    dialog.schema.value.value = 'public';
    await new Promise((r) => setTimeout(r, 40));
    await flush();
    assert.equal(reason(), 'The draft did not answer within 0.02 s');
    assert.equal(buttonNamed('NEXT').disabled, true);
    buttonNamed('CANCEL').click();
    await flush();
    assert.equal(await done, null);
  } finally {
    domains.authoring.BindingWizard.readTimeout = 30000;
  }
});

scoped('opened on Design over a preset whose reads never answer: the draft fails by name on the status line, BACK works', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const never = () => new Promise(() => {});
  grok.dapi.connections.getSchemas = never;
  grok.dapi.domains.draft = never;
  const {BindingDialog} = domains.authoring;
  domains.authoring.BindingWizard.readTimeout = 20;
  try {
    const dialog = new BindingDialog({connection: CONN, schema: 'public', catalog: 'bogus'});
    const done = dialog.open();
    await new Promise((r) => setTimeout(r, 60));
    await flush();
    assert.equal(dialog.wizard.currentStep.value, 'design');
    assert.equal(dialog.schema.value.value, 'public', 'the preset schema is still picked');
    assert.equal(status().textContent, 'The draft did not answer within 0.02 s');
    assert.equal(status().classList.contains('u2-wizard-status-error'), true);
    assert.equal(reason(), 'The draft did not answer within 0.02 s');
    assert.equal(buttonNamed('NEXT').disabled, true);
    dialog.wizard.back();
    await flush();
    assert.equal(dialog.wizard.currentStep.value, 'connection');
    assert.equal(document.querySelector('.u2-binding-facts').textContent,
      'Postgres · The schemas could not be read: The connection did not answer within 0.02 s');
    grok.dapi.domains.draft = async (body) => {
      calls.push(['draft', body]);
      return JSON.parse(JSON.stringify(DRAFT));
    };
    dialog.schema.value.value = null;
    await flush();
    dialog.schema.value.value = 'public';
    await flush();
    assert.equal(status().textContent, '', 'a new read clears the failure');
    assert.equal(buttonNamed('NEXT').disabled, false);
    buttonNamed('CANCEL').click();
    await flush();
    assert.equal(await done, null);
  } finally {
    domains.authoring.BindingWizard.readTimeout = 30000;
  }
});

scoped('a probe the registry did not answer keeps the proposed name and says so, beside the table note', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.domains.schema = () => ({manifest: async () => {
    throw Object.assign(new Error('XMLHttpRequest error.'), {status: 0, code: ''});
  }});
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'audit', table: 'nope'});
  const done = dialog.open();
  await flush();
  assert.equal(dialog.editor.model.name.value, 'postgresnorthwind_audit');
  assert.equal(status().textContent, 'Table nope is not in audit — every table starts included; ' +
    'Registered schemas could not be checked: The connection to the server failed');
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.equal(buttonNamed('NEXT').disabled, false, 'a taken name is refused at CREATE all the same');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

const rail = async (id) => {
  fire(document.querySelector(`.u2-wizard-step[data-id="${id}"]`), 'click');
  await flush();
};
const creates = (calls) => calls.filter((c) => c[0] === 'create' && c[2].dryRun === undefined);

async function validated() {
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(status().textContent, 'Validated');
}

scoped('WO-A5.1 #1 (P0): the validation and the editor belong to the draft picked now; the rail passes every gate on the way', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  const first = dialog.editor;
  await validated();

  await rail('connection');
  const draft = grok.dapi.domains.draft;
  let release;
  grok.dapi.domains.draft = (body) => new Promise((r) => release = () => r(draft(body)));
  dialog.schema.value.value = 'audit';
  await flush();
  assert.equal(status().textContent, '', 'the Validated badge of the other draft is gone');
  await rail('design');
  assert.equal(dialog.wizard.currentStep.value, 'connection', 'the draft is still being read: its gate holds the rail');
  release();
  await flush();
  await rail('review');
  assert.equal(dialog.wizard.currentStep.value, 'connection', 'Design still holds the editor over the other draft');
  assert.equal(dialog._reviewGate(), 'Validate before creating', 'the validation was of the other draft');
  await rail('design');
  assert.equal(dialog.wizard.currentStep.value, 'design');
  assert.notEqual(dialog.editor, first, 'rebuilt over the new draft');
  dialog.editor.model.setSchemaName('nw');
  await rail('review');
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.equal(buttonNamed('CREATE').disabled, true, 'the same name, but never validated over this draft');

  await rail('design');
  dialog.editor.model.includeTables(false);
  await rail('review');
  assert.equal(dialog.wizard.currentStep.value, 'design', 'no table: the Design gate holds the rail');
  dialog.editor.model.includeTables(true);
  dialog.editor.model.setSchemaName('Bad Name');
  await rail('review');
  assert.equal(dialog.wizard.currentStep.value, 'design', 'a bad name: the Design gate holds the rail');
  assert.equal(creates(calls).length, 0);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('WO-A5.1 #2 (P1): findings belong to the payload they were found in; a transport failure is no finding', async () => {
  const calls = [];
  stub(calls, () => {
    throw Object.assign(new Error('Manifest refused'), {code: 'manifest-validation',
      errors: [{path: 'tables.orders.columns.shipcity', code: 'external-column-type', message: 'bad shipcity'}]});
  });
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(document.querySelectorAll('.u2-binding-issue').length, 1);
  dialog.wizard.back();
  dialog.editor.tree.tree.root.querySelector('.u2-list').clientHeight = 800;
  await dialog.editor.select({kind: 'column', table: 'orders', column: 'shipcity'});
  await flush();
  assert.equal(document.querySelectorAll('.u2-manifest-node-problem').length, 1, 'the row is marked');
  assert.equal(document.querySelectorAll('.u2-manifest-panel-diagnostic').length, 1, 'and the panel');
  dialog.editor.model.includeColumn('orders', 'shipcity', false);
  await flush();
  assert.deepEqual(dialog.editor.diagnostics.value, [], 'the edit makes them history');
  assert.equal(document.querySelectorAll('.u2-manifest-node-problem').length, 0, 'no row stays marked');
  assert.equal(document.querySelectorAll('.u2-manifest-panel-diagnostic').length, 0);
  dialog.wizard.next();
  await flush();
  assert.equal(document.querySelectorAll('.u2-binding-issue').length, 0);
  assert.equal(document.querySelector('.u2-binding-issues-title').textContent, 'Not validated');

  grok.dapi.domains.createSchema = async () => {
    throw Object.assign(new Error('Failed to fetch'), {status: 0, code: ''});
  };
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(status().textContent, 'Could not reach the server (Failed to fetch)', 'a domain call nothing answered');
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.deepEqual(dialog.editor.diagnostics.value, [], 'the transport error is on the status line only');
  grok.dapi.domains.createSchema = async () => {
    throw new Error('Internal error');
  };
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(status().textContent, 'Internal error', 'a refusal without the domain shape is worded as it came');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('WO-A5.1 #3 (P2): a change and the change back is still validated, whatever order the keys come back in', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  await validated();
  const model = dialog.editor.model;
  dialog.wizard.back();
  await flush();
  model.setFriendlyName('orders', '');
  model.setFriendlyName('orders', 'Orders');
  model.setNameColumn('orders', 'shipcity');
  model.setNameColumn('orders', 'shipname');
  model.setRequired('orders', 'shipcity', true);
  model.setRequired('orders', 'shipcity', false);
  dialog.wizard.next();
  await flush();
  assert.equal(status().textContent, 'Validated');
  assert.equal(buttonNamed('CREATE').disabled, false);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('WO-A5.1 #4 (P2): validations are numbered — an older answer landing last never overrides a newer one', async () => {
  const calls = [];
  const pending = [];
  stub(calls, () => new Promise((ok, no) => pending.push({ok, no})));
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  const older = dialog._validate();
  const newer = dialog._validate();
  await flush();
  pending[1].ok({status: 'ok', issues: []});
  await newer;
  pending[0].no(Object.assign(new Error('late refusal'), {code: 'manifest-validation'}));
  await older;
  await flush();
  assert.equal(status().textContent, 'Validated', 'the newer answer stands');
  assert.deepEqual(dialog.editor.diagnostics.value, []);
  assert.equal(buttonNamed('CREATE').disabled, false);

  const third = dialog._validate();
  const fourth = dialog._validate();
  await flush();
  pending[3].no(Object.assign(new Error('refused'), {code: 'manifest-validation'}));
  await fourth;
  pending[2].ok({status: 'ok', issues: []});
  await third;
  await flush();
  assert.equal(buttonNamed('CREATE').disabled, true, 'a late older OK does not revive the refused payload');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('WO-A5.1 #5 (P2): the same connection and schema picked again keep the editor; another draft resets it and says so', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await flush();
  const editor = dialog.editor;
  editor.model.includeColumn('orders', 'shipcity', false);
  dialog.wizard.back();
  await flush();
  dialog.schema.value.value = 'audit';
  await flush();
  dialog.schema.value.value = 'public';
  await flush();
  dialog.wizard.next();
  await flush();
  assert.equal(dialog.editor, editor, 'the editor is kept');
  assert.equal(dialog.editor.model.column('orders', 'shipcity').included, false, 'the edit stands');
  assert.equal(status().textContent, '');

  dialog.wizard.back();
  await flush();
  dialog.schema.value.value = 'audit';
  await flush();
  dialog.wizard.next();
  await flush();
  assert.notEqual(dialog.editor, editor);
  assert.equal(dialog.editor.model.column('orders', 'shipcity').included, true);
  assert.equal(status().textContent, 'The design was reset: PostgresNorthwind · audit is a new draft');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

/** A registry the stubbed create writes to, at version 1 like the server's, authored by the
 * signed-in user; `lose` loses the answer of the next create: 'after' commits it first, 'before'
 * commits it only once `land()` is called; `foreign(name, tables)` registers a schema of someone
 * else's under that name, `alter(name, f)` changes a registered one as `f` says, `author(name, id)`
 * says who registered it. */
function registry(calls) {
  const manifests = new Map([['postgresnorthwind_public', {storage: {kind: 'external'}, tables: {}}]]);
  const authors = new Map();
  const names = {has: (name) => manifests.has(name)};
  const state = {lose: null, land: () => {},
    foreign: (name, tables = {t: {columns: {}}}) => {
      manifests.set(name, {version: '1', storage: {kind: 'external', connection: 'Other:Conn', schema: 'x'}, tables});
      authors.set(name, 'u-other');
    },
    alter: (name, f) => manifests.set(name, f(JSON.parse(JSON.stringify(manifests.get(name))))),
    author: (name, id) => authors.set(name, id)};
  grok.dapi.domains.schemas = {filter: (text) => ({list: async () => {
    calls.push(['schemas', text]);
    return [...manifests.keys()].filter((name) => text === `pgSchema = "ext_${name}"`)
      .map((name) => ({name, author: {id: authors.get(name)}}));
  }})};
  const schema = grok.dapi.domains.schema;
  grok.dapi.domains.schema = (name) => ({...schema(name), manifest: async () => {
    calls.push(['registry', name]);
    if (!manifests.has(name))
      throw notFound(name);
    return manifests.get(name);
  }});
  grok.dapi.domains.createSchema = async (name, options) => {
    calls.push(['create', name, options]);
    if (options.dryRun)
      return {status: 'ok', issues: []};
    if (names.has(name))
      throw Object.assign(new Error(`Domain schema "${name}" is already registered`), {code: 'schema-name-taken'});
    const lose = state.lose;
    state.lose = null;
    const add = () => {
      manifests.set(name, {...JSON.parse(JSON.stringify(options.manifest)), version: '1'});
      authors.set(name, grok.shell.user.id);
    };
    if (lose === 'after')
      add();
    if (lose === 'before')
      state.land = add;
    if (lose !== null)
      throw Object.assign(new Error('Gateway Timeout'), {status: 504, code: ''});
    add();
    return {id: 'id-1', name, pgSchema: `ext_${name}`, version: '1'};
  };
  return state;
}

scoped('WO-A5.1 #6 (P0): a create whose answer was lost but which registered goes on to the access and Created', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const server = registry(calls);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  await validated();
  server.lose = 'after';
  buttonNamed('CREATE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'created');
  assert.deepEqual(accessCalls(calls), ['table.grant nw.orders g-sales View', 'table.grant nw.order_details g-sales View']);
  assert.deepEqual(calls.filter((c) => c[0] === 'announce').map((c) => c[1]), ['nw'],
    'the recovered schema is announced to the platform');
  assert.deepEqual(dialog.editor.diagnostics.value, [], 'the timeout is no finding');
  buttonNamed('CLOSE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, 'nw', 'a committed schema is never reported cancelled');
});

scoped('WO-A5.1 #6 (P0): a create the registry does not know yet stays offered; a "taken" answering our retry reads the registry again', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const server = registry(calls);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  await validated();
  server.lose = 'before';
  buttonNamed('CREATE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.equal(status().textContent, 'Gateway Timeout — nw is not registered; CREATE again');
  assert.deepEqual(calls.filter((c) => c[0] === 'announce'), [], 'nothing registered, nothing announced');
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.equal(document.querySelectorAll('.u2-binding-issue').length, 0, 'on the status line, not a finding');
  assert.equal(buttonNamed('CREATE').disabled, false, 'CREATE stays available');

  server.land();
  buttonNamed('CREATE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, 'nw', 'the "taken" was our own create landing late');
  assert.equal(creates(calls).length, 2);
  assert.deepEqual(calls.filter((c) => c[0] === 'announce').map((c) => c[1]), ['nw'],
    'the schema our retry found registered is announced once');
});

scoped('WO-A5.1 #6 (P0): after a lost create, a schema someone else registered under the name is not claimed', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const server = registry(calls);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  await validated();
  server.lose = 'before';
  buttonNamed('CREATE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'review');
  server.foreign('nw');
  buttonNamed('CREATE').click();
  await flush();
  assert.deepEqual(accessCalls(calls), [], 'no grant lands on the other schema');
  assert.equal(dialog.wizard.currentStep.value, 'created', 'the outcome of the lost create is unknown');
  assert.match(document.querySelector('.u2-binding-created-unknown').textContent,
    /nw is registered at version 1 over Other:Conn · x, 1 table \(t\) — not provably this dialog's create/);
  assert.equal(creates(calls).length, 2);
  buttonNamed('CLOSE').click();
  await flush();
  assert.deepEqual(await done, {name: 'nw', access: {applied: [], failed: []}, unknown: true});
});

/** A lost create, the registry then read with `change` applied to what landed, and registered by
 * `author` where one is given. */
async function lostCreate(change, author) {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const server = registry(calls);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  dialog.editor.access.setVisibility('orders', 'freight', [SALES]);
  await validated();
  server.lose = 'after';
  const createSchema = grok.dapi.domains.createSchema;
  grok.dapi.domains.createSchema = async (name, options) => {
    try {
      return await createSchema(name, options);
    } finally {
      if (!options.dryRun) {
        server.alter(name, change);
        if (author !== undefined)
          server.author(name, author);
      }
    }
  };
  buttonNamed('CREATE').click();
  await flush();
  return {calls, dialog, done};
}

scoped('a lost create is ours only at version 1 with every table descriptor as sent: a later apply is an unknown outcome, nothing granted', async () => {
  const {calls, dialog, done} = await lostCreate((m) => ({...m, version: '2'}));
  assert.equal(dialog.wizard.currentStep.value, 'created');
  assert.deepEqual(accessCalls(calls), [], 'no restriction, no grant');
  assert.deepEqual(calls.filter((c) => c[0] === 'announce'), []);
  assert.deepEqual(calls.filter((c) => c[0] === 'altered').map((c) => c[1]), ['nw'],
    'the platform re-reads what is registered under the name');
  assert.equal(document.querySelector('.u2-binding-created-head').textContent, 'Outcome unknownA schema nw is registered');
  assert.deepEqual([...document.querySelectorAll('.u2-binding-created-unknown span')].map((s) => s.textContent), [
    'No answer to the create (Gateway Timeout) — whether it landed is unknown',
    'nw is registered at version 2 over NorthwindBinding:PostgresNorthwind · public, 2 tables (orders, order_details) — not provably this dialog\'s create',
    'No access was applied: grant it from the schema\'s page once you know the binding is yours']);
  assert.equal(status().textContent, 'Outcome unknown — no access was applied');
  assert.equal(buttonNamed('RETRY ACCESS').style.display, 'none');
  assert.equal(buttonNamed('OPEN').style.display, '');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).unknown, true);
});

scoped('a lost create is ours only when the caller registered it: the same content at version 1 by someone else is an unknown outcome', async () => {
  const {calls, dialog, done} = await lostCreate((m) => m, 'u-other');
  assert.equal(dialog.wizard.currentStep.value, 'created');
  assert.deepEqual(accessCalls(calls), [], 'no restriction, no grant');
  assert.ok(document.querySelector('.u2-binding-created-unknown'));
  assert.deepEqual(calls.filter((c) => c[0] === 'schemas').map((c) => c[1]), ['pgSchema = "ext_nw"']);
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).unknown, true);
});

scoped('an unknown outcome over a schema with no tables offers nothing to open', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const server = registry(calls);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  await validated();
  server.lose = 'before';
  buttonNamed('CREATE').click();
  await flush();
  server.foreign('nw', {});
  buttonNamed('CREATE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'created');
  assert.match(document.querySelector('.u2-binding-created-unknown').textContent, /0 tables — not provably/);
  assert.equal(buttonNamed('OPEN').style.display, 'none');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).unknown, true);
});

scoped('a lost create whose registered column differs at version 1 is not claimed; one the registry writes back without a derived friendly name is', async () => {
  let run = await lostCreate((m) => {
    m.tables.orders.columns.freight = {...m.tables.orders.columns.freight, required: true};
    return m;
  });
  assert.deepEqual(accessCalls(run.calls), []);
  assert.ok(document.querySelector('.u2-binding-created-unknown'), 'a changed descriptor: outcome unknown');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await run.done).unknown, true);
  resetDom();

  run = await lostCreate((m) => {
    m.tables.orders.friendlyName = undefined;
    m.tables.orders = JSON.parse(JSON.stringify(m.tables.orders));
    return m;
  });
  assert.equal(document.querySelector('.u2-binding-created-unknown'), null);
  assert.deepEqual(accessCalls(run.calls), ['shareColumn nw.orders freight g-sales View',
    'table.grant nw.orders g-sales View', 'table.grant nw.order_details g-sales View'], 'our own create: the access goes on');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await run.done).unknown, undefined);
});

scoped('WO-A5.1 #7 (P1): a group gone since the pick refuses the create on Review by its label; nothing is created', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.groups.list = async ({filter}) => {
    calls.push(['groups', filter]);
    return [{id: 'g-sales'}];
  };
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  dialog.editor.access.addGroup(ORDERS, DEV);
  dialog.editor.access.setVisibility('orders', 'freight', [DEV, SALES]);
  await validated();
  buttonNamed('CREATE').click();
  await flush();
  assert.deepEqual(calls.filter((c) => c[0] === 'groups').map((c) => c[1]), ['id in ("g-sales", "g-dev")'],
    'every distinct group in one look-up');
  assert.equal(creates(calls).length, 0);
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.equal(status().textContent, 'Group Developers no longer exists — take it out of the access rows');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('WO-A5.1 #8 (P1): OPEN on a registered schema no route opens says so', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const table = domains.table;
  domains.table = async () => {
    throw new Error('Unknown domain table nw.orders');
  };
  try {
    const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
    const done = dialog.open();
    await toReview(dialog, 'nw');
    dialog.editor.access.addGroup(SCHEMA, SALES);
    await validated();
    buttonNamed('CREATE').click();
    await flush();
    buttonNamed('OPEN').click();
    assert.equal((await done).name, 'nw');
    await flush();
    assert.equal(document.querySelector('.u2-notify-error .u2-notify-content').textContent,
      'Domain schema nw is registered but could not be opened');
  } finally {
    domains.table = table;
  }
});

scoped('WO-A5.1 #10 (P1): a column visible only to some is shared with View, and with Edit too once the plan lets anyone edit the table', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.permissions.check = async () => true;
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  const editor = dialog.editor;
  editor.model.setWritable(true);
  editor.access.addGroup(ORDERS, SALES);
  editor.access.addGroup(ORDERS, DEV);
  editor.access.setGrant(ORDERS, 'g-dev', 'edit', true);
  editor.access.setVisibility('orders', 'freight', [SALES, DEV]);
  await validated();
  buttonNamed('CREATE').click();
  await flush();
  assert.deepEqual(accessCalls(calls), ['shareColumn nw.orders freight g-sales View',
    'shareColumn nw.orders freight g-sales Edit', 'shareColumn nw.orders freight g-dev View',
    'shareColumn nw.orders freight g-dev Edit', 'table.grant nw.orders g-sales View', 'table.grant nw.orders g-dev View',
    'table.grant nw.orders g-dev Edit'], 'an Edit share alone hides the column: View always comes with it');
  assert.deepEqual([...document.querySelectorAll('.u2-binding-created-applied span')].slice(1, 5).map((s) => s.textContent),
    ['orders.freight visible to Sales', 'orders.freight editable by Sales (with table Edit)',
      'orders.freight visible to Developers', 'orders.freight editable by Developers (with table Edit)']);
  buttonNamed('CLOSE').click();
  await flush();
  await done;
});

scoped('Review lists the planned access by the Created labels, in apply order, and follows the Design edits', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  grok.dapi.permissions.check = async () => true;
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await flush();
  const editor = dialog.editor;
  editor.model.setSchemaName('nw');
  editor.model.setWritable(true);
  editor.access.addGroup(ORDERS, DEV);
  editor.access.setGrant(ORDERS, 'g-dev', 'edit', true);
  editor.access.setVisibility('orders', 'freight', [SALES]);
  const planned = () => [...document.querySelectorAll('.u2-binding-planned-op')].map((s) => s.textContent);
  const note = () => document.querySelector('.u2-binding-planned-note').textContent;
  dialog.wizard.next();
  await flush();
  assert.equal(document.querySelector('.u2-binding-planned .u2-section-title').textContent, 'Planned access');
  assert.deepEqual(planned(), ['orders.freight visible to Sales', 'orders.freight editable by Sales (with table Edit)',
    'orders: View for Developers', 'orders: Edit for Developers'], 'restrictions first, then the grants — as CREATE applies them');
  assert.equal(note(), 'Applied after CREATE, in this order. Not covered by Validate.');

  dialog.wizard.back();
  await flush();
  editor.access.removeGroup(ORDERS, 'g-dev');
  dialog.wizard.next();
  await flush();
  assert.deepEqual(planned(), ['orders.freight visible to Sales'],
    'the grant is gone from Review, and with it the column Edit share that rode on the table Edit');

  dialog.wizard.back();
  await flush();
  editor.access.setVisibility('orders', 'freight', null);
  dialog.wizard.next();
  await flush();
  assert.deepEqual(planned(), []);
  assert.equal(note(), 'No additional table or column grants planned');
  assert.equal(buttonNamed('CREATE').disabled, true, 'access is no part of the validation');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

/** The access calls of the groups in `failing` fail with `message` while they are listed. */
function failFor(failing, message) {
  const table = grok.dapi.domains.table;
  grok.dapi.domains.table = (address) => {
    const client = table(address);
    const refuse = (group) => {
      if (failing.includes(group))
        throw Object.assign(new Error(message(group)), {code: 'validation'});
    };
    return {...client,
      grant: async (group, permission) => {
        await client.grant(group, permission);
        refuse(group);
      },
      shareColumn: async (column, group, permission) => {
        await client.shareColumn(column, group, permission);
        refuse(group);
      }};
  };
}

scoped('WO-A5.1 #11 (P2): a group gone since the look-up is named and not offered for retry', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  failFor(['g-dev'], (group) => `Unknown group "${group}"`);
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  dialog.editor.access.addGroup({kind: 'table', table: 'order_details'}, DEV);
  dialog.editor.access.setVisibility('orders', 'freight', [DEV]);
  await validated();
  buttonNamed('CREATE').click();
  await flush();
  assert.deepEqual([...document.querySelectorAll('.u2-binding-created-failed span')].map((s) => s.textContent), [
    '3 access steps failed — redo them from the schema\'s page',
    'orders.freight visible to Developers: group Developers no longer exists',
    'orders: View for Sales: withheld — a column restriction of orders failed',
    'order_details: View for Developers: group Developers no longer exists',
  ]);
  assert.equal(buttonNamed('RETRY ACCESS').disabled, true, 'no retry brings a group back');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).access.failed.length, 3);
});

scoped('WO-A5.1 #12 (P2): one restriction step per column and group — a retry repeats no share that went through', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}));
  const failing = ['g-dev'];
  failFor(failing, () => 'column schemas are locked');
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.setVisibility('orders', 'freight', [SALES, DEV]);
  await validated();
  buttonNamed('CREATE').click();
  await flush();
  assert.deepEqual(accessCalls(calls), ['shareColumn nw.orders freight g-sales View', 'shareColumn nw.orders freight g-dev View']);
  failing.length = 0;
  calls.length = 0;
  buttonNamed('RETRY ACCESS').click();
  await flush();
  assert.deepEqual(accessCalls(calls), ['shareColumn nw.orders freight g-dev View'], 'the Sales share is not repeated');
  buttonNamed('CLOSE').click();
  await flush();
  assert.deepEqual((await done).access, {applied: ['orders.freight visible to Sales', 'orders.freight visible to Developers'],
    failed: []});
});

scoped('WO-A5.1 #13 (P2): a stray Enter does not close the Created report; CLOSE does', async () => {
  const calls = [];
  stub(calls, () => ({status: 'ok', issues: []}), {grant: 'refused'});
  const dialog = new domains.authoring.BindingDialog({connection: CONN, schema: 'public'});
  const done = dialog.open();
  await toReview(dialog, 'nw');
  dialog.editor.access.addGroup(SCHEMA, SALES);
  await validated();
  const content = document.querySelector('.u2-wizard-content');
  fire(content, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'created');
  fire(content, 'keydown', {key: 'Enter'});
  await flush();
  assert.notEqual(document.querySelector('.u2-dialog'), null, 'the report with its failures stays');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).access.failed.length, 2);
});
