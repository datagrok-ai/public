/* The "Create domain schema" dialog over stubbed dapi calls: the draft is read once the
   connection and the schema are picked, the Design step is built over it (one table included when
   the caller named one), the review shows the manifest, VALIDATE is the dry run bound to the exact
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
  dataSource: 'Postgres'};
const SALES = {id: 'g-sales', label: 'Sales'};
const DEV = {id: 'g-dev', label: 'Developers'};
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

/** Every dapi call the dialog makes, recorded; the dry run answers what `dryRun` says; a step
 * named in `failing` refuses with its message. */
function stub(calls, dryRun, failing = {}) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.events.fireCustomEvent = (id, args) => calls.push(['event', id, args]);
  grok.dapi.connections = {list: async () => [CONN], getSchemas: async () => ['public', 'audit']};
  grok.dapi.permissions = {check: async (_c, right) => right !== 'DataConnection.RemoveRows'};
  grok.dapi.groups = {filter: () => ({list: async () =>
    [{id: 'g-sales', friendlyName: 'Sales'}, {id: 'g-dev', friendlyName: 'Developers'},
      {id: 'g-dev-2', friendlyName: 'Developers'}]})};
  const refuse = (step) => {
    if (failing[step])
      throw new Error(failing[step]);
  };
  Object.assign(grok.dapi.domains, {
    draft: async (body) => {
      calls.push(['draft', body]);
      return JSON.parse(JSON.stringify(DRAFT));
    },
    createSchema: async (name, options) => {
      calls.push(['create', name, options]);
      return options.dryRun ? dryRun() : {id: 'id-1', name, pgSchema: `ext_${name}`, version: '1'};
    },
    schema: (name) => ({grant: async (group, permission) => calls.push(['schema.grant', name, group, permission])}),
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

/** Connection › Design (the name set) › Review, ready to validate. */
async function toReview(dialog, name) {
  await flush();
  dialog.wizard.next();
  await flush();
  dialog.editor.model.name.value = name;
  dialog.wizard.next();
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
  assert.equal(buttonNamed('NEXT').disabled, false);

  dialog.wizard.next();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'design');
  const editor = dialog.editor;
  assert.ok(editor, 'the editor is built over the draft');
  assert.deepEqual(editor.model.tables.value.filter((t) => t.included).map((t) => t.remote), ['orders'],
    'only the named table starts included');
  const writable = editor.panel.root.querySelector('[data-u2-name="writable"] .u2-input-checkbox');
  assert.equal(writable.disabled, true, 'the Writable switch is off with the reason');
  assert.match(editor.panel.root.textContent, /writes need RemoveRows/);
  assert.equal(buttonNamed('NEXT').disabled, false, 'the drafted name passes');
  editor.model.name.value = 'Bad Name';
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.match(reason(), /^Name: /);
  editor.model.name.value = 'northwind_sales';
  await editor.select({kind: 'schema'});
  await flush();
  const picker = editor.panel.root.querySelector('.u2-access-grid-add');
  assert.deepEqual(picker.querySelectorAll('option').map((o) => o.textContent), ['+ add group or user…', 'Sales',
    'Developers (g-dev)', 'Developers (g-dev-2)'], 'groups are offered by label, a duplicate name disambiguated');
  editor.access.addGroup(SCHEMA, SALES);
  editor.access.addGroup(ORDERS, DEV);
  editor.access.setGrant(ORDERS, 'g-dev', 'edit', true);
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
  // (orders alone here), the table row another — by group id, never through schema.grant
  assert.deepEqual(accessCalls(calls), ['shareColumn northwind_sales.orders freight g-sales View',
    'table.grant northwind_sales.orders g-sales View', 'table.grant northwind_sales.orders g-dev View',
    'table.grant northwind_sales.orders g-dev Edit']);
  assert.equal(dialog.wizard.currentStep.value, 'created', 'access rows: the dialog stays on the report');
  assert.match(document.querySelector('.u2-binding-created').textContent, /northwind_sales is registered/);
  assert.equal(document.querySelectorAll('.u2-binding-created-applied span').length, 5, 'a title and four lines');
  assert.equal(document.querySelector('.u2-binding-created-failed'), null);
  assert.equal(buttonNamed('RETRY ACCESS').disabled, true, 'nothing to retry');
  assert.equal(buttonNamed('CANCEL').style.display, 'none');
  assert.deepEqual(calls.filter((c) => c[0] === 'event').map((c) => [c[1], c[2].name]),
    [['domain-schema-created', 'northwind_sales']]);
  assert.equal(grok.dapi.domains.invalidated > 0, true);
  buttonNamed('CLOSE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, 'northwind_sales');
  assert.deepEqual(result.access, {applied: ['orders.freight visible to Sales', 'orders: View for Sales',
    'orders: View for Developers', 'orders: Edit for Developers'], failed: []});
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
  buttonNamed('CREATE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, 'nw');
  assert.deepEqual(result.access, {applied: [], failed: []});
  assert.equal(document.querySelector('.u2-dialog'), null);
  assert.deepEqual(calls.filter((c) => c[0] === 'event').length, 1);
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

scoped('createBinding refuses while domain databases are off', async () => {
  grok.shell.settings = {enableDomainDatabases: false};
  await assert.rejects(domains.authoring.createBinding({connection: CONN}), /Beta feature/);
  assert.equal(document.querySelector('.u2-dialog'), null);
});
