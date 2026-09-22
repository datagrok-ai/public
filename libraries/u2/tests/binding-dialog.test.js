/* The "Create domain schema" dialog over stubbed dapi calls: the draft is read once the
   connection and the schema are picked, the Design step is built over it (one table included when
   the caller named one), the review shows the manifest, VALIDATE is the dry run (a refusal lands
   on the rows and blocks CREATE), CREATE is the real one followed by the grants and the column
   restriction through the API each takes, and the promise resolves to the name. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {notify} from '../src/components/display/notify.js';

register('./dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {domains} = await import('../src/dg/domain/index.js');

const DRAFT = JSON.parse(readFileSync(new URL('./fixtures/authoring/northwind-draft.json', import.meta.url), 'utf8'));
const CONN = {id: 'c1', nqName: 'NorthwindBinding:PostgresNorthwind', friendlyName: 'PostgresNorthwind',
  dataSource: 'Postgres'};

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

/** Every dapi call the dialog makes, recorded; the dry run answers what `dryRun` says. */
function stub(calls, dryRun) {
  grok.dapi.connections = {list: async () => [CONN], getSchemas: async () => ['public', 'audit']};
  grok.dapi.permissions = {check: async () => true};
  grok.dapi.groups = {filter: () => ({list: async () =>
    [{id: 'g-sales', friendlyName: 'Sales'}, {id: 'g-dev', friendlyName: 'Developers'}]})};
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
      grant: async (group, permission) => calls.push(['table.grant', address, group, permission]),
      shareColumn: async (column, group, permission) => calls.push(['shareColumn', address, column, group, permission]),
      restrictColumn: async (column) => calls.push(['restrictColumn', address, column]),
    }),
  });
}

const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')]
  .find((b) => b.textContent === text);
const reason = () => document.querySelector('.u2-wizard-reason').textContent;
const status = () => document.querySelector('.u2-wizard-status');

scoped('connection › design › review: draft, dry run, create, grants and restrictions, the name', async () => {
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
  assert.equal(document.querySelector('.u2-binding-access').textContent, 'You may introspect and query this connection');
  assert.equal(buttonNamed('NEXT').disabled, false);

  dialog.wizard.next();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'design');
  const editor = dialog.editor;
  assert.ok(editor, 'the editor is built over the draft');
  assert.deepEqual(editor.model.tables.value.filter((t) => t.included).map((t) => t.remote), ['orders'],
    'only the named table starts included');
  assert.equal(buttonNamed('NEXT').disabled, false, 'the drafted name passes');
  editor.model.name.value = 'Bad Name';
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.match(reason(), /^Name: /);
  editor.model.name.value = 'northwind_sales';
  editor.access.addGroup('schema', 'Sales');
  editor.access.setGrant('schema', 'Sales', 'view', true);
  editor.access.addGroup('orders', 'Developers');
  editor.access.setGrant('orders', 'Developers', 'edit', true);
  editor.access.setVisibility('orders', 'freight', ['Sales']);
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
  assert.equal(await done, 'northwind_sales');
  const create = calls.find((c) => c[0] === 'create' && c[2].dryRun === undefined);
  assert.equal(create[1], 'northwind_sales');
  assert.deepEqual(Object.keys(create[2].manifest.tables), ['orders']);
  // a schema-scope row goes to the schema (the server fans it out), a table row to the table, a
  // visibility row to shareColumn — by group ID, never by name
  const access = calls.filter((c) => c[0].endsWith('grant') || c[0] === 'shareColumn').map((c) => c.join(' '));
  for (const expected of ['schema.grant northwind_sales g-sales View', 'table.grant northwind_sales.orders g-dev Edit',
    'shareColumn northwind_sales.orders freight g-sales View'])
    assert.equal(access.includes(expected), true, `${expected} in ${access.join(' | ')}`);
  assert.equal(access.some((c) => c.startsWith('schema.grant') && !c.includes('g-sales')), false);
  assert.equal(document.querySelector('.u2-dialog'), null, 'the dialog is gone');
  assert.equal(grok.dapi.domains.invalidated > 0, true);
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
  await flush();
  dialog.wizard.next();
  await flush();
  dialog.editor.model.name.value = 'x';
  dialog.wizard.next();
  await flush();
  buttonNamed('VALIDATE').click();
  await flush();
  buttonNamed('CREATE').click();
  await flush();
  assert.notEqual(document.querySelector('.u2-dialog'), null, 'still open');
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
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});
