/* The "Edit binding" dialog over stubbed dapi calls: the registry is read before the dialog shows
   (a refused manifest is a refusal by name, the dialog never opens), the other reads settle into
   what the model shows as unknown, the Design step opens over the registered manifest with a no-op
   gated, Review lists the changes in order and VALIDATE annotates them from the dry run bound to
   the exact apply body — what a removal takes with it, what is already so — a destructive plan
   wants the confirmation and SAVE sends it, a version conflict reloads the binding and replays the
   edits with the conflicts and drops reported for a look at Design, a save whose answer was lost
   is confirmed through the registry, and Saved reports the effects as resolved under the lock. */

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

const FIXTURE = JSON.parse(readFileSync(new URL('./fixtures/authoring/northwind-registered.json', import.meta.url), 'utf8'));
const copy = (x) => JSON.parse(JSON.stringify(x));
const NAME = 'northwind_sales';
const INCARNATION = '2026-09-20T10:00:00.000Z';
const ORDERS = {kind: 'table', table: 'orders'};
const SALES = {id: 'g-sales', label: 'Sales'};
const DEV = {id: 'g-dev', label: 'Developers'};
const ENTITY = {name: NAME, friendlyName: 'Northwind sales', description: 'Sales data'};

/** The fixture draft with the warehouse agreeing with everything registered: no blockers. */
function clean() {
  const d = copy(FIXTURE.draft);
  d.manifest.tables.shippers = {businessKey: ['shipperid'], columns: {shipperid: {type: 'int', required: true},
    companyname: {type: 'string'}}};
  d.inventory.tables.push({remote: 'shippers', logical: 'shippers', bindable: true, key: ['shipperid']});
  d.inventory.columns.push({table: 'shippers', remote: 'shipperid', dbType: 'int4', type: 'int'},
    {table: 'shippers', remote: 'companyname', dbType: 'text', type: 'string'},
    {table: 'orders', remote: 'shipaddress', dbType: 'text', type: 'string'});
  d.manifest.tables.orders.columns.shipaddress = {type: 'string'};
  d.manifest.tables.orders.columns.shipcity = {type: 'string'};
  const city = d.inventory.columns.find((c) => c.table === 'orders' && c.remote === 'shipcity');
  city.dbType = 'text';
  city.type = 'string';
  return d;
}

const PLAN = (extra = {}) => ({version: '4', destructive: false, registrationOnly: true,
  creates: {tables: [], columns: {}, uniques: {}, constraints: {}}, drops: {tables: [], columns: []}, typeChanges: [],
  uniqueDrops: [], requiredToggles: [], businessKeyChanges: [], autoNumberChanges: {}, alters: {}, constraintDrops: null,
  violations: [], refusals: [], keyColumnsUnrestricted: [], migration: null, ...extra});

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

/** What a read answers: a value copied, a function called, an Error thrown. */
const answer = (v) => {
  if (v instanceof Error)
    throw v;
  return typeof v === 'function' ? v() : copy(v);
};

/** Every dapi call the dialog makes, recorded; the server's answers are fields a test replaces
 * between steps. */
function stub(calls, fields = {}) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.shell.dart.user = {id: 'u-me', friendlyName: 'askalkin', group: {id: 'g-me'}};
  globalThis.grok_Dapi_Domains_SchemaAltered = (_dart, name) => calls.push(['altered', name]);
  grok.dapi.groups = {getGroupsLookup: async () => [], list: async () => []};
  const server = {manifest: copy(FIXTURE.manifest), snapshot: copy(FIXTURE.snapshot), draft: clean(), entity: ENTITY,
    plan: () => PLAN(), apply: (body) => ({...server.plan(body), applied: true}), ...fields};
  grok.dapi.domains.schemas = {filter: (text) => ({list: async () => {
    calls.push(['schemas', text]);
    const entity = answer(server.entity);
    return entity === null ? [] : [entity];
  }})};
  grok.dapi.domains.schema = (name) => ({
    manifest: async () => {
      calls.push(['manifest', name]);
      return answer(server.manifest);
    },
    access: async () => {
      calls.push(['access', name]);
      return answer(server.snapshot);
    },
    draft: async () => {
      calls.push(['draft', name]);
      return answer(server.draft);
    },
    apply: async (body, options) => {
      calls.push(['apply', name, copy(body), options]);
      return options?.dryRun === true ? server.plan(body) : server.apply(body);
    },
  });
  return server;
}

const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')]
  .find((b) => b.textContent === text);
const reason = () => document.querySelector('.u2-wizard-reason').textContent;
const status = () => document.querySelector('.u2-wizard-status');
const rows = (within = '.u2-binding-changes') => [...document.querySelectorAll(`${within} .u2-binding-change`)]
  .map((r) => [...r.children].map((c) => c.textContent));
const applies = (calls, dryRun) => calls.filter((c) => c[0] === 'apply' && (c[3]?.dryRun === true) === dryRun);

/** The dialog opened and built over the reads. */
async function opened() {
  const dialog = new domains.authoring.BindingEditDialog(NAME);
  const done = dialog.open();
  await flush();
  await flush();
  return {dialog, done};
}

async function validated() {
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(status().textContent, 'Validated');
}

scoped('editBinding refuses while domain databases are off', async () => {
  grok.shell.settings = {enableDomainDatabases: false};
  await assert.rejects(domains.authoring.editBinding(NAME), /Beta feature/);
  assert.equal(document.querySelector('.u2-dialog'), null);
});

scoped('a manifest the registry refuses is a refusal by name; the dialog does not open', async () => {
  const calls = [];
  stub(calls, {manifest: Object.assign(new Error('Domain schema "northwind_sales" not found'), {status: 404, code: 'not-found'})});
  await assert.rejects(domains.authoring.editBinding(NAME),
    /^Error: Binding northwind_sales could not be read: Domain schema "northwind_sales" not found$/);
  await flush();
  assert.equal(document.querySelector('.u2-dialog'), null);
  assert.equal(calls.some((c) => c[0] === 'draft'), false, 'the warehouse is not asked for a schema the registry refused');
});

scoped('a failed draft and snapshot open all the same: the catalog is unknown, access read-only, both said; the facts are read-only', async () => {
  const calls = [];
  stub(calls, {draft: Object.assign(new Error('You don\'t have DataConnection.Query on connection "x"'), {code: 'forbidden'}),
    snapshot: Object.assign(new Error('Internal error'), {status: 500, code: ''})});
  const {dialog, done} = await opened();
  assert.equal(dialog.wizard.currentStep.value, 'design', 'opened on Design, the Connection step a BACK away');
  assert.deepEqual(calls.map((c) => c[0]), ['manifest', 'access', 'schemas', 'draft'],
    'the registry reads at once, the warehouse once the manifest answered');
  const editor = dialog.editor;
  assert.equal(editor.model.catalog, 'unknown');
  assert.equal(editor.access.canEdit(ORDERS), false, 'no snapshot: access is never changed');
  assert.equal(status().textContent, 'The access snapshot could not be read: Internal error; ' +
    'The warehouse catalog could not be read: You don\'t have DataConnection.Query on connection "x"');
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.equal(reason(), 'Nothing changed');
  dialog.wizard.back();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'connection');
  assert.deepEqual([...document.querySelectorAll('.u2-binding-fact')].map((f) => f.textContent),
    ['ConnectionNorthwindBinding:PostgresNorthwind', 'Remote schemapublic', 'Catalog—']);
  assert.equal(document.querySelector('.u2-binding-connection select, .u2-binding-connection input'), null, 'no picker');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('design › review › saved: a no-op is gated, the changes are listed in order, VALIDATE annotates them from the dry run, SAVE is the one apply', async () => {
  const calls = [];
  const server = stub(calls);
  const effects = (revoke) => ({
    grant: [{table: 'orders', group: {id: 'g-dev', friendlyName: 'Developers'}, permission: 'View', effect: 'grant'}],
    revoke: [{table: 'orders', group: {id: 'g-sales', friendlyName: 'Sales'}, permission: 'Edit', effect: revoke}],
    restrict: [{table: 'orders', column: 'shipname', effect: 'restrict', revokes: [],
      grants: [{group: {id: 'g-sales', friendlyName: 'Sales'}, permission: 'View', effect: 'grant'},
        {group: {id: 'g-me', friendlyName: 'askalkin'}, permission: 'View', effect: 'none'},
        {group: {id: 'g-me', friendlyName: 'askalkin'}, permission: 'Edit', effect: 'grant'}]}],
    unrestrict: [],
  });
  server.plan = () => PLAN({creates: {tables: ['employees'], columns: {}, uniques: {}, constraints: {}},
    metadata: {friendlyName: {from: 'Northwind sales', to: 'NW'}}, access: effects('revoke')});
  server.apply = () => ({...PLAN({access: effects('none')}), applied: true});
  const {dialog, done} = await opened();
  assert.equal(status().textContent, '');
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.equal(reason(), 'Nothing changed');
  const editor = dialog.editor;
  assert.deepEqual([editor.model.friendlyName.value, editor.model.description.value], ['Northwind sales', 'Sales data'],
    'the entity\'s caption and description');
  assert.equal(editor.model.version, '3');
  editor.model.includeTable('employees', true);
  editor.model.setSchemaFriendlyName('NW');
  editor.access.addGroup(ORDERS, DEV);
  editor.access.setGrant(ORDERS, 'g-sales', 'edit', false);
  editor.access.setVisibility('orders', 'shipname', [SALES]);
  await flush();
  assert.equal(buttonNamed('NEXT').disabled, false);
  dialog.wizard.next();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.equal(document.querySelector('.u2-binding-changes-title').textContent, 'Changes, not validated');
  assert.deepEqual(rows(), [['Friendly name: "NW"'], ['Table employees added'], ['orders: View for Developers'],
    ['orders: Edit revoked from Sales'], ['orders.shipname restricted — visible to Sales and you; Edit for you alone']],
    'manifest order, then the access ops');
  assert.equal(buttonNamed('SAVE').disabled, true);
  assert.equal(reason(), 'Validate before saving');
  assert.match(document.querySelector('.u2-binding-body .u2-binding-json').textContent, /"ifVersion": "3"/);

  await validated();
  const [dry] = applies(calls, true);
  assert.deepEqual(dry[2], {ifVersion: '3', ifIncarnation: INCARNATION, friendlyName: 'NW',
    tables: {employees: {businessKey: ['employeeid'], columns: {employeeid: {type: 'int', required: true}, lastname: {type: 'string'}}}},
    access: {grant: [{table: 'orders', group: 'g-dev', permission: 'View'}],
      revoke: [{table: 'orders', group: 'g-sales', permission: 'Edit'}],
      restrict: [{table: 'orders', column: 'shipname', from: 'unrestricted', revoke: [], grant: [{group: 'g-sales', permission: 'View'},
        {group: 'g-me', permission: 'View'}, {group: 'g-me', permission: 'Edit'}]}],
      unrestrict: []}}, 'the apply body verbatim, the author kept on the restricted column, the state each column op was made from');
  assert.equal(document.querySelector('.u2-binding-changes-title').textContent, 'Changes, as validated');
  assert.deepEqual(rows(), [['Friendly name: "NW"', 'was "Northwind sales"'], ['Table employees added'],
    ['orders: View for Developers'], ['orders: Edit revoked from Sales'],
    ['orders.shipname restricted — visible to Sales and you; Edit for you alone', 'View for askalkin: already so']]);
  assert.equal(document.querySelector('.u2-binding-confirm').style.display, 'none', 'nothing destructive');
  assert.equal(buttonNamed('SAVE').disabled, false);

  buttonNamed('SAVE').click();
  await flush();
  const [save] = applies(calls, false);
  assert.deepEqual(save[2], dry[2], 'the exact body validated, nothing added');
  assert.equal(dialog.wizard.currentStep.value, 'saved');
  assert.equal(document.querySelector('.u2-binding-created-head').textContent, 'SavedDomain schema northwind_sales is at version 4');
  assert.deepEqual(rows('.u2-binding-saved'), [['Friendly name: "NW"'], ['Table employees added'], ['orders: View for Developers'],
    ['orders: Edit revoked from Sales', 'already so'],
    ['orders.shipname restricted — visible to Sales and you; Edit for you alone', 'View for askalkin: already so']],
    'the effects as resolved under the lock');
  assert.equal(status().textContent, 'Saved');
  assert.equal(grok.dapi.domains.invalidated > 0, true);
  assert.deepEqual(calls.filter((c) => c[0] === 'altered'), [], 'an answered apply needs no announcement');
  assert.equal(buttonNamed('CANCEL').style.display, 'none');
  buttonNamed('CLOSE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, NAME);
  assert.equal(result.applied.applied, true);
  assert.equal(document.querySelector('.u2-dialog'), null);
});

scoped('a validation is bound to the exact apply body: an access edit since is not validated; changed back, it is', async () => {
  const calls = [];
  stub(calls);
  const {dialog, done} = await opened();
  const editor = dialog.editor;
  editor.model.setDescription('Sales data, revised');
  dialog.wizard.next();
  await flush();
  await validated();
  assert.equal(buttonNamed('SAVE').disabled, false);
  dialog.wizard.back();
  await flush();
  editor.access.addGroup(ORDERS, DEV);
  dialog.wizard.next();
  await flush();
  assert.equal(status().textContent, '', 'Review after an edit: not validated');
  assert.equal(buttonNamed('SAVE').disabled, true);
  assert.deepEqual(rows().map((r) => r[0]), ['Description: "Sales data, revised"', 'orders: View for Developers']);
  dialog.wizard.back();
  await flush();
  editor.access.removeGroup(ORDERS, 'g-dev');
  dialog.wizard.next();
  await flush();
  assert.equal(status().textContent, 'Validated', 'the same body again is what was validated');
  assert.equal(buttonNamed('SAVE').disabled, false);
  assert.equal(applies(calls, true).length, 1);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('an apply the server finds already held is not a save: nothing to save, the version unchanged', async () => {
  const calls = [];
  const server = stub(calls);
  server.apply = () => ({...PLAN({version: '3'}), applied: false, noop: true});
  const {dialog, done} = await opened();
  dialog.editor.model.setDescription('Sales data, revised');
  dialog.wizard.next();
  await flush();
  await validated();
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'saved');
  assert.equal(document.querySelector('.u2-binding-created-head').textContent,
    'UnchangedNothing to save — the binding already holds this; domain schema northwind_sales stays at version 3');
  assert.equal(status().textContent, 'Nothing to save — the binding already holds this');
  assert.equal(document.querySelectorAll('.u2-binding-saved .u2-binding-change').length, 0, 'nothing listed as applied');
  buttonNamed('CLOSE').click();
  await flush();
  const result = await done;
  assert.equal(result.applied.noop, true);
  assert.equal(result.applied.version, '3');
});

scoped('a removal is annotated with what goes with it; a destructive plan wants the confirmation, and SAVE sends it', async () => {
  const calls = [];
  const server = stub(calls);
  server.plan = () => PLAN({destructive: true, drops: {tables: [{name: 'shippers'}], columns: [{table: 'orders', column: 'shipcountry'}]},
    lost: {tables: {shippers: {grants: 3, coreSchemaGrants: 0, restrictions: {companyname: 1}, promotedRows: 0, rowGrants: 0,
      savedFilters: 2, affectedFilters: 0}},
    columns: {'orders.shipcountry': {restricted: true, grants: 2, affectedFilters: 1}}, refs: {}}});
  const {dialog, done} = await opened();
  const editor = dialog.editor;
  editor.model.includeTable('shippers', false);
  editor.model.includeColumn('orders', 'shipcountry', false);
  dialog.wizard.next();
  await flush();
  assert.deepEqual(rows(), [['Column orders.shipcountry removed'], ['Table shippers removed']]);
  assert.equal(document.querySelectorAll('.u2-binding-change-removes').length, 2);
  await validated();
  assert.deepEqual(rows(), [
    ['Column orders.shipcountry removed', 'its restriction and all 2 grants on it go with it, whoever holds them; ' +
      '1 saved filter naming it no longer resolves; the warehouse is untouched'],
    ['Table shippers removed', '3 grants, 1 restricted column, 2 saved filters go with it; the warehouse is untouched']]);
  assert.equal(applies(calls, true)[0][2].confirmDestructive, undefined, 'the dry run is the body as it stands');
  const confirm = document.querySelector('.u2-binding-confirm');
  assert.equal(confirm.style.display, '');
  assert.equal(buttonNamed('SAVE').disabled, true);
  assert.equal(reason(), 'Confirm what is removed');
  const box = confirm.querySelector('.u2-input-checkbox');
  box.checked = true;
  fire(box, 'change');
  await flush();
  assert.equal(buttonNamed('SAVE').disabled, false);
  buttonNamed('SAVE').click();
  await flush();
  const [save] = applies(calls, false);
  assert.deepEqual(save[2], {...applies(calls, true)[0][2], confirmDestructive: true});
  assert.equal(dialog.wizard.currentStep.value, 'saved');
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).name, NAME);
});

scoped('a confirmation holds for the plan it was given: an edit and a new validation want it again', async () => {
  const calls = [];
  const server = stub(calls);
  server.plan = () => PLAN({destructive: true, drops: {tables: [{name: 'shippers'}], columns: []},
    lost: {tables: {shippers: {grants: 0, coreSchemaGrants: 0, restrictions: {}, promotedRows: 0, rowGrants: 0,
      savedFilters: 0, affectedFilters: 0}}, columns: {}, refs: {}}});
  const {dialog, done} = await opened();
  dialog.editor.model.includeTable('shippers', false);
  dialog.wizard.next();
  await flush();
  await validated();
  const box = () => document.querySelector('.u2-binding-confirm .u2-input-checkbox');
  box().checked = true;
  fire(box(), 'change');
  await flush();
  assert.equal(buttonNamed('SAVE').disabled, false);
  dialog.wizard.back();
  await flush();
  dialog.editor.model.setDescription('and a note');
  dialog.wizard.next();
  await flush();
  await validated();
  assert.equal(box().checked, false, 'the tick went with the plan it was given');
  assert.equal(buttonNamed('SAVE').disabled, true);
  assert.equal(reason(), 'Confirm what is removed');
  assert.deepEqual(rows()[1], ['Table shippers removed', 'nothing else goes with it; the warehouse is untouched']);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('a version conflict reloads the binding and replays the edits: conflicts and drops want a look at Design, then a new validation', async () => {
  const calls = [];
  const server = stub(calls);
  let conflicts = 1;
  server.apply = () => {
    if (conflicts-- > 0) {
      throw Object.assign(new Error('Domain schema "northwind_sales" is at version 4, not 3'),
        {status: 409, code: 'version-conflict', currentVersion: '4', expectedVersion: '3'});
    }
    return {...PLAN({version: '5'}), applied: true};
  };
  const {dialog, done} = await opened();
  const editor = dialog.editor;
  editor.model.setFriendlyName('orders', 'Sales orders');
  editor.model.includeTable('employees', true);
  editor.access.addGroup({kind: 'table', table: 'shippers'}, DEV);
  dialog.wizard.next();
  await flush();
  await validated();
  // meanwhile: orders was re-captioned, shippers unregistered (still in the warehouse)
  const nb = copy(FIXTURE.manifest);
  nb.version = '4';
  nb.tables.orders.friendlyName = 'Orders (EU)';
  delete nb.tables.shippers;
  const ns = copy(FIXTURE.snapshot);
  ns.version = '4';
  delete ns.tables.shippers;
  for (const key of Object.keys(ns.columns).filter((k) => k.startsWith('shippers.')))
    delete ns.columns[key];
  Object.assign(server, {manifest: nb, snapshot: ns});
  calls.length = 0;
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.deepEqual(calls.map((c) => c[0]), ['apply', 'manifest', 'access', 'schemas', 'draft'], 'one apply, then the reads — no retry');
  assert.equal(status().textContent, 'Reloaded at version 4: 1 edit kept, 1 conflict, 1 dropped — validate again');
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.deepEqual(rows('.u2-binding-reloaded'), [
    ['Table orders: friendly name "Sales orders"', 'the server\'s "Orders (EU)" stands; yours was "Sales orders"'],
    ['shippers: View for Developers', 'no longer there']]);
  assert.deepEqual(rows().map((r) => r[0]), ['Table employees added'], 'the change list is what still applies');
  assert.equal(editor.model.version, '4');
  assert.equal(editor.model.table('orders').friendlyName, 'Orders (EU)');
  assert.equal(buttonNamed('SAVE').disabled, true);
  assert.equal(reason(), 'Look over the reloaded edits on Design');
  buttonNamed('VALIDATE').click();
  await flush();
  assert.equal(reason(), 'Look over the reloaded edits on Design', 'a validation alone is no look');

  dialog.wizard.back();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'design');
  dialog.wizard.next();
  await flush();
  assert.equal(reason(), '', 'looked at, and validated over the reloaded binding');
  assert.equal(status().textContent, 'Validated');
  assert.equal(buttonNamed('SAVE').disabled, false);
  assert.equal(rows('.u2-binding-reloaded').length, 2, 'the report stays until the save');
  assert.equal(applies(calls, true).at(-1)[2].ifVersion, '4');
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'saved');
  assert.equal(applies(calls, false).at(-1)[2].ifVersion, '4');
  assert.match(document.querySelector('.u2-binding-created-head').textContent, /version 5$/);
  buttonNamed('CLOSE').click();
  await flush();
  assert.equal((await done).applied.version, '5');
});

const timeout = () => {
  throw Object.assign(new Error('Gateway Timeout'), {status: 504, code: ''});
};

/** The dialog on Review with the description edited and validated. */
async function edited() {
  const {dialog, done} = await opened();
  dialog.editor.model.setDescription('Moved');
  dialog.wizard.next();
  await flush();
  await validated();
  return {dialog, done};
}

const outcome = () => document.querySelector('.u2-binding-outcome');

scoped('a save whose answer was lost is confirmed through the registry — one version on, the same incarnation, this edit in it — never replayed', async () => {
  const calls = [];
  const server = stub(calls, {apply: timeout});
  const {dialog, done} = await edited();
  Object.assign(server, {manifest: {...copy(FIXTURE.manifest), version: '4'}, entity: {...ENTITY, description: 'Moved'}});
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'saved');
  assert.equal(applies(calls, false).length, 1, 'one apply; nothing is replayed by the dialog');
  assert.deepEqual(calls.filter((c) => c[0] === 'altered'), [['altered', NAME]], 'the confirmed save is announced to the platform once');
  assert.equal(document.querySelector('.u2-binding-created-head').textContent,
    'SavedDomain schema northwind_sales is at version 4 — the answer was lost, the registry confirms it');
  assert.deepEqual(rows('.u2-binding-saved'), [['Description: "Moved"']], 'no effects to report');
  buttonNamed('CLOSE').click();
  await flush();
  const result = await done;
  assert.equal(result.name, NAME);
  assert.equal(result.applied, null);
});

scoped('a lost answer the registry cannot confirm ends the edit: the outcome is named, nothing is replayed, and REOPEN starts over from the server\'s version', async () => {
  const calls = [];
  const server = stub(calls, {apply: timeout});
  const {dialog, done} = await edited();
  // the registry did not move: nobody saw the save land
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(dialog.wizard.currentStep.value, 'review');
  assert.equal(status().textContent, 'Outcome unknown — see Review');
  assert.equal(status().classList.contains('u2-wizard-status-error'), true);
  assert.deepEqual(dialog.editor.diagnostics.value, [], 'the timeout is no finding');
  assert.equal(buttonNamed('SAVE'), undefined, 'SAVE is gone');
  assert.equal(buttonNamed('REOPEN').disabled, false);
  assert.equal(buttonNamed('VALIDATE').style.display, 'none', 'nothing left to validate');
  assert.equal(document.querySelector('.u2-binding-changes-title').textContent, 'Changes');
  assert.equal(document.querySelector('.u2-binding-issues').style.display, 'none');
  assert.deepEqual(rows('.u2-binding-outcome'), [['No answer to the save (Gateway Timeout) — whether it landed is unknown'],
    ['northwind_sales is at version 3 (this edit was made against 3)']]);
  assert.match(outcome().textContent, /^Outcome unknown/);
  assert.match(outcome().textContent, /Nothing was replayed\. REOPEN starts over/);
  assert.deepEqual(calls.filter((c) => c[0] === 'altered'), [['altered', NAME]], 'it may have landed: the platform reads the binding again');
  assert.equal(grok.dapi.domains.invalidated > 0, true);
  dialog.wizard.back();
  await flush();
  assert.match(document.querySelector('.u2-binding-design .u2-binding-ended').textContent, /^Outcome unknown: nothing here will be saved/);
  dialog.wizard.next();
  await flush();
  assert.equal(buttonNamed('REOPEN').disabled, false, 'a look at Design changes nothing: the edit ended');
  assert.equal(buttonNamed('VALIDATE').style.display, 'none');
  assert.equal(document.querySelector('.u2-binding-changes-title').textContent, 'Changes');
  assert.equal(applies(calls, false).length, 1);

  // REOPEN reads the registry first: a failed read leaves this dialog standing, and says so
  calls.length = 0;
  server.manifest = Object.assign(new Error('boom'), {status: 500, code: ''});
  buttonNamed('REOPEN').click();
  await flush();
  await flush();
  assert.equal(document.querySelectorAll('.u2-dialog').length, 1);
  assert.equal(dialog.wizard.currentStep.value, 'review', 'the same dialog');
  assert.deepEqual(calls.map((c) => c[0]), ['manifest', 'access', 'schemas'], 'the registry reads, no warehouse read');
  assert.equal(status().textContent, 'Binding northwind_sales could not be read: boom — REOPEN again');
  assert.equal(buttonNamed('REOPEN').disabled, false);

  server.manifest = copy(FIXTURE.manifest);
  calls.length = 0;
  buttonNamed('REOPEN').click();
  await flush();
  await flush();
  assert.equal(document.querySelectorAll('.u2-dialog').length, 1, 'one dialog: the new one');
  assert.deepEqual(calls.map((c) => c[0]), ['manifest', 'access', 'schemas', 'draft'], 'read afresh');
  assert.equal(reason(), 'Nothing changed', 'nothing carried over');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null, 'the caller gets the reopened dialog\'s outcome');
});

scoped('a lost answer is unknown too when the next version is someone else\'s save, or the name another schema\'s', async () => {
  const calls = [];
  const server = stub(calls, {apply: timeout});
  const first = await edited();
  server.manifest = {...copy(FIXTURE.manifest), version: '4'};
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(status().textContent, 'Outcome unknown — see Review');
  assert.deepEqual(rows('.u2-binding-outcome')[1], ['version 4 of northwind_sales holds someone else\'s save, not this edit']);
  assert.equal(buttonNamed('REOPEN').disabled, false);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await first.done, null);

  server.manifest = copy(FIXTURE.manifest);
  const second = await edited();
  Object.assign(server, {manifest: {...copy(FIXTURE.manifest), version: '4', incarnation: '2026-09-25T00:00:00.000Z'},
    entity: {...ENTITY, description: 'Moved'}});
  buttonNamed('SAVE').click();
  await flush();
  assert.equal(status().textContent, 'Outcome unknown — see Review');
  assert.deepEqual(rows('.u2-binding-outcome')[1],
    ['northwind_sales was deleted and re-created since — another schema is registered under the name now, at version 4']);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await second.done, null);
});

scoped('an access-conflict reloads and replays like a version conflict; a conflict from another incarnation, or a reload that finds one, ends the edit without a replay', async () => {
  const calls = [];
  const server = stub(calls);
  let answer = () => {
    throw Object.assign(new Error('orders.freight is unrestricted now'), {status: 409, code: 'access-conflict'});
  };
  server.apply = () => answer();
  const {dialog, done} = await opened();
  dialog.editor.access.setVisibility('orders', 'freight', [SALES, {id: 'g-me', label: 'askalkin'}, DEV]);
  dialog.wizard.next();
  await flush();
  await validated();
  assert.equal(applies(calls, true)[0][2].access.restrict[0].from, 'restricted');
  const nb = copy(FIXTURE.manifest);
  nb.version = '4';
  const ns = copy(FIXTURE.snapshot);
  ns.version = '4';
  ns.columns['orders.freight'] = {state: 'unrestricted', canShare: true};
  Object.assign(server, {manifest: nb, snapshot: ns});
  calls.length = 0;
  buttonNamed('SAVE').click();
  await flush();
  assert.deepEqual(calls.map((c) => c[0]), ['apply', 'manifest', 'access', 'schemas', 'draft']);
  assert.equal(status().textContent, 'Reloaded at version 4: 0 edits kept, 1 conflict, 0 dropped — validate again');
  assert.deepEqual(rows('.u2-binding-reloaded'),
    [['orders.freight: View for Developers', 'the column was made visible to everyone meanwhile']], 'by name');
  assert.equal(reason(), 'Look over the reloaded edits on Design');

  answer = () => {
    throw Object.assign(new Error('Schema version conflict'), {status: 409, code: 'version-conflict', currentVersion: '1',
      expectedVersion: '4', currentIncarnation: '2026-09-25T00:00:00.000Z', expectedIncarnation: INCARNATION});
  };
  dialog.wizard.back();
  await flush();
  dialog.editor.model.setDescription('again');
  dialog.wizard.next();
  await flush();
  await validated();
  calls.length = 0;
  buttonNamed('SAVE').click();
  await flush();
  assert.deepEqual(calls.map((c) => c[0]), ['apply'], 'no reload, no replay');
  assert.equal(status().textContent, 'The binding you edited is gone — see Review');
  assert.deepEqual(rows('.u2-binding-outcome'),
    [['northwind_sales was deleted and re-created since — another schema is registered under the name now, at version 1']]);
  assert.equal(buttonNamed('REOPEN').disabled, false);
  assert.equal(buttonNamed('SAVE'), undefined);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);

  // a conflict without the incarnations, whose reload finds another one
  answer = () => {
    throw Object.assign(new Error('Schema version conflict'), {status: 409, code: 'version-conflict', currentVersion: '1', expectedVersion: '4'});
  };
  const second = await edited();
  server.manifest = {...nb, version: '1', incarnation: '2026-09-25T00:00:00.000Z'};
  calls.length = 0;
  buttonNamed('SAVE').click();
  await flush();
  assert.deepEqual(calls.map((c) => c[0]), ['apply', 'manifest', 'access', 'schemas', 'draft']);
  assert.equal(status().textContent, 'The binding you edited is gone — see Review');
  assert.equal(second.dialog.editor.model.version, '4', 'the editor was not rebased onto the other schema');
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await second.done, null);
});

scoped('an unrestriction names the state and the ACL revision it was made from', async () => {
  const calls = [];
  stub(calls);
  const {dialog, done} = await opened();
  dialog.editor.access.setVisibility('orders', 'freight', null);
  dialog.wizard.next();
  await flush();
  await validated();
  assert.deepEqual(applies(calls, true)[0][2].access, {grant: [], revoke: [], restrict: [],
    unrestrict: [{table: 'orders', column: 'freight', from: 'restricted', revision: 'rev-freight'}]});
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});

scoped('Escape cancels the dialog even once focus fell out of it onto the body', async () => {
  const calls = [];
  stub(calls);
  const {dialog, done} = await opened();
  dialog.editor.model.setDescription('never saved');
  document.body.focus();
  fire(document.body, 'keydown', {key: 'Escape'});
  await flush();
  assert.equal(document.querySelector('.u2-dialog'), null);
  assert.equal(await done, null);
  assert.deepEqual(applies(calls, false), []);
});

scoped('blockers gate Design and VALIDATE until the items the warehouse no longer binds are taken out', async () => {
  const calls = [];
  stub(calls, {draft: copy(FIXTURE.draft)});
  const {dialog, done} = await opened();
  const editor = dialog.editor;
  assert.equal(editor.model.blockers.value.length, 3, 'shippers missing, shipaddress missing, shipcity drifted');
  editor.model.setDescription('x');
  await flush();
  assert.equal(buttonNamed('NEXT').disabled, true);
  assert.equal(reason(), '3 registered items cannot bind as declared — take them out, or fix the warehouse');
  editor.model.includeTable('shippers', false);
  editor.model.includeColumn('orders', 'shipaddress', false);
  editor.model.includeColumn('orders', 'shipcity', false);
  await flush();
  assert.equal(buttonNamed('NEXT').disabled, false);
  dialog.wizard.next();
  await flush();
  assert.equal(buttonNamed('VALIDATE').disabled, false);
  buttonNamed('CANCEL').click();
  await flush();
  assert.equal(await done, null);
});
