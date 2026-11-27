import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {layout} from '@datagrok-libraries/u2/tests/layout.mjs';
import {scoped} from './domain-fixtures.mjs';
import {byKind, byName, state, user} from '@datagrok-libraries/u2/testing';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');
const {domains} = await import('../src/dg/facade.js');

const copy = (x) => JSON.parse(JSON.stringify(x));
const NAME = 'grit';
const MANIFEST = {name: NAME, version: '5', tables: {issue: {extensible: true,
  columns: {title: {type: 'string'}, mine: {type: 'string', origin: 'user'}}}}};
const PLAN = (extra) => ({version: '5', destructive: false, drops: {tables: [], columns: []}, typeChanges: [],
  uniqueDrops: [], autoNumberChanges: {}, businessKeyChanges: [], constraintDrops: {}, ...extra});

function stub(calls, server) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.dapi.permissions = {check: async () => true};
  DG.Property.propertyOptions = Object.fromEntries(['friendlyName', 'description', 'choices', 'min', 'max', 'format']
    .map((name) => [name, {name, type: name === 'choices' ? 'list' : 'string'}]));
  const tables = () => Object.entries(server.manifest.tables).map(([name, t]) => server.table?.() ?? {name,
    origin: 'package', userColumns: Object.keys(t.columns).filter((c) => t.columns[c].origin === 'user')});
  const named = (entity) => ({filter: (f) => ({list: async () => f === `name = "${NAME}"` ? [entity()] : []})});
  grok.dapi.domains.schemas = {...named(() => ({name: NAME, managedBy: 'package'})),
    include: () => named(() => ({name: NAME, managedBy: 'package', tables: tables()}))};
  grok.dapi.domains.schema = (name) => ({name,
    manifest: async () => JSON.parse(JSON.stringify(server.manifest,
      (key, value) => key === 'origin' ? undefined : value)),
    apply: async (body, options) => {
      calls.push(['apply', copy(body), options?.dryRun === true]);
      if (body.ifVersion !== undefined && body.ifVersion !== String(server.extVersion)) {
        throw Object.assign(new Error('Schema version conflict'),
          {status: 409, code: 'version-conflict', body: {error: 'version-conflict'}});
      }
      return options?.dryRun === true ? PLAN({extVersion: server.extVersion, ...(body.extend && server.plan)}) :
        PLAN({applied: true, extVersion: ++server.extVersion});
    },
  });
}

const extensions = (calls) => calls.filter((c) => c[0] === 'apply' && c[1].extend !== undefined);
const dialog = () => byKind(document, 'dialog')[0];

async function openAndAdd() {
  const done = domains.editProperties(NAME, 'issue');
  await flush();
  await flush();
  const [tree] = byKind(document, 'manifest-tree');
  layout(byKind(tree, 'list')[0], {clientHeight: 800});
  await flush();
  byName(tree, '+ Add property').click();
  await flush();
  return {done};
}

async function save() {
  byName(dialog(), 'SAVE').click();
  await flush();
  await flush();
}

scoped('Edit properties saves against the extension version it was opened at: a concurrent one is a conflict', async () => {
  const calls = [];
  const server = {manifest: MANIFEST, extVersion: 96};
  stub(calls, server);
  const {done} = await openAndAdd();
  server.extVersion = 97;
  await save();
  assert.deepEqual(extensions(calls).map((c) => c[1].ifVersion), ['96']);
  assert.match(byKind(document, 'status-line')[0].textContent, /Someone else changed these properties meanwhile/);
  byName(dialog(), 'CANCEL').click();
  assert.equal(await done, false);
});

scoped('on a server without DomainTable.userColumns Edit properties refuses, and nothing is sent', async () => {
  const calls = [];
  stub(calls, {manifest: MANIFEST, extVersion: 96, table: () => ({name: 'issue'})});
  await assert.rejects(domains.editProperties(NAME, 'issue'), /unavailable on this server version/);
  assert.deepEqual(extensions(calls), []);
});

scoped('the confirmation lists the drops the server plans, not the client diff, and SAVE waits for consent', async () => {
  const calls = [];
  stub(calls, {manifest: MANIFEST, extVersion: 96, plan: {destructive: true,
    drops: {tables: [], columns: [{table: 'issue', column: 'mine', nonNullValues: 3}]},
    typeChanges: [{table: 'issue', column: 'title', to: 'int'}]}});
  const {done} = await openAndAdd();
  await save();
  const confirm = byKind(document, 'dialog').find((d) => state(d).name === 'properties-confirm');
  assert.match(confirm.textContent, /issue\.mine dropped.*3 values deleted/);
  assert.match(confirm.textContent, /issue\.title: type changes to int/);
  const commits = () => extensions(calls).filter((c) => !c[2]);
  byName(confirm, 'SAVE').click();
  await flush();
  assert.deepEqual(commits(), [], 'nothing committed before the deletion is confirmed');
  user.click(byName(byName(confirm, 'confirmDestructive'), 'Delete 3 stored values and change the type of 1 column'));
  await flush();
  byName(confirm, 'SAVE').click();
  await flush();
  await flush();
  assert.equal(commits()[0]?.[1].confirmDestructive, true);
  assert.deepEqual(Object.keys(commits()[0][1].extend.issue.columns).sort(), ['mine', 'property'],
    'the user\'s saved columns go back with it: the extension is full state');
  assert.equal(await done, true);
});
