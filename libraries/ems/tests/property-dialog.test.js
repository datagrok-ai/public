import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {layout} from '@datagrok-libraries/u2/tests/layout.mjs';
import {Scope, notify} from '@datagrok-libraries/u2';
import {byKind, byName, state, user} from '@datagrok-libraries/u2/testing';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');
const {domains} = await import('../src/dg/facade.js');

const copy = (x) => JSON.parse(JSON.stringify(x));
const NAME = 'grit';
// a package schema as manifest() answers it: no extVersion, the server never serves one
const MANIFEST = {name: NAME, version: '5', tables: {issue: {extensible: true,
  columns: {title: {type: 'string'}, mine: {type: 'string', origin: 'user'}}}}};
const PLAN = (extra) => ({version: '5', destructive: false, drops: {tables: [], columns: []}, ...extra});

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

// the extension path as the server runs it: a plan echoes the current ext_version, any other ifVersion is a
// version conflict, a commit moves it on; `server.plan` is what the server plans for an extension beyond that
function stub(calls, server) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.dapi.permissions = {check: async () => true};
  DG.Property.propertyOptions = Object.fromEntries(['friendlyName', 'description', 'choices', 'min', 'max', 'format']
    .map((name) => [name, {name, type: name === 'choices' ? 'list' : 'string'}]));
  // ownership is the registry's word (DomainTable.userColumns): the `origin: 'user'` stamps say which columns it
  // lists, and manifest() answers without them, as the server does; `server.table` lists an older server's table
  const tables = () => Object.entries(server.manifest.tables).map(([name, t]) => server.table?.() ?? {name,
    origin: 'package', userColumns: Object.keys(t.columns).filter((c) => t.columns[c].origin === 'user')});
  const list = async () => (calls.push(['ownership']), [{name: NAME, managedBy: 'package', tables: tables()}]);
  grok.dapi.domains.schemas = {list: async () => [{name: NAME, managedBy: 'package'}], include: () => ({list})};
  grok.dapi.domains.schema = (name) => ({name,
    manifest: async () => (calls.push(['manifest']), JSON.parse(JSON.stringify(server.manifest,
      (key, value) => key === 'origin' ? undefined : value))),
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

scoped('EQ-027: Edit properties saves against the extension version it was opened at: a concurrent one is a conflict', async () => {
  const calls = [];
  const server = {manifest: MANIFEST, extVersion: 96};
  stub(calls, server);
  const {done} = await openAndAdd();
  const opening = calls.filter((c) => c[1]?.extend === undefined);
  server.extVersion = 97; // someone else extends the table meanwhile
  await save();
  assert.deepEqual(extensions(calls).map((c) => c[1].ifVersion), ['96'],
    'every apply of the save names the extension version the dialog was opened at');
  assert.deepEqual(opening, [['manifest'], ['apply', {}, true], ['manifest'], ['ownership']],
    'the version is read before the manifest and the ownership: a change between them is a conflict, not a drop');
  assert.match(byKind(document, 'status-line')[0].textContent, /Someone else changed these properties meanwhile/);
  byName(dialog(), 'CANCEL').click();
  assert.equal(await done, false, 'nothing saved');
});

scoped('EQ-018: a saved user column stays the user\'s, and adding the next column sends it back', async () => {
  const calls = [];
  stub(calls, {manifest: MANIFEST, extVersion: 96});
  const {done} = await openAndAdd();
  const mine = [...document.querySelectorAll('[role="treeitem"]')].find((n) => n.textContent.startsWith('mine'));
  assert.equal(mine.querySelector('.fa-lock'), null, 'the user column is not locked as a plugin column');
  await save();
  const [commit] = extensions(calls).filter((c) => !c[2]);
  assert.deepEqual(Object.keys(commit[1].extend.issue.columns).sort(), ['mine', 'property'],
    'the extension is the user\'s columns as full state: the saved one is not dropped');
  assert.equal(await done, true);
});

// datagrok-api comes from the server: an older one cannot say who owns a column, and a full-state extension sent
// without knowing would drop the user's columns
scoped('EQ-018: on a server without DomainTable.userColumns Edit properties refuses, and nothing is sent', async () => {
  const calls = [];
  stub(calls, {manifest: MANIFEST, extVersion: 96, table: () => ({name: 'issue'})});
  await assert.rejects(domains.editProperties(NAME, 'issue'), /unavailable on this server version/);
  assert.deepEqual(extensions(calls), [], 'no extension sent');
});

scoped('EQ-028: the confirmation lists the drops the server plans, not the client diff, and SAVE waits for consent', async () => {
  const calls = [];
  // the client only adds a column; the server's plan also drops one that holds values
  stub(calls, {manifest: MANIFEST, extVersion: 96, plan: {destructive: true,
    drops: {tables: [], columns: [{table: 'issue', column: 'mine', nonNullValues: 3}]}}});
  const {done} = await openAndAdd();
  await save();
  const confirm = byKind(document, 'dialog').find((d) => state(d).name === 'properties-confirm');
  assert.match(confirm.textContent, /issue\.mine dropped.*3 values deleted/,
    'the confirmation lists the plan\'s drop with its value count');
  const commits = () => extensions(calls).filter((c) => !c[2]);
  byName(confirm, 'SAVE').click();
  await flush();
  assert.deepEqual(commits(), [], 'nothing committed before the deletion is confirmed');
  user.click(byName(byName(confirm, 'confirmDestructive'), 'Delete 3 stored values'));
  await flush();
  byName(confirm, 'SAVE').click();
  await flush();
  await flush();
  assert.equal(commits()[0]?.[1].confirmDestructive, true, 'the confirmed save commits');
  assert.equal(await done, true);
});
