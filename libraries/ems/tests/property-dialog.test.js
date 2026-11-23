import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {layout} from '@datagrok-libraries/u2/tests/layout.mjs';
import {Scope, notify} from '@datagrok-libraries/u2';
import {byKind, byName} from '@datagrok-libraries/u2/testing';

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
// version conflict, a commit moves it on
function stub(calls, server) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.dapi.permissions = {check: async () => true};
  DG.Property.propertyOptions = Object.fromEntries(['friendlyName', 'description', 'choices', 'min', 'max', 'format']
    .map((name) => [name, {name, type: name === 'choices' ? 'list' : 'string'}]));
  grok.dapi.domains.schemas = {list: async () => [{name: NAME, managedBy: 'package'}]};
  grok.dapi.domains.schema = (name) => ({name,
    manifest: async () => (calls.push(['manifest']), copy(server.manifest)),
    apply: async (body, options) => {
      calls.push(['apply', copy(body), options?.dryRun === true]);
      if (body.ifVersion !== undefined && body.ifVersion !== String(server.extVersion)) {
        throw Object.assign(new Error('Schema version conflict'),
          {status: 409, code: 'version-conflict', body: {error: 'version-conflict'}});
      }
      return options?.dryRun === true ? PLAN({extVersion: server.extVersion}) :
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
  assert.deepEqual(opening, [['manifest'], ['apply', {}, true], ['manifest']],
    'the version is read before the manifest: a change between the reads is a conflict, never an overwrite');
  assert.match(byKind(document, 'status-line')[0].textContent, /Someone else changed these properties meanwhile/);
  byName(dialog(), 'CANCEL').click();
  assert.equal(await done, false, 'nothing saved');
});
