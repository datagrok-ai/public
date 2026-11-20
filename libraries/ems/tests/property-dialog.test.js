import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {layout} from '@datagrok-libraries/u2/tests/layout.mjs';
import {Scope, notify} from '@datagrok-libraries/u2';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');
const {domains} = await import('../src/dg/facade.js');

const copy = (x) => JSON.parse(JSON.stringify(x));
const NAME = 'grit';
/** A package schema as manifest() answers it: no extVersion (the server never serves one). */
const MANIFEST = {name: NAME, version: '5', incarnation: '2026-09-20T10:00:00.000Z', tables: {issue: {extensible: true,
  columns: {title: {type: 'string'}, mine: {type: 'string', origin: 'user'}}}}};
const PLAN = (extra = {}) => ({version: '5', destructive: false,
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

/** The extension path as the server runs it: a plan echoes the current ext_version, an `ifVersion` other than it
 * is a version conflict, a commit moves it on. */
function stub(calls, server) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.dapi.permissions = {check: async () => true};
  DG.Property.propertyOptions = Object.fromEntries(['friendlyName', 'description', 'choices', 'min', 'max', 'format']
    .map((name) => [name, {name, type: name === 'choices' ? 'list' : 'string'}]));
  grok.dapi.domains.schemas = {list: async () => [{name: NAME, managedBy: 'package'}]};
  grok.dapi.domains.schema = (name) => ({
    name,
    manifest: async () => {
      calls.push(['manifest']);
      return copy(server.manifest);
    },
    apply: async (body, options) => {
      calls.push(['apply', copy(body), options?.dryRun === true]);
      if (body.ifVersion !== undefined && body.ifVersion !== String(server.extVersion)) {
        throw Object.assign(new Error(`Schema version conflict: expected ${body.ifVersion}, current ${server.extVersion}`),
          {status: 409, code: 'version-conflict', body: {error: 'version-conflict'}});
      }
      if (options?.dryRun === true)
        return PLAN({extVersion: server.extVersion});
      return {...PLAN(), applied: true, extVersion: ++server.extVersion};
    },
  });
}

const commits = (calls) => calls.filter((c) => c[0] === 'apply' && !c[2] && Object.keys(c[1]).length > 0);
const button = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);

async function addProperty() {
  layout(document.querySelector('.u2-manifest-tree .u2-list'), {clientHeight: 800});
  await flush();
  const add = [...document.querySelectorAll('.u2-manifest-tree a, .u2-manifest-tree .u2-link')]
    .find((a) => a.textContent === '+ Add property');
  assert.ok(add, 'the tree offers "+ Add property" on the extensible table');
  add.click();
  await flush();
}

scoped('Edit properties sends the extension version it was opened at: a concurrent extension is a conflict, not a drop', async () => {
  const calls = [];
  const server = {manifest: MANIFEST, extVersion: 96};
  stub(calls, server);
  const done = domains.editProperties(NAME, 'issue');
  await flush();
  await flush();
  await addProperty();

  // someone else extends the table meanwhile
  server.extVersion = 97;
  server.manifest = copy(MANIFEST);
  server.manifest.tables.issue.columns.theirs = {type: 'string', origin: 'user'};

  button('SAVE').click();
  await flush();
  await flush();
  const sent = calls.filter((c) => c[0] === 'apply' && c[1].extend !== undefined);
  assert.ok(sent.length > 0, 'the save reached apply');
  assert.deepEqual(sent.map((c) => c[1].ifVersion), sent.map(() => '96'),
    'every apply of the save names the extension version the dialog was opened at');
  assert.deepEqual(commits(calls).filter((c) => c[1].extend !== undefined), [], 'nothing committed over the newer extension');
  assert.match(document.querySelector('.u2-manifest-properties').textContent,
    /Someone else changed these properties meanwhile/);

  // the reload rebased onto the current version: saving again names it and keeps their column
  calls.length = 0;
  button('SAVE').click();
  await flush();
  await flush();
  const [commit] = commits(calls);
  assert.ok(commit, 'the second save commits');
  assert.equal(commit[1].ifVersion, '97');
  assert.deepEqual(Object.keys(commit[1].extend.issue.columns).sort(), ['mine', 'property', 'theirs']);
  assert.equal(await done, true);
});
