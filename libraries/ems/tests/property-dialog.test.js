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
  // ownership comes from the registry listing (`DomainTable.userColumns`); the `origin: 'user'` stamps below only
  // say which columns it lists — the manifest itself never carries them
  const tables = () => Object.entries(server.manifest.tables).map(([name, t]) => ({name, origin: 'package',
    userColumns: Object.keys(t.columns).filter((c) => t.columns[c].origin === 'user')}));
  const unstamped = () => {
    const m = copy(server.manifest);
    for (const t of Object.values(m.tables))
      for (const c of Object.values(t.columns))
        delete c.origin;
    return m;
  };
  grok.dapi.domains.schemas = {
    list: async () => [{name: NAME, managedBy: 'package'}],
    include: () => ({list: async () => {
      calls.push(['ownership']);
      return [{name: NAME, managedBy: 'package', tables: tables()}];
    }}),
  };
  grok.dapi.domains.schema = (name) => ({
    name,
    manifest: async () => {
      calls.push(['manifest']);
      return unstamped();
    },
    apply: async (body, options) => {
      calls.push(['apply', copy(body), options?.dryRun === true]);
      if (body.ifVersion !== undefined && body.ifVersion !== String(server.extVersion)) {
        throw Object.assign(new Error(`Schema version conflict: expected ${body.ifVersion}, current ${server.extVersion}`),
          {status: 409, code: 'version-conflict', body: {error: 'version-conflict'}});
      }
      if (options?.dryRun === true)
        return PLAN({extVersion: server.extVersion, ...(body.extend !== undefined ? server.plan ?? {} : {})});
      return {...PLAN(), applied: true, extVersion: ++server.extVersion};
    },
  });
}

const commits = (calls) => calls.filter((c) => c[0] === 'apply' && !c[2] && Object.keys(c[1]).length > 0);
const save = () => byName(byKind(document, 'dialog')[0], 'SAVE').click();

async function addProperty() {
  const [tree] = byKind(document, 'manifest-tree');
  layout(byKind(tree, 'list')[0], {clientHeight: 800});
  await flush();
  byName(tree, '+ Add property').click();
  await flush();
}

scoped('Edit properties sends the extension version it was opened at: a concurrent extension is a conflict, not a drop', async () => {
  const calls = [];
  const server = {manifest: MANIFEST, extVersion: 96};
  stub(calls, server);
  const done = domains.editProperties(NAME, 'issue');
  await flush();
  await flush();
  const opening = copy(calls);
  await addProperty();

  // someone else extends the table meanwhile
  server.extVersion = 97;
  server.manifest = copy(MANIFEST);
  server.manifest.tables.issue.columns.theirs = {type: 'string', origin: 'user'};

  save();
  await flush();
  await flush();
  const sent = calls.filter((c) => c[0] === 'apply' && c[1].extend !== undefined);
  assert.ok(sent.length > 0, 'the save reached apply');
  assert.deepEqual(sent.map((c) => c[1].ifVersion), sent.map(() => '96'),
    'every apply of the save names the extension version the dialog was opened at');
  // the first manifest read is the extensible-tables check; the dialog's baseline reads come after it
  assert.deepEqual(opening, [['manifest'], ['apply', {}, true], ['manifest'], ['ownership']],
    'the version is read before the manifest and the ownership: a change between the reads is a conflict, ' +
    'never an overwrite or a drop');
  assert.deepEqual(commits(calls).filter((c) => c[1].extend !== undefined), [], 'nothing committed over the newer extension');
  assert.match(byKind(document, 'status-line')[0].textContent,
    /Someone else changed these properties meanwhile/);

  // the reload rebased onto the current version: saving again names it and keeps their column
  calls.length = 0;
  save();
  await flush();
  await flush();
  const [commit] = commits(calls);
  assert.ok(commit, 'the second save commits');
  assert.equal(commit[1].ifVersion, '97');
  assert.deepEqual(Object.keys(commit[1].extend.issue.columns).sort(), ['mine', 'property', 'theirs']);
  assert.equal(await done, true);
});

scoped('Edit properties confirms the drops the server plans, not the client diff: SAVE waits for the checkbox', async () => {
  const calls = [];
  // the client only adds a column; the server's plan, disagreeing, also drops one that holds values
  const server = {manifest: MANIFEST, extVersion: 96, plan: {destructive: true,
    drops: {tables: [], columns: [{table: 'issue', column: 'mine', nonNullValues: 3}]}}};
  stub(calls, server);
  const done = domains.editProperties(NAME, 'issue');
  await flush();
  await flush();
  await addProperty();
  save();
  await flush();
  await flush();

  const confirm = byKind(document, 'dialog').find((d) => state(d).name === 'properties-confirm');
  assert.ok(confirm, 'a destructive plan asks before saving');
  assert.match(confirm.textContent, /issue\.mine dropped.*3 values deleted/,
    'the confirmation lists the plan\'s drop with its value count');
  const box = byName(confirm, 'confirmDestructive');
  assert.equal(state(box).visible, true, 'the deletion checkbox is shown');
  assert.equal(state(box).label, 'Delete 3 stored values');
  assert.equal(byName(confirm, 'SAVE').disabled, true, 'SAVE waits for the checkbox');
  byName(confirm, 'SAVE').click();
  await flush();
  assert.deepEqual(commits(calls), [], 'nothing committed before the deletion is confirmed');

  user.click(byName(box, 'Delete 3 stored values'));
  await flush();
  assert.equal(byName(confirm, 'SAVE').disabled, false, 'ticking the checkbox enables SAVE');
  byName(confirm, 'SAVE').click();
  await flush();
  await flush();
  const [commit] = commits(calls);
  assert.ok(commit, 'the confirmed save commits');
  assert.equal(commit[1].confirmDestructive, true);
  assert.equal(await done, true);
});
