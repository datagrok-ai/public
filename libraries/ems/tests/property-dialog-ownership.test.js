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

// EQ-018: who owns a column of a plugin schema is the registry's word (`DomainTable.userColumns`); the manifest
// never carries it. A saved user column read as the plugin's was locked in Edit properties, and the next extension
// — the user's columns as FULL state — left it out, which the server plans as dropping it with its values.

const copy = (x) => JSON.parse(JSON.stringify(x));
const NAME = 'grit';
/** A package schema as manifest() answers it: no ownership anywhere. */
const MANIFEST = {name: NAME, version: '5', incarnation: '2026-09-20T10:00:00.000Z', tables: {issue: {extensible: true,
  columns: {title: {type: 'string'}, mine: {type: 'string'}}}}};
const PLAN = (extra = {}) => ({version: '5', destructive: false,
  creates: {tables: [], columns: {}, uniques: {}, constraints: {}}, drops: {tables: [], columns: []}, typeChanges: [],
  uniqueDrops: [], requiredToggles: [], businessKeyChanges: [], autoNumberChanges: {}, alters: {}, constraintDrops: null,
  violations: [], refusals: [], keyColumnsUnrestricted: [], migration: null, extVersion: 7, ...extra});

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

/** The registry listing (`grok.dapi.domains.schemas`) carries ownership: `DomainTable.origin` / `userColumns`. */
function listing(server) {
  const tables = () => [{name: 'issue', origin: 'package', userColumns: [...server.userColumns]}];
  const list = async () => [{name: NAME, managedBy: 'package', tables: tables()}];
  return {list, include: () => ({list})};
}

/** [conflicts]: how many commits answer `version-conflict` first — each one a peer's column landing meanwhile. */
function stub(calls, server) {
  grok.shell.settings = {enableDomainDatabases: true};
  grok.dapi.permissions = {check: async () => true};
  DG.Property.propertyOptions = Object.fromEntries(['friendlyName', 'description', 'choices', 'min', 'max', 'format']
    .map((name) => [name, {name, type: name === 'choices' ? 'list' : 'string'}]));
  grok.dapi.domains.schemas = listing(server);
  grok.dapi.domains.schema = (name) => ({
    name,
    manifest: async () => copy(server.manifest),
    apply: async (body, options) => {
      calls.push(['apply', copy(body), options?.dryRun === true]);
      if (options?.dryRun === true)
        return PLAN();
      if (server.conflicts > 0) {
        server.conflicts--;
        server.manifest.tables.issue.columns.theirs = {type: 'string'};
        server.userColumns.push('theirs');
        throw Object.assign(new Error('Schema version conflict'),
          {status: 409, code: 'version-conflict', body: {error: 'version-conflict'}});
      }
      return {...PLAN(), applied: true};
    },
  });
}

const commits = (calls) => calls.filter((c) => c[0] === 'apply' && !c[2] && c[1].extend !== undefined);
const button = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);
const node = (name) => [...document.querySelectorAll('.u2-manifest-tree .u2-manifest-node')]
  .find((n) => n.querySelector('.u2-manifest-node-name')?.textContent === name);

async function open() {
  const done = domains.editProperties(NAME, 'issue');
  await flush();
  await flush();
  layout(document.querySelector('.u2-manifest-tree .u2-list'), {clientHeight: 800});
  await flush();
  return {done};
}

async function addProperty() {
  const add = [...document.querySelectorAll('.u2-manifest-tree a, .u2-manifest-tree .u2-link')]
    .find((a) => a.textContent === '+ Add property');
  assert.ok(add, 'the tree offers "+ Add property" on the extensible table');
  add.click();
  await flush();
}

async function save() {
  button('SAVE').click();
  await flush();
  await flush();
}

scoped('EQ-018: a saved user column is the user\'s, and adding the next column sends it back', async () => {
  const calls = [];
  const server = {manifest: copy(MANIFEST), userColumns: ['mine'], conflicts: 0};
  stub(calls, server);
  const {done} = await open();
  const mine = node('mine');
  assert.ok(mine, 'the tree lists the saved user column');
  assert.equal(mine.querySelector('.fa-lock'), null,
    'the user column is not locked as a plugin column');
  assert.match(mine.textContent, /custom/, 'the user column is badged as the user\'s own');
  assert.ok(node('title').querySelector('.fa-lock'), 'a plugin column stays locked');
  await addProperty();
  await save();
  const [commit] = commits(calls);
  assert.ok(commit, 'the save commits an extension');
  assert.deepEqual(Object.keys(commit[1].extend.issue.columns).sort(), ['mine', 'property'],
    'the extension is the user\'s columns as full state: the saved one is not dropped');
  assert.equal(await done, true);
});

scoped('EQ-018: after a version conflict the reload re-reads ownership: a peer\'s new column is kept too', async () => {
  const calls = [];
  const server = {manifest: copy(MANIFEST), userColumns: ['mine'], conflicts: 1};
  stub(calls, server);
  const {done} = await open();
  await addProperty();
  await save();
  assert.match(document.querySelector('.u2-manifest-properties').textContent,
    /Someone else changed these properties meanwhile/);
  calls.length = 0;
  await save();
  const [commit] = commits(calls);
  assert.ok(commit, 'the second save commits');
  assert.deepEqual(Object.keys(commit[1].extend.issue.columns).sort(), ['mine', 'property', 'theirs'],
    'every user column the reloaded registry names travels with the extension');
  assert.equal(await done, true);
});
