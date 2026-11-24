import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import {mock} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {layout} from '@datagrok-libraries/u2/tests/layout.mjs';
import {backends, Scope, Registry, SpecContext, renderSpec, registerAll} from '@datagrok-libraries/u2';
import {backend, hierarchyBackend, TREE} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainTree} = await import('../src/dg/tree.js');
const {registerDomainComponents} = await import('../src/dg/registrations.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** The rows the tree has flattened into the list, in order. */
const labels = (t) => [...t.tree.root.querySelectorAll('.u2-tree-label')].map((el) => el.textContent);

/** The test DOM lays nothing out: the scroller inside the tree is the element that needs a viewport. */
function viewport(t, height = 400) {
  layout(t.tree.root.querySelector('.u2-list'), {clientHeight: height});
}

async function tree(options = {}, be = hierarchyBackend()) {
  emsBackends.domain = be;
  const table = await domains.table('stock.location');
  const t = testOwner().build(() => domains.tree(table, options));
  document.body.append(t.root);
  viewport(t);
  await flush();
  return {table, tree: t};
}

/** Opens the branch the label names by clicking its twistie. */
async function open(t, label) {
  const row = [...t.tree.root.querySelectorAll('.u2-tree-row')]
    .find((el) => el.querySelector('.u2-tree-label')?.textContent === label);
  fire(row.querySelector('.u2-tree-twistie'), 'click', {bubbles: true});
  await flush();
  await flush();
}

scoped('the roots are the rows with no parent; a branch reads its children when it is opened', async () => {
  const {tree: t} = await tree();
  assert.equal(t.root.dataset.u2, 'domain-tree');
  assert.equal(t.parentColumn, 'parent_id');
  assert.deepEqual(labels(t), ['Other site', 'Site'], 'roots only, by the name column');
  assert.equal(t.error.value, null);

  await open(t, 'Site');
  assert.deepEqual(labels(t), ['Other site', 'Site', 'Room'], 'one level, not the whole subtree');
  await open(t, 'Room');
  assert.deepEqual(labels(t), ['Other site', 'Site', 'Room', 'Shelf']);
  t.dispose();
});

scoped('selection is a signal over the row behind the node', async () => {
  const {tree: t} = await tree();
  assert.equal(t.selected.value, null);
  fire(t.tree.root.querySelectorAll('.u2-tree-row')[1], 'click', {bubbles: true});
  await flush();
  assert.equal(t.selected.value.id, 'l1');
  assert.equal(t.selected.value.kind, 'site', 'the whole row, not the label');
  t.tree.clearSelection();
  await flush();
  assert.equal(t.selected.value, null);
  t.dispose();
});

scoped('expandTo opens every ancestor of the row and selects it', async () => {
  const {tree: t} = await tree({expandTo: 'l4'});
  assert.deepEqual(labels(t), ['Other site', 'Site', 'Room', 'Shelf', 'Box'],
    'the whole path is open — the ancestors answer excludes the row, so the tree adds it');
  assert.equal(t.selected.value.id, 'l4');
  assert.equal(t.error.value, null);
  t.dispose();
});

scoped('a non-hierarchy table is refused by name; no ancestors, no path', async () => {
  emsBackends.domain = backend();
  const issue = await domains.table('grit.issue');
  assert.throws(() => domains.tree(issue), (e) => e.code === 'filter' &&
    /grit\.issue is not a hierarchy table/.test(e.message));

  // a backend that declares no ancestors keeps the tree, and says why the path did not open
  const be = hierarchyBackend();
  const inner = await be.table('stock.location');
  const stripped = {...inner, query: (spec) => inner.query(spec), access: () => inner.access(),
    count: (f, s, d) => inner.count(f, s, d), transaction: (ops) => inner.transaction(ops),
    frame: (spec) => inner.frame(spec), ancestors: undefined,
    support: {...inner.support, ancestors: false}};
  const {tree: t} = await tree({expandTo: 'l4'}, {table: () => Promise.resolve(stripped)});
  assert.deepEqual(labels(t), ['Other site', 'Site'], 'the roots are still there');
  assert.match(t.error.value, /cannot answer a row's ancestors/);
  t.dispose();
});

scoped('a path truncated at an ancestor out of sight opens nothing, and says so', async () => {
  // the memory backend answers [] for a row whose chain does not resolve — the same shape the
  // server's path has when it stops at an ancestor the caller cannot see
  const orphan = hierarchyBackend({rows: {location: [...TREE.location,
    {id: 'l5', name: 'Orphan', parent_id: 'hidden'}]}});
  const {tree: t} = await tree({expandTo: 'l5'}, orphan);
  assert.equal(t.selected.value, null, 'nothing was reached, so nothing is selected');
  assert.deepEqual(labels(t), ['Other site', 'Site'], 'and the roots are all that opened');
  assert.match(t.error.value, /l5 is not reachable from here — an ancestor is not visible to you/);
  t.dispose();
});

scoped('an address that is not a hierarchy lands in error instead of throwing', async () => {
  emsBackends.domain = backend();
  const t = new DomainTree('grit.issue');
  document.body.append(t.root);
  await flush();
  await flush();
  assert.match(t.error.value, /grit\.issue is not a hierarchy table/);
  assert.equal(t.table, undefined, 'and the handle is never adopted');
  t.dispose();
});

scoped('spec: u2-domain-tree takes the table address', async () => {
  emsBackends.domain = hierarchyBackend();
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const meta = reg.get('u2-domain-tree');
  assert.equal(meta.category, 'Display');
  assert.deepEqual(meta.props.map((p) => p.name).slice(0, 3), ['table', 'expandTo', 'pageSize']);
  assert.match(meta.usage, /under "<id>"/);
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-domain-tree', name: 'places',
    props: {table: 'stock.location'}}}, new SpecContext({data: {}}), reg);
  const t = instance.node('places');
  assert.equal(t instanceof DomainTree, true);
  viewport(t);
  await flush();
  await flush();
  assert.deepEqual(labels(t), ['Other site', 'Site']);
  assert.equal(t.table.address, 'stock.location');
  // `selected` is a declared prop, so a bind path can walk to it — without that it resolves to null
  const bound = instance.resolveBinding('$.places.selected');
  assert.equal(bound.signal.value, null);
  assert.equal(bound.writable, false, 'read-only');
  fire(t.tree.root.querySelectorAll('.u2-tree-row')[1], 'click', {bubbles: true});
  await flush();
  assert.equal(bound.signal.value.id, 'l1', 'and it follows the selection');
  instance.dispose();
});

scoped('the selected node can be given back: clicking it again, Escape, and the all row', async () => {
  const {tree: t} = await tree({allNode: 'All locations'});
  assert.deepEqual(labels(t), ['All locations', 'Other site', 'Site'], 'the all row leads the roots');

  const rowAt = (i) => t.tree.root.querySelectorAll('.u2-tree-row')[i];
  fire(rowAt(2), 'click', {bubbles: true});
  await flush();
  assert.equal(t.selected.value.id, 'l1');

  // a tree is a navigator: "nothing selected" — the whole table — is one of its states
  fire(rowAt(2), 'click', {bubbles: true});
  await flush();
  assert.equal(t.selected.value, null, 'clicking the selected node again gives it back');

  fire(rowAt(2), 'click', {bubbles: true});
  await flush();
  assert.equal(t.selected.value.id, 'l1');
  fire(t.tree.root, 'keydown', {key: 'Escape', bubbles: true});
  await flush();
  assert.equal(t.selected.value, null, 'and so does Escape');

  fire(rowAt(2), 'click', {bubbles: true});
  await flush();
  fire(rowAt(0), 'click', {bubbles: true});
  await flush();
  assert.equal(t.selected.value, null, 'the all row stands for the whole table');
  assert.equal(t.tree.selectedNode.value.id, '~all', 'and is itself selected, so the state is visible');
  t.dispose();
});

scoped('a branch that answers with no children keeps no twistie', async () => {
  const {tree: t} = await tree();
  const twistie = (label) => [...t.tree.root.querySelectorAll('.u2-tree-row')]
    .find((el) => el.querySelector('.u2-tree-label')?.textContent === label)
    .querySelector('.u2-tree-twistie');
  assert.equal(twistie('Other site').className.includes('u2-tree-twistie-hidden'), false,
    'until it has answered, every row may have children');
  await open(t, 'Other site');
  assert.equal(twistie('Other site').className.includes('u2-tree-twistie-hidden'), true,
    'once it has answered with none, it is a leaf');
  assert.equal(twistie('Site').className.includes('u2-tree-twistie-hidden'), false);
  t.dispose();
});

/** A tree that follows a live source of the same table; `tick` is one live poll. */
async function following(options = {}) {
  mock.timers.enable({apis: ['setInterval']});
  const be = emsBackends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const source = testOwner().build(() => table.source({live: true, liveMs: 1000, pageSize: 10}));
  await settle();
  const t = testOwner().build(() => domains.tree(table, {follow: source, ...options}));
  document.body.append(t.root);
  viewport(t);
  await settle();
  const peer = be.tableSync('stock.location');
  const tick = async () => {
    mock.timers.tick(1000);
    await settle();
  };
  const done = () => {
    t.dispose();
    source.dispose();
    mock.timers.reset();
  };
  return {source, tree: t, peer, tick, done};
}

async function settle() {
  for (let i = 0; i < 12; i++)
    await flush();
}

scoped('the tree follows a peer\'s rename, insert and delete of the selected node on the live probe', async () => {
  const f = await following({expandTo: 'l4'});
  try {
    await f.peer.transaction([{op: 'update', table: 'location', id: 'l4', values: {name: 'Box 2'}},
      {op: 'insert', table: 'location', values: {id: 'l0', name: 'Annex', parent_id: 'l1'}}]);
    await f.tick();
    assert.deepEqual(labels(f.tree), ['Other site', 'Site', 'Annex', 'Room', 'Shelf', 'Box 2']);
    assert.equal(f.tree.selected.value.id, 'l4', 'the renamed row stays selected, though a row landed above it');
    await f.peer.transaction([{op: 'delete', table: 'location', id: 'l4'}]);
    await f.tick();
    assert.equal(labels(f.tree).includes('Box 2'), false);
    assert.equal(f.tree.selected.value?.id, 'l3', 'the deleted selection gives way to its parent');
  } finally {
    f.done();
  }
});

scoped('a peer\'s delete of the selected node over pending changes keeps it selected and says so', async () => {
  const f = await following({expandTo: 'l4'});
  try {
    f.source.newRow({name: 'Draft', parent_id: 'l4'});
    await f.peer.transaction([{op: 'delete', table: 'location', id: 'l4'}]);
    await f.tick();
    assert.equal(f.tree.selected.value?.id, 'l4', 'what is bound to the selection is not rebound over pending changes');
    assert.match(f.tree.notice.value, /^Box was deleted/);
    f.source.discard();
    await settle();
    assert.equal(f.tree.selected.value?.id, 'l3', 'once the changes are gone, the parent takes over');
  } finally {
    f.done();
  }
});

scoped('a node the user clicks while a followed tree is refreshing stays selected when the refresh ends', async () => {
  const f = await following({expandTo: 'l4'});
  try {
    const memory = f.peer;
    const real = memory.query;
    let release;
    const gate = new Promise((resolve) => release = resolve);
    memory.query = async (spec) => {
      if (spec.filter?.property === 'parent_id' || spec.filter?.property === 'id')
        await gate;
      return real.call(memory, spec);
    };
    await f.peer.transaction([{op: 'update', table: 'location', id: 'l2', values: {name: 'Room 2'}}]);
    await f.tick();
    const other = [...f.tree.root.querySelectorAll('[role="treeitem"]')].find((n) => n.textContent === 'Other site');
    fire(other, 'click', {bubbles: true});
    await flush();
    release();
    await settle();
    memory.query = real;
    assert.ok(labels(f.tree).includes('Room 2'), 'the refresh ran');
    assert.equal(f.tree.selected.value?.id, 'l9', 'the user\'s pick wins over the row the refresh re-found');
  } finally {
    f.done();
  }
});
