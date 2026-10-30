import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Filters, Registry, registerAll} from '@datagrok-libraries/u2';
import {backend, hierarchyBackend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {SavedSearches} = await import('../src/dg/saved-searches.js');
const {userSettings} = await import('datagrok-api/grok');
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

async function issues(options = {}) {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = testOwner().build(() => table.source({pageSize: 10, ...options}));
  await flush();
  return src;
}

const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);

/** What Enter does in the box: commits the text into the tree. */
function commit(filters, text) {
  filters.input.value.text.value = text;
  return filters.input.value.commit();
}

scoped('a query set in code formats into the box; a committed text writes the query and round-trips to ?q=', async () => {
  const src = await issues();
  const filters = testOwner().build(() => domains.filters(src));
  assert.equal(filters.root.dataset.u2, 'domain-filters');
  await flush();
  const input = filters.input.value;
  assert.equal(input.root.dataset.u2, 'filter-query-input', 'the query box, once the schema is known');
  assert.equal(input.text.value, '');
  src.query.value = 'done=true';
  assert.equal(input.text.value, 'done = true', 'formatted');
  assert.equal(commit(filters, 'title like "a"'), true);
  // the gate answers at once over a clean session, a microtask later
  await flush();
  assert.equal(src.query.value, 'title like "a"', 'text stays text');
  await flush();
  assert.equal(src.total.value, 2, 'Aspirin and Naproxen');
  const path = Filters.queryPath('/apps/T/Issues', src.query.value);
  assert.equal(path, '/apps/T/Issues?q=title%20like%20%22a%22');
  assert.equal(new URLSearchParams(path.slice(path.indexOf('?') + 1)).get('q'), 'title like "a"');
  assert.equal(commit(filters, 'nope = 1'), false, 'an unknown column is a problem, not a query');
  await flush();
  assert.equal(src.query.value, 'title like "a"');
  filters.dispose();
  src.dispose();
});

scoped('a tree query stays a tree; the builder mode shows rows', async () => {
  const src = await issues({query: Filters.group('and', [Filters.cond('done', '=', true)])});
  const filters = testOwner().build(() => domains.filters(src, {mode: 'builder'}));
  await flush();
  const builder = filters.input.value;
  assert.equal(builder.root.dataset.u2, 'filter-builder');
  assert.equal(builder.query.value, 'done = true');
  builder.value = Filters.group('and', [Filters.cond('done', '=', false)]);
  await flush();
  assert.equal(typeof src.query.value, 'object');
  assert.equal(Filters.format(src.query.value), 'done = false');
  await flush();
  assert.equal(src.total.value, 2);
  filters.dispose();
  src.dispose();
});

scoped('a filter from outside opens the builder in the mode its shape needs; the builder\'s own edits keep it', async () => {
  const src = await issues({query: ''});
  const filters = testOwner().build(() => domains.filters(src, {mode: 'builder'}));
  await flush();
  const builder = filters.input.value;
  src.query.value = 'title like "a" or (done = true)';
  await flush();
  assert.equal(builder.mode.value, 'advanced', 'a sub-group needs Advanced');
  src.query.value = 'done = true';
  await flush();
  assert.equal(builder.mode.value, 'simple', 'a flat filter opens Simple');
  builder.setMode('advanced');
  builder.value = Filters.group('and', [Filters.cond('done', '=', false)]);
  await flush();
  assert.equal(src.query.value, 'done = false');
  assert.equal(builder.mode.value, 'advanced', 'an edit in the builder never flips its mode');
  filters.dispose();
  src.dispose();
});

scoped('the builder follows a filter set from outside even when it is the one the builder last sent', async () => {
  const src = await issues({query: ''});
  const filters = testOwner().build(() => domains.filters(src, {mode: 'builder'}));
  await flush();
  const builder = filters.input.value;
  builder.value = Filters.group('and', [Filters.cond('done', '=', true)]);
  await flush();
  src.query.value = 'done = false';
  await flush();
  src.query.value = 'done = true';
  await flush();
  assert.equal(Filters.format(builder.value), 'done = true');
  filters.dispose();
  src.dispose();
});

scoped('the builder keeps a half-written row and sends only the complete conditions', async () => {
  const src = await issues({query: ''});
  const filters = testOwner().build(() => domains.filters(src, {mode: 'builder'}));
  await flush();
  const half = {...Filters.cond('title', 'like', null), property: ''};
  filters.input.value.value = Filters.group('and', [Filters.cond('done', '=', true), half]);
  await flush();
  assert.deepEqual([src.query.value, filters.input.value.value.nodes.length], ['done = true', 2]);
  filters.dispose();
  src.dispose();
});

scoped('a change while the session is dirty asks: cancel puts the filter back, discard applies it', async () => {
  const src = await issues({query: 'done = true'});
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  src.rows.byKey('i1').title = 'Edited';
  assert.equal(src.isDirty.value, true);
  commit(filters, 'done = false');
  await flush();
  assert.notEqual(document.body.querySelector('.u2-dialog'), null, 'the gate asks');
  assert.equal(src.query.value, 'done = true', 'nothing written yet');
  fire(buttonNamed('CANCEL'), 'click');
  await flush();
  assert.equal(filters.input.value.text.value, 'done = true', 'the box shows the filter in force again');
  assert.equal(src.query.value, 'done = true');
  assert.equal(src.isDirty.value, true, 'the edit is kept');
  commit(filters, 'done = false');
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(src.query.value, 'done = false');
  assert.equal(src.isDirty.value, false);
  await flush();
  assert.equal(src.total.value, 2, 're-queried');
  filters.dispose();
  src.dispose();
});

scoped('spec: u2-domain-filters is registered with usage and a mode', () => {
  const reg = new Registry();
  registerAll(reg);
  registerDomainComponents(reg);
  const meta = reg.get('u2-domain-filters');
  assert.equal(meta.usage.length > 0, true);
  assert.deepEqual(meta.props.find((p) => p.name === 'mode').choices, ['query', 'builder']);
});

scoped('a schema the platform refuses leaves the box out, says why through the source, and is tried again', async () => {
  emsBackends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  const filters = testOwner().build(() => domains.filters(src));
  const schema = filters._schema.bind(filters);
  filters._schema = () => Promise.reject(new Error('no filter schema'));
  await flush();
  assert.equal(filters.input.value, null, 'nothing to type into');
  assert.match(String(src.error.value?.message ?? ''), /no filter schema/, 'and the source carries why');
  filters._schema = schema;
  src.query.value = 'done = true';
  await flush();
  await flush();
  assert.notEqual(filters.input.value, null, 'the latch released with the failure: the next state builds the box');
  filters.dispose();
  src.dispose();
});

scoped('"is under" is offered only where the hierarchy it walks exists', async () => {
  emsBackends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  const offers = (name) => filters.input.value.options.schema
    .operators(src.schema.properties.find((p) => p.name === name), Filters.operators.for({name, kind: 'string'}))
    .map((o) => o.id);
  assert.equal(offers('parent_id').includes('under'), true, 'a ref column into a hierarchy table');
  assert.equal(offers('id').includes('under'), true, 'and `id` on the hierarchy itself');
  assert.equal(offers('name').includes('under'), false, 'a plain string column is not a subtree');
  assert.equal(offers('name').includes('like'), true, 'and keeps everything else');
  filters.dispose();
  src.dispose();

  emsBackends.domain = backend();
  const issue = await domains.table('grit.issue');
  const issues = testOwner().build(() => issue.source({pageSize: 10}));
  await flush();
  const plain = testOwner().build(() => domains.filters(issues));
  await flush();
  const prop = issues.schema.properties.find((p) => p.name === 'project_id');
  assert.equal(plain.input.value.options.schema.operators(prop, Filters.operators.for(prop))
    .map((o) => o.id).includes('under'), false, 'the ref target is not a hierarchy');
  plain.dispose();
  issues.dispose();
});

scoped('clearing a refused filter releases it: the box, the status and the stale mark all let go', async () => {
  const src = await issues();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  const input = filters.input.value;
  commit(filters, 'zzz = "1"');
  await flush();
  assert.notEqual(filters.problem.value, null, 'an unknown column is refused');
  assert.equal(input.problems.value.length > 0, true);

  commit(filters, '');
  await flush();
  assert.equal(filters.problem.value, null, 'the box is about what is in it, not about what was');
  assert.deepEqual(input.problems.value, []);
  filters.dispose();
  src.dispose();
});

scoped('an `under` value slot offers the rows of the hierarchy, not the values the column holds', async () => {
  emsBackends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  const schema = filters.input.value.options.schema;
  const prop = src.schema.properties.find((p) => p.name === 'parent_id');
  const signal = new AbortController().signal;

  const under = await schema.values(prop, '', signal, {operator: 'under'});
  assert.deepEqual(under.map((i) => i.label).sort(),
    ['Box', 'Other site', 'Room', 'Shelf', 'Site'], 'every location, including the ones holding nothing');
  assert.deepEqual(under[0].value, {type: 'stock.location', id: under[0].value.id, name: under[0].label},
    'the caption is shown and the id is what lands in the query');
  assert.deepEqual((await schema.values(prop, 'Sh', signal, {operator: 'under'})).map((i) => i.label), ['Shelf'],
    'typed text narrows them');

  // `id under` on the hierarchy itself picks from the same table
  const id = src.schema.properties.find((p) => p.name === 'id');
  assert.equal((await schema.values(id, '', signal, {operator: 'under'})).length, 5);
  filters.dispose();
  src.dispose();
});

scoped('an unquoted `under` value says what to do about it, in the status line as well as the box', async () => {
  emsBackends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  // this one fails in the GRAMMAR, before the schema check that knows what `under` takes —
  // "Expected a value" named neither the operator nor the way out of it
  commit(filters, 'parent_id under Building A');
  await flush();
  assert.match(filters.problem.value, /under takes a location — pick one from the list/);
  assert.equal(filters.problem.value.includes('Expected'), false);

  filters.dispose();
  src.dispose();
});

scoped('over a real hierarchy: a two-word location narrows to it, and Enter does not swap it', async () => {
  emsBackends.domain = hierarchyBackend({rows: {location: [
    {id: 'b1', name: 'Building A', kind: 'site'},
    {id: 'b2', name: 'Building B', kind: 'site'},
    {id: 'l1', name: 'Lab 101', parent_id: 'b1', kind: 'room'}]}});
  const table = await domains.table('stock.location');
  const src = testOwner().build(() => table.source({pageSize: 10}));
  await flush();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  const input = filters.input.value;

  const el = input.root.querySelector('input');
  document.body.append(input.root);
  const suggest = async (text) => {
    el.focus();
    el.value = text;
    fire(el, 'input');
    await flush();
    await flush();
    return [...document.body.querySelectorAll('.u2-fq-option')].map((r) => r.textContent);
  };
  assert.deepEqual(await suggest('parent_id under Building A'), ['Building A'],
    'the whole name is the search — not every location with "A" in it, "Lab 101" first');
  assert.deepEqual(await suggest('parent_id under Building'), ['Building A', 'Building B']);
  filters.dispose();
  src.dispose();
});

scoped('a `basic` backend is offered no under, regex, !like, fuzzy or between, no datetime !=, no bool null test and no bigint ordering', async () => {
  const offersOf = (filters, src) => (name) => {
    const prop = src.schema.properties.find((p) => p.name === name);
    return filters.input.value.options.schema.operators(prop, Filters.operators.for(prop)).map((o) => o.id);
  };
  const memory = backend();
  const issue = memory.tableSync('grit.issue');
  issue.support = {...issue.support, filters: 'basic'};
  emsBackends.domain = memory;
  const issues = await domains.table('grit.issue');
  const src = testOwner().build(() => issues.source({pageSize: 10}));
  await flush();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  const offers = offersOf(filters, src);
  const title = offers('title');
  for (const id of ['matches', '!matches', '!like', 'fuzzy'])
    assert.equal(title.includes(id), false, `${id} is not offered on a string`);
  for (const id of ['=', '!=', 'like', 'starts', 'is null'])
    assert.equal(title.includes(id), true, `${id} stays`);
  assert.equal(offers('due').includes('!='), false, 'no datetime !=');
  assert.equal(offers('due').includes('<'), true, 'the datetime comparisons stay');
  for (const name of ['due', 'number', 'weight'])
    assert.equal(offers(name).includes('between'), false, `no between on ${name}: its payload is a nested group`);
  assert.equal(offers('number').includes('>='), true, 'a flat >= stays, and pairs of them still merge');
  for (const id of ['is null', 'is not null'])
    assert.equal(offers('done').includes(id), false, `no bool ${id}`);
  assert.equal(offers('done').includes('='), true);
  assert.equal(offers('number').includes('is null'), true, 'an int keeps its null test');
  const big = {name: 'molregno', propertyType: 'bigint'};
  const bigOffers = filters.input.value.options.schema.operators(big, Filters.operators.for(big)).map((o) => o.id);
  for (const id of ['>', '>=', '<', '<='])
    assert.equal(bigOffers.includes(id), false, `no bigint ${id}: the connector compares a bigint by equality only`);
  for (const id of ['=', '!=', 'in'])
    assert.equal(bigOffers.includes(id), true, `bigint ${id} stays`);
  filters.dispose();
  src.dispose();

  // `under` goes even where the hierarchy it would walk exists
  const stock = hierarchyBackend();
  const location = stock.tableSync('stock.location');
  location.support = {...location.support, filters: 'basic'};
  emsBackends.domain = stock;
  const locations = await domains.table('stock.location');
  const tree = testOwner().build(() => locations.source({pageSize: 10}));
  await flush();
  const treeFilters = testOwner().build(() => domains.filters(tree));
  await flush();
  assert.equal(offersOf(treeFilters, tree)('parent_id').includes('under'), false, 'no under on a basic backend');
  treeFilters.dispose();
  tree.dispose();
});

scoped('saved and recent searches: the list names the saved search it holds, history keeps 20, newest first', async () => {
  userSettings.values = {};
  const high = {q: 'priority = "high"', search: 'Nap', sort: '!number'};
  SavedSearches.save('grit.issue', 'High', high);
  const src = await issues();
  const filters = testOwner().build(() => domains.filters(src));
  await flush();
  assert.equal(filters.root.querySelectorAll('.u2-domain-filters-tools button').length, 2);
  SavedSearches.apply(src, high);
  await flush();
  assert.equal(SavedSearches.title(src), 'High');
  src.search.value = 'Na';
  await flush();
  assert.equal(SavedSearches.title(src), 'High (edited)');
  assert.deepEqual(SavedSearches.recent('grit.issue').map((r) => r.search), ['Na', 'Nap']);
  for (let i = 0; i < 25; i++)
    SavedSearches.record('grit.issue', {...high, search: `t${i}`});
  assert.equal(SavedSearches.recent('grit.issue').length, 20);
  SavedSearches.remove('grit.issue', 'High');
  assert.equal(SavedSearches.title(src), null);
  filters.dispose();
  src.dispose();
});
