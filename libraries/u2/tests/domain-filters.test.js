/* `domains.filters` (WO 2-7) over the memory backend: the query box appears once the table's
   schema is known, a query set in code formats into it, a committed text writes the source's
   query (text stays text, a tree stays a tree) and round-trips to `?q=`, and a change while the
   session is dirty goes through the gate — cancel puts the previous filter back. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {Filters} from '../src/core/filter/index.js';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';
import {backend, hierarchyBackend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {registerDomainComponents} = await import('../src/dg/domain/registrations.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

async function issues(options = {}) {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10, ...options});
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
  const filters = domains.filters(src);
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
  const filters = domains.filters(src, {mode: 'builder'});
  await flush();
  const builder = filters.input.value;
  assert.equal(builder.root.dataset.u2, 'filter-builder');
  assert.equal(builder.query.value, 'done = true');
  builder.value.value = Filters.group('and', [Filters.cond('done', '=', false)]);
  await flush();
  assert.equal(typeof src.query.value, 'object');
  assert.equal(Filters.format(src.query.value), 'done = false');
  await flush();
  assert.equal(src.total.value, 2);
  filters.dispose();
  src.dispose();
});

scoped('a change while the session is dirty asks: cancel puts the filter back, discard applies it', async () => {
  const src = await issues({query: 'done = true'});
  const filters = domains.filters(src);
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
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const src = table.source({pageSize: 10});
  const filters = domains.filters(src);
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
  backends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const src = table.source({pageSize: 10});
  await flush();
  const filters = domains.filters(src);
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

  backends.domain = backend();
  const issue = await domains.table('grit.issue');
  const issues = issue.source({pageSize: 10});
  await flush();
  const plain = domains.filters(issues);
  await flush();
  const prop = issues.schema.properties.find((p) => p.name === 'project_id');
  assert.equal(plain.input.value.options.schema.operators(prop, Filters.operators.for(prop))
    .map((o) => o.id).includes('under'), false, 'the ref target is not a hierarchy');
  plain.dispose();
  issues.dispose();
});

scoped('clearing a refused filter releases it: the box, the status and the stale mark all let go', async () => {
  const src = await issues();
  const filters = domains.filters(src);
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
  backends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const src = table.source({pageSize: 10});
  await flush();
  const filters = domains.filters(src);
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
  backends.domain = hierarchyBackend();
  const table = await domains.table('stock.location');
  const src = table.source({pageSize: 10});
  await flush();
  const filters = domains.filters(src);
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
  backends.domain = hierarchyBackend({rows: {location: [
    {id: 'b1', name: 'Building A', kind: 'site'},
    {id: 'b2', name: 'Building B', kind: 'site'},
    {id: 'l1', name: 'Lab 101', parent_id: 'b1', kind: 'room'}]}});
  const table = await domains.table('stock.location');
  const src = table.source({pageSize: 10});
  await flush();
  const filters = domains.filters(src);
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
