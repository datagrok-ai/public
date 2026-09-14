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
import {backend} from './domain-fixtures.mjs';

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
  assert.equal(src.query.value, 'title like "a"', 'text stays text');
  await flush();
  assert.equal(src.total.value, 2, 'Aspirin and Naproxen');
  const path = Filters.queryPath('/apps/T/Issues', src.query.value);
  assert.equal(path, '/apps/T/Issues?q=title%20like%20%22a%22');
  assert.equal(new URLSearchParams(path.slice(path.indexOf('?') + 1)).get('q'), 'title like "a"');
  assert.equal(commit(filters, 'nope = 1'), false, 'an unknown column is a problem, not a query');
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
