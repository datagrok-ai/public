/* dapiSource and dapiPager against a fake HttpDataSource: the factory must run per request (the
   real HttpDataSource mutates itself), the default filter is the sanitized free-text query, and
   options pass through. Structural — no datagrok-api involved. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {dapiSource, dapiPager, sanitizeFilterValue} from '../src/dg/dapi-source.js';

function fake() {
  const calls = [];
  let constructed = 0;
  const factory = () => {
    constructed++;
    return {list: (options) => { calls.push(options); return Promise.resolve(['a', 'b']); }};
  };
  return {factory, calls, count: () => constructed};
}

const abort = new AbortController().signal;

test('factory runs per query, so data-source state never accumulates', async () => {
  const {factory, count} = fake();
  const fetch = dapiSource(factory);
  await fetch('one', abort);
  await fetch('two', abort);
  assert.equal(count(), 2);
});

test('default filter: sanitized free-text query; empty query lists unfiltered', async () => {
  const {factory, calls} = fake();
  const fetch = dapiSource(factory);
  await fetch('  ada "x" \\ ', abort);
  assert.deepEqual(calls[0], {pageSize: 20, filter: 'ada x'});
  await fetch('', abort);
  assert.deepEqual(calls[1], {pageSize: 20});
  await fetch('   ', abort);
  assert.deepEqual(calls[2], {pageSize: 20}, 'whitespace-only query is empty');
});

test('custom filter receives the sanitized query and can decline to filter', async () => {
  const {factory, calls} = fake();
  const seen = [];
  const fetch = dapiSource(factory, {
    pageSize: 5,
    filter: (q) => { seen.push(q); return q ? `login like "${q}"` : undefined; },
  });
  await fetch('a"b', abort);
  assert.deepEqual(seen, ['ab']);
  assert.deepEqual(calls[0], {pageSize: 5, filter: 'login like "ab"'});
  await fetch('', abort);
  assert.deepEqual(calls[1], {pageSize: 5});
});

test('sanitizeFilterValue strips quotes and backslashes and trims', () => {
  assert.equal(sanitizeFilterValue(' a"b\\c '), 'abc');
  assert.equal(sanitizeFilterValue('plain'), 'plain');
});

const tick = () => new Promise((resolve) => setTimeout(resolve, 0));

/** Every constructed source records the calls made on it, so "one source per request" and
 * "count is filtered like the pages" are both readable off `sources`. */
function fakePaged() {
  const sources = [];
  const factory = () => {
    const ops = [];
    const source = {
      list: (options) => { ops.push({call: 'list', options}); return Promise.resolve([]); },
      count: () => { ops.push({call: 'count'}); return Promise.resolve(42); },
      filter: (query) => { ops.push({call: 'filter', query}); return source; },
      order: (field, desc) => { ops.push({call: 'order', field, desc}); return source; },
    };
    sources.push(ops);
    return source;
  };
  return {factory, sources};
}

test('dapiPager: a fresh source per request, ordered and filtered for the count too', async () => {
  const {factory, sources} = fakePaged();
  const pager = dapiPager(factory, {pageSize: 3, order: 'createdOn', desc: true, filter: 'isResolved = false'});
  pager.reset();
  await tick();
  assert.equal(sources.length, 2, 'the count and the first page each build their own source');
  assert.deepEqual(sources[0], [
    {call: 'order', field: 'createdOn', desc: true},
    {call: 'filter', query: 'isResolved = false'},
    {call: 'count'},
  ]);
  assert.deepEqual(sources[1], [
    {call: 'order', field: 'createdOn', desc: true},
    {call: 'list', options: {pageSize: 3, pageNumber: 0, filter: 'isResolved = false'}},
  ]);
  assert.equal(pager.total.value, 42);
  pager.dispose();
});

test('dapiPager: the filter thunk is re-read on every reset, and may decline to filter', async () => {
  const {factory, sources} = fakePaged();
  let query = 'ada';
  const pager = dapiPager(factory, {filter: () => query ? `description like "${query}"` : undefined});
  pager.reset();
  await tick();
  assert.deepEqual(sources[1][0].options, {pageSize: 20, pageNumber: 0, filter: 'description like "ada"'});

  query = 'bruno';
  pager.reset();
  await tick();
  assert.deepEqual(sources[3][0].options, {pageSize: 20, pageNumber: 0, filter: 'description like "bruno"'});
  assert.deepEqual(sources[2], [{call: 'filter', query: 'description like "bruno"'}, {call: 'count'}]);

  query = '';
  pager.reset();
  await tick();
  assert.deepEqual(sources[5][0].options, {pageSize: 20, pageNumber: 0}, 'an empty query lists unfiltered');
  assert.deepEqual(sources[4], [{call: 'count'}], 'and counts unfiltered');
  pager.dispose();
});

test('dapiPager: no order option means the source is never reordered', async () => {
  const {factory, sources} = fakePaged();
  const pager = dapiPager(factory);
  pager.reset();
  await tick();
  assert.deepEqual(sources[1], [{call: 'list', options: {pageSize: 20, pageNumber: 0}}]);
  pager.dispose();
});
