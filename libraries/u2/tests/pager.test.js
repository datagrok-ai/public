/* AsyncPager against hand-resolved page promises: pages accumulate, a short page ends the
   collection, reset() invalidates whatever is in flight and re-reads the total, and an error
   keeps the pages already loaded. Platform-free — no DOM, no datagrok-api. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {AsyncPager} from '../src/core/pager.js';
import {Scope} from '../src/core/scope.js';

const tick = () => new Promise((resolve) => setTimeout(resolve, 0));

function fake(options = {}) {
  const pages = [];
  const counts = [];
  const pager = new AsyncPager({
    pageSize: 2,
    fetchPage: (page, abort) => new Promise((resolve, reject) => pages.push({page, abort, resolve, reject})),
    count: (abort) => new Promise((resolve, reject) => counts.push({abort, resolve, reject})),
    ...options,
  });
  return {pager, pages, counts};
}

test('pages accumulate in order, one array signal', async () => {
  const {pager, pages} = fake();
  pager.reset();
  assert.equal(pager.state.value, 'loading');
  assert.equal(pages.length, 1);
  assert.equal(pages[0].page, 0);

  pages[0].resolve(['a', 'b']);
  await tick();
  assert.deepEqual(pager.items.value, ['a', 'b']);
  assert.equal(pager.state.value, 'idle');

  pager.loadMore();
  assert.equal(pages[1].page, 1);
  pages[1].resolve(['c', 'd']);
  await tick();
  assert.deepEqual(pager.items.value, ['a', 'b', 'c', 'd']);
  pager.dispose();
});

test('a short page marks done, and further loadMore() does nothing', async () => {
  const {pager, pages} = fake();
  pager.reset();
  pages[0].resolve(['a']);
  await tick();
  assert.equal(pager.state.value, 'done');
  assert.deepEqual(pager.items.value, ['a']);

  pager.loadMore();
  pager.loadMore();
  assert.equal(pages.length, 1, 'the collection is exhausted');
  assert.equal(pager.total.value, 1, 'an unknown total settles to the loaded count on done');
  pager.dispose();
});

test('a known total is not overwritten by the done fallback', async () => {
  const {pager, pages, counts} = fake();
  pager.reset();
  counts[0].resolve(7);
  pages[0].resolve(['a']);
  await tick();
  assert.equal(pager.state.value, 'done');
  assert.equal(pager.total.value, 7);
  pager.dispose();
});

test('loadMore() while a page is in flight does nothing', async () => {
  const {pager, pages} = fake();
  pager.reset();
  pager.loadMore();
  assert.equal(pages.length, 1);
  pages[0].resolve(['a', 'b']);
  await tick();
  pager.loadMore();
  assert.equal(pages.length, 2, 'the next page starts once the previous one lands');
  pager.dispose();
});

test('reset() drops the stale page and count, and re-reads the total', async () => {
  const {pager, pages, counts} = fake();
  pager.reset();
  counts[0].resolve(120);
  await tick();
  assert.equal(pager.total.value, 120);

  const stale = pages[0];
  pager.reset();
  assert.equal(stale.abort.aborted, true);
  assert.equal(counts.length, 2, 'a new generation counts again');
  assert.equal(pager.total.value, null);
  assert.deepEqual(pager.items.value, []);

  stale.resolve(['x', 'y']);
  counts[0].resolve(999);
  await tick();
  assert.deepEqual(pager.items.value, [], 'the outrun response never lands');
  assert.equal(pager.total.value, null);

  counts[1].resolve(4);
  pages[1].resolve(['a', 'b']);
  await tick();
  assert.deepEqual(pager.items.value, ['a', 'b']);
  assert.equal(pager.total.value, 4);
  pager.dispose();
});

test('error keeps the pages already loaded, and loadMore() retries the failed page', async () => {
  const {pager, pages} = fake();
  pager.reset();
  pages[0].resolve(['a', 'b']);
  await tick();

  pager.loadMore();
  pages[1].reject(new Error('gateway timeout'));
  await tick();
  assert.equal(pager.state.value, 'error');
  assert.deepEqual(pager.items.value, ['a', 'b']);

  pager.loadMore();
  assert.equal(pages[2].page, 1, 'the failed page is the one retried');
  pages[2].resolve(['c']);
  await tick();
  assert.deepEqual(pager.items.value, ['a', 'b', 'c']);
  assert.equal(pager.state.value, 'done');
  pager.dispose();
});

test('a failed count leaves the total unknown and the pages loading', async () => {
  const {pager, pages, counts} = fake();
  pager.reset();
  counts[0].reject(new Error('count unavailable'));
  pages[0].resolve(['a', 'b']);
  await tick();
  assert.equal(pager.total.value, null);
  assert.deepEqual(pager.items.value, ['a', 'b']);
  assert.equal(pager.state.value, 'idle');
  pager.dispose();
});

test('dispose: aborts what is in flight, ignores late responses, owns no scope', async () => {
  const live = Scope.liveCount;
  const {pager, pages, counts} = fake();
  pager.reset();
  pager.dispose();
  assert.equal(pages[0].abort.aborted, true);

  pages[0].resolve(['a', 'b']);
  counts[0].resolve(7);
  await tick();
  assert.deepEqual(pager.items.value, []);
  assert.equal(pager.total.value, null);

  pager.loadMore();
  pager.reset();
  assert.equal(pages.length, 1, 'a disposed pager never fetches again');
  assert.equal(Scope.liveCount, live, 'the pager owns no scope');
});
