const test = require('node:test');
const assert = require('node:assert/strict');
const DG = require('./dg.cjs');

test('LruCache evicts the least recently used entry and reports it', () => {
  const cache = new DG.LruCache(2);
  const evicted = [];
  cache.onItemEvicted = (v) => evicted.push(v);
  cache.set('a', 1);
  cache.set('b', 2);
  cache.get('a');
  cache.set('c', 3);
  assert.equal(cache.has('b'), false);
  assert.equal(cache.get('a'), 1);
  assert.equal(cache.get('c'), 3);
  assert.deepEqual(evicted, [2]);
});

test('StringUtils.levenshteinDistance is normalized to [0, 1]', () => {
  assert.equal(DG.StringUtils.levenshteinDistance('same', 'same'), 0);
  const d = DG.StringUtils.levenshteinDistance('kitten', 'sitting');
  assert.ok(d > 0 && d <= 1, `got ${d}`);
});

test('StringUtils.hashCode is stable and discriminates', () => {
  assert.equal(DG.StringUtils.hashCode('abc'), DG.StringUtils.hashCode('abc'));
  assert.notEqual(DG.StringUtils.hashCode('abc'), DG.StringUtils.hashCode('abd'));
});

test('Point distance', () => {
  assert.equal(new DG.Point(0, 0).distanceTo(new DG.Point(3, 4)), 5);
});
