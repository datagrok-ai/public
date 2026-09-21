import {test} from 'node:test';
import assert from 'node:assert/strict';
import {signal} from '../src/core/signals.js';
import {Filters} from '../src/core/filter/index.js';

test('queryPath: encoded once, an empty query is the base alone', () => {
  assert.equal(Filters.queryPath('/apps/U2Demo/filters', 'age > 30 and sex = "F"'),
    '/apps/U2Demo/filters?q=age%20%3E%2030%20and%20sex%20%3D%20%22F%22');
  assert.equal(Filters.queryPath('/apps/U2Demo/filters', ''), '/apps/U2Demo/filters');
  assert.equal(Filters.queryPath('/x', 'a%20b'), '/x?q=a%2520b', 'a literal percent is not decoded on the way in');
});

test('round trip through every character the grammar uses: URLSearchParams reads the q back', () => {
  for (const q of ['name like "a b" and [Mol Weight] between 200 and 500', 'x = "d\\\'Asie" or y in ("a&b", "c+d")',
    'created > -1w and (status in ("Open") or owner = @current)', 'q = "#tag" && z != 100%', 'ü = "ß"']) {
    const path = Filters.queryPath('/p', q);
    assert.equal(new URLSearchParams(path.slice(path.indexOf('?') + 1)).get('q'), q);
  }
});

test('queryPath over a signal follows it', () => {
  const q = signal('');
  const path = Filters.queryPath('/p', q);
  assert.equal(path.value, '/p');
  q.value = 'a = 1';
  assert.equal(path.value, '/p?q=a%20%3D%201');
});
