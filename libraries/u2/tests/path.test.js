/* The A4 binding-path grammar (WO-11): what parses, what assembles back, where the cycles are,
   and which signals take a write — the pure substrate under the resolver and the editor. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {signal, computed} from '../src/core/signals.js';
import {isWritableSignal} from '../src/core/signals.js';
import {parsePath, segmentsToPath, findBindingCycle} from '../src/spec/path.js';

test('parsePath: identifier chains', () => {
  assert.deepEqual(parsePath('$.a'), ['a']);
  assert.deepEqual(parsePath('$.orders.currentRow'), ['orders', 'currentRow']);
  assert.deepEqual(parsePath('$.a1._b.$c'), ['a1', '_b', '$c']);
  assert.deepEqual(parsePath('$._'), ['_']);
});

test('parsePath: bracket segments with escapes', () => {
  assert.deepEqual(parsePath('$.orders.currentRow[\'Mol Weight\']'), ['orders', 'currentRow', 'Mol Weight']);
  assert.deepEqual(parsePath('$.a[\'it\\\'s\']'), ['a', 'it\'s']);
  assert.deepEqual(parsePath('$.a[\'b\\\\c\']'), ['a', 'b\\c']);
  assert.deepEqual(parsePath('$.a[\'\']'), ['a', '']);
  assert.deepEqual(parsePath('$.a[\'x\'].b[\'y\']'), ['a', 'x', 'b', 'y']);
});

test('parsePath: malformed is null, never a throw', () => {
  const malformed = ['', 'x', '$', '$.', '$$.a', '$.1a', '$.a.', '$.a..b', '$.a b', '$. a',
    '$.a.1x', '$.a[b]', '$.a["b"]', '$.a[\'b\'', '$.a[\'b', '$.a[\'b\']x', '$.a[\'b\\q\']',
    '$.a\'', '$.a]'];
  for (const path of malformed)
    assert.equal(parsePath(path), null, `"${path}" must not parse`);
});

test('segmentsToPath: identifiers dotted, the rest bracketed — the exact inverse of parsePath', () => {
  assert.equal(segmentsToPath(['a']), '$.a');
  assert.equal(segmentsToPath(['orders', 'currentRow', 'Mol Weight']),
    '$.orders.currentRow[\'Mol Weight\']');
  assert.equal(segmentsToPath(['a', 'it\'s']), '$.a[\'it\\\'s\']');
  assert.equal(segmentsToPath(['a', 'b\\c']), '$.a[\'b\\\\c\']');
  for (const segments of [['a'], ['a', 'b'], ['a', 'x y'], ['a', 'it\'s', 'b\\c'], ['a', '', 'z']])
    assert.deepEqual(parsePath(segmentsToPath(segments)), segments);
});

test('findBindingCycle: names the loop, closed, or answers null', () => {
  const spec = (root, components) => ({$schema: 'dg-ui/1', root, components});
  const node = (name, bind) => ({tag: 'u2-x', name, bind});

  assert.equal(findBindingCycle(spec({tag: 'div', children: [
    node('a', {value: '$.b'}), node('b', {value: '$.ctxOnly'}),
  ]})), null, 'a bind out to context data is not an edge');

  assert.deepEqual(findBindingCycle(spec({tag: 'div', children: [
    node('a', {value: '$.b.deep'}), node('b', {value: '$.a'}),
  ]})), ['a', 'b', 'a'], 'first segments make the edges');

  assert.deepEqual(findBindingCycle(spec({tag: 'div', children: [node('a', {value: '$.a'})]})),
    ['a', 'a'], 'a self-loop counts');

  assert.deepEqual(findBindingCycle(spec(
    {tag: 'div', children: [node('daysInput', {value: '$.orders'})]},
    [node('orders', {'params.days': '$.daysInput'})]
  )), ['daysInput', 'orders', 'daysInput'], 'components and dotted sub-bind keys are walked too');

  assert.equal(findBindingCycle(spec({tag: 'div', children: [
    node('a', {value: 'not-a-path'}), node('b', {value: '$.a'}),
  ]})), null, 'a malformed path makes no edge');

  assert.equal(findBindingCycle(spec({tag: 'div'})), null);
});

test('isWritableSignal: a signal takes a write, a computed does not', () => {
  const s = signal(1);
  assert.equal(isWritableSignal(s), true);
  assert.equal(isWritableSignal(computed(() => s.value * 2)), false);
  assert.equal(isWritableSignal(42), false);
  assert.equal(isWritableSignal({value: 1}), false, 'a signal-shaped object is not a signal');
  assert.equal(isWritableSignal(undefined), false);
});
