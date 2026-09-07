/* `Filters.toMask` and the `ColumnEvaluator` leaves (schema-filters WO-4): every core operator per column
   kind over literal raw-data fakes (`INT_NULL`, `FLOAT_NULL`, category indexes, bool bits, bigint
   `get`), per-category string evaluation, group combination with the tail beyond `length`, async
   `bitset` operators, abort. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {BitArray, SEMTYPE, TYPE} from 'datagrok-api/u2core';
import {Filters, FilterError, ColumnEvaluator, KIND} from '../src/core/filter/index.js';
import {INT_NULL, FLOAT_NULL} from '../src/core/filter/evaluate.js';

const NOW = new Date('2026-09-04T12:00:00.000Z');
const ctx = {now: NOW};
const micros = (iso) => Date.parse(iso) * 1000;

function column(name, type, raw, extra = {}) {
  return {name, type, length: extra.length ?? raw.length, getRawData: () => raw, ...extra};
}

function bits(values) {
  const raw = new Uint32Array((values.length + 31) >>> 5);
  values.forEach((v, i) => {
    if (v)
      raw[i >>> 5] |= 1 << (i & 31);
  });
  return raw;
}

const ints = column('n', TYPE.INT, Int32Array.from([1, 5, INT_NULL, 10, -3]));
const floats32 = column('f', TYPE.FLOAT, Float32Array.from([0.1, 2.5, FLOAT_NULL, NaN, 7]));
const floats64 = column('g', TYPE.FLOAT, Float64Array.from([0.1, 2.5, FLOAT_NULL, NaN, 7]));
const dates = column('d', TYPE.DATE_TIME, Float64Array.from([
  micros('2026-09-01T00:00:00Z'), micros('2026-09-03T12:00:00Z'), FLOAT_NULL, micros('2026-09-04T12:00:00Z'),
  micros('2025-01-01T00:00:00Z')]));
const bools = column('b', TYPE.BOOL, bits([true, false, true, false, false]), {length: 5});
const strings = column('s', TYPE.STRING, Int32Array.from([0, 1, 2, 0, 3]),
  {categories: ['Aspirin', 'ibuprofen', '', 'Asp']});
const bigs = column('big', TYPE.BIG_INT, new Int32Array(0),
  {length: 5, get: (i) => [12345678901234567890n, null, 5n, -1n, 7n][i]});
const lists = column('tags', TYPE.LIST, new Int32Array(0),
  {length: 5, get: (i) => [['Aspirin', 'nsaid'], null, [], ['Ibuprofen'], ['x', 'ASP']][i]});
const frame = {
  rowCount: 5,
  column: (name) => [ints, floats32, floats64, dates, bools, strings, bigs, lists].find((c) => c.name === name) ?? null,
};

async function rows(nodes, options = {}) {
  const root = Array.isArray(nodes) ? Filters.group(options.op ?? 'and', nodes, options.not ? {not: true} : undefined) : nodes;
  return Array.from((await Filters.toMask(frame, root, {now: NOW, ...options})).getSelectedIndexes());
}

const one = (property, operator, value, options) => rows([Filters.cond(property, operator, value, options)]);

test('BitArray: fill, invert, and, or keep the tail beyond length clear', () => {
  const all = new BitArray(35, true);
  assert.equal(all.getBuffer().length, 2);
  assert.equal(all.getBuffer()[1], 0b111);
  assert.deepEqual(Array.from(new BitArray(35).getSelectedIndexes()), []);
  const even = BitArray.create(35, (i) => i % 2 === 0);
  const odd = even.clone().invert();
  assert.deepEqual(Array.from(odd.getSelectedIndexes()), [...Array(35).keys()].filter((i) => i % 2 === 1));
  assert.equal(odd.getBuffer()[1] >>> 3, 0, 'invert clears the tail');
  assert.deepEqual(Array.from(even.clone().and(BitArray.create(35, (i) => i < 4)).getSelectedIndexes()), [0, 2]);
  assert.deepEqual(Array.from(BitArray.create(35, (i) => i === 34).or(BitArray.create(35, (i) => i === 1)).getSelectedIndexes()),
    [1, 34]);
  assert.equal(even.get(34), true);
  assert.equal(even.get(33), false);
});

test('int: every operator; INT_NULL matches only is null', async () => {
  assert.deepEqual(await one('n', '=', 5), [1]);
  assert.deepEqual(await one('n', '!=', 5), [0, 3, 4]);
  assert.deepEqual(await one('n', '>', 1), [1, 3]);
  assert.deepEqual(await one('n', '>=', 1), [0, 1, 3]);
  assert.deepEqual(await one('n', '<', 1), [4]);
  assert.deepEqual(await one('n', '<=', 1), [0, 4]);
  assert.deepEqual(await one('n', 'between', [1, 5]), [0, 1]);
  assert.deepEqual(await one('n', 'in', [10, -3, 99]), [3, 4]);
  assert.deepEqual(await one('n', 'not in', [10, -3]), [0, 1]);
  assert.deepEqual(await one('n', 'is null'), [2]);
  assert.deepEqual(await one('n', 'is not null'), [0, 1, 3, 4]);
  assert.deepEqual(await one('n', '=', '5'), [1], 'a string value converts to the column kind');
});

test('float: Float32 compares against the rounded value; FLOAT_NULL and NaN are null in both widths', async () => {
  assert.deepEqual(await one('f', '=', 0.1), [0]);
  assert.deepEqual(await one('g', '=', 0.1), [0]);
  assert.deepEqual(await one('f', '>', 1), [1, 4]);
  assert.deepEqual(await one('g', '<=', 2.5), [0, 1]);
  assert.deepEqual(await one('f', 'between', [0, 3]), [0, 1]);
  assert.deepEqual(await one('f', 'is null'), [2, 3]);
  assert.deepEqual(await one('g', 'is null'), [2, 3]);
  assert.deepEqual(await one('f', 'is not null'), [0, 1, 4]);
  assert.deepEqual(await one('f', '!=', 7), [0, 1]);
});

test('datetime: Date values against µs raw data, spans against the injected now, ISO strings', async () => {
  assert.deepEqual(await one('d', '=', new Date('2026-09-03T12:00:00Z')), [1]);
  assert.deepEqual(await one('d', '=', '2026-09-03T12:00:00.000Z'), [1]);
  assert.deepEqual(await one('d', '>', new Date('2026-09-01T00:00:00Z')), [1, 3]);
  assert.deepEqual(await one('d', '>=', {span: '-1d'}), [1, 3]);
  assert.deepEqual(await one('d', '<', {span: 'now'}), [0, 1, 4]);
  assert.deepEqual(await one('d', '<=', {span: 'now'}), [0, 1, 3, 4]);
  assert.deepEqual(await one('d', 'between', [{span: '-1w'}, new Date('2026-09-03T12:00:00Z')]), [0, 1]);
  assert.deepEqual(await one('d', 'is null'), [2]);
  assert.deepEqual(await one('d', 'is not null'), [0, 1, 3, 4]);
  const later = Array.from((await Filters.toMask(frame,
    Filters.group('and', [Filters.cond('d', '>', {span: '-1d'})]), {now: new Date('2026-09-10T00:00:00Z')})).getSelectedIndexes());
  assert.deepEqual(later, [], 'the span resolves at call time');
});

test('bool: bits LSB-first; is null never matches, is not null always', async () => {
  assert.deepEqual(await one('b', '=', true), [0, 2]);
  assert.deepEqual(await one('b', '=', false), [1, 3, 4]);
  assert.deepEqual(await one('b', '!=', true), [1, 3, 4]);
  assert.deepEqual(await one('b', '!=', false), [0, 2]);
  assert.deepEqual(await one('b', 'is null'), []);
  assert.deepEqual(await one('b', 'is not null'), [0, 1, 2, 3, 4]);
});

test('string: = is exact, like/starts/ends/matches are case-insensitive, "" is null, lists, raw LIKE patterns', async () => {
  assert.deepEqual(await one('s', '=', 'Aspirin'), [0, 3]);
  assert.deepEqual(await one('s', '=', 'aspirin'), []);
  assert.deepEqual(await one('s', '!=', 'Aspirin'), [1, 4], 'the null category never matches');
  assert.deepEqual(await one('s', 'like', 'ASP'), [0, 3, 4]);
  assert.deepEqual(await one('s', '!like', 'asp'), [1]);
  assert.deepEqual(await one('s', 'starts', 'ib'), [1]);
  assert.deepEqual(await one('s', 'ends', 'RIN'), [0, 3]);
  assert.deepEqual(await one('s', 'matches', '^a.*n$'), [0, 3]);
  assert.deepEqual(await one('s', '!matches', '^a'), [1]);
  assert.deepEqual(await one('s', 'in', ['Asp', 'ibuprofen']), [1, 4]);
  assert.deepEqual(await one('s', 'not in', ['Asp', 'ibuprofen']), [0, 3]);
  assert.deepEqual(await one('s', 'is null'), [2]);
  assert.deepEqual(await one('s', 'is not null'), [0, 1, 3, 4]);
  assert.deepEqual(await one('s', 'like', 'a%n', {options: {raw: true}}), [0, 3]);
  assert.deepEqual(await one('s', 'like', 'as_', {options: {raw: true}}), [4]);
  assert.deepEqual(await one('s', '!like', '%i%', {options: {raw: true}}), [4]);
  assert.equal(ColumnEvaluator.likeRegExp('50\\%').test('50%'), true, 'an escaped wildcard is literal');
  assert.equal(ColumnEvaluator.likeRegExp('50\\%').test('500'), false);
  assert.equal(ColumnEvaluator.likeRegExp('a.b%').test('A.BC'), true, 'regex metacharacters are literal');
  assert.equal(ColumnEvaluator.likeRegExp('a.b%').test('axb'), false);
  assert.deepEqual(await one('s', '=', {type: 'User', id: 'Asp'}), [4], 'a ref compares by id');
});

test('string: the predicate runs once per category, not per row', () => {
  let calls = 0;
  const col = column('s', TYPE.STRING, Int32Array.from([0, 1, 0, 1, 0, 1, 2]), {categories: ['a', 'b', '']});
  const mask = ColumnEvaluator.where(col, Filters.cond('s', '=', 'a'), ctx, (cell, value) => (calls++, cell === value));
  assert.equal(calls, 2, 'the null category is skipped');
  assert.deepEqual(Array.from(mask.getSelectedIndexes()), [0, 2, 4]);
});

test('bigint: rows go through get(i); number and digit-string values compare as BigInt', async () => {
  assert.deepEqual(await one('big', '=', '12345678901234567890'), [0]);
  assert.deepEqual(await one('big', '=', 5), [2]);
  assert.deepEqual(await one('big', '>', 0), [0, 2, 4]);
  assert.deepEqual(await one('big', '<=', 5), [2, 3]);
  assert.deepEqual(await one('big', 'between', ['-1', '5']), [2, 3]);
  assert.deepEqual(await one('big', 'in', ['5', -1]), [2, 3]);
  assert.deepEqual(await one('big', 'not in', [5]), [0, 3, 4]);
  assert.deepEqual(await one('big', 'is null'), [1]);
  assert.deepEqual(await one('big', 'is not null'), [0, 2, 3, 4]);
});

test('string list: rows go through get(i); like is "some element contains", !like "none does"; null is null', async () => {
  assert.deepEqual(await one('tags', 'like', 'asp'), [0, 4]);
  assert.deepEqual(await one('tags', '!like', 'asp'), [2, 3], 'an empty list contains nothing; null never passes');
  assert.deepEqual(await one('tags', 'like', 'Aspirin%', {options: {raw: true}}), [0]);
  assert.deepEqual(await one('tags', 'is null'), [1]);
  assert.deepEqual(await one('tags', 'is not null'), [0, 2, 3, 4]);
});

test('groups: and, or, not, nesting; an empty group is all-true; the tail beyond length stays clear', async () => {
  const gt1 = Filters.cond('n', '>', 1);
  const isTrue = Filters.cond('b', '=', true);
  assert.deepEqual(await rows([gt1, isTrue]), []);
  assert.deepEqual(await rows([gt1, isTrue], {op: 'or'}), [0, 1, 2, 3]);
  assert.deepEqual(await rows([gt1], {not: true}), [0, 2, 4], 'not includes the null row');
  assert.deepEqual(await rows([Filters.group('or', [gt1, isTrue]), Filters.cond('s', 'like', 'asp')]), [0, 3]);
  assert.deepEqual(await rows([Filters.group('or', [gt1, isTrue], {not: true}), Filters.cond('n', 'is not null')]), [4]);
  assert.deepEqual(await rows([]), [0, 1, 2, 3, 4]);
  assert.deepEqual(await rows([], {op: 'or'}), [0, 1, 2, 3, 4]);
  assert.deepEqual(await rows([], {not: true}), []);
  assert.deepEqual(await rows([Filters.group('and', []), gt1]), [1, 3], 'an empty sub-group constrains nothing');
  const wide = {rowCount: 40, column: (name) => name === 'x' ?
    column('x', TYPE.INT, Int32Array.from({length: 40}, (_v, i) => i)) : null};
  const mask = await Filters.toMask(wide, Filters.group('and', [Filters.cond('x', '<', 3)], {not: true}));
  assert.equal(mask.length, 40);
  assert.equal(mask.getBuffer().length, 2);
  assert.equal(mask.getBuffer()[1] >>> 8, 0);
  assert.equal(mask.trueCount, 37);
});

test('bitset operators run async between sync ones and see the column and the signal', async () => {
  const seen = [];
  const off = Filters.operators.register({
    id: 'Contains', label: 'Contains', arity: 1, kinds: [KIND.STRING], semType: SEMTYPE.MOLECULE, editor: 'default',
    bitset: async (col, c, signal) => {
      seen.push([col.name, c.value, signal instanceof AbortSignal]);
      await Promise.resolve();
      return BitArray.create(col.length, (i) => (col.categories[col.getRawData()[i]] ?? '').includes(c.value));
    },
  });
  try {
    const mol = {...strings, name: 'smiles', semType: SEMTYPE.MOLECULE};
    const f = {rowCount: 5, column: (name) => name === 'smiles' ? mol : name === 'n' ? ints : null};
    const root = Filters.group('and', [Filters.cond('smiles', 'Contains', 'sp'), Filters.cond('n', '>=', 5)]);
    assert.deepEqual(Array.from((await Filters.toMask(f, root)).getSelectedIndexes()), [3]);
    assert.deepEqual(seen, [['smiles', 'sp', true]]);
    assert.deepEqual(Array.from((await Filters.toMask(f, Filters.group('or', [Filters.cond('smiles', 'like', 'ibu'),
      Filters.cond('smiles', 'Contains', 'sp')]))).getSelectedIndexes()), [0, 1, 3, 4], 'core operators still apply to the semType');
  } finally {
    off();
  }
});

test('errors: an unknown column and an operator without a DataFrame form throw FilterError with the node id', async () => {
  Filters.resetIds('');
  await assert.rejects(rows([Filters.cond('nope', '=', 1)]), (e) => e instanceof FilterError &&
    e.problems[0].code === 'unknown-property' && e.problems[0].nodeId === 'f1');
  Filters.resetIds('');
  await assert.rejects(rows([Filters.cond('s', 'fuzzy', 'x')]), (e) => e instanceof FilterError &&
    e.problems[0].code === 'not-expressible' && e.problems[0].nodeId === 'f1');
  await assert.rejects(rows([Filters.cond('s', 'no-such-op', 'x')]), FilterError);
});

test('errors: a mask of the wrong length is a FilterError on its leaf, not a RangeError from the group', async () => {
  const off = Filters.operators.register({
    id: 'Short', label: 'Short', arity: 1, kinds: [KIND.STRING], semType: SEMTYPE.MOLECULE, editor: 'default',
    bitset: async () => new BitArray(3),
  });
  try {
    const mol = {...strings, name: 'smiles', semType: SEMTYPE.MOLECULE};
    const f = {rowCount: 5, column: (name) => name === 'smiles' ? mol : name === 'n' ? ints : null};
    Filters.resetIds('');
    const root = Filters.group('and', [Filters.cond('n', '>=', 5), Filters.cond('smiles', 'Short', 'x')]);
    await assert.rejects(Filters.toMask(f, root), (e) => e instanceof FilterError &&
      e.problems[0].code === 'evaluation' && e.problems[0].nodeId === 'f2' && /"smiles Short".*3 bits.*5 rows/.test(e.message));
  } finally {
    off();
  }
});

test('abort: an aborted signal rejects before the first leaf and between leaves', async () => {
  const aborted = new AbortController();
  aborted.abort();
  await assert.rejects(rows([Filters.cond('n', '=', 1)], {signal: aborted.signal}), (e) => e.name === 'AbortError');
  const controller = new AbortController();
  let leaves = 0;
  const off = Filters.operators.register({
    id: 'slow', label: 'slow', arity: 1, kinds: [KIND.STRING], semType: 'T', editor: 'default',
    bitset: async (col) => {
      leaves++;
      controller.abort();
      return new BitArray(col.length, true);
    },
  });
  try {
    const t = {...strings, semType: 'T'};
    const f = {rowCount: 5, column: () => t};
    const root = Filters.group('and', [Filters.cond('s', 'slow', 'x'), Filters.cond('s', 'slow', 'y')]);
    await assert.rejects(Filters.toMask(f, root, {signal: controller.signal}), (e) => e.name === 'AbortError');
    assert.equal(leaves, 1, 'the second leaf never ran');
  } finally {
    off();
  }
});
