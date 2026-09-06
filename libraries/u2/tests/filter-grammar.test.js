/* The smart-filter grammar port (schema-filters WO-2), driven by the shared corpus
   `tests/filter-grammar.corpus.json`. The Dart twin of this file is
   `core/shared/grok_shared/test/filter_grammar_corpus_test.dart` (WO-3): both read the SAME JSON,
   so an entry that passes here and fails there (or the reverse) is grammar drift — fix a parser,
   never fork the corpus. Rules both readers implement are in the corpus' `rules`; `canonical`
   and `position` are TS-only. Below the corpus: positions, `negate`, the token under the caret
   and the completion context the query input drives on. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {Filters} from '../src/core/filter/index.js';

const corpus = JSON.parse(readFileSync(fileURLToPath(new URL('./filter-grammar.corpus.json', import.meta.url)), 'utf8'));
const SPAN_TOLERANCE_MS = 5000;

/** `actual` with every Date that a `{$span}` expects and that lies within the tolerance of
 * `now + span` replaced by that `{$span}` object — so a plain deepStrictEqual shows the diff. */
function withSpans(actual, expected, now) {
  if (expected && typeof expected === 'object' && typeof expected.$span === 'string') {
    if (!(actual instanceof Date))
      return actual;
    const target = Filters.resolveSpan(expected.$span, now);
    return Math.abs(actual.getTime() - target.getTime()) <= SPAN_TOLERANCE_MS ? {$span: expected.$span} :
      `Date(${actual.toISOString()}) is not ${expected.$span} of ${now.toISOString()}`;
  }
  if (Array.isArray(expected) && Array.isArray(actual))
    return actual.map((a, i) => withSpans(a, expected[i], now));
  if (expected && typeof expected === 'object' && actual && typeof actual === 'object' && !(actual instanceof Date)) {
    const out = {};
    for (const k of Object.keys(actual))
      out[k] = withSpans(actual[k], expected[k], now);
    return out;
  }
  return actual;
}

function canonicalOf(text) {
  const {tree, errors} = Filters.parseTree(text);
  assert.deepEqual(errors, [], `${JSON.stringify(text)} must parse`);
  return Filters.format(Filters.fromDomainTree(tree));
}

test(`corpus: ${corpus.entries.length} entries, ids unique, every entry has a tree or errors`, () => {
  assert.ok(corpus.entries.length >= 50);
  assert.equal(new Set(corpus.entries.map((e) => e.id)).size, corpus.entries.length);
  for (const e of corpus.entries)
    assert.ok(e.errors === true ? e.tree === undefined : Array.isArray(e.tree), e.id);
});

for (const entry of corpus.entries) {
  test(`corpus ${entry.id}: ${JSON.stringify(entry.input)}`, () => {
    const now = new Date();
    const {tree, errors} = Filters.parseTree(entry.input, {now});
    if (entry.errors === true) {
      assert.ok(errors.length > 0, 'must fail');
      assert.equal(errors[0].code, 'syntax');
      assert.equal(errors[0].nodeId, null);
      assert.deepEqual(tree, []);
      if (entry.position)
        assert.deepEqual(errors[0].position, entry.position);
      return;
    }
    assert.deepEqual(errors, [], `must parse: ${errors[0]?.message}`);
    assert.deepEqual(withSpans(tree, entry.tree, now), entry.tree);
    assert.deepEqual(JSON.parse(JSON.stringify(withSpans(tree, entry.tree, now))), entry.tree, 'JSON sees no span');
    if (entry.canonical !== undefined) {
      assert.equal(Filters.format(Filters.fromDomainTree(tree)), entry.canonical);
      assert.equal(canonicalOf(entry.canonical), entry.canonical, 'the canonical string is a fixed point');
    }
  });
}

test('parseTree: blank text is the empty tree, a syntax problem carries the farthest position', () => {
  assert.deepEqual(Filters.parseTree(''), {tree: [], errors: []});
  assert.deepEqual(Filters.parseTree('  \t '), {tree: [], errors: []});
  const {errors} = Filters.parseTree('name like');
  assert.equal(errors.length, 1);
  assert.equal(errors[0].message, 'Expected a value');
  assert.deepEqual(errors[0].position, {start: 9, end: 9});
  assert.equal(Filters.parseTree('a = 1 or b').errors[0].message, 'Expected an operator');
  assert.deepEqual(Filters.parseTree('a = 1 or b').errors[0].position, {start: 10, end: 10});
  assert.equal(Filters.parseTree('x > 1.5h').errors[0].message, 'Expected a whole number of h, d, w, m or y');
  assert.deepEqual(Filters.parseTree('x > 1.5h').errors[0].position, {start: 4, end: 8});
});

test('parseTree: like interpolates the way Dart does — doubles keep .0, spans print as DateTime.toString()', () => {
  const now = new Date('2026-01-01T12:00:00.000Z');
  const like = (text) => Filters.parseTree(text, {now}).tree[0].value;
  assert.equal(like('x like 5k'), '%5000.0%');
  assert.equal(like('x like 1e3'), '%1000%');
  assert.equal(like('x like 1.5e2'), '%150.0%');
  assert.equal(like('x starts -3'), '-3%');
  assert.equal(like('created like -1h'), '%2026-01-01 11:00:00.000Z%');
  assert.equal(Filters.parseTree('name fuzzy 5.0', {now}).tree[0][2].value, '%5.0%');
  assert.equal(Filters.parseTree('created < 0h', {now}).tree[0].operator, '<', 'a zero span is not after now');
  assert.equal(Filters.parseTree('created < -0h', {now}).tree[0].operator, '<');
});

test('parseTree: spans resolve against the injected now and carry a hidden span', () => {
  const now = new Date('2026-09-04T12:00:00.000Z');
  const [flipped] = Filters.parseTree('created < 1h', {now}).tree;
  assert.equal(flipped.operator, '>');
  assert.equal(flipped.value.toISOString(), '2026-09-04T11:00:00.000Z');
  assert.equal(Filters.spanOf(flipped.value), '-1h');
  assert.deepEqual(Object.keys(flipped.value), []);
  assert.equal(JSON.stringify(flipped.value), '"2026-09-04T11:00:00.000Z"');
  const [now1] = Filters.parseTree('created <= now', {now}).tree;
  assert.equal(now1.value.getTime(), now.getTime());
  assert.equal(Filters.spanOf(now1.value), 'now');
  const [[low, , high]] = Filters.parseTree('created between 1d and -1d', {now}).tree;
  assert.equal(low.value.toISOString(), '2026-09-05T12:00:00.000Z', 'between keeps a positive span literal');
  assert.equal(high.value.toISOString(), '2026-09-03T12:00:00.000Z');
  assert.equal(Filters.spanOf(Filters.parseTree('created = -1m', {now}).tree[0].value), '-1m');
});

test('negate: the D10 table, lists keep their list, the fuzzy pair collapses', () => {
  const tree = [{property: 'a', operator: '=', value: 1}, 'and', {property: 'b', operator: '>', value: 2}, 'or',
    [{property: 'c', operator: '<', value: 3}, 'and', {property: 'd', operator: 'like', value: '%x%'}], 'and',
    {property: 'e', operator: '!=', value: [1, 2]}, 'and', {property: 'f', operator: '~*', value: 'x'}, 'and',
    {property: 'g', operator: 'is', value: null}, 'and',
    [{property: 'h', operator: 'fuzzy', threshold: 0.5, value: 'q'}, 'or', {property: 'h', operator: 'like', value: '%q%'}]];
  assert.deepEqual(Filters.negate(tree), [
    {property: 'a', operator: '!=', value: 1}, 'or', {property: 'b', operator: '<=', value: 2}, 'and',
    [{property: 'c', operator: '>=', value: 3}, 'or', {property: 'd', operator: 'not like', value: '%x%'}], 'or',
    {property: 'e', operator: '=', value: [1, 2]}, 'or', {property: 'f', operator: '!~*', value: 'x'}, 'or',
    {property: 'g', operator: 'is not', value: null}, 'or', {property: 'h', operator: 'not like', value: '%q%'}]);
  assert.deepEqual(Filters.negate(Filters.negate(tree.slice(0, 12))), tree.slice(0, 12), 'an involution but for the pair');
  assert.throws(() => Filters.negate([{property: 'a', operator: 'fuzzy', value: 'x'}]), /Cannot negate "fuzzy"/);
});

const schema = Filters.schema([
  {name: 'name', type: 'string'}, {name: 'age', type: 'int'}, {name: 'created', type: 'datetime'},
  {name: 'status', type: 'string'},
]);

function ctx(textWithCaret) {
  const caret = textWithCaret.indexOf('|');
  return Filters.completionContext(textWithCaret.replace('|', ''), caret, schema);
}

test('token under the caret: names with dots and brackets, strings with either quote, numbers with units, ' +
  'symbols, tags', () => {
  const at = (textWithCaret) => {
    const {prefix, replace} = ctx(textWithCaret);
    return [prefix, replace.start, replace.end];
  };
  assert.deepEqual(at('[Mo|l Weight].x >= 5k && name !like \'ab\' or #tag'), ['[Mo', 0, 14], 'a bracketed dotted name');
  assert.deepEqual(at('[Mol Weight].x >|= 5k && name !like \'ab\' or #tag'), ['>', 15, 17]);
  assert.deepEqual(at('[Mol Weight].x >= 5|k && name !like \'ab\' or #tag'), ['5', 18, 20], 'a number keeps its unit');
  assert.deepEqual(at('[Mol Weight].x >= 5k &|& name !like \'ab\' or #tag'), ['&', 21, 23]);
  assert.deepEqual(at('[Mol Weight].x >= 5k && name !l|ike \'ab\' or #tag'), ['!l', 29, 34]);
  assert.deepEqual(at('[Mol Weight].x >= 5k && name !like \'a|b\' or #tag'), ['a', 35, 39], 'a string without its quote');
  assert.deepEqual(at('[Mol Weight].x >= 5k && name !like \'ab\' or #t|ag'), ['#t', 43, 47]);
  assert.deepEqual(at('(owner=@cur|rent)'), ['@cur', 7, 15]);
  assert.deepEqual(at('a = "unterminated|'), ['unterminated', 4, 17]);
  assert.deepEqual(at('n = -1.5e3|'), ['-1.5e3', 4, 10]);
  assert.deepEqual(at('a *|* b'), ['*', 2, 3], 'an unknown symbol is one character');
  assert.deepEqual(at('a \'x"y|\''), ['x"y', 2, 7], 'only the opening quote closes');
});

test('completionContext: property at the start, after a connector, after "(" and after "not"', () => {
  assert.deepEqual(ctx('|'), {expect: 'property', prefix: '', replace: {start: 0, end: 0}});
  assert.deepEqual(ctx('na|'), {expect: 'property', prefix: 'na', replace: {start: 0, end: 2}});
  assert.deepEqual(ctx('age > 5 and |'), {expect: 'property', prefix: '', replace: {start: 12, end: 12}});
  assert.deepEqual(ctx('age > 5 && st|'), {expect: 'property', prefix: 'st', replace: {start: 11, end: 13}});
  assert.deepEqual(ctx('age > 5 and (|'), {expect: 'property', prefix: '', replace: {start: 13, end: 13}});
  assert.deepEqual(ctx('not (na|'), {expect: 'property', prefix: 'na', replace: {start: 5, end: 7}});
  assert.equal(ctx('not |').expect, 'property');
});

test('completionContext: operator after a property, with the schema property resolved', () => {
  const c = ctx('name |');
  assert.equal(c.expect, 'operator');
  assert.equal(c.property.name, 'name');
  assert.deepEqual(c.replace, {start: 5, end: 5});
  const typed = ctx('name li|');
  assert.equal(typed.expect, 'operator');
  assert.equal(typed.prefix, 'li');
  assert.deepEqual(typed.replace, {start: 5, end: 7});
  assert.equal(ctx('[Mol Weight] |').expect, 'operator');
  assert.equal(ctx('[Mol Weight] |').property, undefined, 'unknown property stays undefined');
  assert.equal(ctx('age > 5 name |').expect, 'operator', 'juxtaposition starts a condition');
});

test('completionContext: value after an operator, inside a string, inside an in-list, between bounds', () => {
  const v = ctx('name like |');
  assert.equal(v.expect, 'value');
  assert.equal(v.operator.id, 'like');
  assert.equal(v.property.name, 'name');
  const s = ctx('name = "as|pirin"');
  assert.equal(s.expect, 'value');
  assert.equal(s.prefix, 'as');
  assert.deepEqual(s.replace, {start: 7, end: 16});
  assert.equal(ctx('name = "aspirin"|').expect, 'connector', 'after the closing quote the value is done');
  assert.equal(ctx('age > 5|').expect, 'value', 'a number is still being typed');
  assert.equal(ctx('age > 5|').prefix, '5');
  const list = ctx('status in ("Open", |');
  assert.equal(list.expect, 'value');
  assert.equal(list.operator.id, 'in');
  assert.equal(ctx('status not in (|').operator.id, 'not in');
  assert.equal(ctx('status !in ("a"|').expect, 'value');
  assert.equal(ctx('status in ("a", "b")|').expect, 'connector');
  assert.equal(ctx('age between 1 and |').expect, 'value');
  assert.equal(ctx('age between 1 and |').operator.id, 'between');
  assert.equal(ctx('age between 1 and 5 |').expect, 'connector');
  assert.equal(ctx('name fuzzy(0.6) |').expect, 'value');
  assert.equal(ctx('name fuzzy(0.6) |').operator.id, 'fuzzy');
  assert.equal(ctx('name ~ |').operator.id, 'matches');
  assert.equal(ctx('created > -1|').expect, 'value');
});

test('completionContext: connector after a value, a tag or ")"', () => {
  assert.deepEqual(ctx('age > 5 |'), {expect: 'connector', prefix: '', replace: {start: 8, end: 8},
    property: schema.properties[1], operator: Filters.operators.get('>')});
  assert.equal(ctx('age > 5 a|').expect, 'connector');
  assert.equal(ctx('age > 5 a|').prefix, 'a');
  assert.equal(ctx('#chem |').expect, 'connector');
  assert.equal(ctx('(age > 5) |').expect, 'connector');
  assert.equal(ctx('not hidden |').expect, 'connector');
  assert.equal(ctx('age > 5 not |').expect, 'property', 'a juxtaposed not starts a condition');
  assert.equal(ctx('age > 5 not na|').prefix, 'na');
  assert.equal(ctx('age > 5 not hidden |').expect, 'connector');
  assert.equal(Filters.completionContext('name = 1 not ', 13, schema).expect, 'property');
});
