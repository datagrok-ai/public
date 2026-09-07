/* `Filters.fromDomainTree` / `toDomainTree` (schema-filters WO-2): the AND-over-OR grouping with
   sticky connectors, every fold and its inverse, the model-level round trip over the whole corpus,
   the hidden span that keeps relative dates live, schema typing, and `Filters.parse`. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {Filters, KIND} from '../src/core/filter/index.js';
import {TYPE, SEMTYPE} from 'datagrok-api/u2core';

const corpus = JSON.parse(readFileSync(fileURLToPath(new URL('./filter-grammar.corpus.json', import.meta.url)), 'utf8'));
const c = (property, operator, value) => value === undefined ? {property, operator} : {property, operator, value};
const shape = (root) => Filters.isGroup(root) ?
  {op: root.op, ...(root.not ? {not: true} : {}), nodes: root.nodes.map(shape)} :
  {property: root.property, operator: root.operator, ...(root.value !== undefined ? {value: root.value} : {}),
    ...(root.options ? {options: root.options} : {})};

test('fromDomainTree: flat and/or lists, AND-over-OR grouping, sticky and omitted connectors', () => {
  assert.deepEqual(shape(Filters.fromDomainTree([c('a', '=', 1), 'and', c('b', '=', 2)])),
    {op: 'and', nodes: [{property: 'a', operator: '=', value: 1}, {property: 'b', operator: '=', value: 2}]});
  assert.deepEqual(shape(Filters.fromDomainTree([c('a', '=', 1), 'or', c('b', '=', 2)])),
    {op: 'or', nodes: [{property: 'a', operator: '=', value: 1}, {property: 'b', operator: '=', value: 2}]});
  const mixed = Filters.fromDomainTree([c('a', '=', 1), 'or', c('b', '=', 2), 'and', c('c', '=', 3), 'or', c('d', '=', 4)]);
  assert.deepEqual(shape(mixed), {op: 'or', nodes: [
    {property: 'a', operator: '=', value: 1},
    {op: 'and', nodes: [{property: 'b', operator: '=', value: 2}, {property: 'c', operator: '=', value: 3}]},
    {property: 'd', operator: '=', value: 4}]});
  assert.equal(shape(Filters.fromDomainTree([c('a', '=', 1), c('b', '=', 2)])).op, 'and', 'omitted = and at first');
  assert.deepEqual(shape(Filters.fromDomainTree([c('a', '=', 1), 'or', c('b', '=', 2), c('c', '=', 3)])).nodes.length, 3,
    'an omitted connector repeats the last one seen');
  assert.deepEqual(shape(Filters.fromDomainTree(c('a', '=', 1))), {op: 'and', nodes: [{property: 'a', operator: '=', value: 1}]});
  assert.deepEqual(shape(Filters.fromDomainTree([])), {op: 'and', nodes: []});
  assert.deepEqual(shape(Filters.fromDomainTree([[c('a', '=', 1), 'or', c('b', '=', 2)]])).op, 'or',
    'a lone sub-list becomes the root');
  assert.deepEqual(shape(Filters.fromDomainTree([[c('a', '=', 1)], 'and', []])),
    {op: 'and', nodes: [{property: 'a', operator: '=', value: 1}]}, 'single-element sub-lists inline, empties drop');
});

test('fromDomainTree: the folds — in, is null, like shapes, regex, fuzzy pair, between, brackets', () => {
  const nodes = Filters.fromDomainTree([
    c('a', '=', [1, 2]), 'and', c('b', '!=', ['x']), 'and', c('d', '=', null), 'and', c('e', '!=', null), 'and',
    c('f', 'is', null), 'and', c('g', 'is not', null), 'and', c('h', 'like', '%x%'), 'and', c('i', 'like', 'x%'), 'and',
    c('j', 'like', '%x'), 'and', c('k', 'not like', '%x%'), 'and', c('l', '~*', 'r'), 'and', c('m', '!~*', 'r'), 'and',
    c('n', 'like', '%a%b%'), 'and', c('o', 'like', '%50\\%%'), 'and', c('p', 'not like', 'x%'), 'and',
    [{property: 'q', operator: 'fuzzy', threshold: 0.6, value: 'z'}, 'or', c('q', 'like', '%z%')], 'and',
    [{property: 'q2', operator: 'fuzzy', threshold: null, value: 'z'}, 'or', c('q2', 'like', '%z%')], 'and',
    [c('r', '>=', 1), 'and', c('r', '<=', 5)], 'and', [c('s', '>=', 1), 'and', c('t', '<=', 5)], 'and',
    c('[Mol Weight]', '>', 1), 'and', c('[a\\]b\\\\c].x', '=', 1), 'and', c('u', 'like', ''), 'and', c('v', 'like', '%'),
  ]).nodes.map(shape);
  assert.deepEqual(nodes, [
    {property: 'a', operator: 'in', value: [1, 2]}, {property: 'b', operator: 'not in', value: ['x']},
    {property: 'd', operator: 'is null'}, {property: 'e', operator: 'is not null'},
    {property: 'f', operator: 'is null'}, {property: 'g', operator: 'is not null'},
    {property: 'h', operator: 'like', value: 'x'}, {property: 'i', operator: 'starts', value: 'x'},
    {property: 'j', operator: 'ends', value: 'x'}, {property: 'k', operator: '!like', value: 'x'},
    {property: 'l', operator: 'matches', value: 'r'}, {property: 'm', operator: '!matches', value: 'r'},
    {property: 'n', operator: 'like', value: '%a%b%', options: {raw: true}},
    {property: 'o', operator: 'like', value: '50%'},
    {property: 'p', operator: '!like', value: 'x%', options: {raw: true}},
    {property: 'q', operator: 'fuzzy', value: 'z', options: {threshold: 0.6}},
    {property: 'q2', operator: 'fuzzy', value: 'z'},
    {property: 'r', operator: 'between', value: [1, 5]},
    {op: 'and', nodes: [{property: 's', operator: '>=', value: 1}, {property: 't', operator: '<=', value: 5}]},
    {property: 'Mol Weight', operator: '>', value: 1}, {property: 'a]b\\c.x', operator: '=', value: 1},
    {property: 'u', operator: 'like', value: '', options: {raw: true}},
    {property: 'v', operator: 'ends', value: ''},
  ]);
});

test('fromDomainTree: a condition that is not a plain object (a Dart map through interop) is a FilterError', () => {
  class DartMap { get property() { return 'a'; } }
  for (const [bad, message] of [[new DartMap(), /not a plain object/], [{operator: '='}, /has no property/]]) {
    assert.throws(() => Filters.fromDomainTree([c('a', '=', 1), 'and', bad]),
      (e) => e.name === 'FilterError' && message.test(e.message));
  }
});

test('fromDomainTree: a root-level >=/<= pair is not a between; ids are fresh', () => {
  Filters.resetIds('');
  const root = Filters.fromDomainTree([c('r', '>=', 1), 'and', c('r', '<=', 5)]);
  assert.deepEqual(root.nodes.map((n) => n.operator), ['>=', '<=']);
  assert.deepEqual(root.nodes.map((n) => n.id), ['f1', 'f2']);
  assert.equal(root.id, 'f3');
});

test('fromDomainTree with a schema: ref ids become FilterRef, ISO strings become Dates, @current stays', () => {
  const schema = Filters.schema([{name: 'owner', type: TYPE.STRING, ref: 'Core.users'}, {name: 'created', type: TYPE.DATE_TIME},
    {name: 'name', type: TYPE.STRING}]);
  const root = Filters.fromDomainTree([c('owner', '=', 'u1'), 'and', c('owner', '=', ['u1', 'u2']), 'and',
    c('owner', '=', '@current'), 'and', c('created', '>', '2026-01-02T00:00:00.000Z'), 'and',
    c('created', '=', ['2026-01-02T00:00:00.000Z']), 'and', c('name', '=', '2026-01-02T00:00:00.000Z'), 'and',
    c('owner.login', '=', 'x')], schema);
  const values = root.nodes.map((n) => n.value);
  assert.deepEqual(values[0], {type: 'Core.users', id: 'u1'});
  assert.deepEqual(values[1], [{type: 'Core.users', id: 'u1'}, {type: 'Core.users', id: 'u2'}]);
  assert.equal(values[2], '@current');
  assert.equal(values[3].getTime(), Date.parse('2026-01-02T00:00:00.000Z'));
  assert.equal(values[4][0].getTime(), Date.parse('2026-01-02T00:00:00.000Z'));
  assert.equal(values[5], '2026-01-02T00:00:00.000Z', 'a string property keeps the string');
  assert.equal(values[6], 'x', 'a dotted path is not typed');
  assert.equal(root.nodes[1].operator, 'in');
});

test('toDomainTree: the inverse of every fold, groups nested, single-node groups inlined, empties dropped', () => {
  const root = Filters.group('and', [
    Filters.cond('a', 'in', [1, 2]), Filters.cond('b', 'not in', ['x']), Filters.cond('d', 'is null'),
    Filters.cond('e', 'is not null'), Filters.cond('h', 'like', 'x'), Filters.cond('i', 'starts', 'x'),
    Filters.cond('j', 'ends', 'x'), Filters.cond('k', '!like', 'x'), Filters.cond('l', 'matches', 'r'),
    Filters.cond('m', '!matches', 'r'), Filters.cond('n', 'like', '%a%b%', {options: {raw: true}}),
    Filters.cond('o', 'like', '50%_\\'), Filters.cond('q', 'fuzzy', 'z', {options: {threshold: 0.6}}),
    Filters.cond('q2', 'fuzzy', 'z'), Filters.cond('r', 'between', [1, 5]), Filters.cond('Mol Weight', '>', 1),
    Filters.group('or', [Filters.cond('s', '=', 1), Filters.cond('t', '=', 2)]),
    Filters.group('or', [Filters.cond('u', '=', 1)]), Filters.group('or', []),
    Filters.cond('w', '=', {type: 'Core.users', id: 'u1', name: 'Alice'}), Filters.cond('x', '=', '@current'),
    Filters.cond('y', '=', new Date('2026-01-02T00:00:00.000Z')),
  ]);
  assert.deepEqual(Filters.toDomainTree(root), [
    c('a', '=', [1, 2]), 'and', c('b', '!=', ['x']), 'and', c('d', '=', null), 'and', c('e', '!=', null), 'and',
    c('h', 'like', '%x%'), 'and', c('i', 'like', 'x%'), 'and', c('j', 'like', '%x'), 'and', c('k', 'not like', '%x%'), 'and',
    c('l', '~*', 'r'), 'and', c('m', '!~*', 'r'), 'and', c('n', 'like', '%a%b%'), 'and', c('o', 'like', '%50\\%\\_\\\\%'), 'and',
    [{property: 'q', operator: 'fuzzy', threshold: 0.6, value: 'z'}, 'or', c('q', 'like', '%z%')], 'and',
    [{property: 'q2', operator: 'fuzzy', threshold: null, value: 'z'}, 'or', c('q2', 'like', '%z%')], 'and',
    [c('r', '>=', 1), 'and', c('r', '<=', 5)], 'and', c('Mol Weight', '>', 1), 'and',
    [c('s', '=', 1), 'or', c('t', '=', 2)], 'and', c('u', '=', 1), 'and', c('w', '=', 'u1'), 'and', c('x', '=', '@current'),
    'and', c('y', '=', '2026-01-02T00:00:00.000Z'),
  ]);
  assert.deepEqual(Filters.toDomainTree(Filters.group('or')), []);
  assert.deepEqual(Filters.toDomainTree(Filters.group('or', [Filters.cond('a', '=', 1), Filters.cond('b', '=', 2)])),
    [c('a', '=', 1), 'or', c('b', '=', 2)]);
});

test('toDomainTree: not pushes down, nested and on a between; an operator without a domain form throws', () => {
  const root = Filters.group('and', [
    Filters.cond('a', '=', 1),
    Filters.group('or', [Filters.cond('b', '>', 2), Filters.cond('c', 'like', 'x'), Filters.cond('d', 'between', [1, 5]),
      Filters.group('and', [Filters.cond('e', 'in', [1]), Filters.cond('f', 'fuzzy', 'q')], {not: true})], {not: true}),
  ]);
  assert.deepEqual(Filters.toDomainTree(root), [c('a', '=', 1), 'and', [
    c('b', '<=', 2), 'and', c('c', 'not like', '%x%'), 'and', [c('d', '<', 1), 'or', c('d', '>', 5)], 'and',
    [c('e', '=', [1]), 'and', c('f', 'like', '%q%')],
  ]], 'a doubly negated fuzzy pair stays the plain like (the pushdown approximation)');
  assert.deepEqual(Filters.toDomainTree(Filters.group('and', [Filters.cond('a', '=', 1)], {not: true})), [c('a', '!=', 1)]);
  const unregister = Filters.operators.register({id: 'Contains', label: 'Contains', arity: 1, kinds: [KIND.STRING],
    semType: SEMTYPE.MOLECULE, editor: 'default', bitset: async () => ({bits: new Uint32Array(0), length: 0})});
  try {
    assert.throws(() => Filters.toDomainTree(Filters.group('and', [Filters.cond('smiles', 'Contains', 'C')])),
      /"Contains" has no domain form/);
    assert.throws(() => Filters.toDomainTree(Filters.group('and', [Filters.cond('a', 'nope', 1)])), /"nope" has no domain form/);
  } finally {
    unregister();
  }
});

test('the hidden span: invisible to JSON and deep equality, printed by format, resolved at call time', () => {
  const now = new Date('2026-09-04T12:00:00.000Z');
  const {tree} = Filters.parseTree('created > -1w and updated < 2d', {now});
  assert.equal(JSON.stringify(tree), JSON.stringify(JSON.parse(JSON.stringify(tree))));
  assert.ok(!JSON.stringify(tree).includes('span'));
  assert.deepEqual(tree[0].value, new Date('2026-08-28T12:00:00.000Z'), 'deepStrictEqual against a plain Date');
  const root = Filters.fromDomainTree(tree);
  assert.deepEqual(root.nodes.map((n) => n.value), [{span: '-1w'}, {span: '-2d'}]);
  assert.equal(Filters.format(root), 'created > -1w and updated > -2d');
  const later = new Date('2026-12-01T00:00:00.000Z');
  assert.deepEqual(Filters.toDomainTree(root, {now: later}), [
    c('created', '>', '2026-11-24T00:00:00.000Z'), 'and', c('updated', '>', '2026-11-29T00:00:00.000Z')]);
  assert.deepEqual(Filters.toDomainTree(Filters.group('and', [Filters.cond('c', '<=', {span: 'now'})]), {now: later}),
    [c('c', '<=', '2026-12-01T00:00:00.000Z')]);
  const plain = Filters.fromDomainTree([c('c', '=', new Date('2026-01-01T00:00:00.000Z'))]);
  assert.ok(plain.nodes[0].value instanceof Date, 'an unmarked Date stays a Date');
});

/** The model with its spans resolved against `now` — what a datetime-typed schema reads back from
 * the ISO strings the tree carries. */
function resolved(node, now) {
  if (Filters.isGroup(node))
    return {...node, nodes: node.nodes.map((n) => resolved(n, now))};
  const date = (v) => Filters.isSpan(v) ? Filters.resolveSpan(v.span, now) : v;
  return node.value === undefined ? node : {...node, value: Array.isArray(node.value) ? node.value.map(date) : date(node.value)};
}

const datetimes = Filters.schema([{name: 'created', type: TYPE.DATE_TIME}, {name: 'updated', type: TYPE.DATE_TIME}]);

for (const entry of corpus.entries.filter((e) => e.errors !== true)) {
  test(`round trip ${entry.id}: the model survives toDomainTree → fromDomainTree`, () => {
    const now = new Date();
    const model = Filters.fromDomainTree(Filters.parseTree(entry.input, {now}).tree, datetimes);
    const again = Filters.fromDomainTree(Filters.toDomainTree(model, {now}), datetimes);
    const expected = resolved(model, now);
    assert.ok(Filters.equals(expected, again), `${Filters.format(expected)} vs ${Filters.format(again)}`);
    assert.equal(Filters.format(again), Filters.format(expected));
  });
}

test('Filters.parse: syntax problems with an empty root, validation problems with a schema, typed values', () => {
  const schema = Filters.schema([{name: 'name', type: TYPE.STRING}, {name: 'age', type: TYPE.INT},
    {name: 'created', type: TYPE.DATE_TIME}, {name: 'owner', type: TYPE.STRING, ref: 'Core.users'}]);
  const bad = Filters.parse('name = ', schema);
  assert.deepEqual(bad.root.nodes, []);
  assert.equal(bad.problems[0].code, 'syntax');
  assert.deepEqual(bad.problems[0].position, {start: 7, end: 7});
  const ok = Filters.parse('age > 30 and owner = "u1" and created > "2026-01-01T00:00:00.000Z" and created < -1w', schema);
  assert.deepEqual(ok.problems, []);
  assert.deepEqual(ok.root.nodes[1].value, {type: 'Core.users', id: 'u1'});
  assert.ok(ok.root.nodes[2].value instanceof Date);
  assert.deepEqual(ok.root.nodes[3].value, {span: '-1w'});
  const problems = Filters.parse('age like "x" and nope = 1 and name fuzzy "q"', schema, 'dataframe').problems;
  assert.deepEqual(problems.map((p) => p.code), ['operator-not-applicable', 'unknown-property', 'not-expressible']);
  assert.deepEqual(Filters.parse('age > 30').problems, [], 'no schema, no validation');
  assert.equal(Filters.parse('').root.nodes.length, 0);
  const now = new Date('2026-09-04T12:00:00.000Z');
  assert.equal(Filters.toDomainTree(Filters.parse('created < 1h', schema, undefined, {now}).root, {now})[0].value,
    '2026-09-04T11:00:00.000Z', 'parse forwards now');
});

test('Filters.parse: a property typed in another case takes the schema spelling, so format prints it', () => {
  const schema = Filters.schema([{name: 'AGE', type: TYPE.INT}, {name: 'Owner', type: TYPE.STRING, ref: 'Core.users'}]);
  const parsed = Filters.parse('age > 1 and OWNER.login = "u1" and owner = "u2"', schema);
  assert.deepEqual(parsed.problems, []);
  assert.deepEqual(parsed.root.nodes.map((n) => n.property), ['AGE', 'Owner.login', 'Owner']);
  assert.deepEqual(parsed.root.nodes[2].value, {type: 'Core.users', id: 'u2'}, 'the canonical property types the value');
  assert.equal(Filters.format(parsed.root), 'AGE > 1 and Owner.login = "u1" and Owner = "u2"');
  assert.equal(Filters.parse('Age > 1', Filters.schema([{name: 'Age', type: TYPE.INT}, {name: 'age', type: TYPE.INT}]))
    .root.nodes[0].property, 'Age', 'an exact name stays');
  assert.equal(Filters.parse('AGE > 1', Filters.schema([{name: 'Age', type: TYPE.INT}, {name: 'age', type: TYPE.INT}]))
    .problems[0].code, 'unknown-property', 'two case variants: neither is meant');
  assert.equal(Filters.fromDomainTree([{property: 'age', operator: '>', value: 1}], schema).nodes[0].property, 'AGE',
    'a tree written by hand canonicalizes the same way');
});

test('a literal % or _ survives model → tree → text → tree → model', () => {
  for (const cond of [Filters.cond('o', 'like', '50%'), Filters.cond('o', 'starts', 'a_b'), Filters.cond('o', 'ends', '%'),
    Filters.cond('o', '!like', '50%_\\')]) {
    const text = Filters.format(Filters.fromDomainTree(Filters.toDomainTree(Filters.group('and', [cond]))));
    const again = Filters.fromDomainTree(Filters.parseTree(text).tree).nodes[0];
    assert.equal(again.operator, cond.operator, text);
    assert.equal(again.value, cond.value, text);
    assert.equal(again.options, undefined, text);
  }
});
