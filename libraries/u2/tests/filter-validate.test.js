import {test} from 'node:test';
import assert from 'node:assert/strict';
import {Filters, KIND} from '../src/core/filter/index.js';
import {TYPE, SEMTYPE} from 'datagrok-api/u2core';

const schema = Filters.schema([
  {name: 'name', type: TYPE.STRING},
  {name: 'status', type: TYPE.STRING, choices: ['Open', 'Closed']},
  {name: 'age', type: TYPE.INT, min: 0, max: 150},
  {name: 'mw', type: TYPE.FLOAT},
  {name: 'big', type: TYPE.BIG_INT},
  {name: 'created', type: TYPE.DATE_TIME},
  {name: 'active', type: TYPE.BOOL},
  {name: 'author', type: TYPE.STRING, ref: 'Core.users'},
  {name: 'tags', type: TYPE.LIST},
]);

function problems(nodes, target, template) {
  return Filters.validate(Filters.group('and', nodes), schema, target, template);
}

function codes(nodes, target) {
  return problems(nodes, target).map((p) => p.code);
}

test('a valid flat tree and an empty group have no problems', () => {
  assert.deepEqual(problems([]), []);
  assert.deepEqual(problems([
    Filters.cond('name', 'like', 'a'), Filters.cond('age', 'between', [1, 2]), Filters.cond('mw', '>', 1.5),
    Filters.cond('big', '=', '123456789012345678901'), Filters.cond('big', '<', 5),
    Filters.cond('created', '>', {span: '-1w'}), Filters.cond('created', '<', new Date()),
    Filters.cond('created', '=', '2026-01-01T00:00:00Z'), Filters.cond('active', '=', true),
    Filters.cond('author', '=', {type: 'User', id: '1'}), Filters.cond('author', 'in', ['1', '@current']),
    Filters.cond('tags', '!like', 'x'), Filters.cond('status', 'in', ['Open']), Filters.cond('name', 'is null'),
    Filters.cond('name', '=', 5),
  ]), []);
});

test('unknown-property', () => {
  Filters.resetIds('');
  const [p] = problems([Filters.cond('nope', '=', 1)]);
  assert.equal(p.code, 'unknown-property');
  assert.equal(p.nodeId, 'f1');
  assert.match(p.message, /nope/);
});

test('operator-not-applicable: unknown ids and ids of another kind', () => {
  assert.deepEqual(codes([Filters.cond('age', 'like', 'x')]), ['operator-not-applicable']);
  assert.deepEqual(codes([Filters.cond('name', '>', 'x')]), ['operator-not-applicable']);
  assert.deepEqual(codes([Filters.cond('name', 'zzz', 'x')]), ['operator-not-applicable']);
  assert.deepEqual(codes([Filters.cond('tags', '=', 'x')]), ['operator-not-applicable']);
});

test('missing-value: arity 1 with undefined/null/"", arity 2 without a pair, n with an empty list', () => {
  assert.deepEqual(codes([Filters.cond('name', '=')]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('name', '=', null)]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('name', '=', '')]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('age', 'between', [1])]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('age', 'between', 1)]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('name', 'in', [])]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('name', 'in', 'x')]), ['missing-value']);
  assert.deepEqual(codes([Filters.cond('name', 'is null', 'ignored')]), [], 'arity 0 takes no value');
});

test('invalid-value: per kind, list elements individually, choices, min/max', () => {
  assert.deepEqual(codes([Filters.cond('age', '=', 1.5)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('age', '=', '5')]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('mw', '=', 'x')]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('mw', '=', Infinity)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('big', '=', '12a')]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('big', '=', 1.5)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('created', '=', 'nope')]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('created', '=', {span: '1x'})]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('created', '=', new Date('x'))]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('active', '=', 'true')]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('author', '=', 5)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('name', '=', true)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('age', 'in', [1, 2.5])]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('age', 'between', [1, 'x'])]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('status', '=', 'Nope')]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('status', 'in', ['Open', 'Nope'])]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('status', 'like', 'Nope')]), [], 'choices bind = and in only');
  assert.deepEqual(codes([Filters.cond('age', '>', -1)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('age', '<', 151)]), ['invalid-value']);
  assert.deepEqual(codes([Filters.cond('age', 'between', [0, 150])]), []);
  const [p] = problems([Filters.cond('age', '=', 200)]);
  assert.match(p.message, /at most 150/);
});

test('dotted paths: the head must exist; a ref hop skips the operator and value checks', () => {
  assert.deepEqual(codes([Filters.cond('author.login', 'like', 5)]), []);
  assert.deepEqual(codes([Filters.cond('name.x', '>', 'y')]), []);
  assert.deepEqual(codes([Filters.cond('nope.login', '=', 'x')]), ['unknown-property']);
  assert.deepEqual(codes([Filters.cond('author.login', '=', 'x')], 'domain'), []);
});

test('target domain: every core operator has a form, a bitset-only semType operator does not', () => {
  const off = Filters.operators.register({id: 'Contains', label: 'Contains', arity: 1, kinds: [KIND.STRING],
    semType: SEMTYPE.MOLECULE, editor: 'default', bitset: async () => ({bits: new Uint32Array(0), length: 0})});
  const mol = Filters.schema([{name: 'smiles', type: TYPE.STRING, semType: SEMTYPE.MOLECULE}]);
  try {
    Filters.resetIds('');
    const root = Filters.group('and', [Filters.cond('smiles', 'Contains', 'c1ccccc1'), Filters.cond('smiles', 'fuzzy', 'x')]);
    assert.deepEqual(Filters.validate(root, mol), []);
    assert.deepEqual(Filters.validate(root, mol, 'domain').map((p) => [p.nodeId, p.code]), [['f1', 'not-expressible']]);
    assert.deepEqual(Filters.validate(root, mol, 'dataframe').map((p) => [p.nodeId, p.code]), [['f2', 'not-expressible']]);
  } finally {
    off();
  }
  assert.deepEqual(codes([Filters.cond('name', 'fuzzy', 'x'), Filters.cond('created', '>', {span: '-1d'})], 'domain'), []);
});

test('target dataframe: every core operator but fuzzy has a mask', () => {
  assert.deepEqual(codes([Filters.cond('name', 'like', 'x'), Filters.cond('age', 'between', [1, 2]),
    Filters.cond('active', '=', true), Filters.cond('created', '>', {span: '-1d'})], 'dataframe'), []);
  assert.deepEqual(codes([Filters.cond('name', 'fuzzy', 'x')], 'dataframe'), ['not-expressible']);
});

test('target dataframe: mask or bitset', () => {
  const off = Filters.operators.register([
    {id: 'm', label: 'm', arity: 1, kinds: [KIND.STRING], semType: 'T', editor: 'default', mask: () => ({bits: new Uint32Array(0), length: 0})},
    {id: 'b', label: 'b', arity: 1, kinds: [KIND.STRING], semType: 'T', editor: 'default', bitset: async () => ({bits: new Uint32Array(0), length: 0})},
    {id: 'd', label: 'd', arity: 1, kinds: [KIND.STRING], semType: 'T', editor: 'default', domain: () => ({property: 'p', operator: '='})},
  ]);
  try {
    const s = Filters.schema([{name: 'p', type: TYPE.STRING, semType: 'T'}]);
    Filters.resetIds('');
    const root = Filters.group('and', [Filters.cond('p', 'm', 'x'), Filters.cond('p', 'b', 'x'), Filters.cond('p', 'd', 'x')]);
    assert.deepEqual(Filters.validate(root, s, 'dataframe').map((p) => p.nodeId), ['f3']);
  } finally {
    off();
  }
});

test('problems come in tree order, one per node, and a template adds its locked problems', () => {
  Filters.resetIds('');
  const template = {root: Filters.group('and', [Filters.cond('name', '=', 'a', {lock: 'all'})])};
  let value = Filters.update(Filters.applyTemplate(template), 'f1', {value: 'b'});
  value = Filters.insert(value, 'f2', 1, Filters.cond('nope', '=', 1));
  value = Filters.insert(value, 'f2', 2, Filters.cond('age', '=', 'x'));
  assert.deepEqual(Filters.validate(value, schema, undefined, template).map((p) => [p.nodeId, p.code]),
    [['f3', 'unknown-property'], ['f4', 'invalid-value'], ['f1', 'locked']]);
});

test('pruneInvalid drops the reported nodes and the groups that leaves empty; the same root when nothing goes', () => {
  Filters.resetIds('');
  const ok1 = Filters.cond('name', 'like', 'a');
  const ok2 = Filters.cond('age', '>', 5);
  const sound = Filters.group('and', [ok1, Filters.group('or', [ok2])]);
  assert.equal(Filters.pruneInvalid(sound, schema), sound, 'nothing to drop: the same root');

  const root = Filters.group('and', [
    ok1,
    Filters.cond('age', '='),
    Filters.group('or', [Filters.cond('nope', '=', 1), Filters.group('and', [Filters.cond('mw', 'between', [1])])]),
    Filters.group('or', [Filters.cond('age', '>', 1000), ok2]),
  ]);
  const pruned = Filters.pruneInvalid(root, schema);
  assert.notEqual(pruned, root);
  assert.equal(pruned.id, root.id);
  assert.deepEqual(pruned.nodes.map((n) => n.id), [ok1.id, root.nodes[3].id], 'the all-invalid group is gone');
  assert.equal(pruned.nodes[0], ok1, 'kept nodes keep their identity');
  assert.deepEqual(pruned.nodes[1].nodes, [ok2]);
  assert.deepEqual(Filters.validate(pruned, schema), []);

  const empty = Filters.pruneInvalid(Filters.group('and', [Filters.cond('age', '=')]), schema);
  assert.deepEqual(empty.nodes, [], 'an emptied root stays a root');

  const target = Filters.group('and', [ok1, Filters.cond('name', 'fuzzy', 'a')]);
  assert.deepEqual(Filters.pruneInvalid(target, schema, 'dataframe').nodes, [ok1], 'the target counts');
  assert.equal(Filters.pruneInvalid(target, schema), target);
});
