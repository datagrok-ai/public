import {test} from 'node:test';
import assert from 'node:assert/strict';
import {Filters} from '../src/core/filter/index.js';

/** and(a > 1, or(b = "x", c = true), d < 5) with predictable ids f1..f6. */
function tree() {
  Filters.resetIds('');
  const a = Filters.cond('a', '>', 1);
  const b = Filters.cond('b', '=', 'x');
  const c = Filters.cond('c', '=', true);
  const inner = Filters.group('or', [b, c]);
  const d = Filters.cond('d', '<', 5);
  const root = Filters.group('and', [a, inner, d]);
  return {root, a, b, c, inner, d};
}

test('newId/resetIds: predictable ids after a reset, distinct otherwise', () => {
  Filters.resetIds('t');
  assert.equal(Filters.newId(), 'ft1');
  assert.equal(Filters.newId(), 'ft2');
  Filters.resetIds();
  assert.equal(Filters.newId(), 'f1');
});

test('group/cond: only the given keys are present', () => {
  Filters.resetIds('');
  assert.deepEqual(Filters.group(), {id: 'f1', op: 'and', nodes: []});
  assert.deepEqual(Filters.cond('a', 'is null'), {id: 'f2', property: 'a', operator: 'is null'});
  assert.deepEqual(Filters.cond('a', '=', 1, {lock: 'value', options: {raw: true}}),
    {id: 'f3', property: 'a', operator: '=', value: 1, lock: 'value', options: {raw: true}});
  assert.deepEqual(Filters.group('or', [], {not: true, lock: 'all'}), {id: 'f4', op: 'or', nodes: [], not: true, lock: 'all'});
  assert.equal(Filters.isGroup(Filters.group()), true);
  assert.equal(Filters.isGroup(Filters.cond('a', '=', 1)), false);
});

test('find/parentOf/walk/count/depth/isFlat', () => {
  const {root, a, b, inner} = tree();
  assert.equal(Filters.find(root, b.id), b);
  assert.equal(Filters.find(root, root.id), root);
  assert.equal(Filters.find(root, 'nope'), null);
  assert.equal(Filters.parentOf(root, b.id), inner);
  assert.equal(Filters.parentOf(root, a.id), root);
  assert.equal(Filters.parentOf(root, root.id), null);
  const seen = [];
  Filters.walk(root, (n, parent, i) => seen.push([n.id, parent?.id ?? null, i]));
  assert.deepEqual(seen, [['f6', null, -1], ['f1', 'f6', 0], ['f4', 'f6', 1], ['f2', 'f4', 0], ['f3', 'f4', 1],
    ['f5', 'f6', 2]]);
  assert.equal(Filters.count(root), 4);
  assert.equal(Filters.isFlat(root), false);
  assert.equal(Filters.isFlat(inner), true);
  assert.equal(Filters.isFlat(Filters.group('and', [], {not: true})), false);
});

test('update: a new root, the untouched sibling subtree keeps identity, root patches apply to the root', () => {
  const {root, a, b, inner, d} = tree();
  const next = Filters.update(root, b.id, {value: 'y'});
  assert.notEqual(next, root);
  assert.equal(Filters.find(next, b.id).value, 'y');
  assert.equal(b.value, 'x', 'the old node is untouched');
  assert.equal(next.nodes[0], a, 'sibling condition keeps identity');
  assert.equal(next.nodes[2], d);
  assert.notEqual(next.nodes[1], inner, 'the path down to the node is copied');
  assert.equal(next.nodes[1].nodes[1], Filters.find(root, 'f3'), 'the sibling inside the group keeps identity');
  assert.equal(Filters.update(root, 'nope', {value: 1}), root, 'unknown id keeps the root');
  const rooted = Filters.update(root, root.id, {op: 'or'});
  assert.equal(rooted.op, 'or');
  assert.equal(rooted.nodes, root.nodes);
});

test('insert: at an index, clamped, into nested groups; unknown parent keeps the root', () => {
  const {root, a, inner, d} = tree();
  const n = Filters.cond('e', '=', 1);
  const next = Filters.insert(root, root.id, 1, n);
  assert.deepEqual(next.nodes.map((x) => x.id), [a.id, n.id, inner.id, d.id]);
  assert.equal(next.nodes[2], inner);
  assert.equal(Filters.insert(root, root.id, 99, n).nodes[3], n, 'index clamps to the end');
  assert.equal(Filters.insert(root, root.id, -5, n).nodes[0], n, 'index clamps to the start');
  const deep = Filters.insert(root, inner.id, 0, n);
  assert.equal(deep.nodes[1].nodes[0], n);
  assert.equal(deep.nodes[0], a);
  assert.equal(Filters.insert(root, 'nope', 0, n), root);
  assert.equal(Filters.insert(root, a.id, 0, n), root, 'a condition is not a parent');
});

test('remove: drops the node, keeps siblings and the other subtree; the root cannot go', () => {
  const {root, a, b, c, inner, d} = tree();
  const next = Filters.remove(root, b.id);
  assert.notEqual(next, root);
  assert.deepEqual(next.nodes[1].nodes, [c]);
  assert.equal(next.nodes[0], a);
  assert.equal(next.nodes[2], d);
  assert.equal(Filters.remove(root, inner.id).nodes.length, 2);
  assert.equal(Filters.remove(root, 'nope'), root);
  assert.equal(Filters.remove(root, root.id), root);
});

test('move: the index is where the node lands after it is detached (designer patch semantics)', () => {
  const {root, a, inner, d} = tree();
  assert.deepEqual(Filters.move(root, a.id, root.id, 2).nodes.map((x) => x.id), [inner.id, d.id, a.id]);
  assert.deepEqual(Filters.move(root, a.id, root.id, 1).nodes.map((x) => x.id), [inner.id, a.id, d.id]);
  assert.deepEqual(Filters.move(root, d.id, root.id, 0).nodes.map((x) => x.id), [d.id, a.id, inner.id]);
  assert.equal(Filters.move(root, a.id, root.id, 0), root, 'landing where it is keeps the root');
  assert.equal(Filters.move(root, d.id, root.id, 2), root);
  assert.equal(Filters.move(root, d.id, root.id, 99).nodes[2], d, 'past the end clamps — same place, same root');
  const across = Filters.move(root, d.id, inner.id, 1);
  assert.deepEqual(across.nodes.map((x) => x.id), [a.id, inner.id]);
  assert.deepEqual(across.nodes[1].nodes.map((x) => x.id), ['f2', d.id, 'f3']);
  assert.equal(across.nodes[0], a);
  assert.equal(Filters.move(root, inner.id, inner.id, 0), root, 'not into itself');
  assert.equal(Filters.move(root, root.id, inner.id, 0), root, 'the root does not move');
  assert.equal(Filters.move(root, a.id, 'nope', 0), root);
  assert.equal(Filters.move(root, a.id, d.id, 0), root, 'a condition is not a parent');
});

test('replace: swaps a node; the root only for a group', () => {
  const {root, a, b, inner} = tree();
  const g = Filters.group('or', [Filters.cond('z', '=', 1)]);
  const next = Filters.replace(root, b.id, g);
  assert.equal(next.nodes[1].nodes[0], g);
  assert.equal(next.nodes[0], a);
  assert.equal(Filters.replace(root, 'nope', g), root);
  assert.equal(Filters.replace(root, root.id, g), g);
  assert.equal(Filters.replace(root, root.id, a), root);
  assert.equal(Filters.find(Filters.replace(root, inner.id, a), inner.id), null);
});

test('flatten: same-op nesting and single-node groups inline; mixed ops or not answer null', () => {
  Filters.resetIds('');
  const a = Filters.cond('a', '=', 1);
  const b = Filters.cond('b', '=', 2);
  const c = Filters.cond('c', '=', 3);
  const flat = Filters.group('and', [a, b]);
  assert.equal(Filters.flatten(flat), flat, 'already flat keeps identity');
  const same = Filters.group('and', [a, Filters.group('and', [b, Filters.group('and', [c])])]);
  assert.deepEqual(Filters.flatten(same).nodes, [a, b, c]);
  assert.equal(Filters.flatten(same).id, same.id);
  const single = Filters.group('and', [a, Filters.group('or', [b])]);
  assert.deepEqual(Filters.flatten(single).nodes, [a, b]);
  assert.deepEqual(Filters.flatten(Filters.group('and', [a, Filters.group('or', [])])).nodes, [a]);
  assert.equal(Filters.flatten(Filters.group('and', [a, Filters.group('or', [b, c])])), null);
  assert.equal(Filters.flatten(Filters.group('and', [a, Filters.group('and', [b], {not: true})])), null);
  assert.equal(Filters.flatten(Filters.group('and', [a, b], {not: true})), null);
});

test('clone: a deep copy, ids kept or fresh, value lists and options not shared', () => {
  const {root, a, b, inner} = tree();
  const listed = Filters.update(root, a.id, {value: [1, 2], options: {threshold: 0.5}});
  const copy = Filters.clone(listed);
  assert.notEqual(copy, listed);
  assert.deepEqual(copy, listed);
  assert.equal(Filters.equals(copy, listed, false), true);
  assert.notEqual(Filters.find(copy, a.id).value, Filters.find(listed, a.id).value);
  assert.notEqual(Filters.find(copy, a.id).options, Filters.find(listed, a.id).options);
  assert.notEqual(Filters.find(copy, inner.id), inner);
  const fresh = Filters.clone(root, true);
  assert.equal(Filters.equals(fresh, root), true);
  assert.equal(Filters.equals(fresh, root, false), false);
  assert.equal(Filters.find(fresh, b.id), null);
  assert.equal(Filters.count(fresh), 4);
});

test('equals: structural, ids and ref display names are noise; dates, spans, lists and options compare by value', () => {
  Filters.resetIds('');
  const t = 1700000000000;
  const x = Filters.group('and', [
    Filters.cond('d', '>', new Date(t)),
    Filters.cond('u', '=', {type: 'User', id: '42', name: 'Alice'}),
    Filters.cond('c', '>', {span: '-1w'}),
    Filters.cond('s', 'in', ['a', 'b']),
    Filters.cond('f', 'fuzzy', 'q', {options: {threshold: 0.6}}),
  ]);
  const y = Filters.group('and', [
    Filters.cond('d', '>', new Date(t)),
    Filters.cond('u', '=', {type: 'User', id: '42'}),
    Filters.cond('c', '>', {span: '-1w'}),
    Filters.cond('s', 'in', ['a', 'b']),
    Filters.cond('f', 'fuzzy', 'q', {options: {threshold: 0.6}}),
  ]);
  assert.equal(Filters.equals(x, y), true);
  assert.equal(Filters.equals(x, y, false), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.nodes[0].id, {value: new Date(t + 1)})), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.nodes[1].id, {value: {type: 'User', id: '43'}})), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.nodes[2].id, {value: {span: '-2w'}})), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.nodes[3].id, {value: ['a']})), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.nodes[4].id, {options: {threshold: 0.7}})), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.id, {op: 'or'})), false);
  assert.equal(Filters.equals(x, Filters.update(y, y.id, {not: true})), false);
  assert.equal(Filters.equals(Filters.group('and', [], {not: false}), Filters.group('and')), true);
  assert.equal(Filters.equals(x.nodes[0], x), false);
});

test('toJson/fromJson: a lossless round trip — semType operator ids, dates, spans, refs, not, locks, options', () => {
  Filters.resetIds('');
  const t = 1700000000000;
  const root = Filters.group('and', [
    Filters.cond('canonical_smiles', 'Contains', 'c1ccccc1'),
    Filters.cond('d', '>', new Date(t)),
    Filters.cond('c', '>', {span: '-1w'}),
    Filters.cond('u', '=', {type: 'Core.users', id: '42', name: 'Alice'}, {lock: 'value'}),
    Filters.group('or', [Filters.cond('s', 'in', ['a', 'b']), Filters.cond('f', 'fuzzy', 'q', {options: {threshold: 0.6}})],
      {not: true, lock: 'all'}),
  ], {lock: 'value'});
  const json = Filters.toJson(root);
  assert.deepEqual(json.nodes[1].value, {date: new Date(t).toISOString()});
  assert.deepEqual(json.nodes[2].value, {span: '-1w'});
  assert.equal(json.nodes[0].operator, 'Contains');
  const back = Filters.fromJson(JSON.parse(JSON.stringify(json)));
  assert.deepEqual(back, root, 'ids and every key come back');
  assert.ok(back.nodes[1].value instanceof Date);
  assert.equal(Filters.equals(back, root, false), true);
});

test('fromJson: missing ids are made up, unknown keys dropped, a Dart map or a bad shape throws', () => {
  Filters.resetIds('');
  const back = Filters.fromJson({op: 'or', nodes: [{property: 'a', operator: '=', value: 1, junk: 1}], junk: 2});
  assert.deepEqual(back, {id: 'f1', op: 'or', nodes: [{id: 'f2', property: 'a', operator: '=', value: 1}]});
  assert.throws(() => Filters.fromJson({property: 'a', operator: '='}), /root is not a group/);
  assert.throws(() => Filters.fromJson({op: 'and', nodes: [{property: 'a'}]}), /no property or operator/);
  assert.throws(() => Filters.fromJson({op: 'xor', nodes: []}), /not and\/or/);
  assert.throws(() => Filters.fromJson({op: 'and', nodes: [], lock: 'some'}), /not none, value or all/);
  assert.throws(() => Filters.fromJson({op: 'and', nodes: [{property: 'a', operator: '=', value: {x: 1}}]}),
    /not a scalar/);
  assert.throws(() => Filters.fromJson(new Map()), /not a plain object/);
  assert.throws(() => Filters.fromJson('a = 1'), /not a plain object/);
});
