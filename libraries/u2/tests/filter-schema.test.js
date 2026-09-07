import {test} from 'node:test';
import assert from 'node:assert/strict';
import {Filters} from '../src/core/filter/index.js';
import {TYPE} from 'datagrok-api/u2core';

const users = Filters.schema([{name: 'login', type: TYPE.STRING}, {name: 'group', type: TYPE.STRING, ref: 'Core.groups'}]);
users.resolveRef = async (prop) => prop.ref === 'Core.groups' ?
  Filters.schema([{name: 'name', type: TYPE.STRING}]) : Promise.reject(new Error('no such ref'));
const issues = Filters.schema([
  {name: 'title', type: TYPE.STRING},
  {name: 'author', type: TYPE.STRING, ref: 'Core.users'},
  {type: TYPE.INT},
], {title: ['Alpha', 'beta', 'Gamma']});
issues.resolveRef = async () => users;

/** and(status in [Active] {value}, created > -30d {all}, or(a = 1, b = 2) {value}) — f1..f6. */
function template(extra = {}) {
  Filters.resetIds('');
  return {
    root: Filters.group('and', [
      Filters.cond('status', 'in', ['Active'], {lock: 'value'}),
      Filters.cond('created', '>', {span: '-30d'}, {lock: 'all'}),
      Filters.group('or', [Filters.cond('a', '=', 1), Filters.cond('b', '=', 2)], {lock: 'value'}),
    ]),
    ...extra,
  };
}

function lockedIds(t, value) {
  const problems = Filters.checkLocks(t, value);
  assert.ok(problems.every((p) => p.code === 'locked'));
  return problems.map((p) => p.nodeId);
}

test('Filters.schema: named properties only; values filtered case-insensitively, absent without a list', async () => {
  assert.deepEqual(issues.properties.map((p) => p.name), ['title', 'author']);
  assert.equal(users.values, undefined);
  assert.deepEqual(await issues.values(issues.properties[0], 'A'),
    [{value: 'Alpha', label: 'Alpha'}, {value: 'beta', label: 'beta'}, {value: 'Gamma', label: 'Gamma'}]);
  assert.deepEqual(await issues.values(issues.properties[0], 'gam'), [{value: 'Gamma', label: 'Gamma'}]);
  assert.deepEqual(await issues.values(issues.properties[1], ''), []);
});

test('property: the head segment of a path', () => {
  assert.equal(Filters.property(issues, 'title').name, 'title');
  assert.equal(Filters.property(issues, 'author.login').name, 'author');
  assert.equal(Filters.property(issues, 'nope'), null);
  assert.equal(Filters.property(issues, 'nope.title'), null);
});

test('property: another case resolves when exactly one property spells the name that way', () => {
  const upper = Filters.schema([{name: 'AGE', type: TYPE.INT}, {name: 'Name', type: TYPE.STRING}]);
  assert.equal(Filters.property(upper, 'age').name, 'AGE');
  assert.equal(Filters.property(upper, 'Age.x').name, 'AGE', 'the head segment only');
  assert.equal(Filters.property(upper, 'NAME').name, 'Name');
  const variants = Filters.schema([{name: 'Age', type: TYPE.INT}, {name: 'age', type: TYPE.INT}]);
  assert.equal(Filters.property(variants, 'age').name, 'age', 'exact wins');
  assert.equal(Filters.property(variants, 'Age').name, 'Age');
  assert.equal(Filters.property(variants, 'AGE'), null, 'two case variants stay unknown');
});

test('lockOf returns on a tree whose nodes never got ids (a spec literal), reporting what it finds', () => {
  const root = {op: 'and', lock: 'value', nodes: [{property: 'a', operator: '=', value: 1},
    {op: 'or', nodes: [{property: 'b', operator: '=', value: 2, lock: 'all'}]}]};
  assert.equal(Filters.lockOf(root, undefined), 'value');
  assert.equal(Filters.lockOf(root, 'nope'), 'none');
  const twins = Filters.group('and', [Filters.cond('a', '=', 1, {lock: 'all'})]);
  twins.nodes[0].id = twins.id;
  assert.equal(Filters.lockOf(twins, twins.id), 'none', 'a repeated id ends the walk too; the root shadows its twin');
});

test('applyTemplate: a copy with the same ids', () => {
  const t = template();
  const value = Filters.applyTemplate(t);
  assert.notEqual(value, t.root);
  assert.deepEqual(value, t.root);
  assert.equal(Filters.equals(value, t.root, false), true);
  assert.deepEqual(Filters.diff(t.root, value), {changed: [], added: [], removed: [], moved: []});
});

test('diff: changed fields by id, top-most additions with their place, top-most removals', () => {
  const t = template();
  let value = Filters.applyTemplate(t);
  value = Filters.update(value, 'f1', {value: ['Active', 'QC']});
  value = Filters.update(value, 'f3', {property: 'c', operator: '!='});
  value = Filters.remove(value, 'f2');
  const added = Filters.cond('x', '=', 1);
  value = Filters.insert(value, value.id, 1, added);
  const group = Filters.group('and', [Filters.cond('y', '=', 2)]);
  value = Filters.insert(value, 'f5', 0, group);
  const d = Filters.diff(t.root, value);
  assert.deepEqual(d.changed, [{id: 'f1', value: ['Active', 'QC']}, {id: 'f3', property: 'c', operator: '!='}]);
  assert.deepEqual(d.added, [{parentId: 'f6', index: 1, node: added}, {parentId: 'f5', index: 0, node: group}]);
  assert.deepEqual(d.removed, ['f2']);
  const gone = Filters.remove(Filters.applyTemplate(t), 'f5');
  assert.deepEqual(Filters.diff(t.root, gone).removed, ['f5'], 'children of a removed group are not listed');
  const swapped = Filters.update(Filters.applyTemplate(t), 'f1', {value: ['Active']});
  assert.deepEqual(Filters.diff(t.root, swapped).changed, [], 'an equal value is no change');
  const optioned = Filters.update(Filters.applyTemplate(t), 'f1', {options: {raw: true}});
  assert.deepEqual(Filters.diff(t.root, optioned).changed, [{id: 'f1', options: {raw: true}}]);
  const cleared = Filters.update(Filters.applyTemplate(t), 'f1', {options: undefined});
  assert.deepEqual(Filters.diff(t.root, cleared).changed, [], 'absent and empty options are the same');
  assert.deepEqual(d.moved, [], 'neither an insertion nor a removal beside a node moves it');
});

test('diff: moved is another parent or another order among the siblings both trees hold', () => {
  const t = template();
  const fresh = () => Filters.applyTemplate(t);
  assert.deepEqual(Filters.diff(t.root, Filters.move(fresh(), 'f2', 'f6', 0)).moved, [
    {id: 'f2', from: {parentId: 'f6', index: 1}, to: {parentId: 'f6', index: 0}},
    {id: 'f1', from: {parentId: 'f6', index: 0}, to: {parentId: 'f6', index: 1}},
  ]);
  const out = Filters.diff(t.root, Filters.move(fresh(), 'f3', 'f6', 2));
  assert.deepEqual(out.moved, [{id: 'f3', from: {parentId: 'f5', index: 0}, to: {parentId: 'f6', index: 2}}],
    'the siblings it left and joined keep their order');
  assert.deepEqual([out.added, out.removed], [[], []]);
  assert.deepEqual(Filters.diff(t.root, Filters.insert(fresh(), 'f6', 0, Filters.cond('x', '=', 1))).moved, []);
  assert.deepEqual(Filters.diff(t.root, Filters.remove(fresh(), 'f1')).moved, []);
});

test('checkLocks: a move out of or into a locked group is locked; allowAdd false pins the parent, ' +
  'not the order', () => {
  const t = template();
  const fresh = () => Filters.applyTemplate(t);
  assert.deepEqual(lockedIds(t, Filters.move(fresh(), 'f4', 'f6', 0)), ['f4'], 'out of the value-locked group');
  assert.deepEqual(lockedIds(t, Filters.move(fresh(), 'f1', 'f5', 0)), ['f1'], 'a locked condition stays');
  Filters.resetIds('');
  const open = Filters.group('and', [Filters.cond('a', '=', 1), Filters.cond('b', '=', 2),
    Filters.group('or', [Filters.cond('c', '=', 3)])]);
  const pinned = {root: open, allowAdd: false};
  assert.deepEqual(lockedIds(pinned, Filters.move(Filters.applyTemplate(pinned), 'f2', 'f5', 0)), []);
  assert.deepEqual(lockedIds(pinned, Filters.move(Filters.applyTemplate(pinned), 'f1', 'f4', 0)), ['f1']);
  assert.match(Filters.checkLocks(pinned, Filters.move(Filters.applyTemplate(pinned), 'f1', 'f4', 0))[0].message,
    /moved/);
  assert.deepEqual(lockedIds({root: open}, Filters.move(Filters.applyTemplate(pinned), 'f1', 'f4', 0)), []);
  const into = Filters.move(Filters.applyTemplate(t), 'f2', 'f5', 0);
  assert.deepEqual(lockedIds(t, into), ['f2'], 'into the locked group: the mover is locked, and so is the target');
});

test('checkLocks: value lock keeps property/operator and removal; all lock keeps everything', () => {
  const t = template();
  const fresh = () => Filters.applyTemplate(t);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f1', {value: ['QC']})), []);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f1', {operator: 'not in'})), ['f1']);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f1', {property: 'owner'})), ['f1']);
  assert.deepEqual(lockedIds(t, Filters.remove(fresh(), 'f1')), ['f1']);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f2', {value: {span: '-1d'}})), ['f2']);
  assert.deepEqual(lockedIds(t, Filters.remove(fresh(), 'f2')), ['f2']);
  assert.match(Filters.checkLocks(t, Filters.update(fresh(), 'f2', {value: 1}))[0].message, /locked/);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f1', {options: {raw: true}})), ['f1'],
    'options are not the value');
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f2', {options: {raw: true}})), ['f2']);
});

test('checkLocks: a group lock covers its subtree, its connector and additions inside it', () => {
  const t = template();
  const fresh = () => Filters.applyTemplate(t);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f3', {value: 5})), []);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f3', {operator: '!='})), ['f3']);
  assert.deepEqual(lockedIds(t, Filters.remove(fresh(), 'f4')), ['f4']);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f5', {op: 'and'})), ['f5']);
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f5', {not: true})), ['f5']);
  assert.deepEqual(lockedIds(t, Filters.remove(fresh(), 'f5')), ['f5']);
  const inside = Filters.cond('z', '=', 1);
  assert.deepEqual(lockedIds(t, Filters.insert(fresh(), 'f5', 0, inside)), [inside.id]);
  const outside = Filters.cond('z', '=', 1);
  assert.deepEqual(lockedIds(t, Filters.insert(fresh(), 'f6', 0, outside)), [], 'the unlocked root takes additions');
  assert.deepEqual(lockedIds(t, Filters.update(fresh(), 'f6', {op: 'or'})), [], 'and the connector change');
});

test('checkLocks: allowAdd, allowedProperties, allowedOperators, allowAdvanced', () => {
  const t = template({allowAdd: false, allowedProperties: ['status', 'created', 'a', 'b', 'owner'],
    allowedOperators: {a: ['=', '!='], owner: ['=']}});
  const fresh = () => Filters.applyTemplate(t);
  const added = Filters.cond('owner', '=', 'me');
  assert.deepEqual(lockedIds(t, Filters.insert(fresh(), 'f6', 0, added)), [added.id]);
  const t2 = {...t, allowAdd: true};
  assert.deepEqual(lockedIds(t2, Filters.insert(fresh(), 'f6', 0, added)), []);
  const wrongProp = Filters.cond('secret', '=', 'x');
  assert.deepEqual(lockedIds(t2, Filters.insert(fresh(), 'f6', 0, wrongProp)), [wrongProp.id]);
  const wrongOp = Filters.cond('owner', 'like', 'x');
  assert.deepEqual(lockedIds(t2, Filters.insert(fresh(), 'f6', 0, wrongOp)), [wrongOp.id]);
  const dotted = Filters.cond('owner.login', '=', 'x');
  assert.deepEqual(lockedIds(t2, Filters.insert(fresh(), 'f6', 0, dotted)), [], 'the head segment is what is allowed');
  assert.deepEqual(lockedIds(t2, Filters.update(fresh(), 'f3', {value: 2})), []);
  assert.deepEqual(lockedIds(t2, Filters.update(fresh(), 'f3', {operator: '>'})), ['f3']);
  Filters.resetIds('');
  const flat = {root: Filters.group('and', [Filters.cond('a', '=', 1)]), allowAdvanced: false};
  const nested = Filters.insert(Filters.applyTemplate(flat), 'f2', 0,
    Filters.group('or', [Filters.cond('b', '=', 'x')]));
  assert.deepEqual(lockedIds(flat, nested), ['f2'], 'allowAdvanced: false refuses nesting at the root');
  assert.deepEqual(lockedIds({...flat, allowAdvanced: true}, nested), []);
  assert.deepEqual(lockedIds(flat, Filters.update(Filters.applyTemplate(flat), 'f1', {value: 2})), []);
});
