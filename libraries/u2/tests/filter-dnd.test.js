/* Where a dragged filter node may land: the pure half of the builder's drag-and-drop. The shim
   lays nothing out, so every test feeds the rects itself, the way tests/dnd.test.js does. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {Filters} from '../src/core/filter/index.js';
import {resolveFilterDrop} from '../src/components/filter/filter-dnd.js';

/** root(and)[a, g(or)[b, c], d, e(or)[]] laid out top to bottom; the group's rows are inset. */
function fixture() {
  Filters.resetIds('');
  const a = Filters.cond('age', '>', 30);
  const b = Filters.cond('sex', '=', 'F');
  const c = Filters.cond('name', 'like', 'a');
  const g = Filters.group('or', [b, c]);
  const d = Filters.cond('mw', '<', 500);
  const e = Filters.group('or');
  const root = Filters.group('and', [a, g, d, e]);
  const rect = (x, y, width, height) => ({x, y, width, height});
  const hits = [
    {id: root.id, rect: rect(0, 0, 200, 200)},
    {id: a.id, rect: rect(0, 0, 200, 40)},
    {id: g.id, rect: rect(0, 40, 200, 80)},
    {id: b.id, rect: rect(10, 50, 190, 30)},
    {id: c.id, rect: rect(10, 90, 190, 30)},
    {id: d.id, rect: rect(0, 120, 200, 40)},
    {id: e.id, rect: rect(0, 160, 200, 30)},
  ];
  return {root, a, b, c, d, e, g, hits};
}

test('a group under the pointer takes the drop into itself, after what it holds; an empty one at 0', () => {
  const {root, d, e, g, hits} = fixture();
  assert.deepEqual(resolveFilterDrop(root, hits, 100, 45, d.id),
    {parentId: g.id, index: 2, kind: 'into', rect: {x: 0, y: 40, width: 200, height: 80}});
  assert.deepEqual(resolveFilterDrop(root, hits, 100, 170, d.id),
    {parentId: e.id, index: 0, kind: 'into', rect: {x: 0, y: 160, width: 200, height: 30}});
});

test('a condition splits its parent at the midpoint: a line before or after it', () => {
  const {root, a, b, d, g, hits} = fixture();
  const before = resolveFilterDrop(root, hits, 100, 125, a.id);
  assert.equal(before.kind, 'line');
  assert.equal(before.parentId, root.id);
  assert.equal(before.index, 1, 'before d, counted with a already gone');
  assert.deepEqual(before.rect, {x: 0, y: 120, width: 200, height: 2});
  const after = resolveFilterDrop(root, hits, 100, 155, a.id);
  assert.equal(after.index, 2, 'after d');
  assert.deepEqual(after.rect, {x: 0, y: 160, width: 200, height: 2});
  assert.deepEqual(resolveFilterDrop(root, hits, 100, 60, d.id),
    {parentId: g.id, index: 0, kind: 'line', rect: {x: 10, y: 50, width: 190, height: 2}},
    'the deepest hit wins: b inside g inside the root');
  assert.equal(resolveFilterDrop(root, hits, 100, 75, d.id).index, 1, 'after b');
  assert.equal(resolveFilterDrop(root, hits, 100, 75, b.id), null, 'right after itself is where b is');
});

test('a move refuses itself, its own subtree and every position it is already in', () => {
  const {root, a, b, c, d, e, g, hits} = fixture();
  assert.equal(resolveFilterDrop(root, hits, 100, 20, a.id), null, 'onto itself');
  assert.equal(resolveFilterDrop(root, hits, 100, 45, g.id), null, 'into itself');
  assert.equal(resolveFilterDrop(root, hits, 100, 60, g.id), null, 'beside its own child');
  assert.equal(resolveFilterDrop(root, hits, 100, 10, a.id), null, 'before itself');
  assert.equal(resolveFilterDrop(root, hits, 100, 35, a.id), null, 'right after itself');
  assert.equal(resolveFilterDrop(root, hits, 100, 45, c.id), null, 'the tail of its own parent is where c is');
  assert.equal(resolveFilterDrop(root, hits, 100, 195, e.id), null, 'the root tail is where the last node is');
  assert.deepEqual(resolveFilterDrop(root, hits, 100, 195, a.id),
    {parentId: root.id, index: 3, kind: 'into', rect: {x: 0, y: 0, width: 200, height: 200}},
    'the root padding appends into the root');
  assert.equal(resolveFilterDrop(root, hits, 100, 45, b.id).index, 1, 'into the parent it already ends');
  assert.equal(resolveFilterDrop(root, hits, 100, 115, b.id).index, 1, 'after c, counted with b gone');
  assert.equal(resolveFilterDrop(root, hits, 100, 115, d.id).index, 2, 'after c, into g');
});

test('nothing under the pointer, unknown ids and an unknown mover resolve to null', () => {
  const {root, a, hits} = fixture();
  assert.equal(resolveFilterDrop(root, hits, 500, 500, a.id), null);
  assert.equal(resolveFilterDrop(root, [], 100, 100, a.id), null);
  assert.equal(resolveFilterDrop(root, hits, 100, 100, 'nope'), null);
  assert.equal(resolveFilterDrop(root, [{id: 'ghost', rect: {x: 0, y: 0, width: 999, height: 999}}], 5, 5, a.id), null,
    'a rect the tree does not know is not a target');
});
