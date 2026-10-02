const test = require('node:test');
const assert = require('node:assert/strict');
const DG = require('./dg.cjs');

const same = (r, x, y, width, height) => assert.deepEqual([r.x, r.y, r.width, r.height], [x, y, width, height]);

test('Rect.fromPoints normalizes the corners', () => {
  same(DG.Rect.fromPoints(10, 20, 0, 5), 0, 5, 10, 15);
});

test('Rect.fromCenterSize is static and centers the rectangle', () => {
  same(DG.Rect.fromCenterSize(50, 50, 20, 10), 40, 45, 20, 10);
});

test('Rect edges and midpoints', () => {
  const r = new DG.Rect(10, 20, 30, 40);
  assert.deepEqual([r.left, r.top, r.right, r.bottom], [10, 20, 40, 60]);
  assert.deepEqual([r.minX, r.maxX, r.minY, r.maxY], [10, 40, 20, 60]);
  assert.deepEqual([r.midX, r.midY], [25, 40]);
});

test('Rect.getGridPart divides both axes by their counts', () => {
  same(new DG.Rect(0, 0, 100, 40).getGridPart(2, 2, 1, 1), 50, 20, 50, 20);
  same(new DG.Rect(0, 0, 90, 30).getGridPart(3, 3, 0, 2), 0, 20, 30, 10);
});

test('Rect cut* and get* slices', () => {
  const r = new DG.Rect(0, 0, 100, 50);
  same(r.cutLeft(10), 10, 0, 90, 50);
  same(r.cutTop(5), 0, 5, 100, 45);
  same(r.cutRight(10), 0, 0, 90, 50);
  same(r.cutBottom(5), 0, 0, 100, 45);
  same(r.getLeft(10), 0, 0, 10, 50);
  same(r.getRight(10), 90, 0, 10, 50);
  same(r.getTop(5), 0, 0, 100, 5);
  same(r.getBottom(5), 0, 45, 100, 5);
});

test('Rect.move and neighbours', () => {
  const r = new DG.Rect(10, 10, 20, 20);
  same(r.move(5, -5), 15, 5, 20, 20);
  same(r.below(10), 10, 30, 20, 10);
  same(r.above(10), 10, 0, 20, 10);
  same(r.toTheLeft(10), 0, 10, 10, 20);
  same(r.toTheRight(10), 30, 10, 10, 20);
});
