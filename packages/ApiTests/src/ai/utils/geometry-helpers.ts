import * as DG from 'datagrok-api/dg';
import {category, expect, test, expectFloat} from '@datagrok-libraries/test/src/test';

// Tests DG.Rect / DG.Point / DG.SemanticValue geometry helpers.
category('AI: Utils: Geometry', () => {
  const expectRect = (r: DG.Rect, x: number, y: number, width: number, height: number): void => {
    expect(r.x, x); expect(r.y, y); expect(r.width, width); expect(r.height, height);
  };
  const expectPoint = (p: DG.Point, x: number, y: number): void => {
    expect(p.x, x); expect(p.y, y);
  };

  test('Rect edges, center, min/max', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expect(r.left, 10);
    expect(r.top, 20);
    expect(r.right, 110);
    expect(r.bottom, 60);
    expect(r.midX, 60);
    expect(r.midY, 40);
    expect(r.minX, 10);
    expect(r.maxX, 110);
    expect(r.minY, 20);
    expect(r.maxY, 60);
  });

  test('Rect side splits and corner', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expectRect(r.getTop(10), 10, 20, 100, 10);
    expectRect(r.getBottom(10), 10, 50, 100, 10);
    expectRect(r.getLeft(25), 10, 20, 25, 40);
    expectRect(r.getRight(25), 85, 20, 25, 40);
    expectRect(r.getTopLeft(30, 15), 10, 20, 30, 15);
  });

  test('Rect cuts', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expectRect(r.cutLeft(20), 30, 20, 80, 40);
    expectRect(r.cutTop(5), 10, 25, 100, 35);
    expectRect(r.cutRight(20), 10, 20, 80, 40);
    expectRect(r.cutBottom(5), 10, 20, 100, 35);
  });

  test('Rect relative placement and move', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expectRect(r.below(15), 10, 60, 100, 15);
    expectRect(r.above(15), 10, 5, 100, 15);
    expectRect(r.toTheLeft(30), -20, 20, 30, 40);
    expectRect(r.toTheRight(30), 110, 20, 30, 40);
    expectRect(r.move(5, 7), 15, 27, 100, 40);
  });

  test('Rect parts, grid part, scaled', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expectRect(r.getTopPart(4, 1), 10, 30, 100, 10);
    expectRect(r.getLeftPart(4, 1), 35, 20, 25, 40);
    // getGridPart has a latent bug (height divides by y index, not yCount); pin as-coded behavior.
    expectRect(r.getGridPart(2, 2, 1, 1), 60, 40, 50, 40);
    expectRect(r.getTopScaled(0.25), 10, 20, 100, 10);
    expectRect(r.getLeftScaled(0.25), 10, 20, 25, 40);
  });

  test('Rect inflate and fit', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expectRect(r.inflate(5, 10), 5, 10, 110, 60);
    expectRect(r.inflateSize(10, 10), 10, 20, 110, 50);
    expectRect(r.inflateRel(2, 2), -90, -20, 300, 120);
    expectRect(r.fitSquare(), 40, 20, 40, 40);
    expectRect(r.fit(10, 10), 40, 20, 40, 40);
    expectRect(r.fit(40, 10), 10, 27.5, 100, 25);
  });

  test('Rect contains and union', async () => {
    const r = new DG.Rect(10, 20, 100, 40);
    expect(r.contains(50, 30), true);
    expect(r.contains(0, 0), false);
    expect(r.contains(10, 20), true);
    expect(r.contains(110, 60), true);
    expect(r.containsPoint(new DG.Point(50, 30)), true);
    expect(r.containsPoint(new DG.Point(200, 200)), false);
    expectRect(r.union(new DG.Rect(200, 200, 10, 10)), 10, 20, 200, 190);
  });

  test('Point distance, fromXY, getBounds', async () => {
    expectFloat(new DG.Point(0, 0).distanceTo(new DG.Point(3, 4)), 5);
    expectPoint(DG.Point.fromXY([7, 9]), 7, 9);
    expectRect(DG.Point.getBounds([]), 0, 0, 1, 1);
    expectRect(DG.Point.getBounds([{x: 1, y: 2}, {x: 5, y: 8}, {x: -3, y: 4}]), -3, 2, 8, 6);
  });

  test('SemanticValue construction and round-trip', async () => {
    const sv = DG.SemanticValue.fromValueType(42.5, DG.SEMTYPE.MOLECULE, 'mg/mL');
    expectFloat(sv.value, 42.5);
    expect(sv.units, 'mg/mL');
    expect(sv.semType, DG.SEMTYPE.MOLECULE);
    expect(typeof sv.dataType === 'string' && sv.dataType.length > 0, true);

    sv.value = 99;
    expectFloat(sv.value, 99);
    sv.semType = DG.SEMTYPE.CONCENTRATION;
    expect(sv.semType, DG.SEMTYPE.CONCENTRATION);
    sv.units = 'uM';
    expect(sv.units, 'uM');

    sv.setMeta('origin', 'unit-test');
    expect(sv.getMeta('origin'), 'unit-test');
  });
}, {owner: 'agolovko@datagrok.ai'});
