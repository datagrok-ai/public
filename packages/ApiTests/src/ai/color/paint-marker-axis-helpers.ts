import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow} from '../helpers';

// Tests DG.Paint marker/axis helpers: markerTypes, markerTypeIndexes, marker, horzAxis, vertAxis.
category('AI: Color: Paint / Marker / Axis helpers', () => {
  const ctx = (): CanvasRenderingContext2D => {
    const c = document.createElement('canvas');
    c.width = 200;
    c.height = 200;
    return c.getContext('2d')!;
  };

  test('markerTypeIndexes returns map with keys matching markerTypes', async () => {
    const types = DG.Paint.markerTypes();
    const indexes = DG.Paint.markerTypeIndexes();
    expect(indexes != null, true);
    const keys = Object.keys(indexes);
    expect(keys.length, types.length);
    for (const t of types) {
      expect(t in indexes, true);
      expect(typeof indexes[t], 'number');
    }
  });

  test('markerTypeIndexes key for circle matches MARKER_TYPE.CIRCLE index', async () => {
    const types = DG.Paint.markerTypes();
    const indexes = DG.Paint.markerTypeIndexes();
    const circle = DG.MARKER_TYPE.CIRCLE;
    expect(types.indexOf(circle) >= 0, true);
    expect(indexes[circle], types.indexOf(circle));
  });

  test('markerTypes() are valid MARKER_TYPE values incl. core circle/square', async () => {
    // markerTypes() returns the user-selectable subset (markerIcons.keys), not the full
    // MARKER_TYPE enum (which also has border variants, outlier, gradient, dot, etc.).
    const types = DG.Paint.markerTypes();
    expect(Array.isArray(types), true);
    expect(types.length > 0, true);
    const enumValues = Object.values(DG.MARKER_TYPE) as string[];
    // core renderable types are always selectable
    expect(types.indexOf(DG.MARKER_TYPE.CIRCLE) >= 0, true);
    expect(types.indexOf(DG.MARKER_TYPE.SQUARE) >= 0, true);
    // every returned type is a non-empty string and a legitimate MARKER_TYPE enum value
    for (const t of types) {
      expect(typeof t, 'string');
      expect(enumValues.indexOf(t) >= 0, true);
    }
    // and the intersection is non-empty
    expect(types.some((t) => enumValues.indexOf(t) >= 0), true);
  });

  test('Paint.marker does not throw for each MARKER_TYPE value (canvas ctx)', async () => {
    const g = ctx();
    for (const v of Object.values(DG.MARKER_TYPE))
      expectNoThrow(() => DG.Paint.marker(g, v, 100, 100, DG.Color.blue, 12));
    // also accepts an html color string
    expectNoThrow(() => DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, 50, 50, '#ff0000', 8));
  });

  test('Paint.horzAxis / vertAxis do not throw (linear, log, inverse)', async () => {
    const g = ctx();
    expectNoThrow(() => DG.Paint.horzAxis(g, 0, 100, 0, 180, 200, 20));
    expectNoThrow(() => DG.Paint.vertAxis(g, 1, 1000, 0, 0, 30, 200, true));
    expectNoThrow(() => DG.Paint.horzAxis(g, 0, 100, 0, 180, 200, 20, false, true));
  });
});
