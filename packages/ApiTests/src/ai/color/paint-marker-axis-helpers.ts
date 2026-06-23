import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Tests DG.Paint marker/axis helpers: markerTypes, markerTypeIndexes.
category('AI: Color: Paint / Marker / Axis helpers', () => {
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
});
