import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectFloat, test} from '@datagrok-libraries/test/src/test';

category('AI: Color: static helpers', () => {
  test('channel extractors', async () => {
    const c = 0x12345678;
    expect(DG.Color.a(c), 0x12);
    expect(DG.Color.r(c), 0x34);
    expect(DG.Color.g(c), 0x56);
    expect(DG.Color.b(c), 0x78);
  });

  test('argb round-trip', async () => {
    const packed = DG.Color.argb(0x12, 0x34, 0x56, 0x78);
    expect(packed, 0x12345678);
    expect(DG.Color.a(packed), 0x12);
    expect(DG.Color.r(packed), 0x34);
    expect(DG.Color.g(packed), 0x56);
    expect(DG.Color.b(packed), 0x78);
  });

  test('setAlpha preserves rgb', async () => {
    const original = 0xFF112233;
    const recoloured = DG.Color.setAlpha(original, 0x80);
    expect(recoloured, 0x80112233);
    expect(DG.Color.a(original), 0xFF);
    expect(DG.Color.r(recoloured), 0x11);
    expect(DG.Color.g(recoloured), 0x22);
    expect(DG.Color.b(recoloured), 0x33);
  });

  test('fromHtml / toHtml round-trip', async () => {
    const i = DG.Color.fromHtml('#00bfff');
    expect((i & 0xFFFFFF) >>> 0, 0x00bfff);
    expect(DG.Color.toHtml(i), '#00bfff');
  });

  test('toRgb formatting', async () => {
    expect(DG.Color.toRgb(0xFF00BFFF), 'rgb(0,191,255)');
  });

  test('scale linear mapping', async () => {
    expectFloat(DG.Color.scale(5, 0, 10), 0.5);
    expectFloat(DG.Color.scale(0, 0, 10), 0);
    expectFloat(DG.Color.scale(10, 0, 10), 1);
    expectFloat(DG.Color.scale(7, 7, 7), 7);
  });

  test('hexToPercentRgb happy path', async () => {
    const res = DG.Color.hexToPercentRgb('#00bfff');
    expect(res != null, true);
    expectArray(res!, [0 / 256, 0xbf / 256, 0xff / 256, 0.3]);
  });

  test('hexToPercentRgb with explicit alpha', async () => {
    const res = DG.Color.hexToPercentRgb('#00bfff80');
    expect(res != null, true);
    expectArray(res!, [0 / 256, 0xbf / 256, 0xff / 256, 0x80 / 256]);
  });

  test('hexToPercentRgb invalid input', async () => {
    expect(DG.Color.hexToPercentRgb('not-a-color'), null);
  });

  test('getCategoricalColor wrap-around', async () => {
    const palette = DG.Color.categoricalPalette;
    expect(Array.isArray(palette), true);
    expect(palette.length > 0, true);
    const n = palette.length;
    const indices = [0, n - 1, n, 2 * n + 3];
    for (var i of indices)
      expect(DG.Color.getCategoricalColor(i), palette[i % n]);
  });
});
