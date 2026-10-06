const test = require('node:test');
const assert = require('node:assert/strict');
const DG = require('./dg.cjs');

test('Color channels round-trip through argb', () => {
  const c = DG.Color.argb(0x80, 0x11, 0x22, 0x33);
  assert.deepEqual([DG.Color.a(c), DG.Color.r(c), DG.Color.g(c), DG.Color.b(c)], [0x80, 0x11, 0x22, 0x33]);
});

test('Color.setAlpha replaces only the alpha channel', () => {
  const c = DG.Color.setAlpha(DG.Color.argb(0xFF, 1, 2, 3), 0x40);
  assert.deepEqual([DG.Color.a(c), DG.Color.r(c), DG.Color.g(c), DG.Color.b(c)], [0x40, 1, 2, 3]);
});

test('Color.toRgb formats the channels', () => {
  assert.equal(DG.Color.toRgb(DG.Color.argb(0xFF, 255, 128, 0)), 'rgb(255,128,0)');
});

test('Color.hexToPercentRgb scales by 255 and rejects junk', () => {
  assert.deepEqual(Array.from(DG.Color.hexToPercentRgb('#ff8000')), [1, 128 / 255, 0, 0.3]);
  assert.deepEqual(Array.from(DG.Color.hexToPercentRgb('#ff800080')), [1, 128 / 255, 0, 128 / 255]);
  assert.equal(DG.Color.hexToPercentRgb('nope'), null);
});
