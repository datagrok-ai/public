/* `DG.U2.BitArray` (the compiled js-api class u2 core reaches through `datagrok-api/u2core`): the
   reviewer repros that failed on the old utils class — stale counts, `fromBytes` byte loss, dirty
   tails, `setRange`/`getRange` bounds, `removeAt`, `findNext`/`findPrev` edges — plus the structured
   clone contract the ml workers rely on. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {BitArray} from 'datagrok-api/u2core';

const throws = (f) => assert.throws(f, RangeError);
const tailClean = (a) => a.length % 32 === 0 || (a.getBuffer()[a.lengthInInts - 1] >>> (a.length & 31)) === 0;
const indexes = (a) => Array.from(a.getSelectedIndexes());

test('trueCount is never stale: set, and direct getBuffer() writes', () => {
  const a = new BitArray(10);
  assert.equal(a.trueCount, 0);
  a.set(3, true);
  assert.equal(a.trueCount, 1);
  const b = new BitArray(10, true);
  assert.equal(b.trueCount, 10);
  b.getBuffer()[0] = 0;
  assert.equal(b.trueCount, 0);
  b.getBuffer()[0] = 0xff;
  assert.equal(b.trueCount, 8);
  assert.equal(b.falseCount, 2);
});

test('fromBytes: 1/2/3 trailing bytes, bitLength, ArrayBuffer and offset views', () => {
  const b = BitArray.fromBytes(Uint8Array.from([0x01, 0x02, 0x04]));
  assert.deepEqual([b.get(0), b.get(9), b.get(18)], [true, true, true]);
  assert.deepEqual([b.trueCount, b.length], [3, 24]);
  assert.equal(BitArray.fromBytes(Uint8Array.from([0x01, 0x02])).trueCount, 2);
  assert.equal(BitArray.fromBytes(Uint8Array.from([0x01])).trueCount, 1);
  assert.equal(BitArray.fromBytes(Uint8Array.from([0x01, 0x02, 0x04, 0x08])).trueCount, 4);
  assert.equal(BitArray.fromBytes(Uint8Array.from([0xff, 0xff, 0xff, 0xff, 0x01, 0x02, 0x04])).trueCount, 35);
  const c = BitArray.fromBytes(Uint8Array.from([0xff, 0xff, 0xff]), 20);
  assert.deepEqual([c.length, c.trueCount, tailClean(c)], [20, 20, true]);
  throws(() => BitArray.fromBytes(new Uint8Array(2), 17));
  assert.deepEqual(indexes(BitArray.fromBytes(Uint8Array.from([1, 0, 1]).buffer)), [0, 16]);
  const off = new Uint8Array([9, 9, 1, 0, 1]);
  assert.deepEqual(indexes(BitArray.fromBytes(off.subarray(2))), [0, 16]);
  const words = new Uint32Array([0b11111101]);
  const d = BitArray.fromBytes(words, 6);
  assert.deepEqual([...d.getBuffer()], [0b111101], 'a Uint32Array is read as words and trimmed');
  assert.notEqual(d.getBuffer(), words, 'a copy');
});

test('the tail beyond length is clear after invert, setAll, setLength, fromUint32Array; equals is a word compare', () => {
  const a = new BitArray(10).invert();
  assert.deepEqual([a.get(10), a.get(31)], [false, false]);
  assert.deepEqual([a.trueCount, tailClean(a)], [10, true]);
  const b = new BitArray(96, true);
  b.setLength(40);
  const c = new BitArray(40, true);
  assert.deepEqual([b.trueCount, c.trueCount, b.equals(c)], [40, 40, true]);
  assert.equal(new BitArray(10, true).getBuffer()[0], 0x3ff);
  assert.equal(b.equals(new BitArray(40).setRange(0, 40, true)), true);
  assert.equal(new BitArray(10).invert().equals(new BitArray(10).setRange(0, 10, true)), true);
  const dirty = new Uint32Array([0xffffffff, 0xffffffff, 0xffffffff]);
  const e = BitArray.fromUint32Array(70, dirty);
  assert.deepEqual([e.trueCount, tailClean(e), dirty[2]], [70, true, 0x3f], 'adopt trims in place');
  const f = new BitArray(10, true);
  f.setLength(20);
  assert.deepEqual([f.trueCount, f.get(10), f.get(19)], [10, false, false], 'grown bits are false');
  f.setLength(5);
  assert.deepEqual([f.trueCount, tailClean(f), f.getBuffer().length], [5, true, 1]);
  f.setLength(100);
  assert.deepEqual([f.trueCount, f.getBuffer().length], [5, 4], 'exact reallocation');
  const g = new BitArray(70, true);
  const h = new BitArray(70, true);
  h.getBuffer()[2] |= 0xffffffc0;
  assert.equal(g.equals(h), false, 'a dirty tail is a caller bug, not masked');
  assert.equal(g.equals(new BitArray(70).setAll(true)), true);
  assert.equal(new BitArray(70).setAll(true).invert().trueCount, 0);
  throws(() => new BitArray(new Uint32Array(1), 33));
  throws(() => new BitArray(-1));
});

test('setRange/getRange: [from, to), to may equal length, against a per-bit oracle', () => {
  const a = new BitArray(8).setRange(2, 5, true);
  assert.equal(a.toString(), '00111000');
  assert.deepEqual([a.getRange(2, 5).length, a.getRange(2, 5).toString()], [3, '111']);
  new BitArray(8).setRange(0, 8, true);
  new BitArray(0).getRange(0, 0);
  const b = new BitArray(8);
  throws(() => b.setRange(5, 2, true));
  assert.equal(b.toString(), '00000000');
  throws(() => b.setRange(0, 9, true));
  let bad = 0;
  for (let t = 0; t < 300; t++) {
    const n = 100;
    const from = Math.floor(Math.random() * (n + 1));
    const to = from + Math.floor(Math.random() * (n + 1 - from));
    const x = Math.random() < 0.5;
    const seed = BitArray.create(n, () => Math.random() < 0.5);
    const oracle = seed.clone();
    for (let i = from; i < to; i++)
      oracle.set(i, x);
    const got = seed.clone().setRange(from, to, x);
    if (!got.equals(oracle) || !tailClean(got))
      bad++;
    const sub = seed.getRange(from, to);
    const subOracle = BitArray.create(to - from, (i) => seed.get(from + i));
    if (!sub.equals(subOracle) || !tailClean(sub))
      bad++;
  }
  assert.equal(bad, 0);
  const c = new BitArray(100, true).setRange(31, 65, false);
  assert.deepEqual([c.trueCount, c.get(30), c.get(31), c.get(64), c.get(65)], [66, true, false, false, true]);
});

test('removeAt shifts the following bits down and reallocates exactly', () => {
  const a = BitArray.fromString('10110');
  a.removeAt(1, 2);
  assert.equal(a.toString(), '110');
  const b = BitArray.fromString('00000');
  b.removeAt(1, 2);
  assert.deepEqual([b.toString(), b.length], ['000', 3]);
  const c = BitArray.fromString('10110');
  c.removeAt(4, 1);
  assert.deepEqual([c.toString(), c.length], ['1011', 4]);
  const d = BitArray.fromString('10110');
  throws(() => d.removeAt(0, -2));
  assert.equal(d.length, 5);
  throws(() => d.removeAt(4, 2));
  const e = new BitArray(70, true);
  e.removeAt(0, 40);
  assert.deepEqual([e.length, e.trueCount, tailClean(e), e.getBuffer().length], [30, 30, true, 1]);
});

test('findNext/findPrev are exclusive of i, -1 means the edge, and cross word boundaries', () => {
  const a = BitArray.fromString('00001');
  assert.deepEqual([a.findNext(-1, true), a.findNext(3, true), a.findNext(4, true)], [4, 4, -1]);
  assert.deepEqual([a.findPrev(-1, true), a.findPrev(4, true), a.findPrev(0, false)], [4, -1, -1]);
  const b = new BitArray(64);
  b.set(63, true);
  assert.deepEqual([b.findNext(-1, true), b.findPrev(-1, true), b.findNext(31, true)], [63, 63, 63]);
  const c = new BitArray(33, true);
  assert.deepEqual([c.findNext(-1, false), c.findPrev(-1, false), c.anyFalse, c.allTrue], [-1, -1, false, true]);
  const d = new BitArray(40, true);
  d.set(39, false);
  assert.deepEqual([d.findNext(-1, false), d.findPrev(-1, false)], [39, 39]);
  const e = new BitArray(1);
  assert.equal(e.findNext(-1, true), -1);
  e.set(0, true);
  assert.deepEqual([e.findNext(-1, true), e.anyTrue], [0, true]);
  const f = new BitArray(2);
  f.set(1, true);
  assert.equal(f.anyTrue, true);
  const g = new BitArray(0);
  assert.deepEqual([g.findNext(-1), g.findPrev(-1), g.anyTrue, g.anyFalse, g.allTrue, g.allFalse],
    [-1, -1, false, false, true, true]);
  const h = new BitArray(100).setRange(32, 64, true);
  assert.deepEqual([h.findPrev(64, true), h.findPrev(32, true), h.findNext(31, false), h.findPrev(70, false),
    h.findPrev(64, false)], [63, -1, 64, 69, 31]);
});

test('andWithCountBits counts without mutating; equals handles length, null, self', () => {
  const a = BitArray.fromString('1100');
  const b = BitArray.fromString('1010');
  assert.deepEqual([a.andWithCountBits(b, true), a.andWithCountBits(b, false), a.andWithCountBits(b)], [1, 3, 1]);
  assert.equal(a.toString(), '1100');
  assert.equal(new BitArray(70, true).andWithCountBits(new BitArray(70, true), true), 70);
  const c = BitArray.fromString('101');
  assert.deepEqual([c.equals(BitArray.fromString('1010')), c.equals(null), c.equals(undefined), c.equals(c)],
    [false, false, false, true]);
  const d = new BitArray(10, true);
  d.setLength(0);
  assert.deepEqual([d.length, d.trueCount], [0, 0]);
  d.setLength(10);
  assert.equal(d.trueCount, 0);
});

test('setLength growing inside an adopted over-long buffer zeroes its stale words', () => {
  const a = new BitArray(new Uint32Array([0x7, 0xffffffff, 0xffffffff]), 3);
  assert.equal(a.toString(), '111');
  a.setLength(40);
  assert.deepEqual([a.toString(), a.trueCount, a.get(35)], ['111' + '0'.repeat(37), 3, false]);
  const b = new BitArray(new Uint32Array([0, 0xdeadbeef, 0]), 5);
  b.setLength(64);
  assert.equal(b.trueCount, 0);
  const c = new BitArray(new Uint32Array(4), 100);
  c.setAll(true);
  c.setLength(3);
  c.setLength(100);
  assert.equal(c.trueCount, 3);
});

test('bit ops, clone, copyFrom, init, strings, indexes, structured clone', () => {
  const a = BitArray.fromString('1100');
  const b = BitArray.fromString('1010');
  assert.deepEqual([a.clone().and(b), a.clone().or(b), a.clone().xor(b), a.clone().andNot(b), a.clone().invert()]
    .map(String), ['1000', '1110', '0110', '0100', '0011']);
  throws(() => a.and(BitArray.fromString('1')));
  throws(() => a.copyFrom(BitArray.fromString('1')));
  assert.equal(new BitArray(4).copyFrom(a).toString(), '1100');
  const d = a.clone();
  d.set(0, false);
  assert.deepEqual([a.get(0), d.get(0)], [true, false]);
  assert.equal(BitArray.create(70, (i) => i % 3 === 0).trueCount, 24);
  assert.equal(new BitArray(5).init((i) => i < 2).toString(), '11000');
  const s = '1'.repeat(33) + '0'.repeat(37);
  assert.equal(BitArray.fromString(s).toString(), s);
  assert.equal(a.toBinaryString(), '1100');
  const e = BitArray.create(100, (i) => i % 7 === 0);
  assert.deepEqual(indexes(e), [0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98]);
  assert.deepEqual(Array.from(e.trueIndexes()), indexes(e));
  assert.deepEqual(indexes(new BitArray(64).setRange(31, 33, true)), [31, 32]);
  assert.deepEqual([e.falseCount, new BitArray(96, true).trueCount], [85, 96]);
  assert.deepEqual([0, 32, 33].map((n) => new BitArray(n).lengthInInts), [0, 1, 2]);
  const f = BitArray.fromString('10110');
  const cloned = structuredClone(f);
  assert.equal(new BitArray(cloned._data, cloned._length).equals(f), true, 'the ml workers rebuild it this way');
  const adopted = new Uint32Array(10);
  adopted[0] = 0b101;
  const g = BitArray.fromUint32Array(3, adopted);
  assert.deepEqual([g.toString(), g.getBuffer() === adopted, g.clone().getBuffer().length], ['101', true, 1]);
});
