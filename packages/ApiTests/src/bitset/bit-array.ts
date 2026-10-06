import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

category('BitArray', () => {
  const pattern = (n: number) => DG.BitArray.create(n, (i) => i % 3 === 0 || i % 7 === 0);
  const patternString = (n: number) => Array.from({length: n}, (_, i) => i % 3 === 0 || i % 7 === 0 ? '1' : '0').join('');

  test('create, fromString, toString', async () => {
    expect(DG.BitArray.create(5).toString(), '00000');
    expect(new DG.BitArray(5, true).toString(), '11111');
    expect(DG.BitArray.fromString('11001').toString(), '11001');
    expect(DG.BitArray.fromString('11001').toBinaryString(), '11001');
    for (const n of [5, 33, 70]) {
      const s = patternString(n);
      expect(pattern(n).toString(), s);
      expect(DG.BitArray.fromString(s).equals(pattern(n)), true);
      expect(DG.BitArray.fromString(s).length, n);
    }
  });

  test('fromBytes', async () => {
    const a = DG.BitArray.fromBytes(new Uint8Array([1, 0, 1]));
    expect(a.length, 24);
    expectArray(Array.from(a.getSelectedIndexes()), [0, 16]);
    const b = DG.BitArray.fromBytes(new Uint8Array([1, 0, 1]), 20);
    expect(b.length, 20);
    expectArray(Array.from(b.getSelectedIndexes()), [0, 16]);
    expectArray(Array.from(DG.BitArray.fromBytes(new Uint8Array([0x80, 0x01]).buffer).getSelectedIndexes()), [7, 8]);
    expectArray(Array.from(DG.BitArray.fromBytes(new Uint8Array([0, 0, 0, 0, 0, 0, 0x80])).getSelectedIndexes()), [55]);
    let thrown = false;
    try {
      DG.BitArray.fromBytes(new Uint8Array(2), 17);
    }
    catch (e) {
      thrown = e instanceof RangeError;
    }
    expect(thrown, true, 'bitLength beyond the bytes must throw');
  });

  test('get, set, trueCount, falseCount', async () => {
    const a = DG.BitArray.fromString('11001');
    expect(a.get(0), true);
    expect(a.get(2), false);
    expect(a.trueCount, 3);
    expect(a.falseCount, 2);
    a.set(2, true);
    expect(a.trueCount, 4);
    a.set(0, false);
    expect(a.toString(), '01101');
    a.getBuffer()[0] = 0;
    expect(a.trueCount, 0, 'no stale count after a direct buffer write');
    expect(a.anyTrue, false);
    expect(a.allFalse, true);
    a.setAll(true);
    expect(a.allTrue, true);
    expect(a.anyFalse, false);
  });

  test('and, or, xor, andNot, invert', async () => {
    const b1 = DG.BitArray.fromString('11100');
    const b2 = DG.BitArray.fromString('00111');
    expect(b1.clone().and(b2).toString(), '00100');
    expect(b1.clone().or(b2).toString(), '11111');
    expect(b1.clone().xor(b2).toString(), '11011');
    expect(b1.clone().andNot(b2).toString(), '11000');
    expect(b1.clone().invert().toString(), '00011');
    expect(b1.andWithCountBits(b2), 1);
    expect(b1.andWithCountBits(b2, false), 4);
    let thrown = false;
    try {
      b1.and(DG.BitArray.create(6));
    }
    catch (e) {
      thrown = e instanceof RangeError;
    }
    expect(thrown, true, 'length mismatch must throw');
  });

  test('tail bits stay clear', async () => {
    const a = new DG.BitArray(70);
    a.invert();
    expect(a.getBuffer()[2] >>> (70 & 31), 0);
    expect(a.trueCount, 70);
    expect(a.equals(new DG.BitArray(70, true)), true);
    const dirty = new DG.BitArray(new Uint32Array([0xffffffff, 0xffffffff, 0xffffffff]), 70);
    expect(dirty.equals(a), true);
    expect(dirty.trueCount, 70);
    a.setLength(75);
    expect(a.trueCount, 70);
    expect(a.get(72), false);
  });

  test('findNext, findPrev', async () => {
    const a = DG.BitArray.fromString('11001');
    expect(a.findNext(2), 4);
    expect(a.findPrev(3, false), 2);
    expect(a.findNext(-1), 0);
    expect(a.findNext(4), -1);
    expect(a.findPrev(-1), 4);
    expect(a.findPrev(0), -1);
    expect(a.findNext(-1, false), 2);
    const b = DG.BitArray.create(70, (i) => i === 40 || i === 69);
    expect(b.findNext(-1), 40);
    expect(b.findNext(40), 69);
    expect(b.findPrev(-1), 69);
    expect(b.findPrev(69), 40);
    expect(b.findPrev(40), -1);
  });

  test('getSelectedIndexes, trueIndexes', async () => {
    const a = DG.BitArray.fromString('11001');
    expectArray(Array.from(a.getSelectedIndexes()), [0, 1, 4]);
    expectArray(Array.from(a.trueIndexes()), [0, 1, 4]);
    expectArray(Array.from(DG.BitArray.create(70, (i) => i === 40).trueIndexes()), [40]);
  });

  test('equals, clone, copyFrom', async () => {
    const a = DG.BitArray.fromString('11001');
    const b = a.clone();
    expect(b.equals(a), true);
    b.set(2, true);
    expect(a.toString(), '11001');
    expect(b.equals(a), false);
    expect(a.equals(DG.BitArray.fromString('110010')), false);
    const c = DG.BitArray.create(5);
    expect(c.copyFrom(b).toString(), '11101');
    let thrown = false;
    try {
      c.copyFrom(DG.BitArray.create(6));
    }
    catch (e) {
      thrown = e instanceof RangeError;
    }
    expect(thrown, true, 'length mismatch must throw');
  });

  test('setRange, getRange', async () => {
    const a = new DG.BitArray(100);
    a.setRange(30, 70, true);
    expect(a.trueCount, 40);
    expect(a.get(29), false);
    expect(a.get(30), true);
    expect(a.get(69), true);
    expect(a.get(70), false);
    expect(a.getRange(28, 72).toString(), '00' + '1'.repeat(40) + '00');
    a.setRange(35, 40, false);
    expect(a.trueCount, 35);
    expect(a.getRange(33, 42).toString(), '110000011');
    expect(a.getRange(0, 0).length, 0);
    expect(a.setRange(0, 100, true).getRange(0, 100).equals(a), true);
    let thrown = false;
    try {
      a.setRange(90, 101, true);
    }
    catch (e) {
      thrown = e instanceof RangeError;
    }
    expect(thrown, true, 'to > length must throw');
  });

  test('setLength, removeAt', async () => {
    const a = DG.BitArray.fromString('11001');
    a.setLength(8);
    expect(a.toString(), '11001000');
    a.setLength(3);
    expect(a.toString(), '110');
    const b = DG.BitArray.fromString('11001');
    b.removeAt(1);
    expect(b.toString(), '1001');
    b.removeAt(0, 2);
    expect(b.toString(), '01');
  });

  test('DG.BitArray is DG.U2.BitArray', async () => {
    expect(DG.BitArray === DG.U2.BitArray, true);
  });

  test('BitSet.fromBitArray, toBitArray', async () => {
    for (const n of [5, 32, 33, 70]) {
      const a = pattern(n);
      const bs = DG.BitSet.fromBitArray(a);
      expect(bs.length, n);
      expect(bs.toBinaryString(), a.toString());
      const back = bs.toBitArray();
      expect(back.equals(a), true);
      back.setAll(true);
      expect(bs.toBinaryString(), a.toString(), 'toBitArray must copy');
    }
    expect(DG.BitSet.fromBitArray(new DG.BitArray(70).invert()).toBinaryString(), '1'.repeat(70));
    const oversized = DG.BitArray.fromUint32Array(5, new Uint32Array([1, 0xffffffff, 0xffffffff]));
    expect(DG.BitSet.fromBitArray(oversized).toBinaryString(), '10000');
  });

  test('DataFrame filter round trip', async () => {
    const df = grok.data.demo.demog(100);
    const every3rd = DG.BitArray.create(100, (i) => i % 3 === 0);
    let changed = 0;
    const sub = df.filter.onChanged.subscribe(() => changed++);
    try {
      df.filter.copyFrom(every3rd);
      expect(df.filter.trueCount, 34);
      expect(changed, 1);
      expect(df.filter.toBinaryString(), every3rd.toString());
      expect(df.filter.toBitArray().equals(every3rd), true);
      df.filter.and(DG.BitSet.fromBitArray(DG.BitArray.create(100, (i) => i < 50)));
      expect(df.filter.trueCount, 17);
      let thrown = false;
      try {
        df.filter.copyFrom(DG.BitArray.create(99));
      }
      catch (e) {
        thrown = e instanceof RangeError;
      }
      expect(thrown, true, 'length mismatch must throw');
      expect(df.filter.trueCount, 17);
    }
    finally {
      sub.unsubscribe();
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
