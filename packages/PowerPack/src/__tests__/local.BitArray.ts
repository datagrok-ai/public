/**
 * @jest-environment jsdom
 */

import * as _package from '../package';
import BitArray from 'utils/src/bit-array'

jest.mock('datagrok-api/dg', () => {
    const originalModule = jest.requireActual('datagrok-api/dg');
    return {
        __esModule: true,
        ...originalModule,
        Func: {find: (o: object) => []},
    };
});

test('BitArray', () => {
    let bitArray = new BitArray(10);
})

test('BitArray', () => {
    let bitArray = new BitArray(10);
    bitArray.setBit(5, true);
    expect(bitArray.getBit(5)).toBe(true);
})

test('BitArray', () => {
    let bitArray = new BitArray(10);
    bitArray.setBit(5, true);
    expect(bitArray.buffer[0]).toBe((1 << 5));
})

test('BitArray', () => {
    let bitArray = new BitArray(2);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    expect(bitArray.buffer[0]).toBe((1 << 5) + (1 << 11));
})

test('BitArray', () => {
    let bitArray = new BitArray(2);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    expect(bitArray.buffer[0]).toBe((1 << 5) + (1 << 11));
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    expect(bitArray.trueCount()).toBe(3);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    expect(bitArray.lengthInInts).toBe(2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    expect(bitArray.lengthInInts).toBe(2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.clear()
    expect(bitArray.lengthInInts).toBe(0);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.clear()
    expect(bitArray.allFalse).toBe(true);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    expect(bitArray.findNext(15)).toBe(35);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    expect(bitArray.findPrev(30)).toBe(11);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let bitArray2 = new BitArray(37);
    bitArray2.setBit(11, true);
    expect(bitArray.andWithCountBits(bitArray2, true)).toBe(1);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let func = (i: number) => {i % 2 == 0};
    expect(bitArray.countWhere(func)).toBe(0);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let func = (i: number) => {return i % 2 == 1};
    expect(bitArray.countWhere(func)).toBe(3);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.removeAt(0, 5);
    expect(bitArray.buffer[0]).toBe(1 + (1 << 6) + (1 << 30));
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.removeAt(0, 5);
    let bitArray2 = new BitArray(32);
    bitArray2.copyFrom(bitArray);
    expect(bitArray2.buffer[0]).toBe(1 + (1 << 6) + (1 << 30));
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let bitArray2 = new BitArray(37);
    bitArray2.setBit(8, true);
    bitArray2.setBit(11, true);
    expect(BitArray.fromAnd(bitArray, bitArray2).trueCount()).toBe(1);
})

test('BitArray', () => {
    expect(BitArray.fromValues(Array.from([false, false, true, false, true])).trueCount()).toBe(2);
})

test('BitArray', () => {
    expect(BitArray.fromValues(Array.from([false, false, true, false, true])).buffer[0]).toBe(20);
})

test('BitArray', () => {
    let func = (i: number) => {return i % 2 == 1};
    expect(BitArray.fromSeq(10, func).trueCount()).toBe(5);
})

test('BitArray', () => {
    let func = (i: number) => {return i % 2 == 1};
    expect(BitArray.fromSeq(10, func).trueCount()).toBe(5);
})

test('BitArray', () => {
    let func = (i: number) => {return i % 2 == 1};
    expect(BitArray.fromString("0011001").trueCount()).toBe(3);
})

test('BitArray', () => {
    let temp = new Uint32Array(10);
    temp[0] = 5;
    temp[1] = 3;
    expect(BitArray.fromUint32Array(60, temp).trueCount()).toBe(4);
})

test('BitArray', () => {
    let temp = new Uint8Array(4);
    temp[0] = 5;
    temp[1] = 3;
    temp[2] = 0;
    temp[3] = 10;
    expect(BitArray.fromBytes(temp).trueCount()).toBe(6);
})

test('BitArray', () => {
    let temp = new Uint8Array(4);
    temp[0] = 5;
    temp[1] = 3;
    temp[2] = 0;
    temp[3] = 10;
    expect(BitArray.fromBytes(temp).toString()).toBe("32 bits, 6 set");
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let bitArray2 = new BitArray(37);
    bitArray2.setBit(5, true);
    bitArray2.setBit(11, true);
    bitArray2.setBit(11, true);
    bitArray2.setBit(35, true);
    bitArray2.setBit(35, false);
    bitArray2.setBit(35, true);
    expect(bitArray.equals(bitArray2)).toBe(true);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let bitArray2 = new BitArray(0);
    bitArray2 = bitArray.clone();
    expect(bitArray.equals(bitArray2)).toBe(true);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.invert();
    expect(bitArray.trueCount()).toBe(37 - 3);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.invert();
    bitArray.setAll(false);
    expect(bitArray.trueCount()).toBe(0);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.invert();
    bitArray.setAll(true);
    expect(bitArray.trueCount()).toBe(37);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.invert();
    bitArray.setAll(true);
    expect(bitArray.trueCount()).toBe(37);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let index = Array.from([0, 5, 8, 36]);
    bitArray.setIndexes(index, true, false);
    expect(bitArray.buffer[0]).toBe(1 + (1 << 5) + (1 << 11) + (1 << 8));
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let index = Array.from([0, 5, 8, 36]);
    bitArray.setIndexes(index, true, false);
    expect(bitArray.everyIndex(index)).toBe(true);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    let index = Array.from([0, 5, 8, 36]);
    expect(bitArray.everyIndex(index)).toBe(false);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setWhere((i: number) => {return i % 2 == 0}, true, false);
    expect(bitArray.trueCount()).toBe(3 + Math.ceil(37 / 2));
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setRange(0, 7, true);
    expect(bitArray.trueCount()).toBe(8 + 2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setRange(0, 7, true);
    let bitArray2 = new BitArray(37);
    bitArray2.or(bitArray);
    expect(bitArray2.trueCount()).toBe(8 + 2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setRange(0, 7, true);
    let bitArray2 = new BitArray(37);
    bitArray2.setBit(12, true);
    bitArray2.setBit(22, true);
    bitArray2.setBit(1, true);
    bitArray2.andNot(bitArray);
    console.log(bitArray2.buffer);
    expect(bitArray2.trueCount()).toBe(2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setRange(0, 7, true);
    let bitArray2 = new BitArray(37);
    bitArray2.setBit(12, true);
    bitArray2.setBit(22, true);
    bitArray2.setBit(1, true);
    bitArray2.andNot(bitArray);
    expect(bitArray2.trueCount()).toBe(2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setRange(0, 7, true);
    bitArray.insertAt(0, 7, false);
    expect(bitArray.trueCount()).toBe(8 + 2);
})

test('BitArray', () => {
    let bitArray = new BitArray(37);
    bitArray.setBit(5, true);
    bitArray.setBit(11, true);
    bitArray.setBit(35, true);
    bitArray.setRange(0, 7, true);
    let bitArray2 = new BitArray(37);
    bitArray2.setBit(12, true);
    bitArray2.setBit(5, true);
    bitArray2.setBit(1, true);
    console.log(bitArray2.buffer);
    bitArray2.removeByMask(bitArray);
    console.log(bitArray2.buffer);
    expect(bitArray2.trueCount()).toBe(1);
})