import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';
import {before, category, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';

category('BitSet', () => {
  let t1: DG.BitSet;
  let t2: DG.BitSet;

  before(async () => {
    t1 = DG.BitSet.create(5);
    t2 = DG.BitSet.fromString('11001');
  });

  test('create', async () => {
    expect(DG.BitSet.create(5).toBinaryString(), '00000');
  });

  test('create fromString', async () => {
    expect(DG.BitSet.fromString('11001').toBinaryString(), '11001');
  });

  test('create fromBytes', async () => {
    const buffer = new ArrayBuffer(1);
    expect(DG.BitSet.fromBytes(buffer, 5).toBinaryString(), '00000');
  });

  test('and method', async () => {
    expect(t1.and(t2).toBinaryString(), '00000');
  });

  test('andNot method', async () => {
    expect(t1.andNot(t2).toBinaryString(), '00000');
  });

  test('clone method', async () => {
    expect(t2.clone().toBinaryString(), '11001');
  });

  test('copyFrom method', async () => {
    expect(t1.copyFrom(t2).toBinaryString(), '11001');
  });

  test('findNext method', async () => {
    expect(t2.findNext(2, true).toString(), '4');
  });

  test('findPrev method', async () => {
    expect(t2.findPrev(3, false).toString(), '2');
  });

  test('get method', async () => {
    expect(t2.get(2).toString(), 'false');
  });

  test('getBuffer method', async () => {
    expectArray(Array.from(t2.getBuffer()), [19]);
  });

  test('getSelectedIndexes method', async () => {
    expect(t2.getSelectedIndexes().toString(), '0,1,4');
  });

  test('invert method', async () => {
    expect(t2.invert().toBinaryString(), '00110');
  });

  test('or method', async () => {
    expect(t1.or(t2).toBinaryString(), '11111');
  });

  test('set method', async () => {
    t2.set(2, true);
    expect(t2.toBinaryString(), '00110');
  });

  test('setAll method', async () => {
    t2.setAll(true);
    expect(t2.toBinaryString(), '11111');
  });

  test('similarityTo method', async () => {
    const t3 = DG.BitSet.fromString('11111');
    expect(t2.similarityTo(t3).toString(), '1');
  });

  test('toBinaryString method', async () => {
    expect(t2.toBinaryString(), '11111');
  });

  test('toString method', async () => {
    expect(t2.toString(), '5 bits, 5 set');
  });

  test('xor method', async () => {
    const t3 = DG.BitSet.fromString('11111');
    expect(t2.xor(t3).toBinaryString(), '00000');
  });
}, {clear: false, owner: 'aparamonov@datagrok.ai'});
