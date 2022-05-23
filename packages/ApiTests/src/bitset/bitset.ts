import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from "@datagrok-libraries/utils/src/test";

category('BitSet', () => {
    test('create', async () => {
        expect(DG.BitSet.create(5).toBinaryString(), '00000');
    });

    test('create fromString', async () => {
        expect(DG.BitSet.fromString('11001').toBinaryString(), '11001');
    });

    test('create fromBytes', async () => {
        let buffer = new ArrayBuffer(1);
        expect(DG.BitSet.fromBytes(buffer, 5).toBinaryString(), '00000');
    })

    let t1 = DG.BitSet.create(5);
    let t2 = DG.BitSet.fromString('11001');

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
        return t2.getBuffer();
    });

    test('getSelectedIndexes method', async () => {
        expect(t2.getSelectedIndexes().toString(), '0,1,4');
    });

    test('invert method', async () => {
        expect(t2.invert().toBinaryString(), '00110');
    });

    test('or method', async () => {
        expect(t1.or(t2).toBinaryString(), '11001');
    });

    test('set method', async () => {
        t2.set(2, true);
        expect(t2.toBinaryString(), '11101');
    });

    test('setAll method', async () => {
        t2.setAll(true);
        expect(t2.toBinaryString(), '11111');
    });

    test('similarityTo method', async () => {
        let t3 = DG.BitSet.fromString('11111');
        expect(t2.similarityTo(t3).toString(), '0.6');
    });

    test('toBinaryString method', async () => {
        expect(t2.toBinaryString(), '11001');
    });

    test('toString method', async () => {
        expect(t2.toString(), '5 bits, 3 set');
    });

    test('xor method', async () => {
        let t3 = DG.BitSet.fromString('11111');
        expect(t2.xor(t3).toBinaryString(), '00110');
    });
});
