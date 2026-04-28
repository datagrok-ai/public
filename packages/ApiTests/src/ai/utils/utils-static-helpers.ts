import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

category('AI: Utils: static helpers', () => {
  test('firstOrNull on empty iterable', async () => {
    expect(DG.Utils.firstOrNull<number>([]), null);
    expect(DG.Utils.firstOrNull<string>(new Set<string>()), null);
  });

  test('firstOrNull on array', async () => {
    expect(DG.Utils.firstOrNull<number>([10, 20, 30]), 10);
  });

  test('firstOrNull on set preserves insertion order', async () => {
    expect(DG.Utils.firstOrNull<string>(new Set(['a', 'b'])), 'a');
  });

  test('identity zero and positive', async () => {
    const zero = DG.Utils.identity(0);
    expect(zero instanceof Uint32Array, true);
    expect(zero.length, 0);
    const five = DG.Utils.identity(5);
    expect(five instanceof Uint32Array, true);
    expectArray(Array.from(five), [0, 1, 2, 3, 4]);
  });

  test('replaceAll replaces every occurrence', async () => {
    expect(DG.Utils.replaceAll('a-b-a-c-a', 'a', 'X'), 'X-b-X-c-X');
    expect(DG.Utils.replaceAll('a-b-a', 'a', ''), '-b-');
    expect(DG.Utils.replaceAll('hello', 'z', 'Q'), 'hello');
  });

  test('isEmpty handles null undefined empty', async () => {
    expect(DG.Utils.isEmpty(undefined), true);
    expect(DG.Utils.isEmpty(null as any), true);
    expect(DG.Utils.isEmpty(''), true);
    expect(DG.Utils.isEmpty(' '), false);
    expect(DG.Utils.isEmpty('x'), false);
  });

  test('nullIfEmpty returns null or original', async () => {
    expect(DG.Utils.nullIfEmpty(''), null);
    expect(DG.Utils.nullIfEmpty(undefined), null);
    expect(DG.Utils.nullIfEmpty(null as any), null);
    expect(DG.Utils.nullIfEmpty(' '), ' ');
    expect(DG.Utils.nullIfEmpty('hi'), 'hi');
  });

  test('getJsonValueType known and falsy and throw', async () => {
    expect(DG.Utils.getJsonValueType('hi'), DG.TYPE.STRING);
    expect(DG.Utils.getJsonValueType(42), DG.TYPE.FLOAT);
    expect(DG.Utils.getJsonValueType(-1.5), DG.TYPE.FLOAT);
    expect(DG.Utils.getJsonValueType(true), DG.TYPE.BOOL);
    for (var falsy of [null, undefined, '', 0, false])
      expect(DG.Utils.getJsonValueType(falsy), null);
    let threw = false;
    try {
      DG.Utils.getJsonValueType({});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, true);
  });

  test('jsonToColumns expands string column', async () => {
    const df = DG.DataFrame.create(2);
    const json = DG.Column.fromStrings('j', ['{"a":1,"b":"x"}', '{"a":2,"b":"y"}']);
    df.columns.add(json);
    DG.Utils.jsonToColumns(json);
    expect(df.columns.contains('a'), true);
    expect(df.columns.contains('b'), true);
    expect(df.col('a')!.type, DG.COLUMN_TYPE.FLOAT);
    expect(df.col('b')!.type, DG.COLUMN_TYPE.STRING);
    expectArray([df.col('a')!.get(0), df.col('a')!.get(1)], [1, 2]);
    expectArray([df.col('b')!.get(0), df.col('b')!.get(1)], ['x', 'y']);
  });

  test('randomString length and alphabet', async () => {
    expect(DG.Utils.randomString(0), '');
    const s = DG.Utils.randomString(20);
    expect(s.length, 20);
    expect(/^[A-Za-z0-9]+$/.test(s), true);
    const t = DG.Utils.randomString(20);
    expect(t.length, 20);
    expect(/^[A-Za-z0-9]+$/.test(t), true);
  });
});
