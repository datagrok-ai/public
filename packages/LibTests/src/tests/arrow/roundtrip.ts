// Per-type Arrow IPC ↔ DG.DataFrame round-trip suite.
//
// Tests go through the Arrow package's registered functions
// (`Arrow:toFeather` / `Arrow:fromFeather`) via `grok.functions.call` to
// exercise the public function-registry path independent of how the lib's
// exports are spelled.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

async function toFeather(table: DG.DataFrame): Promise<Uint8Array> {
  return await grok.functions.call('Arrow:toFeather', {table, asStream: true}) as Uint8Array;
}

async function fromFeather(bytes: Uint8Array): Promise<DG.DataFrame> {
  return await grok.functions.call('Arrow:fromFeather', {bytes}) as DG.DataFrame;
}

async function roundtrip(df: DG.DataFrame): Promise<DG.DataFrame> {
  return fromFeather(await toFeather(df));
}

category('Arrow / IPC round-trip', () => {
  test('roundtrip_int_column', async () => {
    const col = DG.Column.fromInt32Array('vals', new Int32Array([1, 2, 3, -7, 0, 999]));
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    expect(back.rowCount, df.rowCount, 'row count differs');
    expect(back.columns.length, 1, 'column count differs');
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.INT, 'type differs');
    for (let i = 0; i < df.rowCount; ++i)
      expect(c2.get(i), col.get(i), `row ${i} differs`);
  });

  test('roundtrip_float_column', async () => {
    const col = DG.Column.fromFloat64Array('vals', new Float64Array([1.5, -2.25, 3.14159, 0, 1e10]));
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.FLOAT, 'type differs');
    for (let i = 0; i < df.rowCount; ++i)
      expectFloat(c2.get(i) as number, col.get(i) as number, 1e-12);
  });

  test('roundtrip_string_column', async () => {
    const col = DG.Column.fromStrings('vals', ['alpha', 'beta', 'gamma', '', 'alpha', 'beta']);
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.STRING, 'type differs');
    for (let i = 0; i < df.rowCount; ++i)
      expect(c2.get(i), col.get(i), `row ${i} differs`);
  });

  test('roundtrip_bigint_column', async () => {
    // Values straddling the Int32 boundary exercise convertInt64Column's
    // overflow detection branch (api.ts:219-227). BigInt() constructor is
    // used because LibTests' tsconfig targets es6 (no bigint literals).
    // 2^40 (above Int32 max) — built via repeated multiplication because
    // tsconfig target es6 doesn't allow ** on bigints.
    let big = BigInt(1);
    for (let i = 0; i < 40; ++i) big = big * BigInt(2);
    const arr = new BigInt64Array([BigInt(1), big, -big, BigInt(0)]);
    const col = DG.Column.fromBigInt64Array('vals', arr);
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.BIG_INT, `expected bigint, got ${c2.type}`);
    for (let i = 0; i < df.rowCount; ++i)
      expect(String(c2.get(i)), String(col.get(i)), `row ${i} differs`);
  });

  test('roundtrip_datetime_column', async () => {
    // DG.Column.fromList for DATE_TIME accepts Date objects directly; the
    // internal representation is a microsecond-resolution float in raw data.
    const dates = [
      new Date('2024-01-15T12:00:00Z'),
      new Date('2025-06-30T08:30:45Z'),
      new Date('1970-01-01T00:00:00Z'),
    ];
    const col = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'vals', dates as any);
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.DATE_TIME, 'type differs');
    for (let i = 0; i < df.rowCount; ++i) {
      // Both sides serialize to the same UTC ms via valueOf.
      const a = (c2.get(i) as any)?.valueOf();
      const b = (col.get(i) as any)?.valueOf();
      expect(a, b, `row ${i} differs`);
    }
  });

  test('roundtrip_bool_column', async () => {
    const col = DG.Column.fromBitSet('vals', DG.BitSet.create(8, (i) => i % 2 === 0));
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.BOOL, 'type differs');
    expect(c2.length, col.length, 'length differs');
    for (let i = 0; i < col.length; ++i)
      expect(c2.get(i), col.get(i), `row ${i} differs`);
  });

  test('roundtrip_dictionary_strings', async () => {
    // Repeated values exercise the dictionary unpack/reconstruct paths in
    // unpackDictionaryColumn / stringColumnFromDictionary (api.ts:171, 199).
    const values = ['a', 'b', 'a', 'c', 'b', 'a', 'c', 'a'];
    const col = DG.Column.fromStrings('vals', values);
    const df = DG.DataFrame.fromColumns([col]);
    const back = await roundtrip(df);
    const c2 = back.col('vals')!;
    expect(c2.type, DG.COLUMN_TYPE.STRING, 'type differs');
    for (let i = 0; i < values.length; ++i)
      expect(c2.get(i), values[i], `row ${i} differs`);
  });

  test('roundtrip_multi_column', async () => {
    // A mixed-type DataFrame that mirrors the Phase 3 worker payload shape.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('id', new Int32Array([10, 20, 30])),
      DG.Column.fromFloat64Array('value', new Float64Array([1.1, 2.2, 3.3])),
      DG.Column.fromStrings('label', ['x', 'y', 'z']),
    ]);
    const back = await roundtrip(df);
    expect(back.rowCount, df.rowCount, 'row count');
    expect(back.columns.length, df.columns.length, 'column count');
    for (const name of ['id', 'value', 'label'])
      expect(back.col(name) !== null, true, `missing column ${name}`);
    expectDeepEqual(
      back.col('id')!.toList(),
      df.col('id')!.toList(),
    );
    expectDeepEqual(
      Array.from(back.col('value')!.toList()).map((v: any) => Number(v)),
      Array.from(df.col('value')!.toList()).map((v: any) => Number(v)),
      {floatTolerance: 1e-12},
    );
    expectDeepEqual(back.col('label')!.toList(), df.col('label')!.toList());
  });

  test('roundtrip_byte_array_serialization', async () => {
    // Mirrors the existing Arrow-package "serialization" assertion: round-trip
    // through DG.DataFrame.toByteArray / fromByteArray after a Feather pass.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('a', new Int32Array([1, 2, 3])),
      DG.Column.fromStrings('b', ['p', 'q', 'r']),
    ]);
    const after = await roundtrip(df);
    const reconstructed = DG.DataFrame.fromByteArray(after.toByteArray());
    expect(reconstructed.rowCount, df.rowCount);
    expect(reconstructed.columns.length, df.columns.length);
  });
});
