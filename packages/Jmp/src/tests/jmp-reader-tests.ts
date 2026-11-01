import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {jmpSecondsToDate, readHeader, readJmp, readMetadata} from '../formats/jmp/jmp-parser';
import {ByteReader} from '../formats/jmp/jmp-reader';
import {JmpFormatError, JmpTable, SUPPORTED_VERSION, jmpValues} from '../formats/jmp/jmp-types';
import {readColumnProperties} from '../formats/jmp/jmp-properties';
import {compareToGolden, Golden} from './golden';


// Fixtures vendored from JMPReader.jl/test (MIT); see CREDITS.md. files/golden holds what JMPReader.jl
// reads from each (tools/jmp2json.jl).
export const FIXTURES = [
  'bugMWE4', 'byteintegers', 'byteintegers_notcompressed', 'byteintegers_withmissing',
  'compact', 'compact_UInt8', 'compact_UInt16', 'compact_UInt32', 'compressed',
  'currencies', 'date', 'duration', 'example1', 'geographic', 'longcolumnnames',
  'minusfour', 'rowstate', 'singlecolumnsinglerow', 'time',
];

export async function fixtureBytes(name: string): Promise<Uint8Array> {
  return await _package.files.readAsBytes(`jmpreader/${name}.jmp`);
}

export async function fixture(name: string): Promise<JmpTable> {
  return await readJmp(await fixtureBytes(name));
}

function values(table: JmpTable, column: string): (number | string | null)[] {
  const col = table.columns.find((c) => c.info.name === column);
  if (!col)
    throw new Error(`no column "${column}"`);
  return jmpValues(col, table.info.nrows);
}

const isoDate = (v: number | string | null): string | null =>
  v === null ? null : jmpSecondsToDate(v as number).toISOString().slice(0, 10);
const isoDateTime = (v: number | string | null): string | null =>
  v === null ? null : jmpSecondsToDate(v as number).toISOString().slice(0, 19);

category('JMP: header', () => {
  for (const name of FIXTURES) {
    test(name, async () => {
      const header = readHeader(new ByteReader(await fixtureBytes(name)));
      expect(header.version.join(','), SUPPORTED_VERSION.join(','));
    });
  }

  test('rejects non-JMP bytes', async () => {
    await expectExceptionAsync(async () => {
      readHeader(new ByteReader(new TextEncoder().encode('a,b,c\n1,2,3\n')));
    }, (e) => e instanceof JmpFormatError);
  });

  test('rejects old format', async () => {
    const bytes = new Uint8Array(16);
    bytes.set([0xff, 0xff, 0x00, 0x00, 3, 0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 0]);
    await expectExceptionAsync(async () => {
      readHeader(new ByteReader(bytes));
    }, (e) => e instanceof JmpFormatError && /older format/.test(e.message));
    // classic Mac JMP wrote it big-endian
    bytes.set([0x00, 0x00, 0xff, 0xff, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0, 3]);
    await expectExceptionAsync(async () => {
      readHeader(new ByteReader(bytes));
    }, (e) => e instanceof JmpFormatError && /older format/.test(e.message));
  });

  test('rejects truncated header', async () => {
    const bytes = (await fixtureBytes('example1')).subarray(0, 10);
    await expectExceptionAsync(async () => {
      readHeader(new ByteReader(bytes));
    }, (e) => e instanceof JmpFormatError);
  });
}, {node: true});

// Every fixture against the reference output of JMPReader.jl (JMP_PLUGIN.md §9.3).
category('JMP: reference', () => {
  for (const name of FIXTURES) {
    test(name, async () => {
      const golden = JSON.parse(await _package.files.readAsText(`golden/${name}.json`)) as Golden;
      const diffs = compareToGolden(await fixture(name), golden);
      expect(diffs.length, 0, diffs.slice(0, 10).join('; '));
    });
  }
}, {node: true});

// Ported from JMPReader.jl test/runtests.jl, so a failure reads as the case it breaks.
category('JMP: tables', () => {
  test('example1', async () => {
    const t = await fixture('example1');
    expect(t.info.nrows, 4);
    expectArray(t.columns.map((c) => c.info.name), [
      'ints', 'floats', 'charconstwidth', 'time', 'date', 'duration', 'charconstwidth2',
      'charvariable16', 'formula', 'pressures', 'char utf8', 'charvariable8',
    ]);
    expectArray(values(t, 'ints'), [1, 2, 3, 4]);
    expectArray(values(t, 'floats'), [11.1, 22.2, 33.3, 44.4]);
    expectArray(values(t, 'charconstwidth'), ['a', 'b', 'c', 'd']);
    expectArray(values(t, 'time').map(isoDateTime),
      ['1976-04-01T21:12:00', '1984-08-06T23:58:00', '2003-06-02T17:00:00', null]);
    expectArray(values(t, 'date').map(isoDate), ['2024-01-13', '2024-01-14', null, '2032-02-12']);
    expectArray(values(t, 'duration'), [2322, 364, 229, 0]);
    expectArray(values(t, 'charconstwidth2'), ['a', 'bb', 'ccc', 'dddd']);
    expect(values(t, 'charvariable16')[3], 'abcdefghij'.repeat(13));
    expectArray(values(t, 'formula'), ['2', '4', '6', '8']);
    expectArray(values(t, 'pressures'), [101.325, null, 2.6, 4.63309e110]);
    expectArray(values(t, 'char utf8'), ['ꙮꙮꙮ', '🚴💨', 'jäääär', '辛口']);
    expectArray(values(t, 'charvariable8'), ['a', 'bb', 'cc', 'abcdefghijkl']);
  });

  test('compressed', async () => {
    const t = await fixture('compressed');
    expect(t.info.nrows, 10000);
    expect(values(t, 'numeric').every((v) => v === 1), true);
    expect(values(t, 'character130').every((v) => v === 'abcdefghij'.repeat(13)), true);
    expect(values(t, 'min:s').every((v) => v === 196), true);
    expect(isoDateTime(values(t, 'y-m-d h:m:s')[0]), '1904-01-01T00:00:01');
  });

  test('date', async () => {
    expectArray(values(await fixture('date'), 'ddmmyyyy').map(isoDate),
      ['2011-05-25', '1973-05-24', '2027-05-22', '2020-05-01']);
  });

  test('duration', async () => {
    expectArray(values(await fixture('duration'), ':day:hr:m:s'), [88201, 88201, 88201]);
  });

  test('time', async () => {
    const t = await fixture('time');
    for (const c of t.columns.filter((c) => !c.info.name.startsWith(':day'))) {
      expect(c.kind, 'datetime', c.info.name);
      expectArray(jmpValues(c, t.info.nrows).map(isoDateTime), ['1914-04-27T19:54:14', '1978-01-07T06:11:24']);
    }
  });

  test('longcolumnnames', async () => {
    const t = await fixture('longcolumnnames');
    const digits = (n: number) => Array.from({length: n}, (_, i) => String(i + 1).padStart(10, '0')).join('');
    expectArray(t.columns.map((c) => c.info.name), [digits(140), digits(280)]);
  });

  test('singlecolumnsinglerow', async () => {
    const t = await fixture('singlecolumnsinglerow');
    expect(t.info.nrows, 1);
    expectArray(values(t, 'Column 1'), [1]);
  });

  test('byteintegers', async () => {
    const t = await fixture('byteintegers_notcompressed');
    expectArray(values(t, 'onebyte'), [1, 2, -126, 127]);
    expectArray(values(t, 'twobyte'), [32767, null, 0, -32766]);
    expectArray(values(t, 'fourbyte'), [null, 2147483647, -2147483646, null]);
    expectArray(values(t, 'numeric'), [null, 2147483648, null, null]);
    expectArray(values(await fixture('byteintegers_withmissing'), '1-byte integer'), [0, 1, 0, 1, 0, null]);
  });

  test('geographic', async () => {
    const t = await fixture('geographic');
    expect(t.columns.find((c) => c.info.name === 'Latitude_DDD')!.kind, 'latitude');
    expect(t.columns.find((c) => c.info.name === 'Longitude_DDD')!.kind, 'longitude');
    expect(values(t, 'Latitude_DDD')[1], null);
  });

  test('currencies', async () => {
    const t = await fixture('currencies');
    expectArray(values(t, 'COP'), [3.14, 2.78, 1.41]);
    expect(t.columns.every((c) => c.kind === 'currency'), true);
  });

  test('rowstate', async () => {
    const t = await fixture('rowstate');
    expect((values(t, 'rowstate3')[1] as string).split(' ')[0], '▲');
    expect((values(t, 'rowstate3')[2] as string).split(' ')[0], 'ꙮ');
    expect((values(t, 'rowstate2')[2] as string).split(' ')[1], '#C0C0C0');
  });

  test('compact', async () => {
    const t = await fixture('compact');
    expectArray(values(t, 'compactsubtype'), ['aa', 'b', 'ccc', 'dd', 'dd']);
    expectArray(values(t, 'longcompact').map((v) => (v as string).length), [254, 255, 0, 256, 256]);
    expect(values(await fixture('compact_UInt32'), 'A')[65536], 'a65537');
  });

  test('bugMWE4', async () => {
    const t = await fixture('bugMWE4');
    expectArray(values(t, 'c').slice(0, 3), ['cat', 'cat', 'dolor sit amet']);
    expectArray(values(t, 'f').slice(15, 18), ['foo', 'bar', 'hello world']);
  });

  // No fixture has a fixed-decimal column; this follows JMPReader.jl's note that dt5 = 0x5e or 0x63 is fixed
  // decimal with dt3 = width and dt4 = decimals, applied to example1's "floats" column.
  test('fixed decimals', async () => {
    const bytes = (await fixtureBytes('example1')).slice();
    const r = new ByteReader(bytes);
    const info = readMetadata(r, readHeader(r));
    const floats = info.columns.find((c) => c.name === 'floats')!;
    const dt = floats.offset + 2 + new TextEncoder().encode('floats').length;
    bytes.set([0x01, 0x00, 0x0a, 0x02, 0x5e, 0x08], dt);
    const col = (await readJmp(bytes)).columns.find((c) => c.info.name === 'floats')!;
    expect(col.kind, 'float');
    expect(col.fixedDecimals, 2);
    expect((await fixture('example1')).columns.find((c) => c.info.name === 'floats')!.fixedDecimals === undefined,
      true);
  });

  test('currency codes', async () => {
    const t = await fixture('currencies');
    expectArray(t.columns.map((c) => c.currency!), t.columns.map((c) => c.info.name));
  });

  test('hidden-column index', async () => {
    const t = await fixture('byteintegers');
    expectArray(t.info.visibleIndex, [1, 2, 3]);
    expectArray(t.info.hiddenIndex, [0]);
  });

  test('metadata', async () => {
    const t = await fixture('example1');
    expect(t.info.jmpVersion, '17.2.0');
    expect(t.info.charset, 'utf-8');
    expect(t.info.savetime.toISOString(), '2024-03-02T08:06:26.000Z');
  });
}, {node: true});

category('JMP: robustness', () => {
  test('unsupported column comes back empty, the rest loads', async () => {
    const t = await readJmp(await _package.files.readAsBytes('tests/unsupported-column.jmp'));
    const bad = t.columns.filter((c) => c.decodeError);
    expectArray(bad.map((c) => c.info.name), ['formula']);
    expect(bad[0].data.type, 'none');
    expectArray(values(t, 'char utf8'), ['ꙮꙮꙮ', '🚴💨', 'jäääär', '辛口']);
  });

  test('truncated file', async () => {
    await expectExceptionAsync(async () => {
      await readJmp(await _package.files.readAsBytes('tests/truncated.jmp'));
    }, (e) => e instanceof JmpFormatError && /truncated/.test(e.message));
  });

  test('not a JMP file', async () => {
    await expectExceptionAsync(async () => {
      await readJmp(await _package.files.readAsBytes('tests/not-a-jmp.jmp'));
    }, (e) => e instanceof JmpFormatError);
  });

  // Mutated fixtures either read or fail with a JmpFormatError; nothing hangs or throws anything else.
  test('fuzz', async () => {
    let seed = 12345;
    const rnd = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    for (const name of ['example1', 'compressed', 'compact', 'byteintegers_notcompressed', 'rowstate', 'date',
      'compact_UInt8']) {
      const orig = await fixtureBytes(name);
      for (let k = 0; k < 60; k++) {
        let b = orig.slice();
        if (k % 2 === 0) {
          for (let j = 0; j < 1 + Math.floor(rnd() * 8); j++)
            b[Math.floor(rnd() * b.length)] = Math.floor(rnd() * 256);
        } else
          b = b.subarray(0, 16 + Math.floor(rnd() * (b.length - 16)));
        try {
          await readJmp(b);
        } catch (e) {
          if (!(e instanceof JmpFormatError))
            throw new Error(`${name} mutation ${k}: ${e instanceof Error ? e.name + ': ' + e.message : e}`);
        }
      }
    }
  }, {timeout: 120000});
}, {node: true});

const enc = new TextEncoder();

function concat(parts: (Uint8Array | number[])[]): Uint8Array {
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of parts) {
    out.set(p, at);
    at += p.length;
  }
  return out;
}

const u16 = (v: number): number[] => [v & 0xff, v >> 8];
const u32 = (v: number): number[] => [v & 0xff, (v >> 8) & 0xff, (v >> 16) & 0xff, v >>> 24];
const u64 = (v: number): number[] => [...u32(v), 0, 0, 0, 0];
const f64 = (v: number): Uint8Array => new Uint8Array(new Float64Array([v]).buffer);
const nulText = (items: string[]): Uint8Array => enc.encode(items.map((s) => s + '\0').join(''));

/** A column header after the descriptor bytes, laid out as JMP_PLUGIN.md §4.6 describes. */
function columnHeader(gap: number, entries: [number, Uint8Array][]): Uint8Array {
  const body = concat(entries.flatMap(([id, p]) => [u16(id), u64(p.length), p]));
  return concat([new Array(gap).fill(0), u16(entries.length), u64(body.length), body]);
}

function jslProperties(texts: string[]): Uint8Array {
  return concat([...texts.flatMap((t) => [u32(enc.encode(t).length), enc.encode(t)]), u32(0xffffffff)]);
}

function labelsPayload(keys: Uint8Array, labels: string[]): Uint8Array {
  const l = nulText(labels);
  return concat([u32(labels.length), u32(keys.length), keys, u32(labels.length), u32(l.length), l]);
}

category('JMP: column properties', () => {
  test('formula and units in a JMP-made file', async () => {
    const t = await fixture('example1');
    const props = (name: string): string[] =>
      t.columns.find((c) => c.info.name === name)!.properties!.properties.map((p) => p.text);
    expectArray(props('formula'), ['Formula(:ints * 2)']);
    expectArray(props('pressures'), ['Units("kPa")']);
    expectArray(props('ints'), []);
  });

  test('JSL properties, in the order saved', async () => {
    const texts = ['Notes("Was the offer accepted?")', 'Value Order({Custom Order({"Yes", "No"})})',
      '"Value Ordering"({"Yes", "No"})'];
    const p = readColumnProperties(columnHeader(4, [[0x07, jslProperties(texts)]]), 0, 0x02);
    expectArray(p.properties.map((q) => q.name), ['Notes', 'Value Order', 'Value Ordering']);
    expect(p.properties.every((q) => q.node !== undefined), true);
  });

  test('value labels, numeric and text', async () => {
    const numeric = labelsPayload(concat([f64(0), f64(1), f64(2)]), ['Training', 'Validation', 'Test']);
    // a numeric column's header has one byte more before the property list
    const n = readColumnProperties(columnHeader(5, [[0x0e, numeric]]), 0, 0x01);
    expectArray(n.valueLabels!.keys, [0, 1, 2]);
    expectArray(n.valueLabels!.labels, ['Training', 'Validation', 'Test']);
    const text = readColumnProperties(columnHeader(4, [[0x0e, labelsPayload(nulText(['f', 'm']), ['Fri', 'Mon'])]]),
      0, 0x02);
    expectArray(text.valueLabels!.keys, ['f', 'm']);
  });

  test('a header of another shape has none', async () => {
    expectArray(readColumnProperties(new Uint8Array(32), 0, 0x01).properties, []);
    const bad = columnHeader(4, [[0x07, jslProperties(['Notes("x")'])]]);
    bad[6]++; // the total length no longer adds up
    expectArray(readColumnProperties(bad, 0, 0x02).properties, []);
    const torn = labelsPayload(f64(1), ['One']).subarray(0, 20);
    expect(readColumnProperties(columnHeader(5, [[0x0e, torn]]), 0, 0x01).valueLabels === undefined, true);
  });
}, {node: true});
