import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {jmpSecondsToDate, readHeader, readJmp} from '../formats/jmp/jmp-parser';
import {ByteReader} from '../formats/jmp/jmp-reader';
import {JmpFormatError, JmpTable, SUPPORTED_VERSION} from '../formats/jmp/jmp-types';
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
  return col.values;
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
      expectArray(c.values.map(isoDateTime), ['1914-04-27T19:54:14', '1978-01-07T06:11:24']);
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
    expect(bad[0].values.every((v) => v === null), true);
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
