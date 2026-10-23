// Ported from JMPReader.jl (MIT, (c) 2024 Jaakko Ruohio), src/JMPReader.jl, src/metadata.jl
// and src/column.jl, upstream commit 9fd68aa2651a8dde002ecdbd06af4d8d832b18f2
// See ../../../docs/JMP_PLUGIN.md §4 for the format as currently understood.

import {ByteReader} from './jmp-reader';
import {gunzip} from './jmp-decompress';
import {readTableScripts} from './jmp-scripts';
import {
  CURRENCY_FORMATS, DATE_FORMATS, DATETIME_FORMATS, DURATION_FORMATS, GZIP_SECTION_START,
  JMP_EPOCH_OFFSET_S, JMP_INT_MISSING, JmpColumn, JmpColumnInfo, JmpColumnKind, JmpData, JmpFormatError, JmpHeader,
  JmpInfo, JmpTable, LATITUDE_FORMATS, LONGITUDE_FORMATS, MAGIC_JMP, OLD_VERSION, ROWSTATE_COLORS,
  SUPPORTED_VERSION,
} from './jmp-types';

const METADATA_ANCHOR = [0x07, 0x00, 0x08, 0x00, 0x00, 0x00];
const OPTIONAL_BLOCK_TAGS = [0xfc, 0xfd, 0xfe, 0xff];
const PLAIN_NUMERIC_FORMATS = new Set([0x00, 0x03, 0x42, 0x43, 0x44, 0x59, 0x60, 0x63]);
const FIXED_DECIMAL_INPUT_FORMATS = new Set([0x5e, 0x63]);
const MAX_GZIP_RATIO = 1100;
const MAX_FIXED_DECIMALS = 30;
const YIELD_ROWS = 100000;
/** In a currency column's header: a 0x0013 tag, a u32 length of 5, then the ISO code (JMP_PLUGIN.md §4.3). */
const CURRENCY_TAG = [0x13, 0x00, 0x05, 0x00, 0x00, 0x00];
const CURRENCY_SEARCH_BYTES = 256;

const utf8 = new TextDecoder('utf-8');

/** Reads and validates the magic bytes and version triple. */
export function readHeader(r: ByteReader): JmpHeader {
  r.seek(0);
  if (r.length === 0)
    throw new JmpFormatError('The file is empty');
  if (r.length < MAGIC_JMP.length + 12 || MAGIC_JMP.some((b, i) => r.bytes[i] !== b))
    throw new JmpFormatError('Not a JMP data table');
  r.skip(MAGIC_JMP.length);
  const version: [number, number, number] = [r.u32(), r.u32(), r.u32()];
  if (sameVersion(version, OLD_VERSION)) {
    throw new JmpFormatError('This JMP file uses an older format that isn\'t supported yet. ' +
      'Re-save it in a current JMP version.');
  }
  if (!sameVersion(version, SUPPORTED_VERSION))
    throw new JmpFormatError(`Unrecognized JMP file version (${version.map((v) => '0x' + v.toString(16)).join(', ')})`);
  return {version};
}

/** Reads table metadata: row/column counts, build string, save time, column offsets and names. */
export function readMetadata(r: ByteReader, header: JmpHeader): JmpInfo {
  const anchor = r.indexAfter(METADATA_ANCHOR, 0);
  if (anchor < 0)
    throw new JmpFormatError('Could not find the table metadata');
  r.seek(anchor - 38);

  const nrows = r.i64();
  const ncols = r.i32();
  // gzip can't expand more than ~1032:1, so a row count far beyond the file size is corruption
  if (nrows < 0 || ncols < 0 || nrows > r.length * MAX_GZIP_RATIO || ncols > r.length)
    throw new JmpFormatError(`Invalid table size (${nrows} rows, ${ncols} columns)`);
  r.skip(5 * 2);
  const charset = r.lengthPrefixedString(4).replace(/\0+$/, '');
  r.skip(3 * 2);
  const savetime = jmpSecondsToDate(r.f64());
  r.skip(2);
  const buildString = r.lengthPrefixedString(4);
  const propertiesStart = r.pos;
  const m = /Version (.*)$/.exec(buildString);
  if (!m)
    throw new JmpFormatError(`Could not determine JMP version from "${buildString}"`);

  const {nVisible, nHidden} = seekToColumnIndex(r, ncols);
  const visibleIdx: number[] = [];
  for (let i = 0; i < nVisible; i++)
    visibleIdx.push(r.u32());
  const hiddenIdx: number[] = [];
  for (let i = 0; i < nHidden; i++)
    hiddenIdx.push(r.u32());
  const widths: number[] = [];
  for (let i = 0; i < ncols; i++)
    widths.push(r.u16());
  r.skip(7 * 4);

  const offsets = readColumnOffsets(r, ncols);
  const columns: JmpColumnInfo[] = offsets.map((offset, i) => {
    r.seek(offset);
    return {name: r.lengthPrefixedString(2), offset, width: widths[i]};
  });

  const propertiesEnd = offsets.reduce((lo, o) => Math.min(lo, o), r.length);
  return {header, jmpVersion: m[1], buildString, charset, savetime, nrows, ncols, columns,
    visibleIndex: visibleIdx, hiddenIndex: hiddenIdx,
    propertiesRange: [propertiesStart, Math.max(propertiesStart, propertiesEnd)]};
}

/** Brute-force search for the visible/hidden column index (metadata.jl `seek_to_column_data_offsets`). */
function seekToColumnIndex(r: ByteReader, ncols: number): {nVisible: number, nHidden: number} {
  let pos = 2;
  for (;;) {
    pos = r.indexAfter([0xff, 0xff], pos);
    if (pos < 0)
      throw new JmpFormatError('Could not find the column index');
    while (pos < r.length && r.bytes[pos] === 0xff)
      pos++;
    if (pos + 26 > r.length)
      throw new JmpFormatError('Could not find the column index');
    r.seek(pos + 10);
    const nVisible = r.u32();
    const nHidden = r.u32();
    r.skip(8);
    if (nVisible + nHidden === ncols)
      return {nVisible, nHidden};
    pos += 8;
  }
}

/** Skips optional blocks and reads the column offsets (metadata.jl `column_info`). */
function readColumnOffsets(r: ByteReader, ncols: number): number[] {
  // Upstream notes this is a heuristic: a column count could also match these two bytes.
  while (r.pos + 2 <= r.length && r.bytes[r.pos + 1] === 0xff && OPTIONAL_BLOCK_TAGS.includes(r.bytes[r.pos])) {
    r.skip(2);
    const n = r.i64();
    if (n < 0)
      throw new JmpFormatError(`Invalid block length ${n} at ${r.pos - 8}`);
    r.skip(n);
  }
  const ncols2 = r.i32();
  if (ncols2 !== ncols)
    throw new JmpFormatError(`Column counts in the file disagree (${ncols} vs ${ncols2})`);
  const offsets: number[] = [];
  for (let i = 0; i < ncols; i++) {
    const offset = r.i64();
    if (offset < 0 || offset >= r.length)
      throw new JmpFormatError(`Column ${i + 1} starts outside the file; file appears truncated`);
    offsets.push(offset);
  }
  return offsets;
}

/** Decodes the `i`th column (column.jl `column_data`). Never throws for an unsupported layout. */
export async function readColumn(r: ByteReader, info: JmpInfo, i: number): Promise<JmpColumn> {
  const col = info.columns[i];
  r.seek(col.offset);
  r.lengthPrefixedString(2);
  const dt = Array.from(r.bytesAt(6));
  const mark = r.pos;
  const nrows = info.nrows;

  let a: Uint8Array;
  let colEnd = i === info.ncols - 1 ? r.length : info.columns[i + 1].offset;
  if (dt[0] === 0x09 || dt[0] === 0x0a) {
    const start = r.indexAfter(GZIP_SECTION_START, mark);
    if (start < 0)
      throw new JmpFormatError(`Column "${col.name}": compressed section not found`);
    r.seek(start);
    const gzipLen = r.u64();
    const gunzipLen = r.u64();
    a = await gunzip(r.bytesAt(gzipLen), gunzipLen);
  } else {
    if (colEnd < col.offset || colEnd > r.length)
      colEnd = r.length;
    a = r.bytes.subarray(col.offset, colEnd);
  }
  r.seek(mark);

  const done = (kind: JmpColumnKind, data: JmpData): JmpColumn => ({info: col, kind, descriptor: dt, data});
  const [dt1, dt2, dt3, dt4, dt5, dt6] = dt;

  if (dt1 === 0x01 || dt1 === 0x0a) {
    const data = readNumbers(a, nrows, dt6, col.name);
    if ((dt4 === dt5 && PLAIN_NUMERIC_FORMATS.has(dt4)) || FIXED_DECIMAL_INPUT_FORMATS.has(dt5)) {
      const c = done(dt6 === 8 ? 'float' : 'int', data);
      // fixed decimal: dt3 = width, dt4 = decimals (JMPReader.jl column.jl, not confirmed by a fixture)
      if (dt4 !== dt5 && FIXED_DECIMAL_INPUT_FORMATS.has(dt5) && dt4 <= MAX_FIXED_DECIMALS)
        c.fixedDecimals = dt4;
      return c;
    }
    if (dt4 === dt5 && CURRENCY_FORMATS.has(dt4))
      return {...done('currency', data), currency: readCurrencyCode(r.bytes, mark, colEnd)};
    if (dt4 === dt5 && LONGITUDE_FORMATS.has(dt4))
      return done('longitude', data);
    if (dt4 === dt5 && LATITUDE_FORMATS.has(dt4))
      return done('latitude', data);
    if (DATE_FORMATS.has(dt4))
      return done('date', data);
    if (DATETIME_FORMATS.has(dt4))
      return done('datetime', data);
    if (DURATION_FORMATS.has(dt4))
      return done('duration', data);
  }

  if (dt1 === 0xff || dt1 === 0xfe || dt1 === 0xfc)
    return done(dt5 === 8 ? 'float' : 'int', readNumbers(a, nrows, dt5, col.name));

  if (dt1 === 0x09 && dt2 === 0x03)
    return done('rowstate', readRowStates(a, nrows, dt5, col.name));

  if (dt1 === 0x03 && dt2 === 0x03) {
    const v = new DataView(a.buffer, a.byteOffset, a.byteLength);
    const start = tailStart(a, nrows * dt5, col.name);
    const values = new BigInt64Array(nrows);
    for (let k = 0; k < nrows; k++)
      values[k] = v.getBigInt64(start + k * dt5, true);
    return done('bigint', {type: 'bigint', values});
  }

  if ((dt1 === 0x02 || dt1 === 0x09) && (dt2 === 0x01 || dt2 === 0x02)) {
    if ((dt3 === 0 && dt4 === 0 && dt5 > 0) || (dt3 >= 1 && dt3 <= 7 && dt4 === 0))
      return done('string', readFixedWidthStrings(a, nrows, dt5, col.name));
    if (dt3 === 0 && dt4 === 0 && dt5 === 0) {
      if (dt1 === 0x09)
        return done('string', readCompressedStrings(a, nrows, col.name));
      return done('string', readUncompressedStrings(r, a, nrows, col.name));
    }
  }

  const desc = dt.map((b) => '0x' + b.toString(16).padStart(2, '0')).join(',');
  return {
    info: col, kind: 'unsupported', descriptor: dt, data: {type: 'none'},
    decodeError: `Unsupported column layout (${desc})`,
  };
}

/** Parses a .jmp data table. A column that can't be decoded comes back empty with `decodeError` set. */
export async function readJmp(bytes: Uint8Array): Promise<JmpTable> {
  const r = new ByteReader(bytes);
  const header = readHeader(r);
  const info = readMetadata(r, header);
  const columns: JmpColumn[] = [];
  for (let i = 0; i < info.ncols; i++) {
    // let the UI run between columns of a large table, so no stretch blocks it longer than one column
    if (i > 0 && info.nrows >= YIELD_ROWS)
      await new Promise((resolve) => setTimeout(resolve, 0));
    try {
      columns.push(await readColumn(r, info, i));
    } catch (e) {
      columns.push({
        info: info.columns[i], kind: 'unsupported', descriptor: [], data: {type: 'none'},
        decodeError: e instanceof Error ? e.message : String(e),
      });
    }
  }
  let scripts: JmpTable['scripts'] = [];
  try {
    scripts = readTableScripts(bytes, info);
  } catch {
    // the script layout is provisional (jmp-scripts.ts); a surprise there mustn't cost the data
  }
  return {info, columns, scripts};
}

/** The ISO currency code in a currency column's header, or undefined if it isn't where it was in the fixtures. */
function readCurrencyCode(bytes: Uint8Array, from: number, colEnd: number): string | undefined {
  const end = Math.min(colEnd, from + CURRENCY_SEARCH_BYTES, bytes.length);
  outer: for (let p = from; p + CURRENCY_TAG.length + 3 <= end; p++) {
    for (let j = 0; j < CURRENCY_TAG.length; j++) {
      if (bytes[p + j] !== CURRENCY_TAG[j])
        continue outer;
    }
    const code = String.fromCharCode(...bytes.subarray(p + 6, p + 9));
    if (/^[A-Z]{3}$/.test(code))
      return code;
  }
  return undefined;
}

/** Converts seconds since the JMP epoch (1904-01-01) to a Date. */
export function jmpSecondsToDate(s: number): Date {
  return new Date((s - JMP_EPOCH_OFFSET_S) * 1000);
}

function tailStart(a: Uint8Array, n: number, name: string): number {
  if (n > a.length)
    throw new JmpFormatError(`Column "${name}": expected ${n} bytes of data, found ${a.length}`);
  return a.length - n;
}

/** Reads the trailing `nrows` numbers of `width` bytes. Eight-byte numbers come back as doubles (NaN is
 * JMP's missing value); 1-, 2- and 4-byte integers as Int32, with JMP's sentinels made JMP_INT_MISSING. */
function readNumbers(a: Uint8Array, nrows: number, width: number, name: string): JmpData {
  const start = tailStart(a, nrows * width, name);
  const v = new DataView(a.buffer, a.byteOffset, a.byteLength);
  if (width === 8) {
    const out = new Float64Array(nrows);
    if ((a.byteOffset + start) % 8 === 0)
      out.set(new Float64Array(a.buffer, a.byteOffset + start, nrows));
    else {
      for (let k = 0; k < nrows; k++)
        out[k] = v.getFloat64(start + k * 8, true);
    }
    return {type: 'float', values: out};
  }
  const out = new Int32Array(nrows);
  for (let k = 0; k < nrows; k++) {
    const p = start + k * width;
    let x: number;
    switch (width) {
    case 1: x = v.getInt8(p); out[k] = x === -127 ? JMP_INT_MISSING : x; break;
    case 2: x = v.getInt16(p, true); out[k] = x === -32767 ? JMP_INT_MISSING : x; break;
    case 4: x = v.getInt32(p, true); out[k] = x === -2147483647 ? JMP_INT_MISSING : x; break;
    default: throw new JmpFormatError(`Column "${name}": unknown numeric storage width ${width}`);
    }
  }
  return {type: 'int', values: out};
}

/** Row states: marker and colour decoded; the other bytes (0-based 4 and 7) kept as flags. */
function readRowStates(a: Uint8Array, nrows: number, width: number, name: string): JmpData {
  if (width < 7 || width * nrows > a.length)
    throw new JmpFormatError(`Column "${name}": row states don't fit the column`);
  const markers = new Uint32Array(nrows);
  const colors = new Uint32Array(nrows);
  const flags = new Uint16Array(nrows);
  for (let k = 0; k < nrows; k++) {
    const o = width * k;
    markers[k] = (a[o + 6] << 8) | a[o + 5];
    colors[k] = a[o + 3] === 0xff ? (a[o] << 16) | (a[o + 1] << 8) | a[o + 2] :
      parseInt((ROWSTATE_COLORS[a[o]] ?? '#000000').slice(1), 16);
    flags[k] = a[o + 4] | (width > 7 ? a[o + 7] << 8 : 0);
  }
  return {type: 'rowstate', markers, colors, flags};
}

/** Collects a text column's distinct values as it is read; each row gets an index into them. In the browser
 * this and Column.fromIndexes take 0.4 s for 1.9 M unique strings, where Column.fromList takes 7–22 s
 * (JMP_PLUGIN.md §8). */
class Categories {
  readonly values: string[] = [];
  private readonly index = new Map<string, number>();

  of(s: string): number {
    let i = this.index.get(s);
    if (i === undefined) {
      i = this.values.length;
      this.values.push(s);
      this.index.set(s, i);
    }
    return i;
  }
}

function readFixedWidthStrings(a: Uint8Array, nrows: number, width: number, name: string): JmpData {
  const start = tailStart(a, nrows * width, name);
  const cats = new Categories();
  const indexes = new Int32Array(nrows);
  for (let k = 0; k < nrows; k++) {
    const s = start + k * width;
    let e = s;
    while (e < s + width && a[e] !== 0)
      e++;
    indexes[k] = cats.of(utf8.decode(a.subarray(s, e)));
  }
  return {type: 'string', categories: cats.values, indexes};
}

function sliceStrings(data: Uint8Array, widths: Int32Array, name: string): JmpData {
  let total = 0;
  for (const w of widths) {
    if (w < 0)
      throw new JmpFormatError(`Column "${name}": string data runs past the column`);
    total += w;
  }
  if (total > data.length)
    throw new JmpFormatError(`Column "${name}": string data runs past the column`);
  const cats = new Categories();
  const indexes = new Int32Array(widths.length);
  let p = 0;
  // ASCII text: decode once and cut, since byte offsets are character offsets
  if (isAscii(data, total)) {
    const all = utf8.decode(data.subarray(0, total));
    for (let k = 0; k < widths.length; k++) {
      indexes[k] = cats.of(all.substring(p, p + widths[k]));
      p += widths[k];
    }
  } else {
    for (let k = 0; k < widths.length; k++) {
      indexes[k] = cats.of(utf8.decode(data.subarray(p, p + widths[k])));
      p += widths[k];
    }
  }
  return {type: 'string', categories: cats.values, indexes};
}

function isAscii(data: Uint8Array, n: number): boolean {
  for (let i = 0; i < n; i++) {
    if (data[i] >= 0x80)
      return false;
  }
  return true;
}

/** Reads `n` widths of `width` bytes each (1, 2 or 4) from `r`. */
function readWidths(r: ByteReader, n: number, width: number, name: string): Int32Array {
  const out = new Int32Array(n);
  for (let k = 0; k < n; k++)
    out[k] = width === 1 ? r.i8() : width === 2 ? r.i16() : width === 4 ? r.i32() : badWidth(width, name);
  return out;
}

function readCompressedStrings(a: Uint8Array, nrows: number, name: string): JmpData {
  const r = new ByteReader(a);
  if (a.length >= 8 && r.i64() === a.length) {
    // pooled ("compact") storage: per-row indices into a string pool, 0 meaning the empty string
    r.skip(9);
    r.i64();
    const idxWidth = r.i8();
    const idx = new Int32Array(nrows);
    let npools = 0;
    for (let k = 0; k < nrows; k++) {
      idx[k] = idxWidth === 1 ? r.u8() : idxWidth === 2 ? r.u16() : idxWidth === 4 ? r.u32() :
        badWidth(idxWidth, name);
      if (idx[k] > npools)
        npools = idx[k];
    }
    const lenWidth = r.i8();
    // the pool becomes the categories; a duplicate or an empty entry is folded into the first
    const cats = new Categories();
    const remap = new Int32Array(npools + 1);
    remap[0] = cats.of('');
    for (let k = 1; k <= npools; k++) {
      const n = lenWidth === 1 ? r.i8() : lenWidth === 2 ? r.i16() : badWidth(lenWidth, name);
      if (n < 0)
        throw new JmpFormatError(`Column "${name}": negative string length`);
      remap[k] = cats.of(utf8.decode(r.bytesAt(n)));
    }
    for (let k = 0; k < nrows; k++)
      idx[k] = remap[idx[k]];
    return {type: 'string', categories: cats.values, indexes: idx};
  }
  if (a.length < 13 || (a[8] !== 1 && a[8] !== 2))
    badWidth(a[8], name);
  const r2 = new ByteReader(a);
  r2.seek(13);
  const widths = readWidths(r2, nrows, a[8], name);
  return sliceStrings(a.subarray(r2.pos), widths, name);
}

function readUncompressedStrings(r: ByteReader, a: Uint8Array, nrows: number, name: string): JmpData {
  // r is positioned just after dt1..dt6
  r.skip(6);
  r.skip(r.i64());
  r.skip(2);
  r.skip(r.u32() + 8);
  const widthBytes = r.u8();
  r.u32(); // max width
  const widths = readWidths(r, nrows, widthBytes, name);
  let total = 0;
  for (const w of widths)
    total += w;
  if (total < 0 || total > a.length)
    throw new JmpFormatError(`Column "${name}": string data runs past the column`);
  // `a` spans the column block, which ends where the string data ends
  return sliceStrings(a.subarray(a.length - total), widths, name);
}

function badWidth(w: number, name: string): never {
  throw new JmpFormatError(`Column "${name}": unknown width byte ${w}`);
}

function sameVersion(a: number[], b: number[]): boolean {
  return a.length === b.length && a.every((v, i) => v === b[i]);
}
