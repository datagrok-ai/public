// Ported from JMPReader.jl (MIT, (c) 2024 Jaakko Ruohio), src/JMPReader.jl, src/metadata.jl
// and src/column.jl, upstream commit 9fd68aa2651a8dde002ecdbd06af4d8d832b18f2
// See ../../../docs/JMP_PLUGIN.md §4 for the format as currently understood.

import {ByteReader} from './jmp-reader';
import {gunzip} from './jmp-decompress';
import {readTableScripts} from './jmp-scripts';
import {
  CURRENCY_FORMATS, DATE_FORMATS, DATETIME_FORMATS, DURATION_FORMATS, GZIP_SECTION_START,
  JMP_EPOCH_OFFSET_S, JmpColumn, JmpColumnInfo, JmpColumnKind, JmpFormatError, JmpHeader, JmpInfo,
  JmpTable, LATITUDE_FORMATS, LONGITUDE_FORMATS, MAGIC_JMP, OLD_VERSION, ROWSTATE_COLORS,
  ROWSTATE_MARKERS, SUPPORTED_VERSION,
} from './jmp-types';

const METADATA_ANCHOR = [0x07, 0x00, 0x08, 0x00, 0x00, 0x00];
const OPTIONAL_BLOCK_TAGS = [0xfc, 0xfd, 0xfe, 0xff];
const PLAIN_NUMERIC_FORMATS = new Set([0x00, 0x03, 0x42, 0x43, 0x44, 0x59, 0x60, 0x63]);
const FIXED_DECIMAL_INPUT_FORMATS = new Set([0x5e, 0x63]);
const MAX_GZIP_RATIO = 1100;

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
  const hiddenIdx: number[] = [];
  r.skip(nVisible * 4);
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
  return {header, jmpVersion: m[1], buildString, charset, savetime, nrows, ncols, columns, hiddenIndex: hiddenIdx,
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

  const done = (kind: JmpColumnKind, values: (number | string | null)[]): JmpColumn =>
    ({info: col, kind, descriptor: dt, values});
  const [dt1, dt2, dt3, dt4, dt5, dt6] = dt;

  if (dt1 === 0x01 || dt1 === 0x0a) {
    const values = readNumbers(a, nrows, dt6, col.name);
    if ((dt4 === dt5 && PLAIN_NUMERIC_FORMATS.has(dt4)) || FIXED_DECIMAL_INPUT_FORMATS.has(dt5))
      return done(dt6 === 8 ? 'float' : 'int', values);
    if (dt4 === dt5 && CURRENCY_FORMATS.has(dt4))
      return done('currency', values);
    if (dt4 === dt5 && LONGITUDE_FORMATS.has(dt4))
      return done('longitude', values);
    if (dt4 === dt5 && LATITUDE_FORMATS.has(dt4))
      return done('latitude', values);
    if (DATE_FORMATS.has(dt4))
      return done('date', values);
    if (DATETIME_FORMATS.has(dt4))
      return done('datetime', values);
    if (DURATION_FORMATS.has(dt4))
      return done('duration', values);
  }

  if (dt1 === 0xff || dt1 === 0xfe || dt1 === 0xfc)
    return done(dt5 === 8 ? 'float' : 'int', readNumbers(a, nrows, dt5, col.name));

  if (dt1 === 0x09 && dt2 === 0x03)
    return {...done('rowstate', []), ...readRowStates(a, nrows, dt5, col.name)};

  if (dt1 === 0x03 && dt2 === 0x03) {
    const v = new DataView(a.buffer, a.byteOffset, a.byteLength);
    const start = tailStart(a, nrows * dt5, col.name);
    const values: number[] = [];
    for (let k = 0; k < nrows; k++)
      values.push(Number(v.getBigInt64(start + k * dt5, true)));
    return done('bigint', values);
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
    info: col, kind: 'unsupported', descriptor: dt, values: new Array(nrows).fill(null),
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
    try {
      columns.push(await readColumn(r, info, i));
    } catch (e) {
      columns.push({
        info: info.columns[i], kind: 'unsupported', descriptor: [], values: new Array(info.nrows).fill(null),
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

/** Converts seconds since the JMP epoch (1904-01-01) to a Date. */
export function jmpSecondsToDate(s: number): Date {
  return new Date((s - JMP_EPOCH_OFFSET_S) * 1000);
}

function tailStart(a: Uint8Array, n: number, name: string): number {
  if (n > a.length)
    throw new JmpFormatError(`Column "${name}": expected ${n} bytes of data, found ${a.length}`);
  return a.length - n;
}

/** Reads the trailing `nrows` numbers of `width` bytes; JMP's missing-value sentinels become null. */
function readNumbers(a: Uint8Array, nrows: number, width: number, name: string): (number | null)[] {
  const start = tailStart(a, nrows * width, name);
  const v = new DataView(a.buffer, a.byteOffset, a.byteLength);
  const out: (number | null)[] = new Array(nrows);
  for (let k = 0; k < nrows; k++) {
    const p = start + k * width;
    let x: number;
    switch (width) {
    case 1: x = v.getInt8(p); out[k] = x === -127 ? null : x; break;
    case 2: x = v.getInt16(p, true); out[k] = x === -32767 ? null : x; break;
    case 4: x = v.getInt32(p, true); out[k] = x === -2147483647 ? null : x; break;
    case 8: x = v.getFloat64(p, true); out[k] = Number.isNaN(x) ? null : x; break;
    default: throw new JmpFormatError(`Column "${name}": unknown numeric storage width ${width}`);
    }
  }
  return out;
}

/** Row states, as `"<marker> <#RRGGBB>"` values (the form JMPReader.jl reports) and decoded. Only marker
 * and colour are known; the other bytes are kept as flags. */
function readRowStates(a: Uint8Array, nrows: number, width: number, name: string):
  Pick<JmpColumn, 'values' | 'rowStates' | 'rowStateFlags'> {
  if (width < 7 || width * nrows > a.length)
    throw new JmpFormatError(`Column "${name}": row states don't fit the column`);
  const values: string[] = [];
  const rowStates: {marker: number, color: string}[] = [];
  const rowStateFlags: number[] = [];
  for (let k = 0; k < nrows; k++) {
    const o = width * k;
    const markerIdx = (a[o + 6] << 8) | a[o + 5];
    const marker = markerIdx <= 0x1f ? ROWSTATE_MARKERS[markerIdx] : String.fromCodePoint(markerIdx);
    const color = a[o + 3] === 0xff ?
      '#' + [a[o + 2], a[o + 1], a[o]].map((b) => b.toString(16).padStart(2, '0')).join('').toUpperCase() :
      ROWSTATE_COLORS[a[o]] ?? '#000000';
    values.push(`${marker} ${color}`);
    rowStates.push({marker: markerIdx, color});
    rowStateFlags.push(a[o + 4] | (width > 7 ? a[o + 7] << 8 : 0));
  }
  return {values, rowStates, rowStateFlags};
}

function readFixedWidthStrings(a: Uint8Array, nrows: number, width: number, name: string): string[] {
  const start = tailStart(a, nrows * width, name);
  const out: string[] = new Array(nrows);
  for (let k = 0; k < nrows; k++) {
    const s = start + k * width;
    let e = s;
    while (e < s + width && a[e] !== 0)
      e++;
    out[k] = utf8.decode(a.subarray(s, e));
  }
  return out;
}

function sliceStrings(data: Uint8Array, widths: number[], name: string): string[] {
  const out: string[] = new Array(widths.length);
  let p = 0;
  for (let k = 0; k < widths.length; k++) {
    if (widths[k] < 0 || p + widths[k] > data.length)
      throw new JmpFormatError(`Column "${name}": string data runs past the column`);
    out[k] = utf8.decode(data.subarray(p, p + widths[k]));
    p += widths[k];
  }
  return out;
}

function readCompressedStrings(a: Uint8Array, nrows: number, name: string): string[] {
  const r = new ByteReader(a);
  if (a.length >= 8 && r.i64() === a.length) {
    // pooled ("compact") storage: per-row indices into a string pool
    r.skip(9);
    r.i64();
    const idxWidth = r.i8();
    const idx: number[] = [];
    for (let k = 0; k < nrows; k++) {
      idx.push(idxWidth === 1 ? r.u8() : idxWidth === 2 ? r.u16() : idxWidth === 4 ? r.u32() :
        badWidth(idxWidth, name));
    }
    const npools = idx.reduce((m, x) => Math.max(m, x), 0);
    const lenWidth = r.i8();
    const pool: string[] = [];
    for (let k = 0; k < npools; k++) {
      const n = lenWidth === 1 ? r.i8() : lenWidth === 2 ? r.i16() : badWidth(lenWidth, name);
      pool.push(utf8.decode(r.bytesAt(n)));
    }
    return idx.map((j) => j === 0 ? '' : pool[j - 1]);
  }
  const widthBytes = a[8];
  const r2 = new ByteReader(a);
  r2.seek(13);
  const widths: number[] = [];
  for (let k = 0; k < nrows; k++)
    widths.push(widthBytes === 1 ? r2.i8() : widthBytes === 2 ? r2.i16() : badWidth(widthBytes, name));
  return sliceStrings(a.subarray(r2.pos), widths, name);
}

function readUncompressedStrings(r: ByteReader, a: Uint8Array, nrows: number, name: string): string[] {
  // r is positioned just after dt1..dt6
  r.skip(6);
  r.skip(r.i64());
  r.skip(2);
  r.skip(r.u32() + 8);
  const widthBytes = r.u8();
  r.u32(); // max width
  const widths: number[] = [];
  for (let k = 0; k < nrows; k++) {
    widths.push(widthBytes === 1 ? r.i8() : widthBytes === 2 ? r.i16() : widthBytes === 4 ? r.i32() :
      badWidth(widthBytes, name));
  }
  const total = widths.reduce((s, w) => s + w, 0);
  // `a` spans the column block, which ends where the string data ends
  return sliceStrings(a.subarray(a.length - total), widths, name);
}

function badWidth(w: number, name: string): never {
  throw new JmpFormatError(`Column "${name}": unknown width byte ${w}`);
}

function sameVersion(a: number[], b: number[]): boolean {
  return a.length === b.length && a.every((v, i) => v === b[i]);
}
