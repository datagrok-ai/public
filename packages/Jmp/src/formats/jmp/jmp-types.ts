// Ported from JMPReader.jl (MIT, (c) 2024 Jaakko Ruohio), src/types.jl and src/constants.jl,
// upstream commit 9fd68aa2651a8dde002ecdbd06af4d8d832b18f2

import type {JmpScript} from './jmp-scripts';
import type {JmpColumnProperties} from './jmp-properties';

/** First 4 bytes of a .jmp file. */
export const MAGIC_JMP = [0xff, 0xff, 0x00, 0x00];

/** First 4 bytes of a big-endian .jmp file, as classic Mac JMP wrote it (the old format, byte-swapped). */
export const MAGIC_JMP_BIG_ENDIAN = [0x00, 0x00, 0xff, 0xff];

/** Version triple of the format JMPReader.jl supports. */
export const SUPPORTED_VERSION: [number, number, number] = [0x07, 0x04, 0x00010001];

/** Version triple of an older format that JMPReader.jl rejects. */
export const OLD_VERSION: [number, number, number] = [0x03, 0x04, 0x03];

/** A gzip-compressed column section starts with this byte sequence (0xCAFEBEEF, little-endian). */
export const GZIP_SECTION_START = [0xef, 0xbe, 0xfe, 0xca];

/** Seconds between the JMP epoch (1904-01-01) and the Unix epoch (1970-01-01). */
export const JMP_EPOCH_OFFSET_S = 2082844800;

// Column display formats (dt4). 0x6c=":day:hr:m" and 0x6d=":day:hr:m:s" could be time or
// duration; like upstream, we treat them as duration.
export const DATE_FORMATS = new Set([
  0x65, 0x66, 0x67, 0x6e, 0x6f, 0x70, 0x71, 0x72, 0x75, 0x76, 0x7a, 0x7f, 0x88, 0x8b]);
export const DATETIME_FORMATS = new Set([
  0x69, 0x6a, 0x73, 0x74, 0x77, 0x78, 0x79, 0x7b, 0x7c, 0x7d, 0x7e, 0x80, 0x81, 0x82, 0x86, 0x87, 0x89, 0x8a]);
export const DURATION_FORMATS = new Set([0x6c, 0x6d, 0x83, 0x84, 0x85]);
export const CURRENCY_FORMATS = new Set([0x5f]);
export const LONGITUDE_FORMATS = new Set([0x54, 0x55, 0x56]);
export const LATITUDE_FORMATS = new Set([0x51, 0x52, 0x53]);

/** JMP row-state marker glyphs, by index. */
export const ROWSTATE_MARKERS = [
  '•', '+', 'X', '□', '◊', '△', 'Y', 'Z', '◯', '▭', '▯', '*', '⚫', '▬', '▮', '■',
  '⧫', '▽', '◁', '▷', '▲', '▼', '◀', '▶', '∧', '∨', '<', '>', '∣', '─', '/', '\\',
];

/** JMP's row-state colour palette, by index. */
export const ROWSTATE_COLORS = [
  '#000000', '#555555', '#787878', '#C0C0C0', '#FFFFFF',
  '#A00922', '#C91629', '#F03246', '#FF5C76', '#FF98A6',
  '#904700', '#BC5B03', '#E57406', '#FF9138', '#FFB17D',
  '#706F00', '#AFA502', '#DAD109', '#F0E521', '#FFF977',
  '#516A00', '#729400', '#90BF04', '#A2DC06', '#C1FF3D',
  '#00670C', '#11981B', '#21BC2D', '#23E72E', '#6AFF6B',
  '#007254', '#019970', '#04C791', '#06E3AA', '#0FFFBC',
  '#006D71', '#01989C', '#0CBCBC', '#06E2E3', '#67FFF7',
  '#00638F', '#0380B4', '#05A1D2', '#08C5F7', '#75E4FF',
  '#034AB0', '#0557D6', '#2867FD', '#4E9CFF', '#7E89FF',
  '#6A02A7', '#8C05CF', '#AB08FC', '#C170FF', '#D29AFF',
  '#930195', '#B803B9', '#DC06E1', '#FA50FF', '#FD96FF',
  '#9D0170', '#C5048D', '#E906A4', '#FF49BC', '#FF8CCD',
];

/** Thrown for any structural problem in a .jmp file. */
export class JmpFormatError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'JmpFormatError';
  }
}

export interface JmpHeader {
  version: [number, number, number];
}

export interface JmpInfo {
  header: JmpHeader;
  /** JMP version parsed from the build string, e.g. "18.0.1". */
  jmpVersion: string;
  buildString: string;
  charset: string;
  savetime: Date;
  nrows: number;
  ncols: number;
  columns: JmpColumnInfo[];
  /** Raw entries of the visible- and hidden-column indexes; their meaning isn't confirmed yet (JMP_PLUGIN.md §4.4). */
  visibleIndex: number[];
  hiddenIndex: number[];
  /** Byte range between the build string and the first column, where table properties such as saved
   * scripts are kept (JMP_PLUGIN.md §4.5). */
  propertiesRange: [number, number];
}

export interface JmpColumnInfo {
  name: string;
  /** Absolute file offset of the column block. */
  offset: number;
  /** Display width in JMP. */
  width: number;
}

export type JmpColumnKind =
  'float' | 'int' | 'bigint' | 'string' | 'date' | 'datetime' | 'duration' |
  'currency' | 'latitude' | 'longitude' | 'rowstate' | 'unsupported';

/** Missing value of an `int` column: -2^31, which JMP doesn't use (its sentinel is -2^31 + 1) and which
 * is also Datagrok's INT_NULL, so the array can become a column's storage as it is. */
export const JMP_INT_MISSING = -2147483648;

/**
 * A column's values, in typed arrays a Datagrok column can take as its storage (JMP_PLUGIN.md §8).
 * - `float`: every numeric kind but `int` and `bigint`; NaN is missing. Dates and date-times are seconds
 *   since the JMP epoch (1904), durations are seconds.
 * - `int`: 1-, 2- and 4-byte integers; {@link JMP_INT_MISSING} is missing.
 * - `bigint`: 8-byte integers.
 * - `string`: text, as its distinct values and a per-row index into them (Datagrok's own string storage,
 *   and much faster to hand over than one string per row); the empty string is missing.
 * - `rowstate`: per row, the marker index, the colour as 0xRRGGBB, and the bytes not decoded yet.
 * - `none`: a column that couldn't be decoded.
 */
export type JmpData =
  | {type: 'float', values: Float64Array}
  | {type: 'int', values: Int32Array}
  | {type: 'bigint', values: BigInt64Array}
  | {type: 'string', categories: string[], indexes: Int32Array}
  | {type: 'rowstate', markers: Uint32Array, colors: Uint32Array, flags: Uint16Array}
  | {type: 'none'};

export interface JmpColumn {
  info: JmpColumnInfo;
  kind: JmpColumnKind;
  /** Raw descriptor bytes dt1..dt6. */
  descriptor: number[];
  data: JmpData;
  /** Set when the column could not be decoded; `data` is then `none`. */
  decodeError?: string;
  /** Fixed-decimal numeric columns only: the number of decimals JMP shows. */
  fixedDecimals?: number;
  /** Currency columns only: the ISO 4217 code, e.g. "EUR". */
  currency?: string;
  /** Column properties (notes, units, formula, value order and the rest) and value labels (jmp-properties.ts). */
  properties?: JmpColumnProperties;
}

/** The glyph JMP draws for a row-state marker index. */
export function rowStateGlyph(marker: number): string {
  return marker <= 0x1f ? ROWSTATE_MARKERS[marker] : String.fromCodePoint(marker);
}

/** 0xRRGGBB as "#RRGGBB". */
export function hexColor(rgb: number): string {
  return '#' + rgb.toString(16).padStart(6, '0').toUpperCase();
}

/**
 * The value in row `i`, one at a time: null for missing, a number for numeric kinds (in the units of
 * {@link JmpData}), a string for text, and "<glyph> #RRGGBB" for a row state (the form JMPReader.jl
 * reports). For tests and tools; the DataFrame conversion reads the arrays directly.
 */
export function jmpValue(c: JmpColumn, i: number): number | string | null {
  const d = c.data;
  switch (d.type) {
  case 'float': return Number.isNaN(d.values[i]) ? null : d.values[i];
  case 'int': return d.values[i] === JMP_INT_MISSING ? null : d.values[i];
  case 'bigint': return Number(d.values[i]);
  case 'string': return d.categories[d.indexes[i]];
  case 'rowstate': return `${rowStateGlyph(d.markers[i])} ${hexColor(d.colors[i])}`;
  case 'none': return null;
  }
}

/** All of a column's values, as {@link jmpValue} gives them. */
export function jmpValues(c: JmpColumn, nrows: number): (number | string | null)[] {
  return Array.from({length: nrows}, (_, i) => jmpValue(c, i));
}

export interface JmpTable {
  info: JmpInfo;
  columns: JmpColumn[];
  /** Saved table scripts (jmp-scripts.ts); empty when there are none or they couldn't be found. */
  scripts: JmpScript[];
}
