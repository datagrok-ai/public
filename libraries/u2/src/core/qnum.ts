/* Qualified numbers, ported from `ddt/lib/src/utils/qnum.dart:17-93`.

   A qualified number IS a double: the qualifier (`<`, exact, `>`) lives in the two least
   significant mantissa bits, so sorting, aggregation and arithmetic run on it without unpacking,
   at the cost of two bits of precision. Byte order is big-endian on both sides — Dart's
   `ByteData` and the DOM's `DataView` both default to it, which is what puts the low mantissa
   byte at index 7. Dart reads that byte through `getInt8`/`setInt8`; `getUint8`/`setUint8` here
   produce the same byte (the masks only ever touch the low two bits).

   The packing is byte-for-byte the Dart one; the DEFAULT DISPLAY is not. Dart's `toStr` prints
   through `fmt` (`ddt/lib/src/utils/utils.dart:1467-1490`), which gives a non-integral value
   exactly two fraction digits — `<5.2` reads back as `'<5.20'`, and `0.001` as `'0.00'` — while
   {@link QNum.format} below rounds off the packing noise and prints the number as it is
   (`'<5.2'`, `'0.001'`). Anything that must match the platform's text passes its own `format`. */

const bytes = new DataView(new ArrayBuffer(8));

const NUMBER = /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/;

/** How many significant digits survive the packing: the two borrowed bits cost the last one of
 * a double's 17, and 15 is where a packed value rounds back to the number that was packed. */
const DIGITS = 15;

export class QNum {
  static readonly LESS = 1;
  static readonly EXACT = 2;
  static readonly GREATER = 3;

  /** Packs `value` and qualifier `q` into one double (`qnum.dart:49-54`). */
  static create(value: number, q: number = QNum.EXACT): number {
    bytes.setFloat64(0, value);
    bytes.setUint8(7, (bytes.getUint8(7) & 0xfc) | q);
    return bytes.getFloat64(0);
  }

  static less(value: number): number {
    return QNum.create(value, QNum.LESS);
  }

  static greater(value: number): number {
    return QNum.create(value, QNum.GREATER);
  }

  /** {@link LESS}, {@link EXACT} or {@link GREATER} (`qnum.dart:26-30`). */
  static getQ(x: number): number {
    bytes.setFloat64(0, x);
    return bytes.getUint8(7) & 0x03;
  }

  /** The number with its qualifier bits cleared (`qnum.dart:41-46`). */
  static getValue(x: number): number {
    bytes.setFloat64(0, x);
    bytes.setUint8(7, bytes.getUint8(7) & 0xfc);
    return bytes.getFloat64(0);
  }

  /** `''` for exact and for a cleared 0, `<` or `>` otherwise (`qnum.dart:85-90`). */
  static qualifier(q: number): string {
    if (q === 0 || q === QNum.EXACT)
      return '';
    if (q === QNum.LESS)
      return '<';
    if (q === QNum.GREATER)
      return '>';
    throw new Error(`Unknown qualifier type: ${q}`);
  }

  /** `'<5.2'` → the packed double; null for blank text and for text that is not a number
   * (`qnum.dart:62-83`). Dart parses the numeric part through `FloatColumn.parse` (Dart's own
   * `double.parse`, which also takes `Infinity` and `NaN`); this takes the numeric text the rest
   * of u2 takes, so those two do not parse here. */
  static parse(s: string): number | null {
    const text = s.replace(/^ +/, '');
    if (text === '')
      return null;
    const first = text[0];
    const q = first === '<' ? QNum.LESS : first === '>' ? QNum.GREATER : QNum.EXACT;
    const number = (q === QNum.EXACT ? text : text.substring(1)).trim();
    return NUMBER.test(number) ? QNum.create(Number(number), q) : null;
  }

  /** `'<5.2'` — the qualifier glyph in front of the value (`qnum.dart:92`). Without a `format`
   * the value is rounded to the {@link DIGITS} digits that survived the packing, which is what
   * turns `5.199999999999999` back into the `5.2` that was packed. */
  static format(x: number | null, format?: (value: number) => string): string {
    if (x === null)
      return '';
    const value = QNum.getValue(x);
    return QNum.qualifier(QNum.getQ(x)) +
      (format === undefined ? String(Number(value.toPrecision(DIGITS))) : format(value));
  }
}
