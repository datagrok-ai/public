/**
 * Qualified numbers for imprecise measurements.
 * @module dataframe/qnum
 */

import {IDartApi} from "../api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export const QNUM_LESS = 1;
export const QNUM_EXACT = 2;
export const QNUM_GREATER = 3;

let _qnumBuf = new DataView(new ArrayBuffer(8));

/**
 *  A set of static methods for working with qualified numbers.
 *  The internal representation of a qualified number is a regular double precision floating point
 *  number (IEEE 754), except the two least significant bits in mantissa are reserved
 *  for holding the qualifier ([LESS], [EXACT], [GREATER]).
 *
 *  The advantage of that representation is that the standard arithmetic operations could be
 *  performed directly on the number, without unpacking it. This is especially important for batch
 *  operations such as aggregation or sorting. While there is a loss of precision, it is rather
 *  insignificant (50 bits for storing mantissa instead of 52), which makes perfect sense
 *  considering that qualified numbers represent imprecise measurements.
 *
 *  Use [create], [getValue], and [getQ] methods for packing/unpacking.
 * */
export class Qnum {
  /**
   * Extracts the qualifier ({@link QNUM_LESS}, {@link QNUM_EXACT}, {@link QNUM_GREATER}).
   * See also {@link getValue}
   * */
  static getQ(x: number): number {
    _qnumBuf.setFloat64(0, x);
    return _qnumBuf.getInt8(7) & 0x03;
  }

  /**
   * Extracts the value from x, stripping the qualifier .
   * See also {@link getQ}
   * */
  static getValue(x: number): number {
    _qnumBuf.setFloat64(0, x);
    let last = _qnumBuf.getInt8(7) & 0xFC;
    _qnumBuf.setInt8(7, last);
    return _qnumBuf.getFloat64(0);
  }

  /**
   * Creates a QNum value out of the [value] and qualifier [q].
   * */
  static create(value: number, q: number = QNUM_EXACT): number {
    _qnumBuf.setFloat64(0, value);
    let last = _qnumBuf.getInt8(7);
    _qnumBuf.setInt8(7, (last & 0xFC) | q);
    return _qnumBuf.getFloat64(0);
  }

  static exact(x: number): number {
    return Qnum.create(x, QNUM_EXACT)
  };

  static less(x: number): number {
    return Qnum.create(x, QNUM_LESS)
  };

  static greater(x: number): number {
    return Qnum.create(x, QNUM_GREATER);
  }

  /**
   * Parses a string into a qualified number.
   * */
  static parse(s: string): number {
    return api.grok_Qnum_Parse(s);
  }

  /**
   * Converts a qualified number to a string representation.
   * */
  static toString(x: number): string {
    return api.grok_Qnum_ToString(x);
  }

  /**
   * Returns the string representation of the qualifier.
   * */
  static qualifier(x: number): string {
    return api.grok_Qnum_Qualifier(x);
  }
}
