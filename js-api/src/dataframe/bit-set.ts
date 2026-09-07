/**
 * Efficient bit storage and manipulation.
 * @module dataframe/bit-set
 */

import {IndexPredicate, SimilarityMetric} from "../const";
import {SIMILARITY_METRIC} from "../const";
import {observeStream} from "../events";
import {Observable} from "rxjs";
import {IDartApi} from "../api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/**
 * Efficient bit storage and manipulation.
 * See samples: {@link https://public.datagrok.ai/js/samples/data-frame/aggregation/aggregate}
 */
export class BitSet {
  public dart: any;

  /** Creates a {@link BitSet} from the specified Dart object. */
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Creates a {@link BitSet} from the string representing the bitset.
   * @param zerosOnes - A string containing '1' and '0'. */
  static fromString(zerosOnes: string): BitSet {
    return new BitSet(api.grok_BitSet_FromString(zerosOnes));
  }

  /** Creates a {@link BitSet} from the ArrayBuffer representing the bitset.
   * @param buffer - Packed bits, least significant bit first.
   * @param bitLength - count of bits. */
  static fromBytes(buffer: ArrayBuffer, bitLength: number): BitSet {
    if (bitLength == null || !Number.isInteger(bitLength) || bitLength < 0)
      bitLength = buffer.byteLength * 8;
    return new BitSet(api.grok_BitSet_FromBytes(buffer, bitLength));
  }

  /** Creates a {@link BitSet} of the specified length with all bits set to false.
   * @param length - Number of bits.
   * @param f - when specified, Sets all bits by setting i-th bit to the results of f(i) */
  static create(length: number, f?: IndexPredicate | null): BitSet {
    let bitset = new BitSet(api.grok_BitSet(length));
    if (f != null)
      bitset.init(f);
    return bitset;
  }

  /** Returns the underlying storage. Be careful with the
   * direct manipulations, as some statistics (set count, etc) are cached. */
  getBuffer(): Int32Array {
    return api.grok_BitSet_Get_Buffer(this.dart);
  }

  /** A string of `0` and `1` characters, one per bit. */
  toBinaryString(): string {
    return api.grok_BitSet_ToBinaryString(this.dart);
  }

  /** Number of bits in a bitset */
  get length(): number {
    return api.grok_BitSet_Get_Length(this.dart);
  }

  /** Number of set bits */
  get trueCount(): number {
    return api.grok_BitSet_Get_TrueCount(this.dart);
  }

  /** Number of unset bits */
  get falseCount(): number {
    return api.grok_BitSet_Get_FalseCount(this.dart);
  }

  /** Whether any bits are set to true. */
  get anyTrue(): boolean { return this.trueCount > 0; }

  /** Whether any bits are set to false. */
  get anyFalse(): boolean { return this.falseCount > 0; }

  /** Version of the bitset */
  get version(): number { return api.grok_BitSet_Get_Version(this.dart); }

  /** Clones a bitset */
  clone(): BitSet {
    return new BitSet(api.grok_BitSet_Clone(this.dart));
  }

  /** Inverts a bitset. */
  invert(notify: boolean = true): BitSet {
    api.grok_BitSet_Invert(this.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise AND operation against the
   *  specified bitset. Returns this. */
  and(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_And(this.dart, other.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise OR operation against the
   *  specified bitset. Returns this. */
  or(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_Or(this.dart, other.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise XOR operation against the
   *  specified bitset. Returns this. */
  xor(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_Xor(this.dart, other.dart, notify);
    return this;
  }

  /** Modifies this bitset by performing the bitwise AND_NOT operation against the
   *  specified bitset. Returns this. */
  andNot(other: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_AndNot(this.dart, other.dart, notify);
    return this;
  }

  /** Sets all bits to x
   * @param notify - whether BitSet's `changed` event should be fired */
  setAll(x: boolean, notify: boolean = true): BitSet {
    api.grok_BitSet_SetAll(this.dart, x, notify);
    return this;
  }

  /** Finds the first index of value x, going forward from i-th position.
   * @param i - index */
  findNext(i: number, x: boolean): number {
    return api.grok_BitSet_FindNext(this.dart, i, x);
  }

  /** Finds the first index of value x, going backward from i-th position, or -1 if not found.
   * @param i - Index to start searching from.
   * @param x - Value to search for.
   * @returns index of the first bit set to x, or -1 if not found */
  findPrev(i: number, x: boolean): number {
    return api.grok_BitSet_FindPrev(this.dart, i, x);
  }

  /** Gets i-th bit */
  get(i: number): boolean {
    return api.grok_BitSet_GetBit(this.dart, i);
  }

  /** Sets i-th bit to x
   * @param notify - whether BitSet's `changed` event should be fired */
  set(i: number, x: boolean, notify: boolean = true): void {
    api.grok_BitSet_SetBit(this.dart, i, x, notify);
  }

  /** Sets all bits by setting i-th bit to the results of f(i)
   * @param f - function that accepts bit index and returns bit value
   * @param notify - whether BitSet's `changed` event should be fired
   * @returns this */
  init(f: IndexPredicate, notify: boolean = true): BitSet {
    let buf = api.grok_BitSet_Get_Buffer(this.dart);
    let length = this.length;

    for (let i = 0; i < length; i++)
      buf[i] = 0;

    for (let i = 0; i < length; i++) {
      let idx = (i / 0x20) | 0;
      if (f(i))
        buf[idx] |= 1 << (i & 0x1f);
    }

    api.grok_BitSet_SetBuffer(this.dart, buf, notify);
    return this;
  }

  /** Indexes of all set bits. The result is cached.  */
  getSelectedIndexes(): Int32Array {
    return api.grok_BitSet_GetSelectedIndexes(this.dart);
  }

  /** Copies the content from the other {@link BitSet}. */
  copyFrom(b: BitSet, notify: boolean = true): BitSet {
    api.grok_BitSet_CopyFrom(this.dart, b.dart, notify);
    return this;
  }

  /** Fires {@link onChanged} after bits were set without notification. */
  fireChanged(): void {
    api.grok_BitSet_FireChanged(this.dart);
  }

  /** @returns fires when the bitset gets changed. */
  get onChanged(): Observable<any> {
    return observeStream(api.grok_BitSet_Changed(this.dart));
  }

  /** Finds the value of similarity between two BitSets.
   * @param b - second BitSet.
   * @param metric - similarity metric to use.
   * @returns similarity value */
  similarityTo(b: BitSet, metric: SimilarityMetric = SIMILARITY_METRIC.TANIMOTO): number {
    return api.grok_BitSet_SimilarityTo(this.dart, b.dart, metric);
  }

  /** @returns string representation of this bitset, like '0110'. */
  toString(): string {
    return api.grok_Object_ToString(this.dart);
  }

  /** Applies the standard click semantics (Ctrl toggles, Shift extends) to the rows matching [rowIndexPredicate]. */
  handleClick(rowIndexPredicate: IndexPredicate, mouseEvent: MouseEvent, modifiedSelectOnly: boolean = false) {
    api.grok_Utils_SelectRowsWhere(this.dart, rowIndexPredicate, mouseEvent.ctrlKey, mouseEvent.shiftKey, mouseEvent.metaKey, modifiedSelectOnly);
  }
}
