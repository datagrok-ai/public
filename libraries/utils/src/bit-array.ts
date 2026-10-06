import {BitArray as CoreBitArray} from 'datagrok-api/src/u2core/bit-array.js';

/** @deprecated Compatibility wrapper over `DG.BitArray` (`datagrok-api/src/u2core/bit-array`), which new code
 * should use directly. Keeps the old names as delegates; the counts are never cached. */
export default class BitArray {
  // ml workers rebuild structured-cloned instances as `new BitArray(x._data, x._length)`: both fields mirror
  // the core's storage and are re-synced after every operation that may reallocate it.
  private _data: Uint32Array;
  private _length: number;
  private _core: CoreBitArray;
  private _version = 0;

  constructor(data: Uint32Array, length: number);
  constructor(length: number, defaultValue?: boolean);
  constructor(arg: number | Uint32Array, second: boolean | number = false) {
    this._core = typeof arg === 'number' ?
      new CoreBitArray(arg, second as boolean) : new CoreBitArray(arg, second as number);
    this._data = this._core.getBuffer();
    this._length = this._core.length;
  }

  private static _wrap(core: CoreBitArray): BitArray {
    return new BitArray(core.getBuffer(), core.length);
  }

  private _sync(): void {
    this._data = this._core.getBuffer();
    this._length = this._core.length;
    this._version++;
  }

  static fromValues(values: Array<boolean>): BitArray {
    return BitArray.fromSeq(values.length, (i) => values[i]);
  }

  static fromSeq(count: number, flag: (i: number) => boolean): BitArray {
    return BitArray._wrap(CoreBitArray.create(count, flag));
  }

  static fromString(s: string): BitArray {
    return BitArray._wrap(CoreBitArray.fromString(s));
  }

  /** Adopts `data` (no copy). */
  static fromUint32Array(length: number, data: Uint32Array): BitArray {
    return new BitArray(data, length);
  }

  static fromBytes(bytes: Uint8Array): BitArray {
    return BitArray._wrap(CoreBitArray.fromBytes(bytes));
  }

  getRawData(): Uint32Array { return this._data; }

  get buffer(): Uint32Array { return this._data; }

  /** Copies `data` (a short array is zero-padded, as the old aliasing setter read it). */
  set buffer(data: Uint32Array) {
    const words = new Uint32Array(this._core.lengthInInts);
    words.set(data.subarray(0, words.length));
    this._core = new CoreBitArray(words, this._length);
    this._sync();
  }

  get version(): number { return this._version; }

  set version(value: number) { this._version = value; }

  incrementVersion(_notify: boolean = true): void { this._version++; }

  get length(): number { return this._length; }

  get lengthInInts(): number { return this._core.lengthInInts; }

  setLength(value: number): void {
    this._core.setLength(value);
    this._sync();
  }

  removeAt(pos: number, n: number = 1): void {
    this._core.removeAt(pos, n);
    this._sync();
  }

  clear(): void { this.setLength(0); }

  toString(): string {
    return `${this._length} bits, ${this.countBits(true)} set`;
  }

  equals(other: BitArray): boolean { return this._core.equals(other._core); }

  clone(): BitArray {
    const result = BitArray._wrap(this._core.clone());
    result._version = this._version;
    return result;
  }

  copyFrom(other: BitArray): void {
    this._core.copyFrom(other._core);
    this._version++;
  }

  init(flag: (i: number) => boolean, _notify: boolean = true): BitArray {
    this._core.init(flag);
    this._version++;
    return this;
  }

  invert(_notify: boolean = true): void {
    this._core.invert();
    this._version++;
  }

  setAll(value: boolean, _notify: boolean = false): void {
    this._core.setAll(value);
    this._version++;
  }

  /** A copy of the bits in `[from, to)`. */
  getRange(from: number, to: number): BitArray {
    return BitArray._wrap(this._core.getRange(from, to));
  }

  getRangeAsList(from: number, to: number): boolean[] {
    const result: boolean[] = [];
    for (let i = from; i < to; i++)
      result.push(this._core.get(i));
    return result;
  }

  /** Sets the bits in `[from, to]` (inclusive, in either order) to `value`. */
  setRange(from: number, to: number, value: boolean, _notify: boolean = true): BitArray {
    this._core.setRange(Math.min(from, to), Math.max(from, to) + 1, value);
    this._version++;
    return this;
  }

  and(value: BitArray, _notify: boolean = true): BitArray {
    this._core.and(value._core);
    this._version++;
    return this;
  }

  or(value: BitArray, _notify: boolean = true): BitArray {
    this._core.or(value._core);
    this._version++;
    return this;
  }

  xor(value: BitArray, _notify: boolean = true): BitArray {
    this._core.xor(value._core);
    this._version++;
    return this;
  }

  andNot(value: BitArray, _notify: boolean = true): BitArray {
    this._core.andNot(value._core);
    this._version++;
    return this;
  }

  not(_notify: boolean = true): BitArray {
    this._core.invert();
    this._version++;
    return this;
  }

  getBit(pos: number): boolean { return this._core.get(pos); }

  setBit(pos: number, bit: boolean, _notify: boolean = true): void {
    this._core.set(pos, bit);
    this._version++;
  }

  setFast(i: number, value: boolean): void {
    this._core.set(i, value);
    this._version++;
  }

  setTrue(pos: number): void { this.setFast(pos, true); }

  setFalse(pos: number): void { this.setFast(pos, false); }

  trueCount(): number { return this.countBits(true); }

  falseCount(): number { return this.countBits(false); }

  countBits(value: boolean): number {
    return value ? this._core.trueCount : this._core.falseCount;
  }

  andWithCountBits(second: BitArray, value: boolean): number {
    return this._core.andWithCountBits(second._core, value);
  }

  contains(value: boolean): boolean { return this._core.findNext(-1, value) >= 0; }

  get allTrue(): boolean { return this._core.allTrue; }

  get allFalse(): boolean { return this._core.allFalse; }

  get anyTrue(): boolean { return this._core.anyTrue; }

  get anyFalse(): boolean { return this._core.anyFalse; }

  findNext(index: number, value: boolean = true): number { return this._core.findNext(index, value); }

  findPrev(index: number, value: boolean = true): number { return this._core.findPrev(index, value); }
}
