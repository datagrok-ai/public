export type BitPredicate = (i: number) => boolean;

/** A pure-JS bit array with the {@link BitSet} vocabulary: LSB-first `Uint32Array` words, bit `i` at
 * `words[i >>> 5] & (1 << (i & 31))`, and the bits of the last word beyond `length` always zero — so
 * `equals` is a word compare and `DG.BitSet.fromBitArray` reads the words as they are. */
export class BitArray {
  private _data: Uint32Array;
  private _length: number;

  constructor(length: number, fill?: boolean);
  /** Adopts `words` (no copy); they must hold at least `lengthInInts` words. */
  constructor(words: Uint32Array, length: number);
  constructor(arg: number | Uint32Array, second: boolean | number = false) {
    if (typeof arg === 'number') {
      this._length = arg;
      this._data = new Uint32Array(BitArray._wordCount(arg));
      if (second)
        this._data.fill(0xffffffff);
    }
    else {
      this._length = second as number;
      if (arg.length < BitArray._wordCount(this._length))
        throw new RangeError(`${arg.length} words cannot hold ${this._length} bits`);
      this._data = arg;
    }
    this._trim();
  }

  private static _wordCount(length: number): number {
    if (!(Number.isInteger(length) && length >= 0))
      throw new RangeError(`Invalid length ${length}`);
    return (length + 31) >>> 5;
  }

  private _trim(): this {
    const r = this._length & 31;
    if (r !== 0)
      this._data[this._length >>> 5] &= (1 << r) - 1;
    return this;
  }

  private _checkRange(from: number, to: number): void {
    if (from < 0 || to > this._length || from > to)
      throw new RangeError(`Range [${from}, ${to}) is outside [0, ${this._length})`);
  }

  private _checkLength(other: BitArray): void {
    if (other._length !== this._length)
      throw new RangeError(`Lengths differ (${this._length} != ${other._length})`);
  }

  private static _popcount(words: Uint32Array, n: number): number {
    let c = 0;
    for (let i = 0; i < n; i++) {
      let v = words[i];
      v -= (v >>> 1) & 0x55555555;
      v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
      c += Math.imul((v + (v >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24;
    }
    return c;
  }

  /** A bit array of `length` bits, the i-th one `f(i)` when `f` is given. */
  static create(length: number, f?: BitPredicate | null): BitArray {
    const result = new BitArray(length);
    return f == null ? result : result.init(f);
  }

  /** From a string of '0' and '1' characters, like '0110'. */
  static fromString(zerosOnes: string): BitArray {
    return BitArray.create(zerosOnes.length, (i) => zerosOnes.charCodeAt(i) === 49);
  }

  /** A copy of the little-endian bytes; `bitLength` defaults to all of them. */
  static fromBytes(buffer: ArrayBuffer | ArrayBufferView, bitLength?: number): BitArray {
    const bytes = ArrayBuffer.isView(buffer) ?
      new Uint8Array(buffer.buffer, buffer.byteOffset, buffer.byteLength) : new Uint8Array(buffer);
    const length = bitLength ?? bytes.length * 8;
    if (length > bytes.length * 8)
      throw new RangeError(`${bytes.length} bytes cannot hold ${length} bits`);
    const result = new BitArray(length);
    const words = result._data;
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    const full = Math.min(words.length, bytes.length >>> 2);
    for (let i = 0; i < full; i++)
      words[i] = view.getUint32(i * 4, true);
    for (let i = full * 4, end = Math.min(bytes.length, words.length * 4); i < end; i++)
      words[i >>> 2] |= bytes[i] << ((i & 3) * 8);
    return result._trim();
  }

  /** Adopts `words` (no copy), like the constructor. */
  static fromUint32Array(length: number, words: Uint32Array): BitArray {
    return new BitArray(words, length);
  }

  get length(): number {
    return this._length;
  }

  get lengthInInts(): number {
    return BitArray._wordCount(this._length);
  }

  /** The live words; a direct write must keep the bits beyond `length` zero. */
  getBuffer(): Uint32Array {
    return this._data;
  }

  /** Gets the i-th bit; unchecked. */
  get(i: number): boolean {
    return (this._data[i >>> 5] & (1 << (i & 31))) !== 0;
  }

  /** Sets the i-th bit; unchecked. */
  set(i: number, x: boolean): void {
    if (x)
      this._data[i >>> 5] |= 1 << (i & 31);
    else
      this._data[i >>> 5] &= ~(1 << (i & 31));
  }

  get trueCount(): number {
    return BitArray._popcount(this._data, this.lengthInInts);
  }

  get falseCount(): number {
    return this._length - this.trueCount;
  }

  get anyTrue(): boolean { return this.findNext(-1, true) !== -1; }
  get anyFalse(): boolean { return this.findNext(-1, false) !== -1; }
  get allTrue(): boolean { return !this.anyFalse; }
  get allFalse(): boolean { return !this.anyTrue; }

  clone(): BitArray {
    return new BitArray(this._data.slice(0, this.lengthInInts), this._length);
  }

  /** Copies the bits of `other`, which must have the same length. */
  copyFrom(other: BitArray): this {
    this._checkLength(other);
    this._data.set(other._data.subarray(0, this.lengthInInts));
    return this;
  }

  /** Same length and same bits. */
  equals(other: BitArray): boolean {
    if (other === this)
      return true;
    if (other == null || other._length !== this._length)
      return false;
    for (let i = 0, n = this.lengthInInts; i < n; i++) {
      if (this._data[i] !== other._data[i])
        return false;
    }
    return true;
  }

  /** Sets the i-th bit to `f(i)` for every i. */
  init(f: BitPredicate): this {
    this._data.fill(0, 0, this.lengthInInts);
    for (let i = 0; i < this._length; i++) {
      if (f(i))
        this._data[i >>> 5] |= 1 << (i & 31);
    }
    return this;
  }

  setAll(x: boolean): this {
    this._data.fill(x ? 0xffffffff : 0, 0, this.lengthInInts);
    return this._trim();
  }

  invert(): this {
    for (let i = 0, n = this.lengthInInts; i < n; i++)
      this._data[i] = ~this._data[i];
    return this._trim();
  }

  /** Bitwise AND with `other` (same length) in place. */
  and(other: BitArray): this {
    this._checkLength(other);
    for (let i = 0, n = this.lengthInInts; i < n; i++)
      this._data[i] &= other._data[i];
    return this;
  }

  /** Bitwise OR with `other` (same length) in place. */
  or(other: BitArray): this {
    this._checkLength(other);
    for (let i = 0, n = this.lengthInInts; i < n; i++)
      this._data[i] |= other._data[i];
    return this;
  }

  /** Bitwise XOR with `other` (same length) in place. */
  xor(other: BitArray): this {
    this._checkLength(other);
    for (let i = 0, n = this.lengthInInts; i < n; i++)
      this._data[i] ^= other._data[i];
    return this;
  }

  /** Clears the bits set in `other` (same length) in place. */
  andNot(other: BitArray): this {
    this._checkLength(other);
    for (let i = 0, n = this.lengthInInts; i < n; i++)
      this._data[i] &= ~other._data[i];
    return this;
  }

  /** The number of bits equal to `value` in `this AND other`, without modifying either. */
  andWithCountBits(other: BitArray, value: boolean = true): number {
    this._checkLength(other);
    let c = 0;
    for (let i = 0, n = this.lengthInInts; i < n; i++) {
      let v = this._data[i] & other._data[i];
      v -= (v >>> 1) & 0x55555555;
      v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
      c += Math.imul((v + (v >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24;
    }
    return value ? c : this._length - c;
  }

  /** The first index after `i` (pass -1 to start from the beginning) whose bit is `x`, or -1. */
  findNext(i: number, x: boolean = true): number {
    const from = Math.max(i + 1, 0);
    if (from >= this._length)
      return -1;
    const first = from >>> 5;
    for (let w = first, n = this.lengthInInts; w < n; w++) {
      let v = x ? this._data[w] : ~this._data[w];
      if (w === first)
        v &= -1 << (from & 31);
      if (v !== 0) {
        const index = (w << 5) + 31 - Math.clz32(v & -v);
        return index < this._length ? index : -1;
      }
    }
    return -1;
  }

  /** The last index before `i` (pass -1 to start from the end) whose bit is `x`, or -1. */
  findPrev(i: number, x: boolean = true): number {
    const from = i < 0 ? this._length - 1 : Math.min(i, this._length) - 1;
    if (from < 0)
      return -1;
    const last = from >>> 5;
    for (let w = last; w >= 0; w--) {
      let v = x ? this._data[w] : ~this._data[w];
      if (w === last && (from & 31) !== 31)
        v &= (1 << ((from & 31) + 1)) - 1;
      if (v !== 0)
        return (w << 5) + 31 - Math.clz32(v);
    }
    return -1;
  }

  /** Indexes of all set bits. */
  getSelectedIndexes(): Int32Array {
    const result = new Int32Array(this.trueCount);
    let k = 0;
    for (let w = 0, n = this.lengthInInts; w < n; w++) {
      for (let v = this._data[w]; v !== 0; v &= v - 1)
        result[k++] = (w << 5) + 31 - Math.clz32(v & -v);
    }
    return result;
  }

  /** Iterates the indexes of the set bits in ascending order. */
  *trueIndexes(): IterableIterator<number> {
    for (let w = 0, n = this.lengthInInts; w < n; w++) {
      for (let v = this._data[w]; v !== 0; v &= v - 1)
        yield (w << 5) + 31 - Math.clz32(v & -v);
    }
  }

  /** Resizes to `n` bits, the added ones false; reallocates to exactly the words needed. */
  setLength(n: number): void {
    const words = BitArray._wordCount(n);
    if (words !== this._data.length) {
      const data = new Uint32Array(words);
      data.set(this._data.subarray(0, Math.min(words, this.lengthInInts)));
      this._data = data;
    }
    else
      this._data.fill(0, this.lengthInInts, words);
    this._length = n;
    this._trim();
  }

  /** Removes `n` bits starting at `pos`, shifting the following bits down. */
  removeAt(pos: number, n: number = 1): void {
    if (n < 0)
      throw new RangeError(`Invalid count ${n}`);
    this._checkRange(pos, pos + n);
    for (let i = pos, end = this._length - n; i < end; i++)
      this.set(i, this.get(i + n));
    this.setLength(this._length - n);
  }

  /** Sets the bits in `[from, to)` to `x`. */
  setRange(from: number, to: number, x: boolean): this {
    this._checkRange(from, to);
    if (from === to)
      return this;
    const first = from >>> 5;
    const last = (to - 1) >>> 5;
    const head = -1 << (from & 31);
    const tail = (to & 31) === 0 ? -1 : (1 << (to & 31)) - 1;
    const d = this._data;
    if (first === last) {
      if (x)
        d[first] |= head & tail;
      else
        d[first] &= ~(head & tail);
    }
    else if (x) {
      d[first] |= head;
      d.fill(0xffffffff, first + 1, last);
      d[last] |= tail;
    }
    else {
      d[first] &= ~head;
      d.fill(0, first + 1, last);
      d[last] &= ~tail;
    }
    return this;
  }

  /** A copy of the bits in `[from, to)`. */
  getRange(from: number, to: number): BitArray {
    this._checkRange(from, to);
    const result = new BitArray(to - from);
    const shift = from & 31;
    const src = this._data;
    const dst = result._data;
    for (let i = 0, w = from >>> 5, n = this.lengthInInts; i < dst.length; i++, w++) {
      let v = src[w] >>> shift;
      if (shift !== 0 && w + 1 < n)
        v |= src[w + 1] << (32 - shift);
      dst[i] = v;
    }
    return result._trim();
  }

  /** The bits as a string of '0' and '1' characters, like '0110'. */
  toString(): string {
    let s = '';
    for (let i = 0; i < this._length; i++)
      s += this.get(i) ? '1' : '0';
    return s;
  }

  toBinaryString(): string {
    return this.toString();
  }
}
