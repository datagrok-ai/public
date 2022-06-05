export default class BitArray {
  static _onBitCount = Int8Array.from([
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8]);

  static _firstOnBit = Int8Array.from([
    -1, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0]);

  static _lastOnBit = Int8Array.from([
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]);

  private _data: Uint32Array;
  private _length = 0;
  private _version = 0;
  private _updateLevel = 0;
  private _selectedCount = 0;
  private _selectedCountVersion = -1;
  private _selectedIndexesVersion = -1;
  private _versionedName = '';
  private _versionedNameVersion = -1;
  SHRINK_THRESHOLD = 0x100;

  constructor(data: Uint32Array, length: number)
  constructor(length: number, defaultValue?: boolean)
  constructor(arg: number | Uint32Array, defaultValue: boolean | number = false) {
    if (typeof arg === 'number') {
      const length = arg;
      const buff = BitArray._createBuffer(length);
      if (defaultValue) {
        for (let i = 0; i < buff.length; i++)
          buff[i] = -1;
      }
      this._data = buff;
      this._length = length;
    } else if (arg instanceof Uint32Array) {
      this._data = arg as Uint32Array;
      this._length = defaultValue as number;
    } else {
      throw new Error('Invalid constructor');
    }
  }


  getRawData() { return this._data; }

  assureGoez(num: number, argName: String): void {
    if (num < 0) throw new Error(`${argName} should be greater than zero`);
  }

  assureInRange(value: number, min: number, max: number, argName: String): void {
    if ((value < min) || (value > max))
      throw new Error(`Argument ${argName} (${value}) out of range (${min}, ${max})`);
  }

  copy(src: Uint32Array, dst: Uint32Array, count: number): void {
    for (let i = 0; i < count; i++)
      dst[i] = src[i];
  }

  copyFrom(other: BitArray): void {
    if (this._length != other._length)
      throw new Error(`Lengths differ (${this._length} != ${other._length})`);
    this.copy(other._data, this._data, this.lengthInInts);
    this._version++;
  }

  get length() {
    return this._length;
  }

  get buffer() {
    return this._data;
  }

  set buffer(data: Uint32Array) {
    this._data = data;
    this._version++;
  }

  get version() {
    return this._version;
  }

  set version(value: number) {
    this._version = value;
  }

  incrementVersion(notify = true): void {
    this._version++;
  }

  get lengthInInts() {
    return Math.floor((this._length + 0x1f) / 0x20);
  }

  get versionedName() {
    return this._version == this._versionedNameVersion ? this._versionedName : '';
  }

  set versionedName(name: string) {
    this._versionedName = name;
    this._versionedNameVersion = this._version;
  }

  get self() {
    return this;
  }

  setLength(value: number): void {
    if (value < 0)
      throw new Error('should be >= 0');

    if (value == this._length) return;
    const nIntsNeeded = Math.floor((value + 0x1f) / 0x20);
    if ((nIntsNeeded > this._data.length) || ((nIntsNeeded + this.SHRINK_THRESHOLD) < this._data.length)) {
      const newData = new Uint32Array(nIntsNeeded);
      this.copy(this._data, newData, (nIntsNeeded > this._data.length) ? this._data.length : nIntsNeeded);
      this._data = newData;
    }

    if (value > this._length) {
      if (this._length % 0x20 > 0)
        this._data[this.lengthInInts - 1] &= (1 << ((this._length % 0x20) & 0x1f)) - 1;

      this._data.fill(0, this.lengthInInts, nIntsNeeded);
    }
    this._length = value;
    this._version++;
  }

  static fromAnd(set1: BitArray, set2: BitArray): BitArray {
    if (set1._length != set2._length)
      throw new Error(`Lengths differ (${set1._length} != ${set2._length})`);

    const temp = new BitArray(set1._length);
    temp._length = set1._length;
    temp._data = BitArray._createBuffer(temp._length);
    temp._version = 0;

    const len = set1.lengthInInts;
    for (let i = 0; i < len; i++)
      temp._data[i] = set1._data[i] & set2._data[i];

    return temp;
  }

  private static _createBuffer(length: number): Uint32Array {
    return new Uint32Array(Math.floor((length + 0x1f) / 0x20));
  }

  static fromValues(values: Array<boolean>): BitArray {
    const temp = new BitArray(values.length);
    temp._version = 0;

    for (let i = 0; i < temp._length; i++) {
      if (values[i])
        temp._data[Math.floor(i / 0x20)] |= 1 << ((i % 0x20) & 0x1f);
    }
    return temp;
  }

  /// Constructs a [BitSet] of length [count], where idx-th bit is determined by a call to [flag] (idx).
  static fromSeq(count: number, flag: Function): BitArray {
    const temp = new BitArray(count);
    for (let i = 0; i < count; ++i)
      temp.setBit(i, flag(i));

    temp._version = 0;
    return temp;
  }

  /// Constructs a [BitSet] from a string [s] containing '0' or '1'.
  static fromString(s: string): BitArray {
    return BitArray.fromSeq(s.length, (i: number) => s.charAt(i) == '1');
  }

  /// Constructs a [BitSet], based on length [_length] and byte array [_data].
  static fromUint32Array(_length: number, _data: Uint32Array): BitArray {
    const temp = new BitArray(_length);
    temp._data = _data;
    return temp;
  }

  /// Deserializes a [BitSet] from [bytes].
  static fromBytes(bytes: Uint8Array): BitArray {
    const len = bytes.length;
    const temp = new BitArray(len * 8);
    temp._data = new Uint32Array(Math.floor((len + 3) / 4));
    temp._length = len * 8;
    let num1 = 0;
    let num2 = 0;

    while ((len - num2) >= 4) {
      temp._data[num1++] = (
        ((bytes[num2] & 0xff) | ((bytes[num2 + 1] & 0xff) << 8)) |
        ((bytes[num2 + 2] & 0xff) << 0x10)
      ) | ((bytes[num2 + 3] & 0xff) << 0x18);

      num2 += 4;
    }

    if (len - num2 == 3)
      temp._data[num1] = (bytes[num2 + 2] & 0xff) << 0x10;

    if (len - num2 == 2)
      temp._data[num1] |= (bytes[num2 + 1] & 0xff) << 8;

    if (len - num2 == 1)
      temp._data[num1] |= bytes[num2] & 0xff;

    temp._version = 0;
    return temp;
  }

  toString(): string {
    return `${this._length} bits, ${this.countBits(true)} set`;
  }

  /// Performs deep comparison of two bitsets.
  equals(other: BitArray): boolean {
    if (this == other) return true;
    if (other == null) return false;
    if (this._length != other._length) return false;
    if (this._length == 0) return true;

    for (let i = 0; i < this._data.length - 1; i++)
      if (this._data[i] != other._data[i]) return false;

    for (let i = (this._data.length - 1) * 8; i < this._length; i++) {
      if (this.getBit(i) != other.getBit(i))
        return false;
    }
    return true;
  }

  /** Clones a bitset. */
  clone(): BitArray {
    const bitArray = new BitArray(0, false);
    bitArray._data = Uint32Array.from(this._data); // effective length: (lengthInInts)
    bitArray._length = this._length;
    bitArray._version = this._version;
    return bitArray;
  }

  /** Initializes a bitset. */
  init(flag: Function, notify: boolean): BitArray {
    this.setAll(false, false);

    for (let i = 0; i < this._length; i++) {
      if (flag(i))
        this._data[Math.floor(i / 0x20)] |= 1 << ((i % 0x20) & 0x1f);
    }

    this.incrementVersion(notify);
    return this;
  }

  /// Inverts a bitset.
  invert(notify = true): void {
    for (let i = 0; i < this._data.length; i++)
      this._data[i] ^= -1;

    this.incrementVersion(notify);
  }

  /// Sets all bits to [value], optionally suppressing notifications.
  setAll(value: boolean, notify = false): void {
    const flags = value ? -1 : 0;
    const len = this.lengthInInts;

    for (let i = 0; i < len; i++) //todo: optimize
      this._data[i] = flags;

    this.incrementVersion(notify);
  }

  /// Sets bits at [indexes] position to [value].
  /// Clears the bitset if [clear] flag is true.
  /// Change notification is raised when [notify] is true.
  setIndexes(indexes: Array<number>, value = true, clear = true, notify = true): void {
    if (clear)
      this.setAll(!value, false);

    for (const i of indexes)
      this.setFast(i, value);

    this.incrementVersion(notify);
  }

  everyIndex(indexes: Array<number>, value = true): boolean {
    for (const index of indexes) {
      if (this.getBit(index) != value)
        return false;
    }
    return true;
  }

  anyIndex(indexes: Array<number>, value = true): boolean {
    for (const index of indexes) {
      if (this.getBit(index) == value)
        return true;
    }
    return false;
  }

  setWhere(check: Function, value = true, clear = true, notify = true, allowClear = true): void {
    if (clear && allowClear)
      this.setAll(!value, false);

    if (allowClear) {
      for (let i = 0; i < this._length; i++) {
        if (check(i))
          this.setFast(i, value);
      }
    } else {
      for (let i = 0; i < this._length; i++)
        this.setFast(i, check(i) ? value : !value);
    }

    this.incrementVersion(notify);
  }

  setRange(from: number, to: number, value: boolean, notify = true): BitArray {
    this.assureInRange(from, 0, this._length - 1, 'from');
    this.assureInRange(to, 0, this._length - 1, 'to');

    const start = Math.min(from, to);
    const end = Math.max(from, to);

    //todo: optimize
    if (value) {
      for (let i = start; i <= end; i++)
        this.setTrue(i);
    } else {
      for (let i = start; i <= end; i++)
        this.setFalse(i);
    }

    this.incrementVersion(notify);
    return this;
  }

  /// Sets n randomly chosen bits to value, remaining bits to !value.
  setRandom(n: number, value: boolean, notify = true): void {
    if (n < 0 || n > this._length)
      throw new Error('n must be >= 0 && <= Count');

    if (n > this._length / 2)
      this.setRandom(this._length - n, !value);

    this.setAll(!value);

    for (let k = 0; k < n;) {
      const i = Math.floor(Math.random() * this._length);
      if (this.getBit(i) == value) continue;
      this.setFast(i, value);
      k++;
    }

    this.incrementVersion(notify);
  }

  /// Modifies current bitset by performing the bitwise AND operation against the
  /// corresponding elements in the specified bitset.
  and(value: BitArray, notify = true): BitArray {
    if (this._length != value._length)
      throw new Error('Array lengths differ.');

    for (let i = 0, len = this.lengthInInts; i < len; i++)
      this._data[i] &= value._data[i];

    this.incrementVersion(notify);
    return this;
  }

  /// Performs the bitwise AND NOT operation on the elements in the current bitset
  /// against the corresponding elements in the specified bitset.
  andNot(value: BitArray, notify = true): BitArray {
    if (this._length != value._length)
      throw new Error('Array lengths differ.');

    const len = this.lengthInInts;
    for (let num2 = 0; num2 < len; num2++)
      this._data[num2] &= ~value._data[num2];

    this.incrementVersion(notify);
    return this;
  }

  /// Performs the bitwise NOT AND operation on the elements in the current bitset
  /// against the corresponding elements in the specified bitset.
  notAnd(value: BitArray, notify = true): BitArray {
    if (this._length != value._length)
      throw new Error('Array lengths differ.');

    for (let i = 0, len = this.lengthInInts; i < len; i++)
      this._data[i] = (~this._data[i]) & value._data[i];

    this.incrementVersion(notify);
    return this;
  }

  /// Inverts all bit values in the current bitset
  not(notify = true): BitArray {
    for (let i = 0, len = this.lengthInInts; i < len; i++)
      this._data[i] = ~this._data[i];

    this.incrementVersion(notify);
    return this;
  }

  /// Performs the bitwise OR operation on the elements in the current bitset
  /// against the corresponding elements in the specified bitset.
  or(value: BitArray, notify = true) {
    if (this._length != value._length)
      throw new Error('Array lengths differ.');

    for (let i = 0, len = this.lengthInInts; i < len; i++)
      this._data[i] |= value._data[i];

    this.incrementVersion(notify);
    return this;
  }

  /// Performs the bitwise exclusive OR operation on the elements in the current bitset
  /// against the corresponding elements in the specified bitset.
  xor(value: BitArray, notify = true) {
    if (this._length != value._length)
      throw new Error('Array lengths differ.');

    for (let i = 0, len = this.lengthInInts; i < len; i++)
      this._data[i] ^= value._data[i];

    this.incrementVersion(notify);
    return this;
  }

  /// Inserts n 0-bits at position pos, resizing self and shifting bits appropriately.
  insertAt(pos: number, n: number, flag = false): void {
    this.assureInRange(pos, 0, this._length, 'pos');

    if (n == 0) return;

    //TODO: optimize
    //the most primitive implementation, optimize it later!

    // beginUpdate();
    const oldlength = this._length;
    this.setLength(this._length + n);

    //if (!contains(!flag)) return; // nothing to do

    for (let i = oldlength - 1; i >= pos; i--)
      this.setBit(i + n, this.getBit(i));

    for (let i = pos; i < pos + n; i++)
      this.setBit(i, flag);

    // endUpdate();
  }

  /// Deletes n bits beginning at position pos, resizing self and shifting remaining
  /// bits appropriately.
  removeAt(pos: number, n = 1): void {
    // the most primitive implementation, optimize it later!
    if (n < 0)
      throw new Error('n cannot be negative');

    this.assureInRange(pos, 0, this._length - n, 'pos');

    if (this.contains(true)) {
      for (let i = pos; i < this._length - n; i++)
        this.setBit(i, this.getBit(i + n));
    }

    this.setLength(this._length - n);
  }

  removeByMask(mask: BitArray, flag = true): BitArray {
    if (this._length != mask.length)
      throw new Error('length != mask.length');

    if (mask == this) { // no need to iterate
      this.setLength(mask.countBits(!flag));
      this.setAll(!flag);
    } else {
      let dstIdx = 0;

      for (let srcIdx = -1; (srcIdx = mask.findNext(srcIdx, !flag)) != -1;)
        this.setFast(dstIdx++, this.getBit(srcIdx));

      this._length = dstIdx;
      this._version++;
    }

    return this;
  }

  /// Similar to the [] operator.
  getBit(pos: number): boolean {
    return (this._data[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) != 0;
  }

  /// Similar to the [] operator.
  setBit(pos: number, bit: boolean, notify = true) {
    this.setFast(pos, bit);
    if (notify)
      this._version++;
    else
      this._version++;
  }

  /// Sets [i]-th bit to [value], does not check bounds, does not increment version
  setFast(i: number, value: boolean): void {
    if (value)
      this._data[Math.floor(i / 0x20)] |= 1 << (i & 0x1f);
    else
      this._data[Math.floor(i / 0x20)] &= ~(1 << (i & 0x1f));
  }

  setTrue(pos: number): void {
    this._data[Math.floor(pos / 0x20)] |= 1 << (pos & 0x1f);
  }

  setFalse(pos: number) {
    this._data[Math.floor(pos / 0x20)] &= ~(1 << (pos & 0x1f));
  }

  trueCount(): number {
    return this.countBits(true);
  }

  falseCount(): number {
    return this.countBits(false);
  }

  /// Counts bits of the specified value.
  countBits(value: boolean): number {
    if (this._length == 0) return 0;

    if (this._selectedCountVersion != this._version) {
      this._selectedCount = 0;
      const len = this.lengthInInts;
      let i = 0;
      for (; i < len - 1; i++) {
        for (let k = this._data[i]; k != 0; k >>>= 8) { //todo: cast data[i] to uint
          this._selectedCount += BitArray._onBitCount[k & 0xff];
        }
      }

      // The last int.
      let k = this._data[i];
      const remainingBits = this._length & 0x1f;
      if (remainingBits != 0) /* if remainingBits == 0, the last int is fully used and ALL bits should be left as is */
        k &= ~((4294967295) << remainingBits);

      for (; k != 0; k >>>= 8)
        this._selectedCount += BitArray._onBitCount[k & 0xff];

      this._selectedCountVersion = this._version;
    }

    return (value ? this._selectedCount : this._length - this._selectedCount);
  }

  /// Returns a number of set bits where also [check] is true
  countWhere(check: Function): number {
    let result = 0;
    if (this.trueCount() == this._length) {
      for (let i = 0; i < this._length; i++)
        result += check(i) ? 1 : 0;
    } else {
      for (let i = -1; (i = this.findNext(i, true)) != -1;)
        result += check(i) ? 1 : 0;
    }
    return result;
  }

  /// Performs bit "and" and counts bits of the specified value, without bitset modification.
  andWithCountBits(second: BitArray, value: boolean): number {
    if (this._length == 0) return 0;

    let count = 0;
    const len = this.lengthInInts;
    let i = 0;
    for (; i < len - 1; i++) {
      for (let k = this._data[i] & second._data[i]; k != 0; k >>>= 8)
        count += BitArray._onBitCount[k & 0xff];
    }

    // The last int.
    let k = this._data[i] & second._data[i];
    const remainingBits = this._length & 0x1f;
    if (remainingBits != 0)
      k &= ~((4294967295) << remainingBits);
    for (; k != 0; k >>>= 8)
      count += BitArray._onBitCount[k & 0xff];

    return (value ? count : this._length - count);
  }

  clear(): void {
    this.setLength(0);
  }

  contains(value: boolean): boolean {
    return this.findNext(-1, value) >= 0;
  }

  get allTrue() {
    return this.countBits(true) == this._length;
  }

  get allFalse() {
    return this.countBits(false) == this._length;
  }

  get anyTrue() {
    return this.countBits(true) > 0;
  }

  get anyFalse() {
    return this.countBits(false) > 0;
  }

  /// Returns the position of the next bit of the specified value, starting from the specified position.
  /// Returns -1, if there are no such bits.
  findNext(index: number, value = true): number {
    this.assureInRange(index, -1, this._length, 'index');

    if (index >= this._length - 1) return -1;
    index = index < 0 ? 0 : index + 1; // skip start
    let unusedBits = index & 0x1f;
    const numInts = this.lengthInInts;

    for (let i = Math.floor(index / 32); i < numInts; i++) {
      let k = (value ? this._data[i] : ~this._data[i]); // uint cast
      if (unusedBits != 0) {
        k &= ((0xffffffff << unusedBits) & 0xffffffff);
        unusedBits = 0;
      } else if (!value && k == -4294967296) /* looking for false, all bits are set */{
        continue;
      }

      for (let j = 0; k != 0; j += 8, k >>>= 8) {
        const p = BitArray._firstOnBit[k & 0xff];
        if (p >= 0) {
          index = p + (i * 32) + j;
          if (index >= this._length) return -1;
          return index;
        }
      }
    }
    return -1;
  }

  /// Finds previous bit of the specified value in the bitset.
  findPrev(index: number, value = true): number {
    if (index == 0) return -1;
    this.assureInRange(index, -1, this._length, 'index');

    index = index < 0 ? this._length - 1 : index - 1; // skip start

    const lastIntIdx = Math.floor(index / 0x20);
    let remainingBits = (index + 1) & 0x1f;

    for (let i = lastIntIdx; i >= 0; i--) {
      let k = (value ? this._data[i] : ~this._data[i]); // cast
      if (remainingBits != 0) {
        k &= ~((4294967295) << remainingBits);
        remainingBits = 0;
      }
      for (let j = 24; k != 0; j -= 8, k <<= 8) {
        const p = BitArray._lastOnBit[k >>> 0x18];
        if (p >= 0)
          return p + (i * 32) + j;
      }
    }
    return -1;
  }
}
