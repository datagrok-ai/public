export class BitSetFixedArray {
  private _buf: ArrayBuffer;
  private _words: Uint32Array;
  private _wordsPerItem: number;
  private readonly _bitsPerWord = 32;
  private readonly _bytesPerWord = (this._bitsPerWord / 8);
  private _numItems: number | null;
  private _bitsPerItem: number | null;
  constructor(bitsPerItem: number, numItems: number) {
    this._bitsPerItem = bitsPerItem;
    this._wordsPerItem = Math.ceil(bitsPerItem / this._bitsPerWord);
    this._buf = new ArrayBuffer(this._wordsPerItem * this._bytesPerWord * numItems);
    this._words = new Uint32Array(this._buf);
    this._numItems = numItems;
  }

  get length(): number {
    return this._numItems!;
  }

  setTrue(item: number, index: number): void {
    this._words[item * this._wordsPerItem + (index >>> 5)] |= (1 << index);
  }

  setFalse(item: number, index: number): void {
    this._words[item * this._wordsPerItem + (index >>> 5)] &= ~(1 << index);
  }

  private static readonly _onBitCount = [
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
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
  ];

  public count(item: number) {
    let _selectedCount = 0;
    const len = this._wordsPerItem;
    let i = item * len;
    for (; i < item * len + len - 1; i++) {
      for (let k = this._words[i]; k != 0; k >>>= 8) {
        _selectedCount += BitSetFixedArray._onBitCount[k & 0xff];
      }
    }
    let k = this._words[i];
    const remainingBits = this._bitsPerItem! & 0x1f;
    if (remainingBits != 0) {
      k &= ~((4294967295) << remainingBits);
    }
    for (; k != 0; k >>>= 8) {
      _selectedCount += BitSetFixedArray._onBitCount[k & 0xff];
    }
    return _selectedCount;
  }

  public andWithCountBits(item: number, other: BitSetFixedArray) {
    let _selectedCount = 0;
    const len = this._wordsPerItem;
    let i = item * len; let j = 0;
    for (; i < item * len + len - 1; i++, j++) {
      for (let k = this._words[i] & other._words[j]; k != 0; k >>>= 8) {
        _selectedCount += BitSetFixedArray._onBitCount[k & 0xff];
      }
    }
    let k = this._words[i] & other._words[j];
    const remainingBits = this._bitsPerItem! & 0x1f;
    if (remainingBits != 0) {
      k &= ~((4294967295) << remainingBits);
    }
    for (; k != 0; k >>>= 8) {
      _selectedCount += BitSetFixedArray._onBitCount[k & 0xff];
    }
    return _selectedCount;
  }

  public entails(item: number, other: BitSetFixedArray) {
    const len = this._wordsPerItem;
    let i = item * len; let j = 0;
    for (; i < item * len + len; i++, j++) {
      if ((other._words[j] & this._words[i]) !== other._words[j]) {
        return false;
      }
    }
    return true;
  }
};
