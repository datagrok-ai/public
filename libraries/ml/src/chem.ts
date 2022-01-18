import BitArray from '@datagrok-libraries/utils/src/bit-array'

function trueCount(value: any): number {
  if (value._length == 0) return 0;

  if (value._selectedCountVersion != value._version) {
    value._selectedCount = 0;
    const len = Math.floor((value._length + 0x1f) / 0x20);
    let i = 0;
    for (; i < len - 1; i++) {
      for (let k = value._data[i]; k != 0; k >>>= 8) { //todo: cast data[i] to uint
        value._selectedCount += BitArray._onBitCount[k & 0xff];
      }
    }

    // The last int.
    let k = value._data[i];
    const remainingBits = value._length & 0x1f;
    if (remainingBits != 0) {// if remainingBits == 0, the last int is fully used and ALL bits should be left as is.
      k &= ~((4294967295) << remainingBits);
    }

    for (; k != 0; k >>>= 8) {
      value._selectedCount += BitArray._onBitCount[k & 0xff];
    }

    value._selectedCountVersion = value._version;
  }

  return (value ? value._selectedCount : value._length - value._selectedCount);
}

function andWithCountBits(x: any, y: any, value: boolean) {
  if (x._length == 0) return 0;

  let count = 0; const len = Math.floor((x._length + 0x1f) / 0x20); let i = 0;
  for (; i < len - 1; i++) {
    for (var k = x._data[i] & y._data[i]; k != 0; k >>>= 8) {
      count += BitArray._onBitCount[k & 0xff];
    }
  }

  // The last int.
  var k = x._data[i] & y._data[i];
  const remainingBits = x._length & 0x1f;
  if (remainingBits != 0) {
    k &= ~((4294967295) << remainingBits);
  }
  for (; k != 0; k >>>= 8) {
    count += BitArray._onBitCount[k & 0xff];
  }

  return (value ? count : x._length - count);
}

export function tanimoto(x: BitArray, y: BitArray) {
  const total = trueCount(x) + trueCount(y);
  if (total == 0)
    return 1.0;
  const common = andWithCountBits(x, y, true);
  return common / (total - common);
}