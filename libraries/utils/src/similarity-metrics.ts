import BitArray from './bit-array'

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

export function tanimotoSimilarity(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  if (total == 0)
    return 1.0;
  const common = x.andWithCountBits(y, true);
  return common / (total - common);
}

export function diceSimilarity(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  if (total == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return 2 * common / total;
}

export function cosineSimilarity(x: BitArray, y: BitArray) {
  const total = x.trueCount() * y.trueCount();
  if (total == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / Math.sqrt(total);
}

export function sokalSimilarity(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  const common = x.andWithCountBits(y, true);
  return common / (2 * total - 3 * common);
}

export function kulczynskiSimilarity(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  const totalProd = x.trueCount() * y.trueCount();
  if (totalProd == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return (common * total) / (2 * totalProd);
}

export function mcConnaugheySimilarity(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  const totalProd = x.trueCount() * y.trueCount();
  if (totalProd == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return (common * total - totalProd) / totalProd;
}

export function asymmetricSimilarity(x: BitArray, y: BitArray) {
  const min = Math.min(x.trueCount(), y.trueCount());
  if (min == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / min;
}

export function braunBlanquetSimilarity(x: BitArray, y: BitArray) {
  const max = Math.max(x.trueCount(), y.trueCount());
  if (max == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / max;
}

export function russelSimilarity(x: BitArray, y: BitArray) {
  if (x.length == 0)
    return 0.0;
  const common = x.andWithCountBits(y, true);
  return common / x.length;
}

export function rogotGoldbergSimilarity(x: BitArray, y: BitArray) {
  let common = x.andWithCountBits(y, true);
  let total = x.countBits(true) + y.countBits(true);
  let len = x.length;
  let diff = len - total + common;
  if ((common == len) || (diff == len))
    return 1.0;
  else
    return common / total + diff / (2 * len - total);
}