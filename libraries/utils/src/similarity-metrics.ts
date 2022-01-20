import BitArray from './bit-array'

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