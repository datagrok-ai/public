import BitArray from './bit-array'
import {randomInt} from "./random";

export const similarityMetric: {[name: string]: (x: BitArray, y: BitArray) => number} = {
  'Tanimoto': tanimotoSimilarity,
  'Dice': diceSimilarity,
  'Asymmetric': asymmetricSimilarity,
  'Braun-Blanquet': braunBlanquetSimilarity,
  'Cosine': cosineSimilarity,
  'Kulczynski': kulczynskiSimilarity,
  'Mc-Connaughey': mcConnaugheySimilarity,
  'Rogot-Goldberg': rogotGoldbergSimilarity,
  'Russel': russelSimilarity,
  'Sokal': sokalSimilarity,
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
  if ((common == len) || (diff == len)) return 1.0;
  else return common / total + diff / (2 * len - total);
}

export function getDiverseSubset(length: number, n: number, dist: (i1: number, i2: number) => number) {

  function maxBy(values: IterableIterator<number>, orderBy: (i: number) => number) {
    let maxValue = null;
    let maxOrderBy = null;

    for (const element of values) {
      let elementOrderBy = orderBy(element);
      if (maxOrderBy == null || elementOrderBy > maxOrderBy) {
        maxValue = element;
        maxOrderBy = elementOrderBy;
      }
    }
    return maxValue;  
  }

  let subset = [randomInt(length - 1)];
  let complement = new Set();

  for (let i = 0; i < length; ++i) {
    if (!subset.includes(i))
      complement.add(i);
  }

  while (subset.length < n) {
    let idx = maxBy(
      complement.values() as IterableIterator<number>,
      (i) => Math.min.apply(Math, subset.map(function (val, index) {
        return dist(i, val); 
      })));
    if (idx) {
      subset.push(idx);
      complement.delete(idx);
    }
  }
  return subset;
}
