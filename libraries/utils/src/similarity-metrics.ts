import BitArray from './bit-array';
import {randomInt} from './random';

export function getDiverseSubset(fingerprints: BitArray[], length: number, n: number, dist: (i1: number, i2: number) => number) {
  function maxBy(values: IterableIterator<number>, orderBy: (i: number) => number) {
    let maxValue = null;
    let maxOrderBy = null;

    for (const element of values) {
      const elementOrderBy = orderBy(element);
      if (maxOrderBy == null || elementOrderBy > maxOrderBy) {
        maxValue = element;
        maxOrderBy = elementOrderBy;
      }
    }
    return maxValue;
  }

  let randomIdx = randomInt(length - 1);
  while (fingerprints[randomIdx].allFalse)
    randomIdx = randomInt(length - 1);
    
  const subset = [randomIdx];
  const complement = new Set();

  for (let i = 0; i < length; ++i) {
    if (!subset.includes(i) && !fingerprints[i].allFalse)
      complement.add(i);
  }

  while (subset.length < n) {
    const idx = maxBy(
      complement.values() as IterableIterator<number>,
      (i) => Math.min.apply(Math, subset.map(function(val, index) {
        return dist(i, val);
      })));
    if (idx != null) {
      subset.push(idx);
      complement.delete(idx);
    }
  }
  return subset;
}
