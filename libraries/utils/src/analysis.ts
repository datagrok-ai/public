import { randomInt } from "./operations";

export function getDiverseSubset(length: number, n: number, dist: (i1: number, i2: number) => number) {
  function maxBy(
    values: number[], 
    orderBy: (i: number) => number,
  ) {
    let maxValue;
    let maxOrderBy;

    for (const element of values) {
      let elementOrderBy = orderBy(element);
      if (maxOrderBy == null || elementOrderBy <= maxOrderBy) {
        maxValue = element;
        maxOrderBy = elementOrderBy;
      }
    }
    return maxValue;  
  }

  function includeRange(length: number, subset: number[]) {
    let res = [];
    for (let i = 0; i < length; ++i) {
      if (!subset.includes(i))
        res.push(i);
    }
    return res;
  }

  function minId(array: number[], origArray: number[]) {
    let mn = 1e9, id = 0;
    for (let i = 0; i < array.length; ++i) {
      if (mn > array[i]) {
        mn = array[i];
        id = i;
      }
    }
    return origArray[id];
  }

  let subset = [randomInt(length - 1)];
  while (subset.length < n) {
    let idx = maxBy(
      includeRange(length, subset),
      (i) => minId(subset.map(function (val, index) {
        return dist(i, val); 
      }), subset));
    subset.push(idx ? idx : -1);
  }
  return subset;
}