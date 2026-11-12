import {randomInt} from '@datagrok-libraries/utils/src/random';

/**
 * The diverse subset `getDiverseSubset` (@datagrok-libraries/utils similarity-metrics) picks, the same indices in the
 * same order for the same random start: a random first index, then, one at a time, the index farthest from the subset
 * (the largest of its distances' minimum; the lowest index among equals). Each index's minimum distance is kept and
 * updated with the last index added, so a pick costs `length` distances instead of `length` times the subset's size,
 * and nothing is allocated per index: Chem's diversity search over a large column took seconds.
 */
export function diverseSubset(length: number, n: number, dist: (i1: number, i2: number) => number): number[] {
  const subset = [randomInt(length - 1)];
  const inSubset = new Uint8Array(Math.max(length, 0));
  if (subset[0] >= 0 && subset[0] < length)
    inSubset[subset[0]] = 1;
  const minDist = new Float64Array(Math.max(length, 0)).fill(Infinity);
  while (subset.length < n) {
    const last = subset[subset.length - 1];
    let best = -1;
    let bestDist = 0;
    for (let i = 0; i < length; i++) {
      if (inSubset[i])
        continue;
      const d = Math.min(minDist[i], dist(i, last));
      minDist[i] = d;
      if (best === -1 || d > bestDist) {
        best = i;
        bestDist = d;
      }
    }
    if (best === -1) // nothing left to add (`getDiverseSubset` would wait forever)
      break;
    subset.push(best);
    inSubset[best] = 1;
  }
  return subset;
}
