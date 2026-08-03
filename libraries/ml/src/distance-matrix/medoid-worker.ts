import {getAggregationFunction, isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DistanceAggregationMethod} from './types';

/** Computes, for a single cluster, every member's mean distance to the other members. Distances are
 * computed on the fly and reduced to per-element row sums, so the worker never materializes a
 * distance matrix (O(n) memory, O(n^2) compute). `values` are already subset to the cluster members;
 * the returned `meanDistances` array is aligned with that subset (the medoid is its argmin, and the
 * full array lets the caller rank members for Top-N extraction). */
onmessage = (event) => {
  const {values, fnNames, opts, weights, aggregationMethod}: {
    values: any[][], fnNames: KnownMetrics[], opts: {[_: string]: any}[],
    weights: number[], aggregationMethod: DistanceAggregationMethod,
  } = event.data;
  try {
    values.forEach((v, colIdx) => {
      if (isBitArrayMetric(fnNames[colIdx])) {
        for (let i = 0; i < v.length; ++i) {
          if (isNil(v[i])) continue;
          values[colIdx][i] = new BitArray(values[colIdx][i]._data, values[colIdx][i]._length);
        }
      }
    });

    const n = values[0].length;
    const colCount = values.length;
    const distanceFns = new Array(colCount).fill(null)
      .map((_, i) => new Measure(fnNames[i]).getMeasure(opts[i]));
    const rowSums = new Float64Array(n);

    if (colCount === 1) {
      // Single feature: argmin of row sums is invariant under the per-column normalization that
      // multi-column distances apply, so we skip it and accumulate raw distances in one pass.
      const fn = distanceFns[0];
      const col = values[0];
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const value = !isNil(col[i]) && !isNil(col[j]) ? fn(col[i], col[j]) : 1;
          rowSums[i] += value;
          rowSums[j] += value;
        }
      }
    } else {
      // Multi-column: normalize each column's distances (cluster-local min/max) before aggregating,
      // matching the per-column normalize + aggregate that calcDistanceMatrix uses to build the tree.
      const pairCount = n * (n - 1) / 2;
      const perColDist = new Array(colCount).fill(null).map(() => new Float32Array(pairCount));
      const mins = new Array<number>(colCount).fill(Number.POSITIVE_INFINITY);
      const maxs = new Array<number>(colCount).fill(Number.NEGATIVE_INFINITY);
      let p = 0;
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          for (let c = 0; c < colCount; c++) {
            const d = !isNil(values[c][i]) && !isNil(values[c][j]) ? distanceFns[c](values[c][i], values[c][j]) : 1;
            perColDist[c][p] = d;
            if (d < mins[c]) mins[c] = d;
            if (d > maxs[c]) maxs[c] = d;
          }
          p++;
        }
      }
      const aggregate = getAggregationFunction(aggregationMethod, weights);
      const ranges = mins.map((mn, c) => (maxs[c] - mn) || 1);
      const norm = new Array<number>(colCount).fill(0);
      p = 0;
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          for (let c = 0; c < colCount; c++)
            norm[c] = (perColDist[c][p] - mins[c]) / ranges[c];
          const value = aggregate(norm);
          rowSums[i] += value;
          rowSums[j] += value;
          p++;
        }
      }
    }

    const denom = n > 1 ? n - 1 : 1;
    const meanDistances = new Float64Array(n);
    for (let i = 0; i < n; i++)
      meanDistances[i] = rowSums[i] / denom;
    postMessage({meanDistances});
  } catch (e) {
    postMessage({error: e});
  }
};
