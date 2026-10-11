/* Matched-pair swap pooling. A swap is one component value exchanged for another with every other
   component held fixed, so a pool is a grouping on all-but-one component: each component's measured
   pairs come out of one walk, with no rebuild. Pure functions of their arguments. */
import {SarMatrix} from '../sar-matrix-types';
import {RowCell, SUM_ROWS, SWAP_MIN_PAIRS, SWAP_MIN_SERIES, SWAP_ROW_CAP, SwapPool,
  SwapSide} from './sar-matrix-summary-types';

/** Both extremes by value once a group exceeds the cap, so only mid-range pairs are lost. */
function keepExtremes<T>(items: T[], valueOf: (item: T) => number): T[] {
  if (items.length <= SWAP_ROW_CAP)
    return items;
  const half = SWAP_ROW_CAP >> 1;
  const sorted = [...items].sort((a, b) => valueOf(a) - valueOf(b));
  return [...sorted.slice(0, half), ...sorted.slice(sorted.length - half)];
}

/** Measured pairs inside one row, pooled on the unordered fragment pair. Per row, not per core: two rows
 *  can share a core and differ at a folded position. */
export function poolSwaps(pools: Map<string, SwapPool>, matrix: SarMatrix, ri: number, cells: RowCell[],
  dir: number, log: boolean, root: string): void {
  const sampled = cells.length > SWAP_ROW_CAP;
  const kept = keepExtremes(cells, (cell) => cell.value);
  for (let a = 0; a < kept.length; a++) {
    for (let b = a + 1; b < kept.length; b++) {
      const added = addSwapPair(pools,
        {value: matrix.columns[kept[a].ci].substSmiles, activity: kept[a].value, mol: kept[a].molIdx},
        {value: matrix.columns[kept[b].ci].substSmiles, activity: kept[b].value, mol: kept[b].molIdx},
        root, sampled, dir, log);
      if (added === null)
        continue;
      const lower = added.flip ? kept[b] : kept[a];
      const upper = added.flip ? kept[a] : kept[b];
      if (added.pool.best === null || Math.abs(added.delta) > Math.abs(added.pool.best.delta))
        added.pool.best = {matrix, ri, ciFrom: lower.ci, ciTo: upper.ci, delta: added.delta};
    }
  }
}

/** Fold one measured pair into its pool; null when the sides match, no fold exists, or it is counted.
 *  The caller records `best`, since only a pair found through a matrix has cells to name. */
function addSwapPair(pools: Map<string, SwapPool>, a: SwapSide, b: SwapSide,
  root: string, sampled: boolean, dir: number, log: boolean):
  {pool: SwapPool, delta: number, flip: boolean} | null {
  const {value: sa, activity: va, mol: molA} = a;
  const {value: sb, activity: vb, mol: molB} = b;
  if (sa === sb)
    return null;
  const flip = sa > sb;
  const lowerVal = flip ? vb : va;
  const upperVal = flip ? va : vb;
  // A fold multiplies, so the log ratio is pooled; a non-positive raw value has no fold.
  if (!log && (lowerVal <= 0 || upperVal <= 0))
    return null;
  const delta = log ? dir * (upperVal - lowerVal) :
    dir * (Math.log10(upperVal) - Math.log10(lowerVal));
  const key = `${flip ? sb : sa}\0${flip ? sa : sb}`;
  let pool = pools.get(key);
  if (pool === undefined) {
    pool = {from: flip ? sb : sa, to: flip ? sa : sb, n: 0, sum: 0, min: Infinity, max: -Infinity,
      nUp: 0, roots: new Set(), seen: new Set(), best: null, sampled: false};
    pools.set(key, pool);
  }
  const pair = `${Math.min(molA, molB)}:${Math.max(molA, molB)}`;
  if (pool.seen.has(pair))
    return null;
  pool.seen.add(pair);
  pool.n++;
  pool.sum += delta;
  pool.min = Math.min(pool.min, delta);
  pool.max = Math.max(pool.max, delta);
  if (delta > 0)
    pool.nUp++;
  pool.roots.add(root);
  pool.sampled = pool.sampled || sampled;
  return {pool, delta, flip};
}

/** Measured swaps for every component at once: group the measured cells on all components but one, and
 *  each pair in a group that differs in it is a matched pair. No rebuild per component. */
export function poolRoleSwaps(roleNames: string[], roleValues: string[][], roleActivity: number[],
  roleMol: number[], dir: number, log: boolean): Map<string, SwapPool[]> {
  const out = new Map<string, SwapPool[]>();
  for (let r = 0; r < roleNames.length; r++) {
    const buckets = new Map<string, number[]>();
    // One compound recurs in every tier that folded it; its first occurrence forms every pair.
    const placed = new Map<string, Set<number>>();
    for (let k = 0; k < roleActivity.length; k++) {
      let key = '';
      for (let j = 0; j < roleNames.length; j++) {
        if (j !== r)
          key += `${roleValues[j][k]}\u0001`;
      }
      let bucket = buckets.get(key);
      if (bucket === undefined) {
        bucket = [];
        buckets.set(key, bucket);
        placed.set(key, new Set());
      }
      const seen = placed.get(key)!;
      if (seen.has(roleMol[k]))
        continue;
      seen.add(roleMol[k]);
      bucket.push(k);
    }
    const pools = new Map<string, SwapPool>();
    for (const [key, bucket] of buckets) {
      const sampled = bucket.length > SWAP_ROW_CAP;
      const kept = keepExtremes(bucket, (i) => roleActivity[i]);
      for (let a = 0; a < kept.length; a++) {
        for (let b = a + 1; b < kept.length; b++) {
          addSwapPair(pools,
            {value: roleValues[r][kept[a]], activity: roleActivity[kept[a]], mol: roleMol[kept[a]]},
            {value: roleValues[r][kept[b]], activity: roleActivity[kept[b]], mol: roleMol[kept[b]]},
            key, sampled, dir, log);
        }
      }
    }
    out.set(roleNames[r], rankSwaps(pools));
  }
  return out;
}

/** What the swap bought in every pair, in its better direction; a mean would favour the noisiest pool. */
function swapScore(pool: SwapPool): number {
  return Math.max(pool.min, -pool.max);
}

export function rankSwaps(pools: Map<string, SwapPool>): SwapPool[] {
  return [...pools.values()]
    .filter((pool) => pool.n >= SWAP_MIN_PAIRS && pool.roots.size >= SWAP_MIN_SERIES)
    .sort((a, b) => swapScore(b) - swapScore(a) ||
      (a.from < b.from ? -1 : a.from > b.from ? 1 : a.to < b.to ? -1 : 1))
    .slice(0, SUM_ROWS);
}
