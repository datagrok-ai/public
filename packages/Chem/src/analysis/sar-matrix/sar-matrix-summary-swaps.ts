/* Matched-pair swap pooling. A swap is one component value exchanged for another with every other
   component held fixed, so a pool is a grouping on all-but-one component: each component's measured
   pairs come out of one walk, with no rebuild. Pure functions of their arguments. */
import {SarMatrix} from './sar-matrix-types';
import {keepExtremes, RowCell, SUM_ROWS, SWAP_MIN_PAIRS, SWAP_MIN_SERIES, SWAP_ROW_CAP, SwapPool,
  SwapSide} from './sar-matrix-summary-types';

/**
 * Measured pairs inside ONE row index, pooled on the unordered fragment pair.
 *
 * The row index, not the core: a row is keyed on `[coreSmiles, ...foldedValues]`, so two rows can
 * share a core and differ at a folded position — pooling per core would break "everything else
 * identical", which is the whole claim of the card.
 */
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

/**
 * Fold one measured pair into its pool, or decline it. Null when the two sides carry the same value,
 * when a fold cannot be taken, or when this compound pair is already counted.
 *
 * The caller records `best` rather than this: a pair found through a matrix can name the two cells it
 * came from, and a pair found by grouping the component columns cannot.
 */
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
  // A fold multiplies, so its statistics are geometric: pooling the log of the ratio makes the mean,
  // min and max the fold statistics a chemist would quote, and lets one accumulator serve both
  // scales. A non-positive raw value has no fold at all.
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

/**
 * Measured swaps for every component at once, read off the component values of each measured cell
 * rather than off one matrix's columns.
 *
 * Nothing is rebuilt to reach another component. A swap is two compounds alike in every component but
 * one; where the components are given as columns that is a grouping, not a decomposition — group the
 * measured cells on every component except the one being swapped, and each pair inside a group that
 * differs in it is a matched pair. The matrix columns were only ever one route to the same thing, and
 * being one route is what made the other components cost a rebuild.
 */
export function poolRoleSwaps(roleNames: string[], roleValues: string[][], roleActivity: number[],
  roleMol: number[], dir: number, log: boolean): Map<string, SwapPool[]> {
  const out = new Map<string, SwapPool[]>();
  for (let r = 0; r < roleNames.length; r++) {
    const buckets = new Map<string, number[]>();
    // One compound occupies a cell in every tier that folded it, so the same partner tuple arrives
    // several times; the first occurrence can form every pair the later ones could.
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

/** The swap's worth in its better direction: what it bought in EVERY pair we have. Ranking a mean
 *  over many small pools selects the noisiest pool instead. */
function swapScore(pool: SwapPool): number {
  return Math.max(pool.min, -pool.max);
}

/** The single widest measured move anywhere, gate or no gate — so a card with no qualifying pool can
 *  still say what the data does hold rather than only what it does not. */
export function largestSwap(pools: Map<string, SwapPool>): SwapPool | null {
  let best: SwapPool | null = null;
  for (const pool of pools.values()) {
    const reach = Math.max(pool.max, -pool.min);
    if (best === null || reach > Math.max(best.max, -best.min))
      best = pool;
  }
  return best;
}

export function rankSwaps(pools: Map<string, SwapPool>): SwapPool[] {
  return [...pools.values()]
    .filter((pool) => pool.n >= SWAP_MIN_PAIRS && pool.roots.size >= SWAP_MIN_SERIES)
    .sort((a, b) => swapScore(b) - swapScore(a) ||
      (a.from < b.from ? -1 : a.from > b.from ? 1 : a.to < b.to ? -1 : 1))
    .slice(0, SUM_ROWS);
}
