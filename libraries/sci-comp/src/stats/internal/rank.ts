/**
 * Ranking with average-rank tie correction (matches scipy `rankdata(method='average')`).
 */

/**
 * Average ranks: tied values share the mean of the ranks they would occupy.
 *
 * Returns a Float64Array of the same length as the input, with ranks in the
 * range [1, n].
 */
export function rankdata(x: ArrayLike<number>): Float64Array {
  const n = x.length;
  const order = new Int32Array(n);
  for (let i = 0; i < n; i++) order[i] = i;
  // Sort indices by value (stable for ties — irrelevant for average rank)
  const arr = order.slice().sort((a, b) => x[a] - x[b]);

  const ranks = new Float64Array(n);
  let i = 0;
  while (i < n) {
    let j = i;
    // Find end of the tie-run at position i
    while (j + 1 < n && x[arr[j + 1]] === x[arr[i]]) j++;
    // Average rank for positions i..j (1-based ranks)
    const r = (i + j + 2) / 2; // = ((i+1) + (j+1)) / 2
    for (let k = i; k <= j; k++) ranks[arr[k]] = r;
    i = j + 1;
  }
  return ranks;
}

/**
 * Tie-correction term for the Mann-Whitney U normal approximation.
 *
 * Returns Σ (tᵢ³ − tᵢ) over tie groups in the *combined* sample. Used in:
 *   varU = n1·n2·(n+1)/12  −  n1·n2 · T / (12·n·(n−1))
 * where T is this returned sum and n = n1 + n2.
 */
export function tieCorrection(combined: ArrayLike<number>): number {
  const n = combined.length;
  if (n < 2) return 0;
  // Sort a copy
  const arr = new Float64Array(n);
  for (let i = 0; i < n; i++) arr[i] = combined[i];
  arr.sort();
  let T = 0;
  let i = 0;
  while (i < n) {
    let j = i;
    while (j + 1 < n && arr[j + 1] === arr[i]) j++;
    const t = j - i + 1;
    if (t > 1) T += t * t * t - t;
    i = j + 1;
  }
  return T;
}
