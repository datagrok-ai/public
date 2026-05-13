/**
 * Fisher's exact test for a 2×2 contingency table.
 */

import {combinationln} from '../distributions';
import {FisherResult} from '../types';

/**
 * 2×2 Fisher exact test.
 *
 * Table layout:
 *
 * ```
 *        col0  col1 │ row_total
 *  row0:   a     b  │  R0 = a+b
 *  row1:   c     d  │  R1 = c+d
 *  ─────────────────┤
 *          C0    C1 │  N
 * ```
 *
 * Two-sided p-value uses the "minlike" definition (matches scipy default):
 * sum of probabilities of all tables with PMF ≤ that of the observed table.
 *
 * `oddsRatio` is `(a·d)/(b·c)`. Returns `Infinity` when `b·c = 0` and `a·d > 0`,
 * `0` when `a·d = 0` and `b·c > 0`. Throws on negative inputs.
 */
export function fisherExact2x2(table: ReadonlyArray<ReadonlyArray<number>>): FisherResult {
  if (table.length !== 2 || table[0].length !== 2 || table[1].length !== 2)
    throw new Error('fisherExact2x2: expected 2×2 table.');
  const a = table[0][0];
  const b = table[0][1];
  const c = table[1][0];
  const d = table[1][1];
  if (a < 0 || b < 0 || c < 0 || d < 0)
    throw new Error('fisherExact2x2: cells must be non-negative.');

  const R0 = a + b;
  const R1 = c + d;
  const C0 = a + c;
  const N = R0 + R1;

  let oddsRatio: number;
  if (b * c === 0) oddsRatio = a * d > 0 ? Infinity : 0;
  else oddsRatio = (a * d) / (b * c);

  const aMin = Math.max(0, C0 - R1);
  const aMax = Math.min(R0, C0);
  const lnTotal = combinationln(N, C0);
  const pObs = pmf(a, R0, R1, C0, lnTotal);

  let pTwo = 0;
  let pGreater = 0;
  let pLess = 0;
  for (let x = aMin; x <= aMax; x++) {
    const px = pmf(x, R0, R1, C0, lnTotal);
    if (px <= pObs + 1e-15) pTwo += px;
    if (x >= a) pGreater += px;
    if (x <= a) pLess += px;
  }

  return {oddsRatio, pValue: pTwo, pGreater, pLess};
}

/** Hypergeometric PMF: P(X = x) given fixed margins. */
function pmf(x: number, R0: number, R1: number, C0: number, lnTotal: number): number {
  const cVal = C0 - x;
  if (x < 0 || cVal < 0 || x > R0 || cVal > R1) return 0;
  return Math.exp(combinationln(R0, x) + combinationln(R1, cVal) - lnTotal);
}
