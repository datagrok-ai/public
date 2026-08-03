/**
 * Bonferroni multiple-comparison correction.
 */

/**
 * Multiply each p-value by the number of tests, capped at 1.0.
 *
 * - `null` p-values pass through unchanged.
 * - When `nTests` is omitted, it defaults to the count of non-null entries.
 * - When `nTests` is `0`, the input is returned unchanged.
 *
 * @param pValues  Raw p-values (possibly with `null` for missing).
 * @param nTests   Override for the multiplicity factor.
 */
export function bonferroniCorrect(
  pValues: readonly (number | null)[],
  nTests?: number,
): (number | null)[] {
  const n = nTests !== undefined ?
    nTests :
    pValues.reduce((acc: number, p) => p !== null ? acc + 1 : acc, 0);
  if (n === 0)
    return pValues.slice();
  return pValues.map((p) => p === null ? null : Math.min(p * n, 1.0));
}
