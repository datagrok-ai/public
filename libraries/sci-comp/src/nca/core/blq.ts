import type {BlqRule, BlqStrategy, BlqProcessingResult} from './types';

/**
 * Apply a per-phase BLQ pre-processing strategy to one profile.
 *
 * The result is a fresh concentration array (the input is not mutated) plus
 * a list of indices that should be excluded from downstream computation.
 *
 * ## Phases (PKNCA semantics)
 *
 * Phases are determined by position relative to the **first** and **last**
 * measurable observation:
 *
 * - **preFirstMeasurable** — points before the first non-BLQ point.
 * - **embedded**            — points strictly between two non-BLQ points.
 * - **afterLast**           — the first BLQ point after the last non-BLQ.
 * - **consecutiveAfterLast** — second and subsequent BLQ points in the tail.
 *
 * Maps to PKNCA `conc.blq` `first` / `middle` / `last` (last-tail split is
 * an extension this core supports for finer-grained control).
 *
 * Reference: <https://billdenney.github.io/pknca/articles/Selection-of-Calculation-Intervals.html>
 *
 * ## Rules
 *
 * - `set-zero`      — concentration replaced with 0.
 * - `set-half-lloq` — concentration replaced with `lloq[i]/2` (or `lloq/2`
 *                     if `lloq` is scalar).
 * - `exclude`       — the index is added to `excluded`; the concentration
 *                     value is left as-is and the caller drops the row.
 * - `missing`       — concentration replaced with `NaN`; downstream
 *                     functions skip `NaN` entries.
 *
 * @param conc - Original concentrations. Not mutated.
 * @param blqMask - Same length as `conc`; `1` = below LOQ, `0` = measurable.
 * @param lloq - Per-row LLOQ array (same length as `conc`) or a scalar
 *               applied to all rows. Used only by `set-half-lloq`.
 * @param _cmaxIdx - Currently unused. Reserved for future Cmax-relative
 *                   phasing. Pass any value (typically the index of the
 *                   first measurable point or `findCmax` result).
 * @param strategy - The BLQ rule for each of the four phases.
 * @returns The processed concentration array and the list of excluded
 *          indices.
 */
export function applyBlqStrategy(
  conc: Float64Array,
  blqMask: Uint8Array,
  lloq: Float64Array | number,
  _cmaxIdx: number,
  strategy: BlqStrategy,
): BlqProcessingResult {
  const n = conc.length;
  const out = new Float64Array(conc);
  const excluded: number[] = [];

  let firstMeasurable = -1;
  let lastMeasurable = -1;
  for (let i = 0; i < n; i++) {
    if (blqMask[i] === 0) {
      if (firstMeasurable === -1) firstMeasurable = i;
      lastMeasurable = i;
    }
  }

  for (let i = 0; i < n; i++) {
    if (blqMask[i] === 0) continue;

    let rule: BlqRule;
    if (firstMeasurable === -1 || i < firstMeasurable)
      rule = strategy.preFirstMeasurable;
    else if (i > lastMeasurable) {
      rule = (i === lastMeasurable + 1) ?
        strategy.afterLast :
        strategy.consecutiveAfterLast;
    } else
      rule = strategy.embedded;

    applyRule(out, excluded, i, rule, lloq);
  }

  return {
    conc: out,
    excluded: Int32Array.from(excluded),
  };
}

function applyRule(
  out: Float64Array,
  excluded: number[],
  i: number,
  rule: BlqRule,
  lloq: Float64Array | number,
): void {
  switch (rule) {
  case 'set-zero':
    out[i] = 0;
    break;
  case 'set-half-lloq':
    out[i] = (typeof lloq === 'number' ? lloq : lloq[i]) / 2;
    break;
  case 'exclude':
    excluded.push(i);
    break;
  case 'missing':
    out[i] = NaN;
    break;
  }
}
