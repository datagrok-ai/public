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
 * This kernel only WRITES the substitution. What each rule then means for
 * AUC/AUMC, the λz regression, Cmax and Tlag is the contract on {@link BlqRule}
 * and is implemented by `augmentProfile` through two masks: the substituted
 * values (`set-zero`, `set-half-lloq`) reach the integrator and the λz
 * eligibility filter as the values written here; `exclude`d and `NaN` points
 * form the drop set; Cmax/Tmax and Tlag read the BLQ mask under every rule.
 * (Before GROK-20960 the integrator read the BLQ mask too, so every rule
 * integrated as `exclude`.)
 *
 * `blqMask` marks below-LOQ points with `1`. `lloq` is a per-row array or a
 * scalar applied to all rows, and is read only by `set-half-lloq`.
 *
 * @param _cmaxIdx - Unused; kept for signature stability. Pass `0`.
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
