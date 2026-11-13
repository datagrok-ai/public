/**
 * `aucInterval` — the AUC of one profile over an arbitrary interval
 * `[start, end]`, with PKNCA `aucint.last` semantics (`pk.calc.aucint`,
 * `auc.type = "AUClast"`).
 *
 * It integrates EXACTLY the profile `computeNca` integrates for AUClast — the
 * same `augmentProfile` (BLQ rules, route-aware dose-time row, IV-bolus c0) and
 * the same measurable-point collection (drop set skipped, trailing
 * non-positive run trimmed) — and the same trapezoid kernel for
 * `rules.aucMethod` / `rules.compensatedSummation`. So
 * `aucInterval(inputs, rules, 0, Infinity).value === computeNca(...).values.aucLast`
 * bit for bit. Nothing here is new numerics except the bound interpolation.
 *
 * Bound rules (PKNCA 0.12.1 `pk.calc.aucint` → `interp.extrap.conc[.dose]` →
 * `interpolate.conc` / `extrapolate.conc`):
 *
 * - A bound ON a sample time uses that sample.
 * - A bound strictly between two samples is INTERPOLATED with the method the
 *   integration uses on that interval (`interpolate.conc` picks it from
 *   `choose_interval_method` on the whole profile): log
 *   (`exp(ln c1 + f·(ln c2 − ln c1))`, PKNCA `interpolate_conc_log`) on a log
 *   interval, linear (`c1 + f·(c2 − c1)`) otherwise. Under
 *   `linear-up-log-down` a log interval is a strictly descending one between
 *   two positive concentrations; under `linear` none is.
 * - The sub-profile `[start, bound points…, end]` is then integrated by the
 *   method kernel, which — like PKNCA's `choose_interval_method` on the
 *   interpolated times — decides each sub-interval afresh.
 * - **No extrapolation beyond tlast.** Past the last measurable point the
 *   concentration is 0 (`extrapolate.conc` with `auc.type = "AUClast"`), so an
 *   `end` past tlast integrates to tlast only (`truncatedAtTlast`). `end` may
 *   be `Infinity`.
 * - A `start` before the first point of the integrated profile adds no area.
 *   That first point is the dose time `t = 0` whenever the profile was
 *   augmented (always, unless the input carries times before 0), so the
 *   pre-dose concentration is 0 — PKNCA's `conc.origin = 0`.
 *
 * `log-linear` has no PKNCA counterpart (PKNCA's `lin-log` is a different
 * scheme); it follows the same rule — the bound is log-interpolated wherever
 * the kernel integrates log. Validated against PKNCA for `linear` and
 * `linear-up-log-down` (`fixtures/08_aucint.json`).
 */

import type {ProfileInputs, NcaRules, AucMethod} from './types';
import {isKnownRouteCode} from './types';
import {augmentProfile} from './augment';
import {collectMeasurable, pickAucFn} from './compute-nca';

/** Outcome of {@link aucInterval}. */
export interface AucIntervalResult {
  /**
   * AUC over `[startUsed, endUsed]`. `0` when the interval lies wholly outside
   * the measured range. `NaN` when the profile has no measurable point
   * (`computeNca` status `'failed'`), `inputs.route` is not a known `RouteCode`,
   * or the interval is invalid (`start` not finite, `end` NaN, or `end < start`).
   */
  readonly value: number;
  /** `start` clipped to the integrated profile `[tFirst, tLast]`. */
  readonly startUsed: number;
  /** `end` clipped to `[tFirst, tLast]` — what was actually integrated to. */
  readonly endUsed: number;
  /** Last measurable time of the profile (AUClast's end); `NaN` on failure. */
  readonly tLast: number;
  /**
   * `true` when `end > tLast`: the area stops at tlast (no extrapolation), so
   * the value is an AUC to tlast, not to `end`. A caller presenting
   * "AUC0–24" for a profile sampled to 12 h must say so.
   */
  readonly truncatedAtTlast: boolean;
}

const FAILED: AucIntervalResult = Object.freeze({
  value: NaN, startUsed: NaN, endUsed: NaN, tLast: NaN, truncatedAtTlast: false,
});

/**
 * AUC of one profile over `[start, end]` (PKNCA `aucint.last`). Stateless;
 * never throws. See the module header for the exact contract.
 *
 * Reads `rules.aucMethod`, `rules.blq` and `rules.compensatedSummation`; the
 * λz settings are irrelevant (no extrapolation).
 */
export function aucInterval(
  inputs: ProfileInputs, rules: NcaRules, start: number, end: number,
): AucIntervalResult {
  if (!Number.isFinite(start) || Number.isNaN(end) || end < start) return FAILED;
  if (!isKnownRouteCode(inputs.route)) return FAILED; // as computeNca: never guessed extravascular

  const aug = augmentProfile(inputs, rules.blq);
  if (aug === null) return FAILED;
  const {time, conc} = collectMeasurable(aug.time, aug.conc, aug.dropMask);
  const n = time.length;
  if (n === 0) return FAILED;

  const tFirst = time[0];
  const tLast = time[n - 1];
  const lo = Math.min(Math.max(start, tFirst), tLast);
  const hi = Math.min(Math.max(end, tFirst), tLast);
  const truncatedAtTlast = end > tLast;
  if (hi <= lo)
    return {value: 0, startUsed: lo, endUsed: hi, tLast, truncatedAtTlast};

  // Sub-profile: the bound at `lo`, every sample strictly inside, the bound at `hi`.
  const subT: number[] = [lo];
  const subC: number[] = [concAt(time, conc, lo, rules.aucMethod)];
  for (let i = 0; i < n; i++) {
    if (time[i] > lo && time[i] < hi) {
      subT.push(time[i]);
      subC.push(conc[i]);
    }
  }
  subT.push(hi);
  subC.push(concAt(time, conc, hi, rules.aucMethod));

  const aucFn = pickAucFn(rules.aucMethod, rules.compensatedSummation);
  const value = aucFn(Float64Array.from(subT), Float64Array.from(subC), 0, subT.length - 1);
  return {value, startUsed: lo, endUsed: hi, tLast, truncatedAtTlast};
}

/**
 * Concentration at `t` within `[time[0], time[n−1]]`: the sample when `t` is a
 * sample time, otherwise interpolated on the bracketing interval with the
 * method the AUC kernel integrates that interval with (the per-interval rule of
 * `auc.ts`, which equals PKNCA's `choose_interval_method` for `linear` and
 * `lin up/log down`).
 */
function concAt(
  time: Float64Array, conc: Float64Array, t: number, method: AucMethod,
): number {
  let j = 0;
  while (j < time.length - 1 && time[j + 1] <= t) j++;
  if (time[j] === t) return conc[j];
  const t1 = time[j];
  const t2 = time[j + 1];
  const c1 = conc[j];
  const c2 = conc[j + 1];
  const f = (t - t1) / (t2 - t1);
  const logInterval =
    method === 'linear-up-log-down' ? (c2 < c1 && c1 > 0 && c2 > 0) :
      method === 'log-linear' ? (c1 > 0 && c2 > 0 && c1 !== c2) :
        false;
  // PKNCA interpolate_conc_log / interpolate_conc_linear, verbatim.
  return logInterval ?
    Math.exp(Math.log(c1) + f * (Math.log(c2) - Math.log(c1))) :
    c1 + f * (c2 - c1);
}
