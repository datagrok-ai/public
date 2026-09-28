/**
 * `augmentProfile` — the ONE place that turns a profile's raw inputs into the
 * array space `computeNca` integrates and fits λz on (pipeline Steps 1–3):
 *
 * 1. BLQ pre-processing (`applyBlqStrategy`) → substituted concentrations and
 *    the effective BLQ mask (input mask ∪ `exclude`d points).
 * 2. Observed Cmax/Tmax on the un-augmented post-BLQ profile.
 * 3. Route-aware dose-time augmentation:
 *    - **IV bolus** — a positive measured dose-time sample is used as-is
 *      (`method: 'observed'`). Otherwise the row is ABSENT: with no dose-time
 *      row, `(0, c0)` is prepended; with a dose-time row that is BLQ-flagged or
 *      non-positive/non-finite, that row is REMOVED and `(0, c0)` takes its
 *      slot (`replacedDoseTimeRow`). c0 comes from `insertC0` (the PKNCA
 *      `c0 → logslope → c1 → cmin → set0` chain on the post-dose points). A
 *      pre-dose sample can never supply the back-extrapolation anchor, whatever
 *      the substitution rule — the assay fact governs c0, and PKNCA's own `c0`
 *      behaves the same way (measured 2026-09-22).
 *    - **IV infusion / extravascular** — a dose-time row that survives the
 *      drop set is kept (a pre-dose 0 is a valid observation); otherwise
 *      `(0, 0)` is prepended by convention.
 *
 * Exported because two callers need exactly this index space: `computeNca`
 * (Steps 4–7) and a consumer that maps the engine's `pointsUsed` / manual λz
 * indices back to its own rows (nca-studio's λz editor). A mirror of the
 * orchestration cannot track REPLACE positionally, so the kernel reports
 * `sourceIndex` explicitly. Only index 0 is ever a dose-time row: a second
 * `t = 0` row is an ordinary observation (same shape as the pre-kernel gate).
 *
 * The kernel is stateless and never throws; the ONLY numerical primitives it
 * calls are the existing `applyBlqStrategy`, `findCmax` and `insertC0`.
 */

import type {
  ProfileInputs, BlqStrategy, BlqProcessingResult, CmaxResult, C0Estimate,
} from './types';
import {ROUTE_IV_BOLUS} from './types';
import {applyBlqStrategy} from './blq';
import {findCmax} from './cmax';
import {insertC0} from './c0';

/** The augmented profile — the array space `computeNca` works in. */
export interface AugmentedProfile {
  /** Times of the augmented profile (sorted ascending; `time[0] === 0` after
   *  any augmentation). */
  readonly time: Float64Array;
  /** Post-BLQ (substituted) concentrations; the synthetic `(0, c0)` / `(0, 0)`
   *  sits in slot 0 when inserted. */
  readonly conc: Float64Array;
  /** EFFECTIVE BLQ mask — input `blqMask` ∪ `exclude`d points. What the
   *  OBSERVED quantities read: Cmax/Tmax skip it, Tlag treats it as 0. */
  readonly blqMask: Uint8Array;
  /** What integration (AUC/AUMC) and the λz regression SKIP: `exclude`d and
   *  non-finite (`missing`) points. Substituted (`set-zero` / `set-half-lloq`)
   *  values are NOT in it — they reach the integrator as the values the rule
   *  wrote. The λz regression additionally drops `conc <= 0` itself. */
  readonly dropMask: Uint8Array;
  /** Input-array index of augmented element `k`, or `-1` for a synthetic
   *  point. No augmentation → `[0..n-1]`; prepend → `[-1, 0..n-1]`; REPLACE →
   *  `[-1, 1..n-1]` (input row 0 absent). */
  readonly sourceIndex: Int32Array;
  /** λz window anchor in augmented space (0 when c0 was inserted and is the
   *  peak of the augmented profile). */
  readonly cmaxIdx: number;
  /** Step 2 — the OBSERVED peak on the un-augmented post-BLQ profile. This is
   *  what `computeNca` reports as Cmax/Tmax; an inserted c0 is never the
   *  reported peak. */
  readonly observedCmax: CmaxResult;
  /** Step 1 — BLQ pre-processing on the raw inputs (`excluded` indexes the
   *  INPUT arrays). */
  readonly blqApplied: BlqProcessingResult;
  /** IV bolus only — how the dose-time concentration was obtained; `null`
   *  for every other route. */
  readonly c0: C0Estimate | null;
}

/**
 * Build the augmented profile for one set of inputs. Returns `null` when the
 * profile has no measurable point at all (every concentration BLQ / excluded /
 * non-finite) — `computeNca` reports that as `status: 'failed'`.
 */
export function augmentProfile(
  inputs: ProfileInputs, blq: BlqStrategy,
): AugmentedProfile | null {
  const n = inputs.time.length;

  // Step 1: BLQ pre-processing on the raw concentrations.
  const blqApplied = applyBlqStrategy(inputs.conc, inputs.blqMask, inputs.lloq, 0, blq);
  const procConc = blqApplied.conc;
  const effBlq = new Uint8Array(inputs.blqMask);
  for (let k = 0; k < blqApplied.excluded.length; k++)
    effBlq[blqApplied.excluded[k]] = 1;
  const dropMask = buildDropMask(procConc, blqApplied);

  // Step 2: observed Cmax/Tmax on the un-augmented post-BLQ profile.
  const observedCmax = findCmax(inputs.time, procConc, effBlq);
  if (observedCmax === null) return null;

  const identity = new Int32Array(n);
  for (let i = 0; i < n; i++) identity[i] = i;
  const unaugmented: AugmentedProfile = {
    time: inputs.time, conc: procConc, blqMask: effBlq, dropMask,
    sourceIndex: identity, cmaxIdx: observedCmax.cmaxIdx,
    observedCmax, blqApplied, c0: null,
  };
  const hasDoseTimeRow = n > 0 && inputs.time[0] === 0;

  // Step 3: route-aware dose-time augmentation.
  if (inputs.route === ROUTE_IV_BOLUS) {
    // A positive, unflagged, finite measured dose-time value is the c0.
    if (hasDoseTimeRow && effBlq[0] === 0 && procConc[0] > 0 && Number.isFinite(procConc[0])) {
      return {
        ...unaugmented,
        c0: {value: procConc[0], method: 'observed', replacedDoseTimeRow: false},
      };
    }
    // Otherwise the dose-time row is ABSENT (missing, or present-but-BLQ /
    // non-positive) — strip it if present and estimate c0 on the post-dose
    // points. The gate reads the EFFECTIVE BLQ mask + the post-BLQ value, never
    // `dropMask`, so a `set-zero` BLQ at t=0 is "absent" here even though the
    // same point is a kept 0 for extravascular integration.
    const s = hasDoseTimeRow ? 1 : 0;
    const aug = insertC0(
      inputs.time.subarray(s), procConc.subarray(s), effBlq.subarray(s));
    if (aug === null) {
      // Defensive — unreachable from computeNca: Step 2 found a measurable
      // point, and with the default chain `estimateC0` ends in `set0`, so
      // `insertC0` can only fail on a caller-supplied chain. Preserves the
      // pre-kernel behaviour (profile left un-augmented).
      return unaugmented;
    }
    const sourceIndex = new Int32Array(aug.time.length);
    sourceIndex[0] = -1;
    for (let k = 1; k < sourceIndex.length; k++) sourceIndex[k] = k - 1 + s;
    return {
      time: aug.time, conc: aug.conc, blqMask: aug.blqMask,
      dropMask: prependByte(dropMask.subarray(s), 0),
      sourceIndex, cmaxIdx: aug.cmaxIdx, observedCmax, blqApplied,
      c0: {value: aug.c0, method: aug.method, replacedDoseTimeRow: s === 1},
    };
  }

  // IV infusion / extravascular: an existing dose-time row that is NOT dropped
  // counts even when conc(0) = 0 — a valid pre-dose observation, not a missing
  // value. Otherwise the pre-dose concentration is 0 by convention.
  const hasT0 = hasDoseTimeRow && dropMask[0] === 0 && Number.isFinite(procConc[0]);
  if (hasT0) return unaugmented;
  const sourceIndex = new Int32Array(n + 1);
  sourceIndex[0] = -1;
  for (let k = 1; k <= n; k++) sourceIndex[k] = k - 1;
  return {
    time: prependScalar(inputs.time, 0),
    conc: prependScalar(procConc, 0),
    blqMask: prependByte(effBlq, 0),
    dropMask: prependByte(dropMask, 0),
    sourceIndex,
    cmaxIdx: observedCmax.cmaxIdx + 1,
    observedCmax, blqApplied, c0: null,
  };
}

/**
 * The drop set — what integration and λz skip: `exclude`d points and
 * non-finite (`missing`) values. Deliberately a separate mask from the
 * effective BLQ mask: the observed quantities (Cmax, Tlag, BLQ fraction) keep
 * reading the BLQ mask, while AUC/AUMC/λz must honour the values the
 * substitution rules wrote — `set-zero` integrates as 0 (and the λz
 * regression's own `conc <= 0` guard keeps it out of the fit), `set-half-lloq`
 * integrates as LLOQ/2 and IS λz-eligible. Before GROK-20960 the integrator and
 * the λz filter read the effective BLQ mask, so every rule was numerically
 * `exclude` and a declared `set-zero` was a no-op (DATA-GAP-42).
 *
 * `missing` ≡ `exclude` for every parameter: both land here; they differ only
 * in whether the concentration was overwritten with NaN.
 */
function buildDropMask(
  conc: Float64Array, blqApplied: BlqProcessingResult,
): Uint8Array {
  const dropMask = new Uint8Array(conc.length);
  for (let k = 0; k < blqApplied.excluded.length; k++)
    dropMask[blqApplied.excluded[k]] = 1;
  for (let i = 0; i < conc.length; i++)
    if (!Number.isFinite(conc[i])) dropMask[i] = 1;
  return dropMask;
}

function prependScalar(src: Float64Array, value: number): Float64Array {
  const out = new Float64Array(src.length + 1);
  out[0] = value;
  out.set(src, 1);
  return out;
}

function prependByte(src: Uint8Array, value: number): Uint8Array {
  const out = new Uint8Array(src.length + 1);
  out[0] = value;
  out.set(src, 1);
  return out;
}
