/**
 * Default NCA rules for UC-01 single-compound analysis, plus small
 * builder helpers for nested-partial overrides used by the (Phase 4)
 * controller.
 *
 * Every type here is imported
 * from `@datagrok-libraries/sci-comp`. The adapter never redefines the rules
 * shape, only supplies values.
 */

import {nca} from '@datagrok-libraries/sci-comp';

/**
 * BLQ-handling strategy used by UC-01 default rules.
 *
 * `set-zero` for pre-first / embedded / after-last BLQ points; `exclude` for
 * consecutive BLQ points after the last measurable observation.
 *
 * **This is the ICH M13A house choice, NOT PKNCA's defaults** (the header used
 * to claim the latter, which is false for `embedded`): PKNCA's own `conc.blq`
 * defaults are `first = "keep", middle = "drop", last = "keep"` — it DROPS an
 * embedded BLQ, we integrate it as 0, which is what M13A §2.2.2.2 prescribes
 * ("all BLQ → 0 for AUC"). The PKNCA option that corresponds to this strategy
 * is `"keep"` on all three positions, and that is what the validation oracles
 * now pass explicitly.
 *
 * The choice only started deciding numbers in sci-comp 0.11.0 — before it, the
 * integrator read the BLQ mask and every rule behaved as `exclude` (DATA-GAP-42).
 * The trailing `exclude` still costs nothing: the engine trims the trailing run
 * of non-positive values regardless (BUG-05), so `set-zero` and `exclude`
 * converge there — asserted upstream by the `06_blq_rules` R-E ≡ R-A block.
 */
export const DEFAULT_BLQ_STRATEGY: nca.BlqStrategy = {
  preFirstMeasurable: 'set-zero',
  embedded: 'set-zero',
  afterLast: 'set-zero',
  consecutiveAfterLast: 'exclude',
};

/**
 * Default NCA rules baked into UC-01 (spec §3.3). No UI exposes these for
 * editing in this slice; per-profile overrides are introduced with the
 * interactive λz override (later UC).
 *
 * Note on `adjRSquaredFactor: 1e-4` — without it the auto-best-fit selector
 * diverges from PKNCA on Theoph subject 6 (interface Addendum A2). Do not
 * change this value without re-running the reference suite (UC-01 spec §8).
 */
export const DEFAULT_RULES: nca.NcaRules = {
  aucMethod: 'linear-up-log-down',
  compensatedSummation: false,
  blq: DEFAULT_BLQ_STRATEGY,
  lambdaZ: {
    mode: 'auto-best-fit',
    minPoints: 3,
    minRSquared: 0.85,
    excludeCmax: true,
    adjRSquaredFactor: 1e-4,
  },
  extrapWarnPct: 20,
  extrapErrorPct: 50,
  // %AUMCextrap soft warning (FR-200 moment params). House heuristic — no
  // regulatory cutoff exists; mirrors extrapWarnPct (sci-comp 0.7.0). REGRESSION_PIN.
  extrapWarnPctAumc: 20,
};

/**
 * Default λz adj-R² ACCEPTANCE threshold (`R_acc`) — the value compared against a
 * FINISHED λz fit during profile acceptance (`adjR² < R_acc` → Flagged).
 *
 * **nca-studio SAP / presentation layer, NOT a sci-comp compute parameter (thesis #1).**
 * It is deliberately absent from `nca.NcaRules`: sci-comp never reads it, and the value
 * the user enters is applied VERBATIM — no derived recomputation (the removed `goodR`
 * formula). This decouples the two roles the single `λz min adj-R²` field used to conflate
 * (lambdaz-adjr2-threshold-spec-v2 §3): `minRSquared` stays the sci-comp point-selection
 * gate, `R_acc` is the acceptance criterion.
 *
 * `0.90` = PKNCA's `min.hl.r.squared` default and the spec's strict-context recommendation
 * (§3.3). It is `≥` the default core gate `minRSquared = 0.85`, so the consistency
 * constraint `R_acc ≥ minR` (§3.2) holds at the defaults with no change to the gate.
 */
export const DEFAULT_ACCEPTANCE_RSQUARED = 0.90;
/** Allowed range for the acceptance-threshold input (spec §3.3). The lower bound
 *  reflects the industry spread (SAPs use 0.70–0.90+); the effective floor on any given
 *  study is `max(ACCEPTANCE_RSQUARED_MIN, minRSquared)` via the `R_acc ≥ minR` guard. */
export const ACCEPTANCE_RSQUARED_MIN = 0.70;
export const ACCEPTANCE_RSQUARED_MAX = 0.95;

/**
 * Return a new `NcaRules` with `blq` shallow-merged from the supplied partial.
 * Other top-level fields are preserved by reference; `lambdaZ` is untouched.
 *
 * Used by the (Phase 4) controller when a per-profile BLQ strategy is needed
 * without disturbing the rest of the rules object.
 */
export function withBlq(base: nca.NcaRules, partial: Partial<nca.BlqStrategy>): nca.NcaRules {
  return {
    ...base,
    blq: {...base.blq, ...partial},
  };
}

/**
 * Return a new `NcaRules` with `lambdaZ` shallow-merged from the supplied
 * partial. Other top-level fields are preserved by reference; `blq` is
 * untouched.
 *
 * Used by the (Phase 4) controller for per-profile λz overrides — in
 * particular for the interactive λz override (`mode: 'manual-points'` +
 * `manualPoints`) introduced after UC-01.
 */
export function withLambdaZ(
  base: nca.NcaRules, partial: Partial<nca.LambdaZStrategy>,
): nca.NcaRules {
  return {
    ...base,
    lambdaZ: {...base.lambdaZ, ...partial},
  };
}
