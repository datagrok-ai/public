/**
 * Recipe and review-decision fingerprints (nca-workflow-stages F-S1.4, research P3/P6; D13, D17).
 *
 * - {@link computeRecipeFingerprint} answers "has Setup changed since the last Confirm?" (P6). It
 *   hashes every Setup input that can change a fit, and nothing else.
 * - {@link computeDecisionFingerprint} answers "is what the analyst saw when deciding still what
 *   is on screen?" (D13). A mismatch resets the decision to `To review`.
 *
 * Both must be identical wherever they are computed — at Confirm, at the Confirm dry run (AD-7),
 * at save, and at reopen — so they live here, once.
 *
 * Reuses `fnv1a` (synchronous, deterministic; `study-schema.ts`). Pure: no `grok` / `ui`.
 */

import {fnv1a} from './study-schema';
import type {AnalysisRecipe} from './study-schema';
import {DEFAULT_ACCEPTANCE_RSQUARED} from './defaults';

/** Deterministic JSON: object keys sorted at every level; non-finite numbers written as tokens
 *  (`JSON.stringify` would turn NaN / ±Infinity into `null`, conflating them with a real null). */
function stableStringify(v: unknown): string {
  if (v === null || v === undefined) return 'null';
  if (typeof v === 'number') return Number.isFinite(v) ? String(v) : `"#${String(v)}"`;
  if (typeof v === 'string') return JSON.stringify(v);
  if (typeof v === 'boolean') return v ? 'true' : 'false';
  if (Array.isArray(v)) return `[${v.map(stableStringify).join(',')}]`;
  if (typeof v === 'object') {
    const o = v as Record<string, unknown>;
    const keys = Object.keys(o).filter((k) => o[k] !== undefined).sort();
    return `{${keys.map((k) => `${JSON.stringify(k)}:${stableStringify(o[k])}`).join(',')}}`;
  }
  return JSON.stringify(String(v));
}

/**
 * Fingerprint of the Setup inputs of a recipe: `columns` (mapping + context + sex map), `modes`,
 * `ncaRules`, `units` and `acceptanceRSquared` (absent → the shipped default, so a recipe that
 * never set it and one that set the default fingerprint alike).
 *
 * **Not hashed:** `groupByCols` and `passthroughCols` — Report / add-to-study inputs that change
 * no fit. Grouping on Report must never lock Review.
 *
 * A toggle-and-untoggle yields the original fingerprint (AC-S1.4.2): the hash reads values, not
 * an edit history.
 */
export function computeRecipeFingerprint(recipe: AnalysisRecipe): string {
  return fnv1a(stableStringify({
    columns: recipe.columns,
    modes: recipe.modes,
    ncaRules: recipe.ncaRules,
    units: recipe.units,
    acceptanceRSquared: recipe.acceptanceRSquared ?? DEFAULT_ACCEPTANCE_RSQUARED,
  }));
}

/** Version of {@link computeRecipeHash}'s input projection (`recipeOf`) and canonical form. Bump it
 *  whenever either changes: a stored hash of another version cannot be compared, and the pipeline
 *  then treats the entry as touched (fail-safe — it saves beside, never over). */
export const RECIPE_HASH_VERSION = 1;

/**
 * Hash of the WHOLE recipe — every field, `groupByCols` and `passthroughCols` included — for the
 * pipeline's "is this analysis still exactly what I wrote?" test (nca-pipeline-save-rework F4).
 * Unlike {@link computeRecipeFingerprint} (Setup inputs only), a Report-stage Group-by or an
 * annotation pick changes it: either is the scientist's work, and the pipeline must not overwrite
 * it. Pass `recipeOf(entry)` so a stored entry and a freshly built recipe hash the same projection.
 */
export function computeRecipeHash(recipe: AnalysisRecipe): string {
  return fnv1a(stableStringify(recipe));
}

/** A numeric fit summary at 6 significant digits (the F-S1.4 precision contract). */
function sig6(x: number): string {
  return Number.isFinite(x) ? x.toPrecision(6) : 'NA';
}

/** Everything a review decision was made on, for one review unit. Strings are the canonical
 *  `getString` forms (CSV round-trip safe) — never re-formatted numbers. */
export interface DecisionFingerprintInput {
  /** The unit's observations as canonical `[time, conc]` pairs, in source order. For a composite
   *  unit, the MEMBER observations (AC-S1.4.3 — a cell-reassignment changes the hash even when the
   *  cell id string does not). */
  readonly points: ReadonlyArray<readonly [string, string]>;
  /** The per-profile dose and route as displayed (canonical strings). */
  readonly dose: string;
  readonly route: string;
  /** The times (canonical strings) of the points the reported λz fit used. */
  readonly pointsUsedTimes: ReadonlyArray<string>;
  readonly lambdaZ: number;
  readonly adjR2: number;
  readonly aucInf: number;
  readonly pctExtrap: number;
  /** The acceptance verdict and its violated-criterion codes (stable ids, never display text). */
  readonly verdict: string;
  readonly flagCodes: ReadonlyArray<string>;
  /** The row's warnings as `code + severity` (never the message — R2/F7). */
  readonly warnings: ReadonlyArray<{readonly code: string; readonly severity: string}>;
}

/** Warning codes about the body-weight normalisation of CL / Vz / Vss only (nca-cl-bw-normalised-units
 *  AD-1 / AD-5b), left out of every decision fingerprint — dense and composite. */
export const BASIS_NOTE_CODES: ReadonlySet<string> = new Set([
  'PER_KG_NOT_FORMED', 'BW_VARIES_WITHIN_PROFILE', 'SPARSE_PER_KG_WEIGHT_BASED', 'COMPOSITE_BW_PARTIAL',
  'COMPOSITE_NOT_WEIGHT_BASED',
]);

/**
 * Fingerprint of what the analyst saw when deciding (D13; R1/A1, R1/A2, R2/F7).
 *
 * - **Copy-edit invariant** (AC-S1.4.1): warnings are hashed by code + severity and acceptance by
 *   flag code, so rewording a message or a reason never moves a fingerprint.
 * - **Rules are not hashed**: a rules change that leaves this profile's fit unchanged keeps the
 *   decision; one that changes the fit changes the fit summary below.
 * - **Precision**: fit numbers at 6 significant digits. A sub-6th-digit change never resets a
 *   decision (same contract as the P6 no-op); one that crosses a 6th-digit rounding boundary may.
 *   The failure direction is a spurious "please re-look", never a silently kept stale decision.
 * - Every hashed value is closed-form and non-stochastic, so the Confirm dry run (AD-7) predicts
 *   exactly. A bootstrap-derived field must never be added here without re-verifying AC-S3.5.1.
 * - **Reporting-basis notes are not hashed** ({@link BASIS_NOTE_CODES}): they say how CL / Vz /
 *   Vss were normalised by body weight, which no hashed value depends on — so the per-kg notes
 *   arriving with nca-cl-bw-normalised-units never reset a decision already taken on a fit.
 */
export function computeDecisionFingerprint(u: DecisionFingerprintInput): string {
  const warnings = u.warnings.filter((w) => !BASIS_NOTE_CODES.has(w.code))
    .map((w) => `${w.code}:${w.severity}`).sort();
  return fnv1a(stableStringify({
    points: u.points.map(([t, c]) => `${t}|${c}`),
    dose: u.dose,
    route: u.route,
    used: [...u.pointsUsedTimes],
    fit: [sig6(u.lambdaZ), sig6(u.adjR2), sig6(u.aucInf), sig6(u.pctExtrap)],
    verdict: u.verdict,
    flags: [...u.flagCodes],
    warnings,
  }));
}
