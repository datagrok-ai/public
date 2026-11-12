/**
 * NcaController — single-entry orchestrator for the UC-01 pipeline:
 *   validate → index → extract per profile → computeNca → audit → write results.
 *
 * @see libraries/sci-comp/src/nca (the `computeNca` contract)
 *
 * **Pull-based, single function.** UC-01 has no manual λz override and no
 * recompute trigger graph — pull is the right shape; the rich `NcaController`
 * class shape from interface §9.3 is deferred to the interactive-λz UC. The
 * function is reusable from the column-mapping dialog "Compute" button (Phase 5)
 * and from the Datagrok grok-function entry point (Phase 6).
 *
 * **Error policy (spec §6.1).**
 * - Validator errors (≥ 1 `ValidationError`) → the function returns
 *   `{status:'validation-failed', errors}` without computing anything. The UI
 *   renders the error panel; compute is blocked.
 * - Extractor / engine throws on one subject → that subject's row is marked
 *   `status:'failed'` with all values NaN, an `ENGINE_ERROR` warning carrying
 *   the thrown message, a `console.error` log, and an audit-log entry; the
 *   loop continues for the remaining subjects.
 * - Extractor warnings (e.g. `FLOAT32_WIDENED` per Phase 3) are merged into
 *   the corresponding `ComputeResult.provenance.warnings` so they surface in
 *   the per-profile `WARNINGS_JSON` column. Per [[NOTE-01]] / Phase 3 review:
 *   no dedup — each profile's provenance is self-contained for traceability.
 *
 * **Audit-log policy (DoD §9 row 14).** Every `runNcaPipeline` call writes one
 * entry per profile to *the package's* `AuditLog`, with `action:'compute'`, the
 * 5 PRD IR-12 fields (timestamp, subjectKey, ruleHash, status — plus action),
 * and an opaque `payload` carrying the full `ProfileProvenance` snapshot. The
 * controller defaults to the package-level singleton (`getAuditLog()` from
 * `audit-log.ts`) so DoD row 14 is satisfied unconditionally — callers cannot
 * silently skip audit. Tests override the singleton via `options.log` for
 * isolation.
 *
 * **rule-hash policy.** Stable JSON.stringify with explicit field ordering
 * (insensitive to object-literal key order); UC-01 has no rule UI so the only
 * realistic differences come from multiple runs against different presets in
 * one session. The hashing concern lives here (not in audit-log.ts) per
 * audit-log.ts's "hashing algorithm is owned by the (Phase 4) controller".
 */

import type * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {nca} from '@datagrok-libraries/sci-comp';

import type {ColumnContract, ValidationError} from '../adapter/column-contract';
import type {ProfileIndex} from '../adapter/profile-index';
import type {ProfileScalars} from '../adapter/profile-inputs';
import {validateObservations, validateUnits} from '../adapter/validators';
import {buildProfileIndices} from '../adapter/profile-index';
import {extractProfileInputs} from '../adapter/profile-inputs';
import {writeResultsToDataFrame, patchResultRow} from '../adapter/results-writer';
import type {AuditLog} from '../adapter/audit-log';
import {getAuditLog} from '../adapter/audit-log';
import {impliedClPerKgWarning} from '../adapter/cl-plausibility';
import {NO_BODY_WEIGHT, normaliseClVolumes} from '../adapter/unit-conversion';
import type {BodyWeightContext} from '../adapter/unit-conversion';
import {withLambdaZ} from '../adapter/defaults';
import {
  applyMinPointsFloor, buildAugmentedProfile, buildManualPoints, pointsUsedToRows,
  subjectHasManualOverride, isLambdaZNonContiguous, LAMBDAZ_NONCONTIGUOUS,
} from '../adapter/lambdaz-manual-points';
import type {LzAction} from '../adapter/lambdaz-manual-points';
import type {
  LzOverride, ReviewAction, ReviewLedger, ReviewReason, ReviewDecisionRecord,
} from '../adapter/study-schema';
import {evaluateAcceptance} from '../adapter/acceptance';
import type {AcceptanceThresholds} from '../adapter/acceptance';
import {numAt} from '../adapter/adapter-utils';
import {classifyReviewUnits, profileMeasurability} from '../adapter/review-units';
import type {CompositeUnitInput, ProfileUnitInput, ReviewUnit} from '../adapter/review-units';
import {computeDecisionFingerprint} from '../adapter/fingerprint';
import type {DecisionFingerprintInput} from '../adapter/fingerprint';
import {
  ledgerMapOf, ledgerListOf, recordDecision, supersedeDecision, reconcileLedgers,
} from '../adapter/review-ledger';

/**
 * Per-profile slice of the pipeline outcome — kept on the result so tests
 * and downstream views can inspect individual subjects without re-parsing
 * the results DataFrame.
 */
export interface PerProfileResult {
  readonly index: ProfileIndex;
  readonly route: nca.RouteCode;
  readonly result: nca.ComputeResult;
  /**
   * Per-profile ENGINE dose basis (nca-cl-vz-convention AD-2): `'absolute'`
   * (the engine ran on mg) or `'per_kg'` (`/kg` dose with no body weight — the
   * engine ran per kg). Threaded from `extractProfileInputs` into the writer's
   * `DOSE_NORM` column, exactly like `route`; the reported CL / Vz / Vss are per
   * kg either way (`normaliseClVolumes`).
   */
  readonly doseNormalization: 'absolute' | 'per_kg';
  /** The profile's body weight for the per-kg rule (nca-cl-bw-normalised-units AD-2) — threaded
   *  to the writer beside `doseNormalization`, so the full build and the λz patch form CL / Vz /
   *  Vss per kg from the same inputs. */
  readonly bodyWeight: BodyWeightContext;
  /**
   * Extractor warnings already merged into `result.provenance.warnings` —
   * exposed separately so test code can assert the merge contract without
   * re-filtering the merged list.
   */
  readonly extractorWarnings: ReadonlyArray<nca.ParameterWarning>;
  /**
   * Source-df row indices ACTUALLY used in this profile's λz fit (VR-11, F6) —
   * translated from `provenance.lambdaZ.pointsUsed` (augmented indices) back to
   * original rows via the {@link buildAugmentedProfile} `augRowId`. Empty when the
   * fit failed / was suppressed. This is the SINGLE source for the facet's used-point
   * markers AND the click-toggle decision, replacing the time-window proxy so a
   * marker cannot disagree with the reported fit (the window rule mislabels a point
   * that sits inside `[tFirst,tLast]` but is not in a gappy manual set).
   */
  readonly fitRows: ReadonlySet<number>;
  /**
   * Whether this profile's λz fit is NON-CONTIGUOUS (decision 2) — a measurable
   * point sits inside the fitted window but is excluded (only possible from a manual
   * selection; auto-best-fit is contiguous by construction). Drives the
   * `LAMBDAZ_NONCONTIGUOUS` warning + the facet caption note.
   */
  readonly nonContiguous: boolean;
}

/** Shared empty membership set — no per-call allocation for the no-fit path. */
const EMPTY_ROW_SET: ReadonlySet<number> = new Set<number>();

/** The non-contiguous-λz advisory (decision 2), prepended to a profile's warnings
 *  when its manual selection skips an interior measurable point. */
const NONCONTIGUOUS_WARNING: nca.ParameterWarning = {
  code: LAMBDAZ_NONCONTIGUOUS,
  severity: 'warning',
  message:
    'λz selection is non-contiguous — a measurable point inside the fitted terminal ' +
    'window is excluded. By convention λz is fit to consecutive terminal points; ' +
    'confirm the excluded point is a justified outlier.',
};

/**
 * Outcome discriminator. Tests, views, and the eventual grok-function entry
 * point pattern-match on `status` and consume the matching variant.
 */
export type NcaRunOutcome =
  | {readonly status: 'validation-failed'; readonly errors: ReadonlyArray<ValidationError>}
  | {
      readonly status: 'completed';
      readonly results: DG.DataFrame;
      /** Read-only to callers; `NcaController.recomputeProfile` replaces an entry IN PLACE so every
       *  holder of this array (the view's saved-summary input) sees the post-override result. */
      readonly perProfile: ReadonlyArray<PerProfileResult>;
    };

/**
 * Optional knobs.
 *
 * `log` overrides the package-level singleton (tests only — production code
 * must never supply this so that DoD §9 row 14 is enforced uniformly via the
 * singleton). `now` is injected for deterministic timestamp assertions in
 * tests; defaults to `Date.now`.
 */
export interface RunNcaPipelineOptions {
  readonly log?: AuditLog;
  readonly now?: () => number;
  /**
   * Per-subject route-normalization warnings (per-row-route AD-1/3/5), keyed by
   * stringified subject id. Produced by `normalizeRouteColumn` in the mapping
   * pane (Option A: encoding happens pane-side, extractor route read unchanged)
   * and merged into each profile's `provenance.warnings` here so they surface in
   * `WARNINGS_JSON` — without coupling the extractor to route strings.
   */
  readonly routeWarningsBySubject?: ReadonlyMap<string, ReadonlyArray<nca.ParameterWarning>>;
  /**
   * Interactive λz manual-override hook (Chunk B). For each profile, returns the
   * `manualPoints` (augmented indices) to force the λz fit onto, or `null` to
   * leave that subject on the base rules' `auto-best-fit`. Supplied by
   * {@link NcaController} (which resolves its keyed override store); absent for
   * the plain UC-01 path, so a frame with no overrides is byte-identical to before
   * (AC-3.2). When it returns non-null, the strict ≥`minPoints` floor (AC-4.3) is
   * applied to the resulting fit.
   */
  readonly lambdaZResolver?: (index: ProfileIndex, inputs: nca.ProfileInputs) => Int32Array | null;
}

/**
 * Run the full UC-01 pipeline on one observations DataFrame.
 *
 * Pure function with respect to its inputs: no globals read, no mutations of
 * the input DataFrame or contract. Side effects are limited to the supplied
 * `AuditLog` (one append per profile) and `console.error` on engine throws.
 */
export function runNcaPipeline(
  df: DG.DataFrame,
  contract: ColumnContract,
  scalars: ProfileScalars,
  rules: nca.NcaRules,
  options: RunNcaPipelineOptions = {},
): NcaRunOutcome {
  const errors = validateObservations(df, contract);
  if (errors.length > 0)
    return {status: 'validation-failed', errors};

  const profiles = buildProfileIndices(df, contract);

  // Unit-load-bearing validators (nca-unit-semantics — E6 / AD-4 / AD-8 /
  // AD-9). Runs AFTER the column-level validators (so molar-unit errors
  // don't pile on top of a MISSING_COLUMN) and AFTER profile indices are
  // built (the LLOQ sanity validator iterates per-profile). Non-empty ⇒
  // validation-failed before any extract / compute.
  const unitErrors = validateUnits(df, contract, scalars, profiles);
  if (unitErrors.length > 0)
    return {status: 'validation-failed', errors: unitErrors};

  const ruleHash = computeRuleHash(rules);
  const now = options.now ?? Date.now;
  const log = options.log ?? getAuditLog();

  const perProfile: PerProfileResult[] = new Array(profiles.length);
  const routes: nca.RouteCode[] = new Array(profiles.length);
  const doseNorms: Array<'absolute' | 'per_kg'> = new Array(profiles.length);
  const bodyWeights: BodyWeightContext[] = new Array(profiles.length);
  const results: nca.ComputeResult[] = new Array(profiles.length);
  // C0 back-extrapolation result per profile (Ed obs #4). NaN-init: only
  // populated for IV-bolus profiles where `nca.estimateC0` returns a value.
  // Independent of `computeNca` success — c0 is a function of inputs only,
  // so it surfaces even when downstream computation reports `failed`.
  const c0s = new Float64Array(profiles.length).fill(NaN);

  for (let i = 0; i < profiles.length; i++) {
    const index = profiles[i];
    const single = computeSingleProfile(df, index, contract, scalars, rules, options);
    perProfile[i] = single.perProfile;
    routes[i] = single.perProfile.route;
    doseNorms[i] = single.perProfile.doseNormalization;
    bodyWeights[i] = single.perProfile.bodyWeight;
    results[i] = single.perProfile.result;
    c0s[i] = single.c0;

    log.append({
      timestamp: now(),
      subjectKey: index.subjectKey,
      action: 'compute',
      ruleHash,
      status: single.perProfile.result.status,
      payload: single.perProfile.result.provenance,
    });
  }

  const resultsDf = writeResultsToDataFrame(profiles, routes, results, c0s, doseNorms, bodyWeights);
  return {status: 'completed', results: resultsDf, perProfile};
}

/**
 * Compute ONE profile: extract → (optional λz manual-override resolution) →
 * `computeNca` → (strict floor if manual) → warning merge. The per-profile body
 * shared by {@link runNcaPipeline}'s loop and {@link NcaController.recomputeProfile}
 * so the auto and single-subject paths cannot drift. Does NOT append to the audit
 * log — each caller appends with its own `action` (`'compute'` vs `'override-lambdaz'`).
 */
function computeSingleProfile(
  df: DG.DataFrame,
  index: ProfileIndex,
  contract: ColumnContract,
  scalars: ProfileScalars,
  rules: nca.NcaRules,
  options: RunNcaPipelineOptions,
): {perProfile: PerProfileResult; c0: number} {
  const routeWarnings = options.routeWarningsBySubject?.get(index.subjectKey) ?? [];
  let extractorWarnings: ReadonlyArray<nca.ParameterWarning> = routeWarnings;
  let route: nca.RouteCode = scalars.route;
  let doseNormalization: 'absolute' | 'per_kg' = 'absolute';
  let bodyWeight: BodyWeightContext = NO_BODY_WEIGHT;
  let result: nca.ComputeResult;
  let c0 = NaN;
  // Real λz fit membership (F6) — the source rows the reported fit actually used.
  // Filled inside the try once the final result + its inputs are known; stays empty
  // for the failed/threw path.
  let fitRows: ReadonlySet<number> = EMPTY_ROW_SET;
  let nonContiguous = false;

  try {
    const extracted = extractProfileInputs(df, index, contract, scalars);
    extractorWarnings = [...routeWarnings, ...extracted.warnings];
    route = extracted.inputs.route;
    doseNormalization = extracted.doseNormalization;
    bodyWeight = extracted.bodyWeight;
    // estimateC0 — see the long note that used to live in the loop: IV-bolus
    // only, NaN means "not IV bolus / extraction threw", the default chain never
    // returns null.
    if (route === nca.ROUTE_IV_BOLUS) {
      const c0v = nca.estimateC0(
        extracted.inputs.time, extracted.inputs.conc, extracted.inputs.blqMask);
      if (c0v !== null) c0 = c0v;
    }
    // λz manual override (Chunk B): resolve this subject's forced points, if any.
    const manualPoints = options.lambdaZResolver?.(index, extracted.inputs) ?? null;
    const effectiveRules = manualPoints !== null ?
      withLambdaZ(rules, {mode: 'manual-points', manualPoints}) : rules;
    let raw = normalizeFailedValues(nca.computeNca(extracted.inputs, effectiveRules));
    // Strict ≥minPoints floor for a manual fit (AC-4.3) — auto fits are already
    // ≥minPoints, so this only bites a manual selection.
    if (manualPoints !== null)
      raw = applyMinPointsFloor(raw, effectiveRules.lambdaZ.minPoints);
    // The plausibility check reads the ENGINE-basis CL (architect B1), never the reported one.
    const clWarning = impliedClPerKgWarning(
      raw.values.cl, extracted.inputs.bodyWeight, doseNormalization);
    // The per-kg rule's own notes (PER_KG_NOT_FORMED / BW_VARIES_WITHIN_PROFILE) join the
    // profile's warning stream here; the writer forms the values from the same inputs.
    const perKgNotes = normaliseClVolumes(raw.values, doseNormalization, bodyWeight).warnings;
    const postExtras = [...extractorWarnings, ...(clWarning === null ? [] : [clWarning]), ...perKgNotes];
    result = postExtras.length === 0 ? raw : mergeWarnings(raw, postExtras);
    // Translate the reported λz points (augmented indices) back to source rows via
    // the SAME augmented profile the fit was built on (BLQ unchanged by an override,
    // so `rules.blq` == `effectiveRules.blq`). This is the honest membership marker
    // source — never the time-window range.
    if (result.provenance.lambdaZ !== null) {
      const aug = buildAugmentedProfile(extracted.inputs, index.observationRows, rules.blq);
      fitRows = new Set<number>(pointsUsedToRows(result.provenance.lambdaZ.pointsUsed, aug.rowId));
      // Decision 2: flag a non-contiguous manual selection (auto is always contiguous,
      // so this only fires on a hand edit that skips an interior measurable point).
      nonContiguous = isLambdaZNonContiguous(result.provenance.lambdaZ.pointsUsed, aug);
      if (nonContiguous) result = mergeWarnings(result, [NONCONTIGUOUS_WARNING]);
    }
  } catch (e) {
    const message = e instanceof Error ? e.message : String(e);
    console.error(`runNcaPipeline: subject '${index.subjectKey}' failed: ${message}`);
    const engineWarning: nca.ParameterWarning = {code: 'ENGINE_ERROR', severity: 'error', message};
    result = makeFailedResult(rules, [engineWarning, ...extractorWarnings]);
  }

  return {
    perProfile: {index, route, result, doseNormalization, bodyWeight, extractorWarnings, fitRows, nonContiguous},
    c0,
  };
}

/**
 * Enforce spec §6.1 row 3 at the controller boundary: when the engine returns
 * `status:'failed'`, the subject row must carry NaN in every parameter slot
 * regardless of what the engine put in `values`. Sci-comp's `computeNca`
 * contract documents this case as "Cmax not available" — in practice the
 * engine already sets every value to NaN — but the spec places the invariant
 * on the controller (`Subject row written with all values = NaN`), so we
 * enforce it here in case sci-comp ever drifts. Provenance (warnings,
 * blqApplied, lambdaZ) passes through untouched.
 */
function normalizeFailedValues(result: nca.ComputeResult): nca.ComputeResult {
  if (result.status !== 'failed') return result;
  const v = result.values;
  const allNaN =
    Number.isNaN(v.cmax) && Number.isNaN(v.tmax) &&
    Number.isNaN(v.aucLast) && Number.isNaN(v.aucInf) &&
    Number.isNaN(v.pctExtrap) && Number.isNaN(v.lambdaZ) &&
    Number.isNaN(v.halfLife) && Number.isNaN(v.cl) && Number.isNaN(v.vz) &&
    // FR-200 moment / lag params (sci-comp 0.7.0) — drift-guard must cover
    // every column or the §6.1 "failed → all NaN" invariant silently leaks.
    Number.isNaN(v.aumcLast) && Number.isNaN(v.aumcInf) &&
    Number.isNaN(v.mrt) && Number.isNaN(v.vss) &&
    Number.isNaN(v.tlag) && Number.isNaN(v.pctExtrapAumc);
  if (allNaN) return result;
  return {
    ...result,
    values: {
      cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
      lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN,
      aumcLast: NaN, aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN,
      pctExtrapAumc: NaN,
    },
  };
}

/**
 * Prepend extractor warnings onto `provenance.warnings`. Fresh provenance and
 * fresh warnings array — the input `ComputeResult` is never mutated.
 */
function mergeWarnings(
  result: nca.ComputeResult,
  extras: ReadonlyArray<nca.ParameterWarning>,
): nca.ComputeResult {
  const warnings: nca.ParameterWarning[] = [...extras, ...result.provenance.warnings];
  return {
    ...result,
    provenance: {...result.provenance, warnings},
  };
}

/**
 * Synthesize a `ComputeResult` for a subject whose pipeline threw. All values
 * NaN, `lambdaZ` provenance null, empty BLQ-applied buffers, status `'failed'`.
 * The supplied warnings (engine-error + extractor warnings) are attached so
 * WARNINGS_JSON tells the auditor what happened.
 */
function makeFailedResult(
  rules: nca.NcaRules,
  warnings: ReadonlyArray<nca.ParameterWarning>,
): nca.ComputeResult {
  const values: nca.ParameterValues = {
    cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
    lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN,
    aumcLast: NaN, aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN,
    pctExtrapAumc: NaN,
  };
  const provenance: nca.ProfileProvenance = {
    lambdaZ: null,
    blqApplied: {conc: new Float64Array(0), excluded: new Int32Array(0)},
    aucMethod: rules.aucMethod,
    compensated: rules.compensatedSummation,
    warnings,
  };
  return {status: 'failed', values, provenance};
}

/**
 * Stable string fingerprint of `NcaRules`. Serializes fields in a fixed order
 * so the hash is invariant under object-literal key reordering. Treats
 * `adjRSquaredFactor: undefined` as `null` for canonicalization, and serializes
 * `manualPoints` (Int32Array) as a plain number array.
 *
 * This is a fingerprint, not a cryptographic hash — Phase 3 P2 (per AD-04 /
 * 21 CFR Part 11 retrofit) replaces it with a record-hash chain when the
 * GxP-track trigger fires.
 */
function computeRuleHash(rules: nca.NcaRules): string {
  // Explicit field ordering — JSON.stringify respects insertion order, which
  // differs between {aucMethod, blq, lambdaZ} and {lambdaZ, aucMethod, blq}
  // depending on how the rules object was constructed. Pinning order here
  // makes the fingerprint invariant to caller-side key ordering.
  const blq = rules.blq;
  const lz = rules.lambdaZ;
  return JSON.stringify({
    aucMethod: rules.aucMethod,
    compensatedSummation: rules.compensatedSummation,
    extrapWarnPct: rules.extrapWarnPct,
    extrapErrorPct: rules.extrapErrorPct,
    // FR-200: a run differing only in the %AUMCextrap knob must produce a
    // distinct rule-hash (AC-A5) — else audit fingerprints collide.
    extrapWarnPctAumc: rules.extrapWarnPctAumc,
    blq: {
      preFirstMeasurable: blq.preFirstMeasurable,
      embedded: blq.embedded,
      afterLast: blq.afterLast,
      consecutiveAfterLast: blq.consecutiveAfterLast,
    },
    lambdaZ: {
      mode: lz.mode,
      minPoints: lz.minPoints,
      minRSquared: lz.minRSquared,
      excludeCmax: lz.excludeCmax,
      adjRSquaredFactor: lz.adjRSquaredFactor ?? null,
      manualPoints: lz.manualPoints !== undefined ? Array.from(lz.manualPoints) : null,
      manualRange: lz.manualRange ?? null,
    },
  });
}

// === Interactive λz manual override — NcaController (interface §9.5) ============

/** One per-point λz-membership change requested by the UI (F4). */
export interface OverrideAction {
  /** Source-df row index of the point. */
  readonly rowIdx: number;
  /** `force-in` / `force-out` to override, `auto` to reset that point (F-4). */
  readonly action: LzAction;
}

/**
 * Stateful compute orchestrator for the interactive λz feature (interface §9.5).
 * Holds `{df, contract, scalars, rules}`; wraps {@link runNcaPipeline} for the
 * whole-cohort path ({@link computeAll}) and {@link computeSingleProfile} +
 * {@link patchResultRow} for the fast single-subject path
 * ({@link recomputeProfile}) the toggle uses.
 *
 * **Override state = one keyed store** (nca-workflow-stages AD-2, superseding
 * manual-point-selection D5 "stored on the source frame"). Each override is keyed by
 * `(profileKey, timeKey, ordinal)` — the results-row `SUBJECT`, the observation time in
 * canonical `getString` form, and its index among that profile's same-time rows — so its
 * identity is independent of any frame's row index. That is what lets an override survive
 * a full Compute and a Save → reopen in EVERY mode (D22): in NONMEM / composite / sex-remap
 * modes the frame this controller computes on is rebuilt every Compute, and a row index
 * does not survive that. The source frame carries no `_lz_*` column. The store is the ONE
 * representation; every reader (resolver, facet markers, ✎ flag, points grid) goes through
 * the methods below.
 *
 * Lifecycle (AC-1.5): {@link setInputs} refreshes the held inputs on each Compute
 * (the pane rebuilds contract/scalars/rules from live UI); `recomputeProfile`
 * between Computes reuses the last-held inputs and never re-reads the UI.
 */
export class NcaController {
  private readonly log: AuditLog;
  private readonly nowFn: () => number;
  private outcome: NcaRunOutcome | null = null;
  private readonly profileByKey = new Map<string, ProfileIndex>();
  /** Per-subject real λz fit membership (source rows), the facet marker + toggle
   *  source. Populated whole-cohort by {@link computeAll}, patched one-subject by
   *  {@link recomputeProfile} so a toggle keeps every other subject's set intact.
   *  `recomputeProfile` patches `outcome.results` and replaces the subject's entry of
   *  `outcome.perProfile` in place; this cache is the keyed read of the same membership. */
  private readonly fitRowsCache = new Map<string, ReadonlySet<number>>();
  /** Subjects whose λz selection is non-contiguous (decision 2) — the facet caption
   *  note source. Maintained alongside {@link fitRowsCache}. */
  private readonly nonContiguousSet = new Set<string>();
  private routeWarnings: ReadonlyMap<string, ReadonlyArray<nca.ParameterWarning>> | undefined;
  /** The keyed λz override store (AD-2): profileKey → its overrides. */
  private readonly lzStore = new Map<string, LzOverride[]>();
  /** Per-row key cache for the CURRENT frame: row → `{profileKey, pointKey}`. Rebuilt lazily
   *  whenever the frame / profiles change ({@link invalidateRowKeys}). */
  private rowKeys: Map<number, {profileKey: string; pointKey: string}> | null = null;
  /** Overrides dropped by the last {@link computeAll} because their point no longer exists
   *  or its concentration changed (AC-S2.1.3). Read-and-cleared by {@link consumeReconciled}. */
  private reconciled = 0;
  /** Review-decision ledgers by unit key (F-S2.2). */
  private ledgers = new Map<string, ReviewLedger>();

  constructor(
    private df: DG.DataFrame,
    private contract: ColumnContract,
    private scalars: ProfileScalars,
    private rules: nca.NcaRules,
    options: {
      log?: AuditLog; now?: () => number;
      routeWarnings?: ReadonlyMap<string, ReadonlyArray<nca.ParameterWarning>>;
    } = {},
  ) {
    this.log = options.log ?? getAuditLog();
    this.nowFn = options.now ?? Date.now;
    this.routeWarnings = options.routeWarnings;
  }

  /** Refresh the held inputs from live UI state (AC-1.5, called on each Compute). The
   *  override store is untouched — its keys are frame-independent (AD-2). */
  setInputs(
    df: DG.DataFrame, contract: ColumnContract, scalars: ProfileScalars, rules: nca.NcaRules,
    routeWarnings?: ReadonlyMap<string, ReadonlyArray<nca.ParameterWarning>>,
  ): void {
    this.df = df; this.contract = contract; this.scalars = scalars; this.rules = rules;
    this.routeWarnings = routeWarnings;
    this.rowKeys = null;
  }

  /** The frame this controller computes on (the source frame, or the derived NONMEM /
   *  composite / sex-mapped frame) — the frame whose row indices `OverrideAction.rowIdx`,
   *  {@link currentLzAction} and {@link fitRowsBySubject} refer to. */
  get frame(): DG.DataFrame {return this.df;}
  get heldContract(): ColumnContract {return this.contract;}
  get heldRules(): nca.NcaRules {return this.rules;}

  // ── Keyed override store (AD-2) ──────────────────────────────────────────

  /** Every override in the store, for persistence on the analysis entry (v15
   *  `lzOverrides`). Ordered by profile then insertion. */
  getLzOverrides(): LzOverride[] {
    return [...this.lzStore.values()].flat();
  }

  /** Replace the store (reopen / Revert). The next {@link computeAll} applies them and
   *  drops any that no longer resolve (AC-S2.1.3). */
  setLzOverrides(list: ReadonlyArray<LzOverride>): void {
    this.lzStore.clear();
    for (const o of list) {
      let arr = this.lzStore.get(o.profileKey);
      if (arr === undefined) {arr = []; this.lzStore.set(o.profileKey, arr);}
      arr.push(o);
    }
  }

  /** Overrides dropped by the last compute because their point vanished or its
   *  concentration changed; returns the count and clears it (AC-S2.1.3). */
  consumeReconciled(): number {
    const n = this.reconciled;
    this.reconciled = 0;
    return n;
  }

  /** True when `profileKey` has at least one override in the store. */
  hasOverrides(profileKey: string): boolean {
    return (this.lzStore.get(profileKey)?.length ?? 0) > 0;
  }

  private invalidateRowKeys(): void {this.rowKeys = null;}

  /** Build (once per frame) the row → `(profileKey, timeKey#ordinal)` map from the held
   *  profiles. `ordinal` counts earlier rows of the same profile at the same time key. */
  private rowKeyMap(): Map<number, {profileKey: string; pointKey: string}> {
    if (this.rowKeys !== null) return this.rowKeys;
    const m = new Map<number, {profileKey: string; pointKey: string}>();
    const timeCol = this.df.col(this.contract.observations.time.name);
    for (const [profileKey, index] of this.profileByKey) {
      const seen = new Map<string, number>();
      const rows = index.observationRows;
      for (let i = 0; i < rows.length; i++) {
        const r = rows[i];
        const t = timeCol !== null && !timeCol.isNone(r) ? timeCol.getString(r) : '';
        const ord = seen.get(t) ?? 0;
        seen.set(t, ord + 1);
        m.set(r, {profileKey, pointKey: `${t}#${ord}`});
      }
    }
    this.rowKeys = m;
    return m;
  }

  /** Canonical concentration string of a row (the CSV-round-trip-safe `getString` form). */
  private concKeyAt(row: number): string {
    const c = this.df.col(this.contract.observations.concentration.name);
    return c === null || c.isNone(row) ? '' : c.getString(row);
  }

  /** Row → action for one profile, from the store. Only overrides whose point exists on
   *  the current frame AND whose concentration still matches apply; the rest are ignored
   *  here and pruned after the compute (AC-S2.1.3). */
  private actionsForIndex(index: ProfileIndex): Map<number, LzAction> {
    const out = new Map<number, LzAction>();
    const list = this.lzStore.get(index.subjectKey);
    if (list === undefined || list.length === 0) return out;
    const byPoint = new Map<string, LzOverride>();
    for (const o of list) byPoint.set(`${o.timeKey}#${o.ordinal}`, o);
    const timeCol = this.df.col(this.contract.observations.time.name);
    const seen = new Map<string, number>();
    const rows = index.observationRows;
    for (let i = 0; i < rows.length; i++) {
      const r = rows[i];
      const t = timeCol !== null && !timeCol.isNone(r) ? timeCol.getString(r) : '';
      const ord = seen.get(t) ?? 0;
      seen.set(t, ord + 1);
      const o = byPoint.get(`${t}#${ord}`);
      if (o !== undefined && o.concKey === this.concKeyAt(r)) out.set(r, o.action);
    }
    return out;
  }

  /** The stored override for one row, or null (not overridden / no longer matching). */
  private overrideAtRow(rowIdx: number): LzOverride | null {
    const k = this.rowKeyMap().get(rowIdx);
    if (k === undefined) return null;
    const [timeKey, ordStr] = splitPointKey(k.pointKey);
    const o = this.lzStore.get(k.profileKey)?.find(
      (x) => x.timeKey === timeKey && x.ordinal === ordStr);
    return o !== undefined && o.concKey === this.concKeyAt(rowIdx) ? o : null;
  }

  /** Drop every stored override that does not resolve on the current frame (its profile,
   *  time/ordinal point, or concentration changed). Returns how many were dropped. */
  private pruneUnresolved(): number {
    let dropped = 0;
    const keys = this.rowKeyMap();
    const live = new Set<string>();
    for (const [row, k] of keys) live.add(`${k.profileKey}\u0001${k.pointKey}\u0001${this.concKeyAt(row)}`);
    for (const [profileKey, list] of this.lzStore) {
      const kept = list.filter((o) =>
        live.has(`${profileKey}\u0001${o.timeKey}#${o.ordinal}\u0001${o.concKey}`));
      dropped += list.length - kept.length;
      if (kept.length === 0) this.lzStore.delete(profileKey);
      else this.lzStore.set(profileKey, kept);
    }
    return dropped;
  }

  private runOptions(): RunNcaPipelineOptions {
    return {
      log: this.log, now: this.nowFn, lambdaZResolver: this.makeResolver(),
      ...(this.routeWarnings ? {routeWarningsBySubject: this.routeWarnings} : {}),
    };
  }

  /** The current results frame (same instance across recomputes — §9.4), or null. */
  get results(): DG.DataFrame | null {
    return this.outcome?.status === 'completed' ? this.outcome.results : null;
  }
  get lastOutcome(): NcaRunOutcome | null {return this.outcome;}

  /** Whole-cohort compute; per-subject overrides are re-applied from the keyed store.
   *  Afterwards any stored override that did not resolve on this frame is pruned and
   *  counted (AC-S2.1.3) — it was already ignored by the resolver, so pruning changes no
   *  number. A validation failure prunes nothing (there are no profiles to resolve on). */
  computeAll(options: {log?: AuditLog} = {}): NcaRunOutcome {
    const run = options.log !== undefined ? {...this.runOptions(), log: options.log} : this.runOptions();
    const outcome = runNcaPipeline(this.df, this.contract, this.scalars, this.rules, run);
    this.outcome = outcome;
    this.profileByKey.clear();
    this.fitRowsCache.clear();
    this.nonContiguousSet.clear();
    this.invalidateRowKeys();
    if (outcome.status === 'completed') {
      for (const pp of outcome.perProfile) {
        this.profileByKey.set(pp.index.subjectKey, pp.index);
        this.fitRowsCache.set(pp.index.subjectKey, pp.fitRows);
        if (pp.nonContiguous) this.nonContiguousSet.add(pp.index.subjectKey);
      }
      this.reconciled += this.pruneUnresolved();
    }
    return outcome;
  }

  /** The held profile index of one results-row subject, or undefined. */
  profileIndex(subjectKey: string): ProfileIndex | undefined {
    return this.profileByKey.get(subjectKey);
  }

  /**
   * AUC over [start, end] (hours, time after dose) of one subject's profile — sci-comp
   * `nca.aucInterval` (PKNCA `aucint.last`) on the SAME extracted, unit-canonical inputs and
   * BLQ rules the profile's NCA ran on (guardrail "one function, every caller"). NaN when the
   * subject is unknown or extraction fails. tk-study-day-stratification F4 (AUC0–τ past tlast).
   */
  intervalAuc(subjectKey: string, start: number, end: number): number {
    const index = this.profileByKey.get(subjectKey);
    if (index === undefined) return NaN;
    try {
      const {inputs} = extractProfileInputs(this.df, index, this.contract, this.scalars);
      return nca.aucInterval(inputs, this.rules, start, end).value;
    } catch (e) {
      // NaN reads as "AUC0–τ not available" downstream (the AUClast fallback label); keep the cause visible.
      console.warn(`NcaController.intervalAuc(${subjectKey}): ${e instanceof Error ? e.message : String(e)}`);
      return NaN;
    }
  }

  /** Recompute ONE subject in place (AC-1.2) — patches its results row, preserving
   *  frame identity; every other row is untouched. `auditAction` distinguishes a
   *  forward override (`'override-lambdaz'`) from a reset (`'reset-lambdaz'`, F-4). */
  recomputeProfile(subjectKey: string, auditAction = 'override-lambdaz'): void {
    if (this.outcome === null || this.outcome.status !== 'completed')
      throw new Error('NcaController.recomputeProfile: call computeAll() first.');
    const index = this.profileByKey.get(subjectKey);
    if (index === undefined)
      throw new Error(`NcaController.recomputeProfile: unknown subject '${subjectKey}'.`);
    const single = computeSingleProfile(
      this.df, index, this.contract, this.scalars, this.rules, this.runOptions());
    patchResultRow(this.outcome.results, index, single.perProfile.route,
      single.perProfile.result, single.c0, single.perProfile.doseNormalization, single.perProfile.bodyWeight);
    // The per-profile entry too, IN PLACE: the view holds this same array and builds the saved
    // summary (warnings, status) from it — a stale entry saved the pre-override warnings, and
    // the next open, a full recompute, rewrote them.
    const all = this.outcome.perProfile as PerProfileResult[];
    const at = all.findIndex((p) => p.index.subjectKey === subjectKey);
    if (at >= 0) all[at] = single.perProfile;
    // Keep the marker/toggle membership in sync with the patched fit — this ONE
    // subject changed; every other subject's cached set is untouched.
    this.fitRowsCache.set(subjectKey, single.perProfile.fitRows);
    if (single.perProfile.nonContiguous) this.nonContiguousSet.add(subjectKey);
    else this.nonContiguousSet.delete(subjectKey);
    this.log.append({
      timestamp: this.nowFn(), subjectKey, action: auditAction,
      ruleHash: computeRuleHash(this.rules), status: single.perProfile.result.status,
      payload: single.perProfile.result.provenance,
    });
  }

  /** Replace the global rules and recompute the cohort (AC-1.4). Sticky overrides
   *  are preserved: they re-apply on top of the new rules via `_lz_action`. */
  setRules(rules: nca.NcaRules): NcaRunOutcome {
    this.rules = rules;
    return this.computeAll();
  }

  /** Record per-point overrides in the keyed store (AC-1.3, F3; AD-2), then recompute that
   *  subject (AC-4.4: a reset logs `reset-lambdaz`). A `force-*` action stamps who/when; an
   *  `auto` action removes that point's override. */
  overrideLambdaZPoints(subjectKey: string, actions: ReadonlyArray<OverrideAction>): void {
    if (!this.profileByKey.has(subjectKey))
      throw new Error(`NcaController.overrideLambdaZPoints: unknown subject '${subjectKey}'.`);
    const keys = this.rowKeyMap();
    const {user, userId} = currentUser();
    const at = new Date(this.nowFn()).toISOString();
    let list = this.lzStore.get(subjectKey) ?? [];
    for (const a of actions) {
      const k = keys.get(a.rowIdx);
      if (k === undefined || k.profileKey !== subjectKey) continue;
      const [timeKey, ordinal] = splitPointKey(k.pointKey);
      list = list.filter((o) => !(o.timeKey === timeKey && o.ordinal === ordinal));
      if (a.action !== 'auto') {
        list.push({profileKey: subjectKey, timeKey, ordinal, concKey: this.concKeyAt(a.rowIdx),
          action: a.action, user, ...(userId !== null ? {userId} : {}), at});
      }
    }
    if (list.length === 0) this.lzStore.delete(subjectKey);
    else this.lzStore.set(subjectKey, list);
    // A batch that only resets points is an undo, not a forward override (F-4).
    const isReset = actions.length > 0 && actions.every((a) => a.action === 'auto');
    this.recomputeProfile(subjectKey, isReset ? 'reset-lambdaz' : 'override-lambdaz');
    // D13 "replaced": a point edit after Accept as is / Exclude λz replaces the decision — the
    // analyst is now looking at a different fit. The old decision moves to history; D11(1) then
    // applies to the new fit (Accepted → Not needed; still Flagged → a new decision is needed).
    this.supersede(subjectKey, 'replaced');
  }

  // ── Review decisions (F-S2.2) ────────────────────────────────────────────

  /** The persisted form of every ledger (v15 `reviewDecisions`). */
  getReviewLedgers(): ReviewLedger[] {return ledgerListOf(this.ledgers);}

  /** Replace the ledgers (reopen / Revert). Call {@link reconcileDecisions} after the next
   *  compute to drop decisions the fit no longer supports. */
  setReviewLedgers(list: ReadonlyArray<ReviewLedger>): void {this.ledgers = ledgerMapOf(list);}

  /** The current decision of one unit (valid after the last reconcile), or null. */
  decisionFor(unitKey: string): ReviewDecisionRecord | null {
    return this.ledgers.get(unitKey)?.current ?? null;
  }

  /** The ledger of one unit (current + history), or undefined. */
  ledgerFor(unitKey: string): ReviewLedger | undefined {return this.ledgers.get(unitKey);}

  /** Results-row keys whose current decision is Exclude λz — the {@link reportableResults}
   *  input (AD-6). Valid decisions only: {@link reconcileDecisions} has already moved every
   *  mismatched one to history. */
  excludedKeys(): Set<string> {
    const out = new Set<string>();
    for (const [k, l] of this.ledgers) if (l.current?.action === 'exclude') out.add(k);
    return out;
  }

  /** Record a decision (D11/D12). `fingerprint` is what the analyst saw (the unit's current
   *  decision fingerprint from {@link reviewSnapshot}). Audited as `review-<action>`. */
  decide(unitKey: string, action: ReviewAction, reasons: ReadonlyArray<ReviewReason>,
    fingerprint: string): void {
    const {user, userId} = currentUser();
    const rec: ReviewDecisionRecord = {
      action, reasons: [...reasons], user, ...(userId !== null ? {userId} : {}),
      at: new Date(this.nowFn()).toISOString(), fingerprint,
    };
    this.ledgers = recordDecision(this.ledgers, unitKey, rec);
    this.auditReview(unitKey, `review-${action}`, rec);
  }

  /** The analyst's Reset: the current decision moves to history; the unit is `To review`. */
  resetDecision(unitKey: string): void {this.supersede(unitKey, 'reset');}

  private supersede(unitKey: string, by: 'reset' | 'replaced'): void {
    const at = new Date(this.nowFn()).toISOString();
    const prev = this.ledgers.get(unitKey)?.current ?? null;
    const {map, changed} = supersedeDecision(this.ledgers, unitKey, by, at);
    this.ledgers = map;
    if (changed) {
      this.auditReview(unitKey, by === 'reset' ? 'review-reset' : 'review-auto-reset',
        {supersededBy: by, decision: prev});
    }
  }

  /** D13 auto-reset: every decision whose fingerprint no longer matches moves to history.
   *  Returns the reset count (the "N review decisions reset" notice, F-S3.6). `inScope` limits
   *  the pass to a subset of units (the dense profiles reconcile before the aggregates read the
   *  exclusion set; the sparse composites reconcile after their panel computes). */
  reconcileDecisions(
    fingerprints: ReadonlyMap<string, string>, inScope: (unitKey: string) => boolean = () => true,
  ): number {
    const at = new Date(this.nowFn()).toISOString();
    const before = this.ledgers;
    const {map, resetKeys} = reconcileLedgers(this.ledgers, fingerprints, at, inScope);
    this.ledgers = map;
    for (const k of resetKeys)
      this.auditReview(k, 'review-auto-reset', {supersededBy: 'auto-reset', decision: before.get(k)?.current});
    return resetKeys.length;
  }

  private auditReview(unitKey: string, action: string, payload: unknown): void {
    this.log.append({
      timestamp: this.nowFn(), subjectKey: unitKey, action,
      ruleHash: computeRuleHash(this.rules), status: 'ok', payload,
    });
  }

  /**
   * The review state of the held results (F-S1.3 / F-S1.4): every review unit, classified,
   * with its CURRENT decision fingerprint. The single reader for the Review column, the stage
   * strip, `Next to review`, the save-time `stageCache`, and the Confirm dry run — none of them
   * re-derives a unit or a fingerprint.
   *
   * @param t           the live acceptance thresholds (`R_acc` included).
   * @param opts.viable `denseProfilesViable` of the frame the results were computed on.
   * @param opts.composites  UC-04 sparse group composites (+ their fingerprints), if any.
   * @param opts.memberHashByProfile  NONMEM composite cells: the member-observation hash per
   *                    results-row key (AC-S1.4.3).
   */
  reviewSnapshot(t: AcceptanceThresholds, opts: {
    viable: boolean;
    composites?: ReadonlyArray<CompositeUnitInput & {fingerprint: string}>;
    memberHashByProfile?: ReadonlyMap<string, string>;
  }): {units: ReviewUnit[]; fingerprints: Map<string, string>} {
    const fingerprints = new Map<string, string>();
    const results = this.results;
    const profiles: ProfileUnitInput[] = [];
    if (results !== null) {
      const obs = this.contract.observations;
      const meas = profileMeasurability(this.df, obs.subjectId.name, obs.time.name,
        obs.concentration.name, obs.blqFlag?.name ?? null, obs.lloq?.name ?? null, NaN);
      const col = (n: string): DG.Column | null => results.col(n);
      const subj = col('SUBJECT');
      const status = col('STATUS');
      const warn = col('WARNINGS_JSON');
      const route = col('ROUTE');
      const num = (n: string, i: number): number => {const c = col(n); return c === null ? NaN : numAt(c, i);};
      const timeCol = this.df.col(obs.time.name);
      const concCol = this.df.col(obs.concentration.name);
      const doseCol = obs.dose !== undefined ? this.df.col(obs.dose.name) : null;
      const str = (c: DG.Column | null, r: number): string => (c === null || c.isNone(r) ? '' : c.getString(r));
      for (let i = 0; i < (subj?.length ?? 0); i++) {
        const key = subj!.getString(i);
        const acceptance = evaluateAcceptance({
          spanRatio: num('SPAN_RATIO', i), adjR2: num('LAMBDAZ_ADJ_R2', i),
          nPoints: num('LAMBDAZ_NPOINTS', i), pctExtrap: num('PCT_EXTRAP', i),
        }, t);
        const m = meas.get(key);
        profiles.push({
          key, status: status?.getString(i) ?? '', acceptance,
          measurableCount: m?.measurable ?? 0, blqTimes: m?.blqTimes ?? [],
        });
        const index = this.profileByKey.get(key);
        const rows = index !== undefined ? Array.from(index.observationRows) : [];
        const fitRows = this.fitRowsCache.get(key) ?? EMPTY_ROW_SET;
        let warnings: {code: string; severity: string}[] = [];
        try {
          const parsed = JSON.parse(str(warn, i) || '[]') as Array<{code?: string; severity?: string}>;
          warnings = parsed.map((w) => ({code: String(w.code ?? ''), severity: String(w.severity ?? '')}));
        } catch {/* malformed JSON → no warning component (the verdict still is) */}
        const member = opts.memberHashByProfile?.get(key);
        const fpIn: DecisionFingerprintInput = {
          points: [
            ...rows.map((r) => [str(timeCol, r), str(concCol, r)] as const),
            ...(member !== undefined ? [['members', member] as const] : []),
          ],
          dose: rows.length > 0 ? str(doseCol, rows[0]) || String(this.scalars.dose) : String(this.scalars.dose),
          route: str(route, i),
          pointsUsedTimes: rows.filter((r) => fitRows.has(r)).map((r) => str(timeCol, r)),
          lambdaZ: num('LAMBDA_Z', i), adjR2: num('LAMBDAZ_ADJ_R2', i),
          aucInf: num('AUC_INF', i), pctExtrap: num('PCT_EXTRAP', i),
          verdict: acceptance.verdict, flagCodes: acceptance.flagCodes, warnings,
        };
        fingerprints.set(key, computeDecisionFingerprint(fpIn));
      }
    }
    for (const c of opts.composites ?? []) fingerprints.set(`sparse:${c.groupKey}`, c.fingerprint);
    const units = classifyReviewUnits({
      profiles, viable: opts.viable, composites: opts.composites ?? [],
      ledgers: this.ledgers, fingerprints,
    });
    // Only units that NEED review may hold a decision; pass just those to reconcile.
    const reviewFps = new Map<string, string>();
    for (const u of units) if (u.need === 'review') reviewFps.set(u.key, fingerprints.get(u.key) ?? '');
    return {units, fingerprints: reviewFps};
  }

  /** Per-subject REAL λz fit membership (source rows) — the single source for the
   *  facet used-point markers AND the click-toggle's in/out decision (F6). Replaces
   *  the time-window proxy so a marker cannot disagree with the reported fit. Empty
   *  set for a subject with no fit; empty map before the first Compute. */
  fitRowsBySubject(): ReadonlyMap<string, ReadonlySet<number>> {
    return this.fitRowsCache;
  }

  /** Subjects whose λz selection is NON-CONTIGUOUS (decision 2) — the facet caption
   *  note source. Empty before the first Compute; a subject leaves the set when a
   *  later edit / reset makes its selection contiguous again. */
  nonContiguousSubjects(): ReadonlySet<string> {
    return this.nonContiguousSet;
  }

  /** The stored λz action for one row of {@link frame} (`auto` when unset) — the reliable
   *  basis for a toggle decision (with {@link fitRowsBySubject} for membership). */
  currentLzAction(rowIdx: number): LzAction {
    return this.overrideAtRow(rowIdx)?.action ?? 'auto';
  }

  /** For each subject currently under manual λz override, the analyst label from
   *  `_lz_user` (F7 / AC-1.3) — the source for the results grid's ✏️ indicator.
   *  Empty map when nothing is overridden. Reads the SAME keyed store as
   *  {@link makeResolver}, so the indicator cannot disagree with what actually
   *  forced the fit. Value may be `''` (override made where no shell user was
   *  resolvable); map MEMBERSHIP, not a non-empty value, is the overridden test. */
  overriddenSubjects(): Map<string, string> {
    const out = new Map<string, string>();
    for (const [key, index] of this.profileByKey) {
      const user = this.overrideUserForIndex(index);
      if (user !== null) out.set(key, user);
    }
    return out;
  }

  /** The analyst label for ONE subject's override, or `null` when it is on
   *  auto-best-fit — the single-subject counterpart to {@link overriddenSubjects}
   *  for the fast post-toggle flag update (F7). A non-null `''` still means
   *  overridden (unknown user). */
  overrideUserFor(subjectKey: string): string | null {
    const index = this.profileByKey.get(subjectKey);
    return index === undefined ? null : this.overrideUserForIndex(index);
  }

  /** The `_lz_time` stamp for ONE subject's override, or `null` when it is on
   *  auto-best-fit — the timestamp counterpart to {@link overrideUserFor} for the
   *  fast post-toggle results-tooltip update. Same representative point as the
   *  analyst (both from {@link overrideProvenanceForIndex}). Raw ISO; the display
   *  surface formats it via `formatOverrideTimestamp`. */
  overrideTimeFor(subjectKey: string): string | null {
    const index = this.profileByKey.get(subjectKey);
    return index === undefined ? null : (this.overrideProvenanceForIndex(index)?.time ?? null);
  }

  /** For each subject currently under manual λz override, the `_lz_time` stamp of the
   *  SAME representative point {@link overriddenSubjects} takes the analyst from — the
   *  full-pass source for the results-grid ✎ tooltip's who/when line. Map MEMBERSHIP
   *  matches {@link overriddenSubjects} (an override with a blank clock still appears,
   *  value `''`). */
  overriddenSubjectTimes(): Map<string, string> {
    const out = new Map<string, string>();
    for (const [key, index] of this.profileByKey) {
      const prov = this.overrideProvenanceForIndex(index);
      if (prov !== null) out.set(key, prov.time);
    }
    return out;
  }

  /** Who + when for ONE source observation row, or `null` when that row is `auto`
   *  (not overridden) — the row-keyed provenance accessor the λz-points context-panel
   *  grid reads for its `Analyst` / `Recorded` columns (TODO.md AC-5.1). Parallel to
   *  {@link currentLzAction}: a null result is the per-row not-overridden test, so the
   *  grid can decide column visibility from it. `user`/`time` may be `''`; `time` is
   *  raw ISO (the caller formats it via `formatOverrideTimestamp`). */
  lzProvenanceForRow(rowIdx: number): {user: string; time: string} | null {
    const o = this.overrideAtRow(rowIdx);
    return o === null ? null : {user: o.user, time: o.at};
  }

  /** The `_lz_user` of the first force-in/out point in a subject's rows, or `null`
   *  when none is forced. Delegates to {@link overrideProvenanceForIndex} so the
   *  analyst and the timestamp always come from the SAME representative point. */
  private overrideUserForIndex(index: ProfileIndex): string | null {
    return this.overrideProvenanceForIndex(index)?.user ?? null;
  }

  /** Who + when of the first force-in/out point in a subject's rows, or `null` when
   *  none is forced — the ONE scan behind {@link overrideUserFor} / {@link overrideTimeFor}
   *  / {@link overriddenSubjects} / {@link overriddenSubjectTimes}, so who and when
   *  can never refer to different points. `user`/`time` may be `''` (override made
   *  where no shell user / clock was resolvable); a non-null result is the overridden
   *  test, exactly as before. */
  private overrideProvenanceForIndex(index: ProfileIndex): {user: string; time: string} | null {
    if (!this.hasOverrides(index.subjectKey)) return null;
    // First overridden point in ROW order (the representative point), so who and when
    // always refer to the same point.
    const rows = index.observationRows;
    for (let i = 0; i < rows.length; i++) {
      const o = this.overrideAtRow(rows[i]);
      if (o !== null) return {user: o.user, time: o.at};
    }
    return null;
  }

  private makeResolver(): (index: ProfileIndex, inputs: nca.ProfileInputs) => Int32Array | null {
    const blq = this.rules.blq;
    const autoRules = this.rules.lambdaZ.mode === 'auto-best-fit' ?
      this.rules : withLambdaZ(this.rules, {mode: 'auto-best-fit', manualPoints: undefined});
    return (index, inputs) =>
      resolveManualPoints(this.actionsForIndex(index), index, inputs, blq, autoRules);
  }
}

/** Split a `timeKey#ordinal` point key. The ordinal is the LAST `#` segment, so a time
 *  string that itself contains `#` still round-trips. */
function splitPointKey(pointKey: string): [string, number] {
  const i = pointKey.lastIndexOf('#');
  return [pointKey.slice(0, i), Number(pointKey.slice(i + 1))];
}

/**
 * Resolve one subject's `manualPoints` from its row → action map (from the keyed store),
 * or `null` when the subject is not under manual override (stays auto-best-fit). Composes
 * the verified core: build the augmented profile + `augRowId`, read the current
 * auto-window rows (a fresh auto `computeNca`), then map the force-in/out/auto
 * actions onto augmented indices.
 */
function resolveManualPoints(
  actionByRow: Map<number, LzAction>, index: ProfileIndex, inputs: nca.ProfileInputs,
  blqStrategy: nca.BlqStrategy, autoRules: nca.NcaRules,
): Int32Array | null {
  if (!subjectHasManualOverride(actionByRow)) return null;
  const aug = buildAugmentedProfile(inputs, index.observationRows, blqStrategy);
  const autoRes = nca.computeNca(inputs, autoRules);
  const autoInRows = new Set<number>(
    autoRes.provenance.lambdaZ !== null ?
      pointsUsedToRows(autoRes.provenance.lambdaZ.pointsUsed, aug.rowId) : []);
  return buildManualPoints(aug.rowId, autoInRows, actionByRow);
}

/** The current platform user, stamped on an override (AC-1.3) and a decision: the display label
 *  the UI shows, and the user id for attribution (IR-12). Guarded — the headless/test context may
 *  have no shell user (label '', id null). */
function currentUser(): {user: string; userId: string | null} {
  try {
    const u = grok.shell.user;
    return u ? {user: u.friendlyName ?? u.id ?? '', userId: u.id ?? null} : {user: '', userId: null};
  } catch {
    return {user: '', userId: null};
  }
}
