/**
 * P-B1c — index translation between source-observation rows and the sci-comp
 * engine's internal AUGMENTED index space, for the interactive λz manual
 * point-selection feature (VR-11 / FR-224). Overrides are keyed by source-row identity.
 *
 * **Why this exists.** `computeNca` does not fit λz on the caller's `ProfileInputs`
 * arrays — it fits on an AUGMENTED array (post-BLQ, with a synthetic t=0 point
 * inserted for IV-bolus C0 / extravascular 0 unless a kept dose-time row already
 * exists; sci-comp `augment.ts`). So `NcaRules.lambdaZ.manualPoints` are indices
 * into THAT array, and `LambdaZResult.pointsUsed` comes back in the same space —
 * NOT indices of the source observation rows. A naive `index ± offset` translation
 * is the documented DOC-FIX trap (and is WRONG since the route-aware gate — see
 * below), and a `by time value` translation (an earlier draft) is ambiguous for
 * duplicate-time rows (replicate / QC re-draws) and fragile under unit-converted
 * float times.
 *
 * **The mechanism (AC-2.1): a stable `augRowId`.** The augmented profile is the
 * engine's OWN — `nca.augmentProfile` (sci-comp ≥ 0.11.0) is the exported kernel
 * `computeNca` itself runs for pipeline Steps 1–3 (BLQ → observed Cmax →
 * route-aware dose-time augmentation), so this module re-derives NO orchestration
 * (rule 5). The kernel reports `sourceIndex[k]`: the INPUT-array index of
 * augmented element `k`, or `-1` for a synthetic point; `augRowId[k]` is simply
 * `observationRows[sourceIndex[k]]` (or `-1`). Every translation keys on this row
 * id, never on float time, so duplicate-time rows get distinct ids and no
 * float-equality tolerance is needed (N-1/N-2, AC-2.3).
 *
 * **Why `sourceIndex`, not a positional offset.** The pre-0.11.0 mirror assumed
 * augmented element `k` maps to `observationRows[k − offset]`. That broke with the
 * route-aware gate (GROK-20960): an IV-bolus dose-time row that is BLQ-flagged or
 * non-positive is REPLACED — index 0 becomes the synthetic `(0, c0)` and the
 * input's row 0 is gone — so the augmented array has the input's length while
 * input row 0 appears nowhere in it. Only an explicit map can say that; the
 * kernel is the one place that knows.
 */

import {nca} from '@datagrok-libraries/sci-comp';

/** Per-point λz-membership action stored in the `_lz_action` source column (D5). */
export type LzAction = 'auto' | 'force-in' | 'force-out';

// --- Service columns on the SOURCE observations df (F3, D5) ---------------------
// The override state lives ONLY here (single source of truth, D4); these columns
// are `_lz`-prefixed and hidden from every user-facing surface via the shared
// {@link isLzServiceColumn} predicate (AC-3.1: passthrough picker, mapping-pane
// role pickers, computeContentHash, raw source grid). Kept in this DG-free module
// so all four exclusion sites + the controller import one spelling.

/** `_lz_action` (string enum) — the ONLY service column read by compute. */
export const LZ_ACTION_COL = '_lz_action';
/** `_lz_user` — platform user who made the override (metadata). */
export const LZ_USER_COL = '_lz_user';
/** `_lz_time` — timestamp of the override (metadata). */
export const LZ_TIME_COL = '_lz_time';
/** All λz service columns, in creation order. (`_lz_srchash` retired with the keyed override
 *  store — F-S2.1; a legacy frame carrying it is still caught by the `_lz` prefix rule.) */
export const LZ_SERVICE_COLUMNS: readonly string[] = [
  LZ_ACTION_COL, LZ_USER_COL, LZ_TIME_COL,
];

/**
 * Is `name` a λz override service column? The single predicate the four
 * exclusion surfaces share (AC-3.1) — a `_lz`-prefix test, so a future
 * `_lz_*` column is covered with no edit at the call sites.
 */
export function isLzServiceColumn(name: string): boolean {
  return name.startsWith('_lz');
}

/**
 * Format a stored `_lz_time` ISO stamp (`new Date(now).toISOString()`, UTC) into a
 * compact, scientist-readable `YYYY-MM-DD HH:MM` in the viewer's LOCAL time — the
 * single formatter shared by the results-grid ✎ tooltip (who/when line) and the
 * λz-points context-panel grid's `Recorded` column, so the two cannot drift.
 *
 * DG-free (pure string/Date), lives here beside the `_lz_time` column contract.
 * Local (not UTC) because the scientist reading "when did I change this" expects
 * their own clock; the surfaces that show it label it "(local time)". A blank input
 * returns blank (no override); an unparseable stamp returns the raw string rather
 * than an empty cell, so a corrupted value is visible, not silently dropped.
 */
export function formatOverrideTimestamp(iso: string): string {
  if (iso === '') return '';
  const d = new Date(iso);
  if (Number.isNaN(d.getTime())) return iso;
  const p = (n: number): string => String(n).padStart(2, '0');
  return `${d.getFullYear()}-${p(d.getMonth() + 1)}-${p(d.getDate())} ` +
    `${p(d.getHours())}:${p(d.getMinutes())}`;
}

/**
 * The engine's augmented profile plus the parallel source-row id map. `time` /
 * `conc` / `blqMask` are index-identical to what `computeNca` fits λz on
 * (`blqMask` is the kernel's EFFECTIVE BLQ mask — input mask ∪ excluded — the
 * Cmax/Tlag semantics; the drop set the integrator uses is `dropMask`);
 * `rowId[k]` is the source-df row index of augmented element `k`, or `-1` for a
 * synthetic point (a prepended C0 / extravascular-0, or the `(0, c0)` that
 * REPLACED a BLQ / non-positive IV-bolus dose-time row).
 */
export interface AugmentedProfile {
  readonly time: Float64Array;
  readonly conc: Float64Array;
  readonly blqMask: Uint8Array;
  readonly rowId: Int32Array;
}

/** An empty augmented profile — the kernel's `null` (no measurable point). */
const EMPTY_AUGMENTED: AugmentedProfile = {
  time: new Float64Array(0), conc: new Float64Array(0),
  blqMask: new Uint8Array(0), rowId: new Int32Array(0),
};

/**
 * The engine's augmented profile + `augRowId` for one subject — a thin wrapper
 * over `nca.augmentProfile`, the SAME kernel `computeNca` runs (rule 5: no
 * orchestration re-derived here).
 *
 * `observationRows` are the source-df row indices in the SAME order they were
 * copied into `inputs` by `extractProfileInputs` (i.e. `index.observationRows`);
 * `rowId[k] === observationRows[aug.sourceIndex[k]]` for every non-synthetic
 * element and `-1` where `sourceIndex[k] < 0`. A profile with no measurable
 * point (the kernel returns `null`) yields an empty profile — such a subject is
 * `'failed'`, has no λz and nothing to translate.
 */
export function buildAugmentedProfile(
  inputs: nca.ProfileInputs,
  observationRows: Int32Array,
  blqStrategy: nca.BlqStrategy,
): AugmentedProfile {
  const aug = nca.augmentProfile(inputs, blqStrategy);
  if (aug === null) return EMPTY_AUGMENTED;
  const rowId = new Int32Array(aug.sourceIndex.length);
  for (let k = 0; k < rowId.length; k++) {
    const src = aug.sourceIndex[k];
    rowId[k] = src < 0 ? -1 : observationRows[src];
  }
  return {time: aug.time, conc: aug.conc, blqMask: aug.blqMask, rowId};
}

/**
 * Map an engine-returned `pointsUsed` (augmented-relative) back to the source
 * observation rows that were actually used in the λz fit, dropping the synthetic
 * prepended point (`rowId === -1`). Used to render used-point markers and to read
 * the auto-best-fit window as source rows (AC-2.2).
 */
export function pointsUsedToRows(pointsUsed: Int32Array, augRowId: Int32Array): number[] {
  const out: number[] = [];
  for (let k = 0; k < pointsUsed.length; k++) {
    const augIdx = pointsUsed[k];
    if (augIdx < 0 || augIdx >= augRowId.length) continue;
    const row = augRowId[augIdx];
    if (row >= 0) out.push(row);
  }
  return out;
}

/** Warning code minted when a manual λz selection is NON-CONTIGUOUS — a measurable
 *  point sits strictly inside the fitted time window but is excluded from the fit.
 *  Surfaced in `provenance.warnings` (WARNINGS_JSON) + the facet caption. */
export const LAMBDAZ_NONCONTIGUOUS = 'LAMBDAZ_NONCONTIGUOUS';

/**
 * Is the reported λz fit **non-contiguous** — i.e. some MEASURABLE point whose time
 * lies strictly between the fit's first and last used times is NOT in the fit
 * (decision 2, 2026-09-01)? By accepted NCA practice λz is fit to CONSECUTIVE
 * terminal points; the auto-best-fit is always contiguous by construction, so this
 * can only fire on a manual selection that skips an interior point (a legitimate
 * outlier drop, but non-standard — worth flagging, not blocking).
 *
 * Works in TIME (not augmented-index) space, so it is robust to any array ordering:
 * measurability is `conc > 0 && finite` — the SAME eligibility a λz fit uses, so a
 * BLQ / ≤0 interior point is not miscounted as a gap.
 */
export function isLambdaZNonContiguous(pointsUsed: Int32Array, aug: AugmentedProfile): boolean {
  if (pointsUsed.length < 2) return false; // a 0/1-point fit cannot have an interior gap
  const used = new Set<number>(Array.from(pointsUsed));
  let tFirst = Infinity;
  let tLast = -Infinity;
  for (let k = 0; k < pointsUsed.length; k++) {
    const t = aug.time[pointsUsed[k]];
    if (t < tFirst) tFirst = t;
    if (t > tLast) tLast = t;
  }
  for (let i = 0; i < aug.time.length; i++) {
    if (used.has(i)) continue;
    const t = aug.time[i];
    const c = aug.conc[i];
    if (t > tFirst && t < tLast && c > 0 && Number.isFinite(c)) return true;
  }
  return false;
}

/**
 * Build `NcaRules.lambdaZ.manualPoints` (augmented indices) for the effective
 * "in the λz fit" set of a subject in manual mode (AC-2.1).
 *
 * The effective set starts from the auto-best-fit window membership at click time
 * (`autoInRows` — the source rows the engine's auto fit used, via
 * {@link pointsUsedToRows}) and applies the per-row deltas: `force-in` adds a row,
 * `force-out` removes it, `auto` follows the auto window. Rows are mapped to their
 * augmented index via `augRowId`; the synthetic prepended point (id `-1`) is never
 * selectable. Result is sorted ascending, as the engine expects.
 */
export function buildManualPoints(
  augRowId: Int32Array,
  autoInRows: ReadonlySet<number>,
  actionByRow: ReadonlyMap<number, LzAction>,
): Int32Array {
  const out: number[] = [];
  for (let augIdx = 0; augIdx < augRowId.length; augIdx++) {
    const row = augRowId[augIdx];
    if (row < 0) continue; // synthetic prepend — never manually selectable
    const action = actionByRow.get(row) ?? 'auto';
    const isIn =
      action === 'force-in' ? true :
        action === 'force-out' ? false :
          autoInRows.has(row);
    if (isIn) out.push(augIdx);
  }
  return Int32Array.from(out);
}

/**
 * Is this subject in manual λz mode? True iff any of its rows carries a
 * non-`auto` action. A fully-auto subject stays on `auto-best-fit` (the caller
 * must NOT build `manualPoints` for it — AC-3.2 clean-frame identity).
 */
export function subjectHasManualOverride(actionByRow: ReadonlyMap<number, LzAction>): boolean {
  for (const action of actionByRow.values())
    if (action === 'force-in' || action === 'force-out') return true;
  return false;
}

/** Warning code emitted when a manual λz selection falls below the ≥minPoints floor. */
export const LAMBDAZ_BELOW_MIN_POINTS = 'LAMBDAZ_BELOW_MIN_POINTS';

/**
 * Strict ≥`minPoints` floor for a MANUAL λz fit (AC-4.3, user decision 2026-08-31).
 *
 * `sci-comp`'s `lambdaZManual` applies NO minimum-points floor — it fits any set
 * of ≥2 points — but the user chose the universal PKNCA/Phoenix convention
 * (`min.hl.points` = 3): a fit using fewer than `minPoints` points is not a
 * low-confidence estimate, it is `Not estimated`. This is nca-studio's acceptance
 * layer, applied AFTER `computeNca` (the gate is the caller's, not sci-comp's —
 * rule 5). Auto-best-fit already rejects sub-`minPoints` windows, so this only
 * ever bites a manual selection.
 *
 * Fires on EITHER limb of AC-4.3: (a) fewer than `minPoints` points, OR (b) a
 * **non-decaying fit (`λz ≤ 0`, or non-finite)** — an analyst who force-includes
 * rising/absorption-phase points can make `lambdaZManual` return a non-positive
 * slope, which sci-comp deliberately leaves for the caller to reject
 * (`lambda-z.ts` "the returned slope can have any sign"). Either way λz and every
 * λz-derived parameter (AUCinf, t½, CL, Vz, AUMCinf, MRT, Vss, %extrap, %AUMCextrap)
 * are suppressed, `provenance.lambdaZ` is nulled, status is downgraded to `'partial'`
 * (Cmax/AUClast survive — λz-independent), and a `LAMBDAZ_BELOW_MIN_POINTS` warning
 * is attached. The result is then indistinguishable from a genuine no-λz `'partial'`,
 * so the grid renders `Not estimated`. Returns the input unchanged only when λz was
 * already absent, OR the fit has ≥`minPoints` points AND a positive slope.
 */
export function applyMinPointsFloor(
  result: nca.ComputeResult, minPoints: number,
): nca.ComputeResult {
  const lz = result.provenance.lambdaZ;
  if (lz === null) return result;
  const belowPoints = lz.pointsUsed.length < minPoints;
  const nonDecaying = !(lz.lambdaZ > 0); // true for λz ≤ 0 AND for a NaN slope
  if (!belowPoints && !nonDecaying) return result;

  const v = result.values;
  const reason = belowPoints ?
    `used ${lz.pointsUsed.length} point(s), below the minimum of ${minPoints} ` +
      `(PKNCA/Phoenix convention)` :
    `gave a non-decaying fit (λz = ${lz.lambdaZ} ≤ 0)`;
  return {
    status: 'partial',
    values: {
      // λz-independent — preserved.
      cmax: v.cmax, tmax: v.tmax, aucLast: v.aucLast, aumcLast: v.aumcLast, tlag: v.tlag,
      // λz-derived — suppressed.
      aucInf: NaN, pctExtrap: NaN, lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN,
      aumcInf: NaN, mrt: NaN, vss: NaN, pctExtrapAumc: NaN,
    },
    provenance: {
      ...result.provenance,
      lambdaZ: null,
      warnings: [
        ...result.provenance.warnings,
        {
          code: LAMBDAZ_BELOW_MIN_POINTS,
          severity: 'warning',
          message: `Manual λz selection ${reason} — λz not estimated.`,
        },
      ],
    },
  };
}

