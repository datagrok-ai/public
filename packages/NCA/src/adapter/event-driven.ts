/**
 * Event-driven (NONMEM) ingestion — `DataFrame → DataFrame` normalization stage.
 *
 * **AD-1 — event ingestion is a pure normalization stage.** A NONMEM-layout
 * observations table (EVID / AMT / CMT / RATE / MDV [+ ADDL / II / SS]) is
 * transformed into the *canonical* observations DataFrame the existing,
 * PKNCA-validated pipeline already consumes. The NCA math, BLQ phasing, λz,
 * CL/Vz/Vss, sparse composites and TK accumulation are untouched — this module
 * only changes the *input shape*:
 *
 *  - row classification by EVID/MDV          (G-1, P-5, AD-6) — dose/obs/reset
 *    rows are split into separate streams, which gives the "dose precedes a
 *    co-timed obs" guarantee (G-9, P-5) structurally: an obs anchors via binary
 *    search over the sorted dose-ladder, so no global `(ID,TIME,EVID desc)` row
 *    sort is needed (dose events ARE time-sorted before the ladder is built)
 *  - ADDL/II → dose ladder                   (G-4, P-3, AD-3) — both off-by-one guards
 *  - TAD = TIME − latest-dose ≤ TIME         (G-5, P-3) — inclusion `0 ≤ TAD ≤ II`
 *  - AMT as the per-interval dose of record  (G-2, P-2, AD-4)
 *  - route from ROUTE col / RATE / CMT       (G-3, P-2, AD-5) — inference flagged for confirm
 *  - covariate hoist onto every obs row      (G-6)
 *  - occasion derivation + AD-8 pseudo-key   (U-6, AD-8)
 *
 * The emitted frame carries the canonical observations contract columns
 * ({@link SUBJECT_COL} = TAD-keyed pseudo-subject, {@link TAD_COL} as `time`,
 * {@link CONC_COL} as `concentration`, {@link DOSE_COL} = AMT, the derived
 * `~RouteCode` column, plus hoisted covariates), so it flows unchanged into
 * `runNcaPipeline` (dense / serial designs) **or** `computeSparseNca` (the
 * destructive-composite TK path) — the consumer's choice, not this module's.
 *
 * Deterministic and pure: reads `df`, returns a fresh frame. Unit-testable
 * against the two PKNCA NONMEM cards without a Datagrok UI runtime.
 *
 * **Scope / honest-uncertainty (CLAUDE.md §4, DATA-GAP-27/28).** RATE>0
 * (infusion), SS≠0, EVID∈{2,3,4} and per-day-reset clocks are **built** but are
 * NOT exercised by the two oracle fixtures; each is gated behind a named,
 * severity-classified warning (AD-6/AD-7) or a blocking `ValidationError`
 * (SS reject) — never silent best-effort.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';
import {nca} from '@datagrok-libraries/sci-comp';

import type {ValidationError} from './column-contract';
import {stringifyKey} from './adapter-utils';
import {normalizeRouteString, DERIVED_ROUTE_CODE_COL} from './route-normalize';

// ── Emitted (canonical) column names ───────────────────────────────────────
/** Pseudo-subject key column (AD-8 — occasion-scoped, e.g. `"R001_Day28"`). */
export const SUBJECT_COL = 'Subject';
/** Time-after-dose column — becomes the canonical `time` (G-5, P-3). */
export const TAD_COL = 'TAD';
/** Concentration column (DV of usable observation rows). */
export const CONC_COL = 'Conc';
/** Per-interval dose-of-record column (AMT, absolute mg — AD-4 / AC-A6). */
export const DOSE_COL = 'Dose_mg';
/** Derived occasion / study-day column (DAY when present, else ladder index). */
export const OCCASION_COL = 'Occasion';
/** Per-subject infusion duration (AMT/RATE for RATE>0 IV-infusion, else null →
 *  NaN). Emitted so the per-subject duration reaches the pipeline the same way
 *  per-subject route does via `~RouteCode` (Slice A2 / BD-2): `runNcaPipeline`
 *  takes ONE global `scalars`, so a global `infusionDuration` would apply one
 *  subject's duration to every subject — wrong for multi-subject infusion
 *  datasets. `extractProfileInputs` reads this column per profile, falling back
 *  to `scalars.infusionDuration` when the cell is null. */
export const DERIVED_INFUSION_DURATION_COL = '~InfusionDuration';

// ── New ParameterWarning codes (open `nca.ParameterWarning.code` union — no
//    sci-comp change; AC-D4). Severity is curated per AD-6 (info vs warning). ──
/** A subject's first sorted row is not a dose → covariate forward-fill cannot
 *  seed; prevents silent NA fill (G-6). */
export const NONMEM_SUBJECT_FIRST_ROW_NOT_DOSE = 'NONMEM_SUBJECT_FIRST_ROW_NOT_DOSE';
/** Two distinct subjects collapsed to one occasion-scoped key (S-6, AD-8). */
export const NONMEM_OCCASION_KEY_COLLISION = 'NONMEM_OCCASION_KEY_COLLISION';
/** RATE>0 infusion: duration = AMT/RATE inferred (no explicit duration column);
 *  the scientist should confirm the route. Informational, not a math-validity
 *  claim — the RATE>0 MRT/Vss path is PKNCA-parity-validated (fixture 1, AC-1.3)
 *  and threaded per subject (Slice A2). BD-1 (renamed from
 *  `NONMEM_INFUSION_RATE_UNVALIDATED`, 2026-06-19): "UNVALIDATED" became a
 *  misnomer once the oracle landed; severity downgraded warning→info to cut
 *  warning-stream noise (gap-18 lesson). `routeConfirmRequired` still fires. */
export const NONMEM_INFUSION_RATE_INFERRED = 'NONMEM_INFUSION_RATE_INFERRED';
/** Route was inferred from RATE/CMT (no explicit ROUTE column) — requires a
 *  mandatory scientist confirmation before Compute (AD-5 / FM-4). */
export const NONMEM_ROUTE_INFERRED = 'NONMEM_ROUTE_INFERRED';
/** EVID∈{3,4} reset → a new profile episode opens for that ID (AD-6). */
export const EVID_RESET_NEW_EPISODE = 'EVID_RESET_NEW_EPISODE';
/** EVID=2 (other event) logged + dropped from the conc stream (AD-6, info). */
export const EVID_OTHER_IGNORED = 'EVID_OTHER_IGNORED';
/** SS=1 honored as a steady-state interval (τ=II) — UNVALIDATED (AD-7). */
export const NONMEM_STEADY_STATE_ASSUMED = 'NONMEM_STEADY_STATE_ASSUMED';

/**
 * NONMEM role → source-column-name binding the transform consumes. Distinct
 * from `ColumnContract` (which describes the *output* frame): some roles are
 * legitimately absent (`03` has no `addl`/`ii`/`ss`). The mapping pane bridges
 * the all-optional `EventDrivenDetection` (detector output) → this
 * required-present contract at the confirm step (architect SIMPLIFY-2).
 */
export interface EventContract {
  /** Subject / animal id column (NONMEM `ID`, or a CRO `AnimalID`). Required. */
  readonly id: string;
  /** Absolute clock time column (NONMEM `TIME`). Required. */
  readonly time: string;
  /** Dependent-variable (concentration) column (NONMEM `DV`). Required. */
  readonly dv: string;
  /** Dose-amount column (NONMEM `AMT`). Required (the dose of record, AD-4). */
  readonly amt: string;
  /** Event-id column (NONMEM `EVID`). Required (the classification key). */
  readonly evid: string;
  /** Missing-dependent-variable flag column (NONMEM `MDV`). Optional. */
  readonly mdv?: string;
  /** Compartment column (NONMEM `CMT`). Optional (route inference, AD-5). */
  readonly cmt?: string;
  /** Infusion-rate column (NONMEM `RATE`). Optional (RATE>0 ⇒ infusion, AD-5). */
  readonly rate?: string;
  /** Additional-doses count column (NONMEM `ADDL`). Optional (AD-3). */
  readonly addl?: string;
  /** Inter-dose interval column (NONMEM `II`). Optional (AD-3). */
  readonly ii?: string;
  /** Steady-state flag column (NONMEM `SS`). Optional (AD-7). */
  readonly ss?: string;
  /** Explicit route column (e.g. `ROUTE` = IV/PO/SC). Optional; when present it
   *  takes precedence over RATE/CMT inference (AD-5). */
  readonly route?: string;
  /** Study-day / occasion column (e.g. `DAY`). Optional; drives occasion when
   *  present, else occasion is derived from the dose-ladder index (U-6). */
  readonly day?: string;
  /** Body-weight column (e.g. `BW_kg`). Optional; hoisted + used downstream. */
  readonly bodyWeight?: string;
  /** Lower-limit-of-quantification column. Optional; hoisted. */
  readonly lloq?: string;
  /** BLQ-flag column. Optional; hoisted. */
  readonly blqFlag?: string;
  /** Extra covariate columns to forward-fill onto every obs row (Sex, DOSE_MGKG,
   *  NOAEL_mg_kg, Human*…). Optional. */
  readonly covariates?: ReadonlyArray<string>;
}

/** Options for {@link normalizeEventDrivenObservations}. */
export interface NormalizeEventOptions {
  /**
   * When route is *inferred* (no explicit ROUTE column), a scientist-confirmed
   * route applied uniformly to every inferred subject. `null`/absent ⇒ the
   * transform reports {@link NormalizeEventResult.inferredRoutes} and emits a
   * {@link NONMEM_ROUTE_INFERRED} warning; the mapping pane blocks Compute until
   * the route is confirmed (AD-5 mandatory-confirm gate).
   */
  readonly confirmedRoute?: nca.RouteCode | null;
}

/** Per-subject inferred-route record (route-confirm step, AD-5). */
export interface InferredRoute {
  readonly subjectKey: string;
  readonly code: nca.RouteCode;
  /** True when the code came from RATE/CMT inference (needs confirm), false when
   *  resolved from an explicit ROUTE column. */
  readonly inferred: boolean;
  /** Infusion duration `AMT/RATE` when RATE>0 (AD-5), else null. */
  readonly infusionDuration: number | null;
}

/** Row-classification tally (returned for the pane's confirm panel + tests). */
export interface ClassificationCounts {
  readonly dose: number;
  readonly obs: number;
  readonly reset: number;
  readonly other: number;
  readonly missingDv: number;
}

/** Result of {@link normalizeEventDrivenObservations}. Ingestion provenance is
 *  flattened onto the return (architect SIMPLIFY-1 — single confirmed consumer:
 *  the pane's confirm step). */
export interface NormalizeEventResult {
  /** The canonical observations frame the existing pipeline consumes. */
  readonly observations: DG.DataFrame;
  /** Per-subject ingestion warnings (open-union codes), keyed by pseudo-subject. */
  readonly warningsBySubject: ReadonlyMap<string, ReadonlyArray<nca.ParameterWarning>>;
  /** Blocking errors (non-empty ⇒ the pane shows the error panel, no compute). */
  readonly errors: ReadonlyArray<ValidationError>;
  /** Expanded dose ladder per source subject id (AD-3) — for the confirm panel. */
  readonly doseLadderBySubject: ReadonlyMap<string, ReadonlyArray<number>>;
  /** Inferred / resolved routes per pseudo-subject (AD-5 confirm gate). */
  readonly inferredRoutes: ReadonlyArray<InferredRoute>;
  /** True iff ≥1 route was inferred (no explicit ROUTE column) and not confirmed
   *  — the mapping pane must block Compute (AD-5). */
  readonly routeConfirmRequired: boolean;
  /** Classification tally over the whole input. */
  readonly classificationCounts: ClassificationCounts;
}

/** EVID values classified as a dose event (`1` = dose, `4` = reset + dose). */
const DOSE_EVIDS: ReadonlySet<number> = new Set([1, 4]);

/**
 * Build the ADDL/II dose ladder for one explicit dose row (AD-3).
 *
 * Ladder = `{t, t+II, …, t+ADDL·II}` — **exactly `1 + ADDL` doses, starting at
 * the dose row's own time** (Bauer 2019: ADDL is *additional* doses, so the
 * total is `1 + ADDL`). Both off-by-one directions are guarded (R2 N-1): a
 * phantom dose at `(ADDL+1)·II` and a missing initial dose (ladder starting at
 * `t+II`) each corrupt the Day-28 anchor. The function asserts
 * `length === 1 + ADDL` and `ladder[0] === t`.
 *
 * @throws if the constructed ladder violates either invariant (defensive — the
 *   construction below cannot, but a future edit that does must fail loudly;
 *   AC-B4 also calls {@link assertDoseLadder} directly on hand-built ladders).
 */
export function expandDoseLadder(tDose: number, addl: number, ii: number): number[] {
  if (!Number.isFinite(addl) || addl <= 0) return [tDose];
  if (!Number.isFinite(ii) || ii <= 0) return [tDose];
  const n = 1 + Math.round(addl);
  const ladder = new Array<number>(n);
  for (let k = 0; k < n; k++) ladder[k] = tDose + k * ii;
  assertDoseLadder(ladder, tDose, Math.round(addl));
  return ladder;
}

/**
 * Both-direction off-by-one guard for an ADDL/II ladder (AC-B4). Throws when the
 * ladder does not have exactly `1 + addl` entries starting at `tDose`.
 */
export function assertDoseLadder(ladder: ReadonlyArray<number>, tDose: number, addl: number): void {
  if (ladder.length !== 1 + addl) {
    throw new Error(
      `ADDL/II ladder length ${ladder.length} ≠ 1 + ADDL (${1 + addl}); ` +
      `off-by-one (phantom or missing dose) would corrupt the dosing schedule.`);
  }
  if (ladder[0] !== tDose) {
    throw new Error(
      `ADDL/II ladder must start at the dose time ${tDose}; got ${ladder[0]} ` +
      `(a missing-initial-dose ladder corrupts the repeat-dose anchor).`);
  }
}

/** A classified dose event for one subject (one explicit dose row pre-expansion). */
interface DoseEvent {
  readonly time: number;
  readonly amt: number;
  readonly addl: number;
  readonly ii: number;
  readonly rate: number;
  readonly cmt: number;
  readonly ss: number;
  /** True when this dose row also resets the episode clock (EVID=4). */
  readonly resets: boolean;
  /** Source-row index (for covariate hoist seed + route read). */
  readonly row: number;
}

/** A usable observation row for one subject. */
interface ObsEvent {
  readonly time: number;
  readonly dv: number;
  readonly row: number;
}

/** Read a numeric cell as a JS number, mapping platform nulls → NaN. */
function readNum(col: DG.Column, row: number): number {
  const raw = col.getRawData();
  const v = raw[row];
  if (raw instanceof Float64Array || raw instanceof Float32Array)
    return v === DGNS.FLOAT_NULL ? NaN : v;
  if (raw instanceof Int32Array)
    return v === DGNS.INT_NULL ? NaN : v;
  return typeof v === 'number' ? v : NaN;
}

/** Read an integer-valued cell (EVID/MDV/CMT/SS), NaN when null/non-numeric. */
function readInt(col: DG.Column | null, row: number): number {
  if (col === null) return NaN;
  const v = readNum(col, row);
  return Number.isFinite(v) ? Math.round(v) : NaN;
}

/**
 * Normalize a NONMEM event-layout observations DataFrame into the canonical
 * observations frame the validated pipeline consumes (AD-1). See module header.
 */
export function normalizeEventDrivenObservations(
  df: DG.DataFrame,
  ec: EventContract,
  opts: NormalizeEventOptions = {},
): NormalizeEventResult {
  const errors: ValidationError[] = [];
  const warningsBySubject = new Map<string, nca.ParameterWarning[]>();
  const doseLadderBySubject = new Map<string, number[]>();
  const inferredRoutes: InferredRoute[] = [];
  let cDose = 0; let cObs = 0; let cReset = 0; let cOther = 0; let cMissingDv = 0;

  const idCol = mustCol(df, ec.id);
  const timeCol = mustCol(df, ec.time);
  const dvCol = mustCol(df, ec.dv);
  const amtCol = mustCol(df, ec.amt);
  const evidCol = mustCol(df, ec.evid);
  const mdvCol = ec.mdv !== undefined ? df.col(ec.mdv) : null;
  const cmtCol = ec.cmt !== undefined ? df.col(ec.cmt) : null;
  const rateCol = ec.rate !== undefined ? df.col(ec.rate) : null;
  const addlCol = ec.addl !== undefined ? df.col(ec.addl) : null;
  const iiCol = ec.ii !== undefined ? df.col(ec.ii) : null;
  const ssCol = ec.ss !== undefined ? df.col(ec.ss) : null;
  const routeCol = ec.route !== undefined ? df.col(ec.route) : null;
  const dayCol = ec.day !== undefined ? df.col(ec.day) : null;

  const pushWarn = (key: string, w: nca.ParameterWarning): void => {
    let list = warningsBySubject.get(key);
    if (list === undefined) {list = []; warningsBySubject.set(key, list);}
    if (!list.some((e) => e.code === w.code)) list.push(w);
  };

  // ── Group source rows by subject id, preserving file order. ───────────────
  const bySubject = new Map<string, number[]>();
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    const sid = stringifyKey(idCol.get(i));
    let rows = bySubject.get(sid);
    if (rows === undefined) {rows = []; bySubject.set(sid, rows);}
    rows.push(i);
  }

  // Is this a multi-occasion design? (AD-8 occasion-scoped pseudo-key gate.)
  const multiOccasion = isMultiOccasion(dayCol, n);

  // ── Per-subject normalization → output row buffers. ───────────────────────
  const outSubject: string[] = [];
  const outTad: number[] = [];
  const outConc: number[] = [];
  const outDose: number[] = [];
  const outRoute: number[] = [];
  const outInfusionDur: number[] = [];
  const outOccasion: number[] = [];
  /** Observation source row per output row — covariates are hoisted from here
   *  first (NONMEM repeats covariates on every row, incl. per-obs BLQ/LLOQ),
   *  falling back to the anchor dose row for any null cell (G-6 forward-fill). */
  const outObsRow: number[] = [];
  /** Anchor dose-row per output row — covariate forward-fill fallback. */
  const outDoseRow: number[] = [];
  const seenKeys = new Map<string, string>(); // pseudo-key → source subject id

  for (const [sid, rows] of bySubject) {
    const doses: DoseEvent[] = [];
    const obs: ObsEvent[] = [];

    for (const r of rows) {
      const evid = readInt(evidCol, r);
      if (DOSE_EVIDS.has(evid)) {
        doses.push({
          time: readNum(timeCol, r),
          amt: readNum(amtCol, r),
          addl: addlCol !== null ? readInt(addlCol, r) : NaN,
          ii: iiCol !== null ? readNum(iiCol, r) : NaN,
          rate: rateCol !== null ? readNum(rateCol, r) : NaN,
          cmt: readInt(cmtCol, r),
          ss: readInt(ssCol, r),
          resets: evid === 4,
          row: r,
        });
        cDose++;
        if (evid === 4) cReset++;
      } else if (evid === 0 || Number.isNaN(evid)) {
        // Observation row (EVID=0, or absent EVID treated as obs). MDV=1 ⇒
        // planned-missing — never enters the conc stream (AD-6).
        const mdv = readInt(mdvCol, r);
        if (mdv === 1) {cMissingDv++; continue;}
        obs.push({time: readNum(timeCol, r), dv: readNum(dvCol, r), row: r});
        cObs++;
      } else if (evid === 3) {
        cReset++;
        // Pure reset (no dose) — handled as an episode boundary below.
        doses.push({time: readNum(timeCol, r), amt: NaN, addl: NaN, ii: NaN,
          rate: NaN, cmt: NaN, ss: NaN, resets: true, row: r});
      } else if (evid === 2) {
        cOther++;
        pushWarn(sid, {code: EVID_OTHER_IGNORED, severity: 'info',
          message: `EVID=2 (other event) at row ${r} logged and dropped from the ` +
            `concentration stream (no AUC effect).`});
      }
    }

    // SS handling (AD-7) — reject the unsupported branches with a blocking error
    // before building a wrong AUCτ. Honor SS=1 only when II is present.
    const ssError = checkSteadyState(sid, doses, errors, pushWarn);
    if (ssError) continue; // subject rejected; its rows are not emitted

    // Sort dose events ascending by time (the file is usually sorted, but the
    // anchor logic must not depend on it). Stable enough for the ladder.
    doses.sort((a, b) => a.time - b.time);

    // Build the per-subject anchor set: every explicit dose time ∪ each dose's
    // expanded ADDL/II ladder (AD-3). Union — explicit multi-dose rows are NOT
    // re-expanded into each other (prevents the phantom-dose blow-up, FM-2).
    const anchors: number[] = [];
    const realDoses = doses.filter((d) => Number.isFinite(d.amt));
    for (const d of realDoses) {
      const ladder = expandDoseLadder(d.time, d.addl, d.ii);
      for (const t of ladder) anchors.push(t);
    }
    // Sorted, de-duplicated anchor times (merge of all ladders, AD-3).
    const ladderTimes = Array.from(new Set(anchors)).sort((a, b) => a - b);
    doseLadderBySubject.set(sid, ladderTimes);
    // Explicit dose-row times (k=0). Distinguished from IMPLIED (expanded)
    // ladder times for the closed-interval TAD convention below.
    const explicitDoseTimes = new Set(realDoses.map((d) => d.time));

    // Episode boundaries (AD-6): each EVID∈{3,4} (reset / reset+dose) opens a new
    // profile episode. An obs's episode = #(reset times ≤ its time) — a reset
    // sorts before a co-timed obs, so the obs falls in the NEW episode. Pre/post
    // -reset AUCs are not pooled (the episode suffix splits them into profiles).
    const resetTimes = doses.filter((d) => d.resets).map((d) => d.time).sort((a, b) => a - b);
    const multiEpisode = resetTimes.length > 0;

    // Occasion-scope this subject's key (AD-8) when it is a multi-occasion
    // DESIGN — either the DAY column carries ≥2 distinct values (global), OR a
    // dose row carries an ADDL/II repeat schedule (`addl > 0`). The latter is the
    // no-DAY path (spec Slice B step 7 / AC-B9 / FM-NEW-2): a repeat-dose study
    // with no DAY column must STILL split Day-1 vs Day-N into separate profiles —
    // the occasion is derived from the dose-ladder index, so the KEY (not just
    // the OCCASION column) must carry it, else `buildProfileIndices` merges all
    // occasions of one animal into a single (wrong) dense profile. Explicit
    // multi-dose rows WITHOUT ADDL (loading+maintenance, or an EVID=4 reset) are
    // handled by episode-scoping below, not occasion-scoping — so we key on the
    // ADDL-expansion signal, not raw ladder length.
    const hasAddlExpansion = realDoses.some((d) => Number.isFinite(d.addl) && d.addl > 0);
    const occasionScoped = multiOccasion || hasAddlExpansion;
    if (multiEpisode) {
      pushWarn(sid, {code: EVID_RESET_NEW_EPISODE, severity: 'warning',
        message: `Subject '${sid}': ${resetTimes.length} reset event(s) (EVID 3/4) ` +
          `split the record into ${resetTimes.length + 1} episodes; pre/post-reset ` +
          `exposures are analysed as separate profiles.`});
    }

    // Map each anchor time → the DoseEvent that produced it (for AMT-of-record
    // and covariate-hoist reads). For an expanded ladder, every implied dose
    // inherits the originating row's AMT / covariates.
    const anchorOwner = buildAnchorOwner(realDoses);

    // Warn if the subject's first row (file order) is not a dose — covariate
    // forward-fill cannot be seeded from a dose record (G-6).
    if (realDoses.length > 0 && rows.length > 0) {
      const firstEvid = readInt(evidCol, rows[0]);
      if (!DOSE_EVIDS.has(firstEvid)) {
        pushWarn(sid, {code: NONMEM_SUBJECT_FIRST_ROW_NOT_DOSE, severity: 'warning',
          message: `Subject '${sid}' first row (EVID=${firstEvid}) is not a dose; ` +
            `covariate forward-fill uses the nearest preceding dose record.`});
      }
    }

    // Resolve route for this subject (AD-5): explicit ROUTE column wins; else
    // infer from RATE/CMT and flag for confirm. Uniform per subject (v1).
    const routeRes = resolveRoute(sid, realDoses, routeCol, opts.confirmedRoute ?? null, pushWarn);
    inferredRoutes.push({subjectKey: sid, code: routeRes.code, inferred: routeRes.inferred,
      infusionDuration: routeRes.infusionDuration});

    // Emit one output row per usable observation (anchored to its interval).
    for (const o of obs) {
      if (!Number.isFinite(o.time)) continue;
      let anchorIdx = latestAnchorIdx(ladderTimes, o.time);
      if (anchorIdx < 0) continue; // obs before any dose → cannot assign TAD; drop
      // Closed-interval [0, II] convention (AC-B3): an obs coinciding with an
      // IMPLIED (expanded) dose is the END of the previous interval (TAD = II),
      // not the start of this one (TAD = 0) — e.g. a Day-1 24 h trough lands on
      // the 2nd implied dose but belongs to interval 1. Only an EXPLICIT dose row
      // co-timed with an obs yields TAD = 0 (a genuine co-dose / predose sample,
      // e.g. the 03 t=0 first observation).
      if (ladderTimes[anchorIdx] === o.time && !explicitDoseTimes.has(o.time) && anchorIdx > 0)
        anchorIdx -= 1;
      const anchorTime = ladderTimes[anchorIdx];
      const tad = o.time - anchorTime;
      if (tad < 0) continue; // defensive (anchoring guarantees ≥ 0)

      // Occasion (U-6): DAY when present, else the 1-based ladder index.
      const occasion = dayCol !== null && !dayCol.isNone(o.row) ?
        readInt(dayCol, o.row) : anchorIdx + 1;

      const ownerRow = anchorOwner.get(anchorTime) ?? (realDoses[0]?.row ?? o.row);
      const amt = anchorOwner.size > 0 ?
        readNum(amtCol, ownerRow) : NaN;

      // Occasion-scoped pseudo-subject key (AD-8): DAY/ladder occasion + reset
      // episode suffixes split a multi-occasion / multi-episode animal into
      // distinct profiles; a single-occasion design (03, single dose) keeps the
      // bare subject id.
      const episode = multiEpisode ? countLE(resetTimes, o.time) : 0;
      const key = sid +
        (occasionScoped ? `_Day${occasion}` : '') +
        (multiEpisode ? `_ep${episode}` : '');
      registerKey(seenKeys, key, sid, pushWarn);

      outSubject.push(key);
      outTad.push(tad);
      outConc.push(o.dv);
      outDose.push(Number.isFinite(amt) ? amt : NaN);
      outRoute.push(routeRes.code);
      // Per-subject infusion duration (constant within a subject; null → NaN so
      // the extractor falls back to the scalar — Slice A2 / BD-2).
      outInfusionDur.push(routeRes.infusionDuration ?? NaN);
      outOccasion.push(occasion);
      outObsRow.push(o.row);
      outDoseRow.push(ownerRow);
    }
  }

  // ── Assemble the canonical observations frame. ────────────────────────────
  const observations = assembleFrame(
    df, ec, outSubject, outTad, outConc, outDose, outRoute, outInfusionDur,
    outOccasion, outObsRow, outDoseRow);

  const routeConfirmRequired = (opts.confirmedRoute ?? null) === null &&
    inferredRoutes.some((r) => r.inferred);

  return {
    observations,
    warningsBySubject,
    errors,
    doseLadderBySubject,
    inferredRoutes,
    routeConfirmRequired,
    classificationCounts: {dose: cDose, obs: cObs, reset: cReset, other: cOther, missingDv: cMissingDv},
  };
}

/**
 * The dose column and its unit the app's sparse composite and its dose decorators read (BUG-19
 * AD-6) — ONE expression for both mapping-pane call sites. Standard mode: the mapped dose column in
 * its declared unit. NONMEM, non-composite: {@link DOSE_COL} (AMT, the absolute mg dose of record) —
 * CL on that dose, never AMT × BW. NONMEM composite: the per-kg dose-group column (`mg/kg`,
 * unchanged), or {@link DOSE_COL} as absolute mg when no such column exists.
 */
export function sparseDoseSource(a: {
  mode: 'standard' | 'nonmem' | 'nonmem-composite';
  /** Standard mode: the mapped dose column and its declared unit. */
  doseColName: string | null;
  doseUnit: nca.DoseUnit;
  /** NONMEM composite: the per-kg dose-group column, when one was found. */
  compositeDoseGroupCol: string | null;
}): {doseColName: string | null; doseUnit: nca.DoseUnit} {
  switch (a.mode) {
  case 'standard': return {doseColName: a.doseColName, doseUnit: a.doseUnit};
  case 'nonmem': return {doseColName: DOSE_COL, doseUnit: 'mg'};
  case 'nonmem-composite':
    return a.compositeDoseGroupCol !== null ? {doseColName: a.compositeDoseGroupCol, doseUnit: 'mg/kg'} :
      {doseColName: DOSE_COL, doseUnit: 'mg'};
  }
}

// ── Helpers ────────────────────────────────────────────────────────────────

function mustCol(df: DG.DataFrame, name: string): DG.Column {
  const col = df.col(name);
  if (col === null)
    throw new Error(`normalizeEventDrivenObservations: required column '${name}' missing.`);
  return col;
}

/** Multi-occasion iff a DAY column is present with ≥ 2 distinct non-null values. */
function isMultiOccasion(dayCol: DG.Column | null, n: number): boolean {
  if (dayCol === null) return false;
  const seen = new Set<number>();
  for (let i = 0; i < n; i++) {
    if (dayCol.isNone(i)) continue;
    seen.add(readInt(dayCol, i));
    if (seen.size >= 2) return true;
  }
  return false;
}

/** Index of the latest anchor time ≤ `t` (binary search on sorted `ladder`). */
function latestAnchorIdx(ladder: ReadonlyArray<number>, t: number): number {
  let lo = 0; let hi = ladder.length - 1; let ans = -1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (ladder[mid] <= t) {ans = mid; lo = mid + 1;} else hi = mid - 1;
  }
  return ans;
}

/** Count of values in sorted `times` that are ≤ `t`. */
function countLE(times: ReadonlyArray<number>, t: number): number {
  let c = 0;
  for (const x of times) if (x <= t) c++; else break;
  return c;
}

/** Map each anchor time → the source dose-row index that owns it (for AMT +
 *  covariate reads). An expanded ladder inherits its originating row. */
function buildAnchorOwner(realDoses: ReadonlyArray<DoseEvent>): Map<number, number> {
  const owner = new Map<number, number>();
  for (const d of realDoses) {
    const ladder = expandDoseLadder(d.time, d.addl, d.ii);
    for (const t of ladder) if (!owner.has(t)) owner.set(t, d.row);
  }
  return owner;
}

/** Register a pseudo-subject key; warn on a collision between two distinct ids. */
function registerKey(
  seen: Map<string, string>, key: string, sid: string,
  pushWarn: (key: string, w: nca.ParameterWarning) => void,
): void {
  const prev = seen.get(key);
  if (prev === undefined) seen.set(key, sid);
  else if (prev !== sid) {
    pushWarn(key, {code: NONMEM_OCCASION_KEY_COLLISION, severity: 'warning',
      message: `Occasion-scoped key '${key}' is shared by source subjects ` +
        `'${prev}' and '${sid}'; verify the subject-id / DAY mapping.`});
  }
}

/** SS handling (AD-7). Returns true when the subject must be rejected. */
function checkSteadyState(
  sid: string, doses: ReadonlyArray<DoseEvent>, errors: ValidationError[],
  pushWarn: (key: string, w: nca.ParameterWarning) => void,
): boolean {
  let reject = false;
  for (const d of doses) {
    if (!Number.isFinite(d.ss) || d.ss === 0) continue;
    if (d.ss === 2) {
      errors.push({code: 'NONMEM_SS_SUPERPOSITION_UNSUPPORTED',
        message: `Subject '${sid}': SS=2 (steady state by superposition over a ` +
          `non-uniform prior schedule) is not supported — routing it to the SS=1 ` +
          `honor path would compute a wrong AUCτ. Provide an explicit dosing history.`});
      reject = true;
    } else if (d.ss === 1) {
      if (Number.isFinite(d.ii) && d.ii > 0) {
        pushWarn(sid, {code: NONMEM_STEADY_STATE_ASSUMED, severity: 'warning',
          message: `Subject '${sid}': SS=1 treated as steady-state dosing (τ=II=${d.ii}). ` +
            `NOTE: the AUCτ window is NOT mechanically enforced — all of the subject's ` +
            `observations are included (no [0,II] filter); verify the record holds one ` +
            `interval. Unvalidated by the oracle (DATA-GAP-27).`});
      } else {
        errors.push({code: 'NONMEM_SS_UNSUPPORTED',
          message: `Subject '${sid}': SS=1 without an inter-dose interval (II) is ` +
            `ambiguous — rejected rather than emit a wrong AUCτ. Supply II.`});
        reject = true;
      }
    } else {
      errors.push({code: 'NONMEM_SS_UNSUPPORTED',
        message: `Subject '${sid}': SS=${d.ss} is not a supported steady-state flag.`});
      reject = true;
    }
  }
  return reject;
}

interface RouteResolution {
  readonly code: nca.RouteCode;
  readonly inferred: boolean;
  readonly infusionDuration: number | null;
}

/**
 * Resolve a subject's route (AD-5): explicit ROUTE column > RATE/CMT inference.
 * CMT → route is NEVER hard-coded to a fixed map — inference recovers only
 * `{IV-bolus, IV-infusion, extravascular}`; PO vs SC is not recoverable, so an
 * extravascular depot maps to `ROUTE_OTHER` and the result is flagged
 * `inferred` for the mandatory confirm gate.
 */
function resolveRoute(
  sid: string, realDoses: ReadonlyArray<DoseEvent>, routeCol: DG.Column | null,
  confirmedRoute: nca.RouteCode | null,
  pushWarn: (key: string, w: nca.ParameterWarning) => void,
): RouteResolution {
  const firstDoseRow = realDoses.length > 0 ? realDoses[0].row : -1;

  // Explicit ROUTE column wins (precedence). Tolerant string normalization.
  if (routeCol !== null && firstDoseRow >= 0 && !routeCol.isNone(firstDoseRow)) {
    const rate = realDoses[0].rate;
    const infusionPositive = Number.isFinite(rate) && rate > 0;
    const res = normalizeRouteString(String(routeCol.get(firstDoseRow)), infusionPositive);
    if (res.status === 'mapped')
      return {code: res.code, inferred: false, infusionDuration: infusionDurFor(realDoses[0])};
    if (res.status === 'other')
      return {code: res.code, inferred: false, infusionDuration: null};
    // block / empty → an explicit ROUTE value exists but is unmappable; fall
    // through to RATE/CMT inference, but warn so the confirm gate does not
    // falsely claim "no ROUTE column" (the scientist must know their value was
    // unrecognised, not absent).
    pushWarn(sid, {code: NONMEM_ROUTE_INFERRED, severity: 'warning',
      message: `Subject '${sid}': ROUTE value '${routeCol.get(firstDoseRow)}' is not a ` +
        `recognised route; falling back to RATE/CMT inference. Confirm the route before Compute.`});
  }

  // RATE/CMT inference (no explicit route). RATE>0 ⇒ infusion (unvalidated).
  const d0 = realDoses.length > 0 ? realDoses[0] : null;
  const rate = d0 !== null ? d0.rate : NaN;
  const cmt = d0 !== null ? d0.cmt : NaN;
  let code: nca.RouteCode;
  let infusionDuration: number | null = null;
  if (Number.isFinite(rate) && rate > 0) {
    code = nca.ROUTE_IV_INFUSION;
    infusionDuration = infusionDurFor(d0!);
    pushWarn(sid, {code: NONMEM_INFUSION_RATE_INFERRED, severity: 'info',
      message: `Subject '${sid}': RATE>0 ⇒ IV-infusion (duration = AMT/RATE = ` +
        `${infusionDuration} h inferred); confirm the route.`});
  } else if (cmt === 1)
    code = nca.ROUTE_IV_BOLUS; // central compartment, no rate → IV bolus
  else
    code = nca.ROUTE_OTHER; // depot / unknown compartment → extravascular


  if (confirmedRoute !== null)
    return {code: confirmedRoute, inferred: true, infusionDuration};

  pushWarn(sid, {code: NONMEM_ROUTE_INFERRED, severity: 'warning',
    message: `Subject '${sid}': route inferred from RATE/CMT (no explicit ROUTE ` +
      `column). Confirm the route before Compute — a wrong IV/extravascular label ` +
      `mislabels CL vs CL/F (~3-fold error at F≈30%).`});
  return {code, inferred: true, infusionDuration};
}

function infusionDurFor(d: DoseEvent): number | null {
  if (Number.isFinite(d.rate) && d.rate > 0 && Number.isFinite(d.amt) && d.amt > 0)
    return d.amt / d.rate;
  return null;
}

/**
 * Assemble the canonical observations frame from the per-output-row buffers.
 * Hoisted covariates are read from each output row's *anchor dose row* (G-6) so
 * the constant-per-subject metadata (BW, Sex, DAY, NOAEL, route string, LLOQ,
 * BLQ, …) lands as real columns that `appendGroupingColumns` carries unchanged.
 */
function assembleFrame(
  df: DG.DataFrame, ec: EventContract,
  subject: string[], tad: number[], conc: number[], dose: number[],
  route: number[], infusionDur: number[], occasion: number[],
  obsRow: number[], doseRow: number[],
): DG.DataFrame {
  const m = subject.length;
  const cols: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, SUBJECT_COL, subject),
    DGNS.Column.fromFloat64Array(TAD_COL, Float64Array.from(tad)),
    DGNS.Column.fromFloat64Array(CONC_COL, Float64Array.from(conc)),
    DGNS.Column.fromFloat64Array(DOSE_COL, Float64Array.from(dose)),
    DGNS.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, Int32Array.from(route)),
    // Per-subject infusion duration (Slice A2 / BD-2); FLOAT_NULL where null.
    DGNS.Column.fromFloat64Array(DERIVED_INFUSION_DURATION_COL, Float64Array.from(infusionDur)),
    DGNS.Column.fromInt32Array(OCCASION_COL, Int32Array.from(occasion)),
  ];

  // Hoisted covariate columns: read at each output row's *observation* row,
  // falling back to the anchor dose row for any null cell (G-6 forward-fill).
  const covNames = collectCovariateNames(ec);
  for (const name of covNames) {
    const src = df.col(name);
    if (src === null) continue;
    cols.push(buildHoistedColumn(src, obsRow, doseRow, m));
  }

  return DGNS.DataFrame.fromColumns(cols);
}

/** Covariate columns to hoist: explicit BW/LLOQ/BLQ/DAY/route + extras, minus
 *  the event-mechanics columns the output frame does not carry. */
function collectCovariateNames(ec: EventContract): string[] {
  const names = new Set<string>();
  const add = (s?: string): void => {if (s !== undefined) names.add(s);};
  add(ec.bodyWeight);
  add(ec.lloq);
  add(ec.blqFlag);
  add(ec.day);
  add(ec.route);
  if (ec.covariates !== undefined) for (const c of ec.covariates) names.add(c);
  // Never hoist the event-mechanics columns (consumed, not carried).
  for (const mech of [ec.id, ec.time, ec.dv, ec.amt, ec.evid, ec.mdv, ec.cmt,
    ec.rate, ec.addl, ec.ii, ec.ss])
    if (mech !== undefined) names.delete(mech);
  return Array.from(names);
}

/** Build one hoisted column by reading `src` at each output row's observation
 *  row, falling back to the anchor dose row for null cells (G-6). Preserves the
 *  source column type. */
function buildHoistedColumn(
  src: DG.Column, obsRow: ReadonlyArray<number>, doseRow: ReadonlyArray<number>, m: number,
): DG.Column {
  /** Resolved source row for output row i: obs row if non-null there, else the
   *  anchor dose row (forward-fill). */
  const pick = (i: number): number => {
    const o = obsRow[i];
    if (!src.isNone(o)) return o;
    return doseRow[i];
  };
  switch (src.type) {
  case DGNS.COLUMN_TYPE.INT: {
    const out = new Int32Array(m);
    for (let i = 0; i < m; i++) {
      const r = pick(i);
      out[i] = src.isNone(r) ? DGNS.INT_NULL : (src.get(r) as number);
    }
    return DGNS.Column.fromInt32Array(src.name, out);
  }
  case DGNS.COLUMN_TYPE.FLOAT: {
    const out = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      const r = pick(i);
      out[i] = src.isNone(r) ? DGNS.FLOAT_NULL : (src.get(r) as number);
    }
    return DGNS.Column.fromFloat64Array(src.name, out);
  }
  case DGNS.COLUMN_TYPE.BOOL: {
    const out = new Array<boolean | null>(m);
    for (let i = 0; i < m; i++) {
      const r = pick(i);
      out[i] = src.isNone(r) ? null : (src.get(r) as boolean);
    }
    return DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, src.name, out);
  }
  default: {
    const out = new Array<string | null>(m);
    for (let i = 0; i < m; i++) {
      const r = pick(i);
      out[i] = src.isNone(r) ? null : src.getString(r);
    }
    return DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, src.name, out);
  }
  }
}
