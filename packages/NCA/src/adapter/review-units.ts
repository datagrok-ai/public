/**
 * Review-unit classifier (nca-workflow-stages F-S1.3, research P2; D10, D11, SF-3).
 *
 * A **review unit** is the thing a λz review decision is recorded against — the reporting unit
 * of the analysis mode:
 *   - dense / NONMEM occasion / NONMEM composite cell → the results row (`SUBJECT` string, AD-1);
 *   - UC-04 sparse destructive / batch → the GROUP COMPOSITE (`sparse:<groupKey>`); its per-animal
 *     rows are `Not applicable` whatever their measurability (SF-3, user-ratified 2026-09-25).
 *
 * Rules, in precedence order, for a unit that is its own reporting unit:
 *   (d) zero measurable concentrations → `no-exposure`, needs review, Acknowledge only
 *       (R2/F6: takes precedence over Not estimated — there is no fit to accept or exclude);
 *   (c) Not estimated (incl. STATUS `failed`) → `not-estimated`, needs review, Acknowledge only;
 *   (b) Flagged → `flagged`, needs review, Accept as is / Exclude λz;
 *   (a) Accepted → `not-needed` (including a profile an override made Accepted — D11(1)).
 * And (e): in a non-viable design every per-animal row is `not-applicable`. A per-animal row with
 * zero measurable concentrations keeps `blqOnly` + its BLQ times so the Review tooltip discloses
 * it (AC-S1.3.1a) — disclosure, not a counted review item.
 *
 * Pure: reads plain values (the caller reads the frames). The one DG-reading helper,
 * {@link profileMeasurability}, reads columns only — no `grok` / `ui`.
 */

import type * as DG from 'datagrok-api/dg';

import {ACCEPT_ACCEPTED, ACCEPT_FLAGGED, ACCEPT_NA} from './acceptance';
import {DEFAULT_RULES} from './defaults';
import type {AcceptanceFlagCode, AcceptanceReading} from './acceptance';
import {ACTIONS_BY_CATEGORY} from './review-reasons';
import type {ReviewCategory} from './review-reasons';
import type {ReviewAction, ReviewDecisionRecord, ReviewLedger} from './study-schema';

/** Unit-key prefix of a UC-04 sparse group composite (AD-1). */
export const SPARSE_UNIT_PREFIX = 'sparse:';

/** The UC-04 composite review-unit key for a sparse dose group. */
export function sparseUnitKey(groupKey: string): string {
  return `${SPARSE_UNIT_PREFIX}${groupKey}`;
}

export type ReviewNeed = 'review' | 'not-needed' | 'not-applicable';
export type UnitCategory = ReviewCategory | 'accepted';

/** One results-row profile, as the classifier needs it. */
export interface ProfileUnitInput {
  /** The results-row `SUBJECT` string (AD-1). */
  readonly key: string;
  /** The results-row `STATUS` (`ok` / `partial` / `failed`). */
  readonly status: string;
  readonly acceptance: AcceptanceReading;
  /** Measurable (positive, finite, not BLQ) concentrations of this profile. */
  readonly measurableCount: number;
  /** Canonical time strings of this profile's BLQ samples (AC-S1.3.1a disclosure). */
  readonly blqTimes: ReadonlyArray<string>;
}

/** One UC-04 sparse group composite, as the classifier needs it. */
export interface CompositeUnitInput {
  /** The sparse group key (becomes `sparse:<groupKey>`). */
  readonly groupKey: string;
  readonly acceptance: AcceptanceReading;
  /** Nominal timepoints of the composite with a measurable (positive, finite) mean. */
  readonly measurableTimepoints: number;
}

export interface ClassifyInput {
  readonly profiles: ReadonlyArray<ProfileUnitInput>;
  /** `denseProfilesViable` of the frame the results were computed on — passed in, never
   *  re-derived here. `false` = destructive / batch: per-animal rows are Not applicable. */
  readonly viable: boolean;
  /** UC-04 group composites (empty unless the sparse panel ran). */
  readonly composites: ReadonlyArray<CompositeUnitInput>;
  /** Persisted decision ledgers by unit key. */
  readonly ledgers: ReadonlyMap<string, ReviewLedger>;
  /** The CURRENT decision fingerprint by unit key. A ledger's current decision is valid only
   *  when its fingerprint equals this (D13). A unit with no entry here cannot validate a decision. */
  readonly fingerprints: ReadonlyMap<string, string>;
}

export interface ReviewUnit {
  readonly key: string;
  readonly kind: 'profile' | 'composite';
  readonly need: ReviewNeed;
  readonly category: UnitCategory;
  /** The current decision, when it is still valid (fingerprint match). */
  readonly decision: ReviewDecisionRecord | null;
  /** A current decision that no longer matches (the fit changed — D13). The controller moves it
   *  to history on the next reconcile; the classifier only reports it. */
  readonly staleDecision: ReviewDecisionRecord | null;
  /** The actions the unit allows (empty unless it needs review). */
  readonly allowedActions: ReadonlyArray<ReviewAction>;
  /** Violated-criterion codes (Flagged only) — the criterion reasons offered and the bulk-review
   *  intersection (F-S4.4). */
  readonly flagCodes: ReadonlyArray<AcceptanceFlagCode>;
  /** Not-applicable per-animal row with zero measurable concentrations (AC-S1.3.1a). */
  readonly blqOnly: boolean;
  /** Its BLQ sample times, for the disclosure tooltip. */
  readonly blqTimes: ReadonlyArray<string>;
}

/** Category of a unit that is its own reporting unit (rules (a)–(d)). */
function categorise(verdict: string, status: string, measurable: number): UnitCategory {
  if (measurable === 0) return 'no-exposure';
  if (status === 'failed' || verdict === ACCEPT_NA) return 'not-estimated';
  if (verdict === ACCEPT_FLAGGED) return 'flagged';
  // Defensive: an unknown verdict string is never silently Accepted (thesis #4).
  return verdict === ACCEPT_ACCEPTED ? 'accepted' : 'not-estimated';
}

function decisionState(
  key: string, need: ReviewNeed, input: ClassifyInput,
): {decision: ReviewDecisionRecord | null; staleDecision: ReviewDecisionRecord | null} {
  const cur = input.ledgers.get(key)?.current ?? null;
  if (cur === null || need !== 'review') return {decision: null, staleDecision: cur};
  const fp = input.fingerprints.get(key);
  return fp !== undefined && fp === cur.fingerprint ?
    {decision: cur, staleDecision: null} : {decision: null, staleDecision: cur};
}

function unitOf(
  key: string, kind: 'profile' | 'composite', category: UnitCategory,
  acceptance: AcceptanceReading, input: ClassifyInput,
): ReviewUnit {
  const need: ReviewNeed = category === 'accepted' ? 'not-needed' : 'review';
  const {decision, staleDecision} = decisionState(key, need, input);
  return {
    key, kind, need, category, decision, staleDecision,
    allowedActions: category === 'accepted' ? [] : ACTIONS_BY_CATEGORY[category],
    flagCodes: category === 'flagged' ? acceptance.flagCodes : [],
    blqOnly: false, blqTimes: [],
  };
}

/** Classify every review unit of an analysis (per-animal / per-subject rows first, in input
 *  order, then the sparse composites). */
export function classifyReviewUnits(input: ClassifyInput): ReviewUnit[] {
  const out: ReviewUnit[] = [];
  for (const p of input.profiles) {
    if (!input.viable) {
      // (e) SF-3: the group composite is the unit; the per-animal row is disclosed, not counted.
      const blqOnly = p.measurableCount === 0;
      out.push({
        key: p.key, kind: 'profile', need: 'not-applicable',
        category: categorise(p.acceptance.verdict, p.status, p.measurableCount),
        decision: null, staleDecision: input.ledgers.get(p.key)?.current ?? null,
        allowedActions: [], flagCodes: [],
        blqOnly, blqTimes: blqOnly ? p.blqTimes : [],
      });
      continue;
    }
    out.push(unitOf(p.key, 'profile',
      categorise(p.acceptance.verdict, p.status, p.measurableCount), p.acceptance, input));
  }
  for (const c of input.composites) {
    const key = sparseUnitKey(c.groupKey);
    out.push(unitOf(key, 'composite',
      categorise(c.acceptance.verdict, 'ok', c.measurableTimepoints), c.acceptance, input));
  }
  return out;
}

/** Units still `To review` (need review, no valid decision) — the "N left" count's members. */
export function pendingUnits(units: ReadonlyArray<ReviewUnit>): ReviewUnit[] {
  return units.filter((u) => u.need === 'review' && u.decision === null);
}

/** Per-profile measurability, read once from the frame the results were computed on. */
export interface Measurability {
  measurable: number;
  readonly blqTimes: string[];
}

/**
 * Count each profile's measurable concentrations and collect its BLQ sample times.
 *
 * Measurable = a finite, positive concentration that is not BLQ. BLQ = the BLQ flag is set, or
 * the concentration is below the LLOQ (column value when present, else the scalar), or it is ≤ 0.
 * A missing concentration is neither (it was not measured). Keys are the subject column's
 * `getString` — the same string the results writer puts in `SUBJECT`.
 */
export function profileMeasurability(
  df: DG.DataFrame, subjectName: string, timeName: string | null, concName: string,
  blqFlagName: string | null = null, lloqName: string | null = null, lloqScalar: number = NaN,
): Map<string, Measurability> {
  const out = new Map<string, Measurability>();
  const subj = df.col(subjectName);
  const conc = df.col(concName);
  if (subj === null || conc === null) return out;
  const time = timeName !== null ? df.col(timeName) : null;
  const blq = blqFlagName !== null ? df.col(blqFlagName) : null;
  const lloq = lloqName !== null ? df.col(lloqName) : null;
  for (let i = 0; i < df.rowCount; i++) {
    const key = subj.isNone(i) ? '' : subj.getString(i);
    let m = out.get(key);
    if (m === undefined) {m = {measurable: 0, blqTimes: []}; out.set(key, m);}
    if (conc.isNone(i)) continue;
    const c = Number(conc.get(i));
    if (!Number.isFinite(c)) continue;
    const flag = blq !== null && !blq.isNone(i) && (blq.get(i) === true || Number(blq.get(i)) === 1);
    const l = lloq !== null && !lloq.isNone(i) ? Number(lloq.get(i)) : lloqScalar;
    const below = Number.isFinite(l) && c < l;
    if (flag || below || c <= 0) {
      if (time !== null) m.blqTimes.push(time.isNone(i) ? '' : time.getString(i));
    } else
      m.measurable++;
  }
  return out;
}

/**
 * Are the per-subject profiles dense enough for the dense decorators (F%,
 * dose-proportionality, TK)? They all consume a per-subject AUCinf, which needs
 * a fittable terminal phase (≥ `minPoints` samples). A destructive / batch sparse
 * design has ~1-2 samples per subject ⇒ every per-subject NCA is degenerate
 * (AUCinf NA), making those decorators inapplicable. Returns false in that case.
 *
 * Deliberately a PER-SUBJECT degeneracy test (median samples/subject), NOT the
 * `isSparseDesign` pooling test: a *serial* design pools many animals per nominal
 * time (pooling = true) yet each animal has a full curve (dense = viable), so the
 * pooling test would wrongly suppress the decorators for serial studies.
 */
export function denseProfilesViable(df: DG.DataFrame, subjectName: string | null): boolean {
  if (subjectName === null) return true;
  const col = df.col(subjectName);
  if (col === null) return true;
  const counts = new Map<string, number>();
  for (let i = 0; i < df.rowCount; i++) {
    const k = col.isNone(i) ? '' : col.getString(i);
    counts.set(k, (counts.get(k) ?? 0) + 1);
  }
  if (counts.size === 0) return true;
  const sorted = [...counts.values()].sort((a, b) => a - b);
  const median = sorted[Math.floor(sorted.length / 2)];
  return median >= DEFAULT_RULES.lambdaZ.minPoints; // < minPoints ⇒ degenerate per-subject NCA
}

/** Intersection of the units' violated-criterion codes, in the first unit's order — the criterion
 *  reasons a bulk decision pre-checks (F-S4.4). One rule for the Review pane and the demo seed. */
export function sharedFlagCodes(units: ReadonlyArray<ReviewUnit>): AcceptanceFlagCode[] {
  if (units.length === 0) return [];
  return units[0].flagCodes.filter((c) => units.every((u) => u.flagCodes.includes(c)));
}
