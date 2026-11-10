/**
 * λz review reason vocabulary (nca-workflow-stages F-S1.6, research P1; D12).
 *
 * A review decision (Accept as is / Acknowledge / Exclude λz) needs at least one reason. The
 * picker offers three kinds, in this order:
 *   (a) **criterion reasons** — generated from the profile's OWN violated acceptance flags, in the
 *       shipped short form (`adj. R² < 0.9` …), pre-checked: the analyst confirms rather than
 *       retypes what the engine already knows (PKNCA `exclude` criterion-string form);
 *   (b) **judgement reasons** — what the engine cannot know, keyed by `(action, category)`
 *       (AC-S1.6.0), so `Acknowledge × no-exposure` and `Acknowledge × not-estimated` are
 *       disjoint and the all-BLQ reason exists only in the former (R2/F6);
 *   (c) `Other` — free text.
 *
 * Codes are STABLE persisted identifiers, decoupled from any export vocabulary (R1/F3): the
 * internal NC/ND category (NC = not calculated — insufficiency, Phoenix `Flag_N_Samples`-class;
 * ND = not determinable — criteria failed, Phoenix `Not_Accepted`-class) is carried for a future
 * CDISC export to map to PPSTAT / PPREASND / PKSUMXFL. A code, once shipped, is never renamed.
 *
 * Leaf module: no imports beyond the pure acceptance predicate.
 */

import type {AcceptanceFlagCode, AcceptanceThresholds} from './acceptance';
import {criterionReasonText} from './acceptance';
import type {ReviewAction, ReviewReason} from './study-schema';

/** The review categories a unit can be in when it needs a decision (from `classifyReviewUnits`). */
export type ReviewCategory = 'flagged' | 'not-estimated' | 'no-exposure';

/** Internal export category of a reason (R1/F3). `null` = no NC/ND semantics (Accept as is). */
export type ReasonExportCategory = 'NC' | 'ND' | null;

export interface ReasonDef {
  readonly code: string;
  /** Short picker label (user-directed 2026-09-30). Display only — never persisted. */
  readonly label: string;
  /** The full reason: the picker tooltip AND the persisted decision text (the audit record). */
  readonly text: string;
  readonly exportCategory: ReasonExportCategory;
}

/** The free-text reason's code. Its `text` is whatever the analyst typed. */
export const OTHER_REASON_CODE = 'OTHER';

/** Prefix of a criterion reason's code: `CRIT_<flagCode>`, e.g. `CRIT_ADJR2_LT`. */
export const CRITERION_REASON_PREFIX = 'CRIT_';

/** Which actions a category allows (D11): a Flagged fit can be accepted or excluded; a profile
 *  with no fit (Not estimated) or no exposure can only be acknowledged — there is no fit to
 *  accept or exclude (R2/F6). The single source for the classifier's `allowedActions`. */
export const ACTIONS_BY_CATEGORY: Readonly<Record<ReviewCategory, ReadonlyArray<ReviewAction>>> = {
  'flagged': ['accept', 'exclude'],
  'not-estimated': ['acknowledge'],
  'no-exposure': ['acknowledge'],
};

type Key = `${ReviewAction}:${ReviewCategory}`;

/**
 * Judgement reasons, keyed by `(action, category)` (AC-S1.6.0). SAP-grounded drafts (research P1:
 * PKNCA `exclude_nca.R`; MMV367 / TBI-223 / AZP-3601 SAPs; EMA BE guideline). In `no-exposure` the
 * "true finding" reason is FIRST (R2/F8): a genuine below-LLOQ result is a PK finding, and a list
 * offering only error causes would force it to be mislabelled.
 */
const JUDGEMENT_REASONS: Readonly<Partial<Record<Key, ReadonlyArray<ReasonDef>>>> = {
  'accept:flagged': [
    {code: 'ACC_VISUAL_OK', label: 'Visually adequate', exportCategory: null,
      text: 'Terminal phase adequately characterised on visual inspection'},
    {code: 'ACC_SPAN_SCHEDULE', label: 'Sampling-limited span', exportCategory: null,
      text: 'Span constrained by sampling schedule — slope consistent with other subjects in group'},
    {code: 'ACC_CRITERION_NA', label: 'Criterion not applicable', exportCategory: null,
      text: 'Criterion not meaningful for this profile (e.g. flip-flop / slow absorption)'},
  ],
  'exclude:flagged': [
    {code: 'EXC_NOT_REACHED', label: 'Terminal phase not reached', exportCategory: 'ND',
      text: 'Terminal phase not reached within sampling window'},
    {code: 'EXC_FEW_POINTS', label: 'Too few points', exportCategory: 'NC',
      text: 'Insufficient post-Cmax points for a reliable slope'},
    {code: 'EXC_ERRATIC', label: 'Erratic terminal phase', exportCategory: 'ND',
      text: 'Non-monotonic / erratic terminal concentrations'},
    {code: 'EXC_NEAR_LLOQ', label: 'Near LLOQ', exportCategory: 'ND',
      text: 'Terminal concentrations near or below LLOQ'},
  ],
  'acknowledge:not-estimated': [
    {code: 'ACK_FEW_POINTS', label: 'Fewer than 3 points', exportCategory: 'NC',
      text: 'Fewer than 3 measurable post-Cmax points'},
    {code: 'ACK_NO_DECLINE', label: 'No terminal decline', exportCategory: 'ND',
      text: 'No declining terminal phase'},
    {code: 'ACK_COMPUTE_FAILED', label: 'Computation failed', exportCategory: 'NC',
      text: 'Computation failed — data issue recorded'},
  ],
  'acknowledge:no-exposure': [
    {code: 'NOEXP_TRUE_FINDING', label: 'True finding (BLQ)', exportCategory: 'NC',
      text: 'Dosed as planned — concentrations below LLOQ throughout (true finding)'},
    {code: 'NOEXP_DOSING_FAILURE', label: 'Dosing failure', exportCategory: 'NC',
      text: 'Possible dosing failure'},
    {code: 'NOEXP_ASSAY', label: 'Sample / assay issue', exportCategory: 'NC',
      text: 'Sample handling / assay issue'},
    {code: 'NOEXP_EXPECTED', label: 'Expected (control)', exportCategory: 'NC',
      text: 'Expected (pre-dose / vehicle / control)'},
  ],
};

/** The judgement reasons for one `(action, category)`; empty when the pair is not allowed. */
export function judgementReasons(action: ReviewAction, category: ReviewCategory): ReadonlyArray<ReasonDef> {
  return JUDGEMENT_REASONS[`${action}:${category}`] ?? [];
}

/** Criterion reasons for a profile's violated flags, in the shipped short form. Offered for the
 *  two Flagged actions only — a Not-estimated / No-exposure unit has no evaluated criterion. */
export function criterionReasons(
  flagCodes: ReadonlyArray<AcceptanceFlagCode>, t: AcceptanceThresholds,
): ReadonlyArray<ReasonDef> {
  return flagCodes.map((c) => {
    const text = criterionReasonText(c, t); // already short — the label IS the text
    return {code: `${CRITERION_REASON_PREFIX}${c}`, label: text, text, exportCategory: 'ND' as const};
  });
}

/** Display form of a decision's reasons (AC-S1.6.1): joined with `"; "` (PKNCA convention). */
export function formatReasons(reasons: ReadonlyArray<ReviewReason>): string {
  return reasons.map((r) => r.text).join('; ');
}

/** D12 / AC-S1.6.2: a decision can be saved only with ≥ 1 reason whose text is non-blank
 *  (an `Other` with nothing typed does not count). */
export function reasonsValid(reasons: ReadonlyArray<ReviewReason>): boolean {
  return reasons.some((r) => r.text.trim().length > 0);
}

/** The Review STATE a current decision puts a unit in — the Review-column value and the Review
 *  pane title (D15). ONE map for both surfaces (review: was defined twice). */
export const REVIEW_STATE_LABEL: Readonly<Record<ReviewAction, string>> = {
  accept: '✓ Accepted*', acknowledge: 'Acknowledged', exclude: 'Excluded',
};

/** Human label of an action — the button / tooltip wording (D11/D15). */
export const REVIEW_ACTION_LABEL: Readonly<Record<ReviewAction, string>> = {
  accept: 'Accept as is',
  acknowledge: 'Acknowledge',
  exclude: 'Exclude λz',
};
