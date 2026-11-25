/**
 * The Hub demo studies (nca-hub-demo-studies PD-8): five simulated studies of one discovery program
 * (Series X, lead DGK-0007), seeded once per server by `ncaSeedDemoStudies` and shared View to
 * `All users`. Pure data — read by the seed (`pipeline/demo-seed.ts`) and by the tests, so the
 * tests assert against the SAME constants the seed applies.
 *
 * **All data are simulated.** The SMILES are real molecules carrying fictitious PK; NOAEL, the human
 * dose and the exposure margin are design constants, not a determined NOAEL. Generator, seeds, model
 * and the planted matrix: `files/REGEN.md`.
 *
 * **`REGRESSION_PIN`** marks a value measured once on the shipped engine and pinned so the demo stays
 * stable (R1/F-7) — it proves the demo did not drift, not that the verdict is scientifically right.
 */

// Type-only: erased at build, so this module stays a dependency-free leaf.
import type {DerivedChoices} from './occasion';

/** The structure column of every demo table (SMILES). Declared, not detected: the seed types it
 *  `Molecule` itself, so the demo installs without the Chem package (Chem then draws it). */
export const DEMO_STRUCTURE_COL = 'SMILES';

/** Where the seed reads its tables (served from `files/demo-studies/`). */
export const DEMO_STUDIES_DIR = 'demo-studies';

/**
 * The acceptance thresholds the planted windows were designed for (R1/F-5), written explicitly into
 * every demo recipe so a later change of the app defaults neither moves a demo verdict nor stales the
 * planted matrix. The span threshold is the code constant `SPAN_RATIO_POOR` — the tests assert it is
 * still 2.
 */
export const DEMO_THRESHOLDS = {
  minRSquared: 0.85,
  acceptanceRSquared: 0.90,
  minPoints: 3,
  extrapWarnPct: 20,
  spanRatio: 2,
} as const;

/** The rules JSON of every demo recipe (`parseRules` input): the pinned thresholds, nothing else. */
export const DEMO_RULES_JSON = JSON.stringify({
  lambdaZ: {minRSquared: DEMO_THRESHOLDS.minRSquared, minPoints: DEMO_THRESHOLDS.minPoints},
  extrapWarnPct: DEMO_THRESHOLDS.extrapWarnPct,
  acceptanceRSquared: DEMO_THRESHOLDS.acceptanceRSquared,
});

/** A pre-recorded λz point override (CF-2: the edit is its own justification; the reason is in the
 *  source `Comment` of the overridden rows). */
export interface DemoOverride {
  readonly subject: string;
  /** Observation times (h) of the points. */
  readonly times: ReadonlyArray<number>;
  readonly action: 'force-in' | 'force-out';
}

/** A pre-recorded review decision on named units (one batch: same reasons, same time — what the UI
 *  bulk action persists). Criterion reasons are added for a Flagged unit exactly as the Review pane
 *  pre-checks them (the units' shared flag codes). */
export interface DemoDecision {
  readonly ref: string;
  /** Results-row `SUBJECT` keys; `sparse:*` for a composite. Empty + `compositeVerdict` = every
   *  composite of the analysis. */
  readonly units: ReadonlyArray<string>;
  readonly action: 'accept' | 'exclude' | 'acknowledge';
  /** Registered judgement reason codes (`review-reasons.ts`). */
  readonly codes: ReadonlyArray<string>;
  /** The free-text `Other` reason. */
  readonly other: string;
}

/**
 * Every unplanned profile that the engine flags ONLY by `N points ≤ 3` gets one pre-recorded bulk
 * Accept as is per analysis (user decision at build, 2026-11-04): measured on the shipped engine,
 * the λz best-fit's 1e-4 adj. R² tie-break keeps a 3-point window on 30–41 % of ordinary profiles
 * (44 of 108 on the app's own rat / batch fixtures). The seed asserts each such fit is good
 * (adj. R² ≥ {@link N_ONLY_MIN_ADJ_R2}) before it records the decision.
 */
export const N_ONLY_DECISION = {
  action: 'accept' as const,
  codes: ['ACC_VISUAL_OK'],
  other: '3-point terminal phase (the auto-fit kept the last three points); visually log-linear, ' +
    'adj. R² ≥ 0.98 — accepted as is.',
};
export const N_ONLY_MIN_ADJ_R2 = 0.98;

export interface DemoAnalysisSpec {
  /** Analysis label (= the compound id; Study 2 adds ` · F` / ` · M`). */
  readonly label: string;
  /** Rows of the study table this analysis holds: a column = value filter, or all rows. */
  readonly rows?: {readonly col: string; readonly value: string};
  readonly overrides?: ReadonlyArray<DemoOverride>;
  readonly decisions?: ReadonlyArray<DemoDecision>;
  /** Profiles left `To review` (C7, C8) — never bulk-decided. */
  readonly leaveToReview?: ReadonlyArray<string>;
  /** Expected `To review` count after seeding (AC-2.4). */
  readonly expectedPending: number;
}

export interface DemoStudySpec {
  readonly ref: 1 | 2 | 3 | 4 | 5;
  readonly label: string;
  /** File under {@link DEMO_STUDIES_DIR}. */
  readonly file: string;
  /** The source frame name every analysis carries (the table name the view and exports show). */
  readonly frameName: string;
  /** UC-04 destructive design. */
  readonly sparse: boolean;
  /** Source columns to Group by (as the Group-by ribbon names them; asserted against its candidates). */
  readonly groupBy: ReadonlyArray<string>;
  /** Annotation (passthrough) columns. */
  readonly passthrough: ReadonlyArray<string>;
  /** Study 5 (PD-5): the infusion-duration column the seed computes with, and the uniform duration
   *  (h) the recipe stores as a Value. The recipe is NOT column-mode: the analysis view ignores a
   *  per-row infusion-duration column (a known gap — MRT / Vss would reopen
   *  without −T_inf/2), so a column recipe would reopen to different numbers. Every row of the
   *  column must equal the value (the seed asserts it), so both compute the same infusion. */
  readonly infusion?: {readonly col: string; readonly hours: number};
  /** The derived-analysis choices the recipe stores (tk-study-day-stratification): Study 4
   *  declares τ = 24 h (once daily, time after dose) so its dose proportionality and
   *  accumulation report AUC0–τ rather than the labelled AUClast fallback. */
  readonly derived?: DerivedChoices;
  readonly analyses: ReadonlyArray<DemoAnalysisSpec>;
}

const SERIES = ['DGK-0001', 'DGK-0002', 'DGK-0003', 'DGK-0004', 'DGK-0005', 'DGK-0006', 'DGK-0007', 'DGK-0008'];
const PASSTHROUGH = ['SMILES', 'CompoundID'];

/** C1–C3 (research §C/§D with CF-1): Study 1 is shown fully reviewed. */
function study1Analysis(cid: string): DemoAnalysisSpec {
  const decisions: DemoDecision[] = [];
  if (cid === 'DGK-0008') {
    decisions.push({ref: 'C1', units: ['R801', 'R802', 'R803', 'R804', 'R805', 'R806'], action: 'accept',
      codes: ['ACC_SPAN_SCHEDULE'],
      other: 'Sampling covers < 2 half-lives; t½ and AUC∞ indicative — extend sampling to 48–72 h in the next study.'});
  }
  if (cid === 'DGK-0001') {
    decisions.push({ref: 'C2', units: ['R104', 'R105', 'R106'], action: 'accept', codes: [],
      other: 'High-clearance compound below LLOQ from 6 h; 3-point terminal phase, t½ indicative only.'});
  }
  if (cid === 'DGK-0006') {
    decisions.push({ref: 'C3', units: ['R604'], action: 'acknowledge', codes: [],
      other: 'Double peak (enterohepatic recirculation or dual absorption site); no monoexponential terminal phase.'});
  }
  return {label: cid, rows: {col: 'CompoundID', value: cid}, decisions, expectedPending: 0};
}

export const DEMO_STUDIES: ReadonlyArray<DemoStudySpec> = [
  {
    ref: 1, label: 'Demo · rat IV/PO screen', file: 'demo_study1_series_x_rat_iv_po.csv',
    frameName: 'Rat IV/PO screen — Series X (simulated)', sparse: false,
    groupBy: ['ROUTE'], passthrough: PASSTHROUGH, analyses: SERIES.map(study1Analysis),
  },
  {
    ref: 2, label: 'Demo · mouse sparse PO', file: 'demo_study2_dgk0007_mouse_sparse_po.csv',
    frameName: 'Mouse sparse PO — DGK-0007 (simulated)', sparse: true,
    groupBy: [], passthrough: PASSTHROUGH,
    // C11: decided from the pinned composite verdict (`DEMO_SPARSE_VERDICTS`).
    analyses: [
      {label: 'DGK-0007 · F', rows: {col: 'Sex', value: 'F'}, expectedPending: 0},
      {label: 'DGK-0007 · M', rows: {col: 'Sex', value: 'M'}, expectedPending: 0},
    ],
  },
  {
    ref: 3, label: 'Demo · rat dose escalation', file: 'demo_study3_dgk0007_rat_dose_escalation.csv',
    frameName: 'Rat dose escalation — DGK-0007 (simulated)', sparse: false,
    groupBy: ['Dose_mg_kg', 'ROUTE'], passthrough: PASSTHROUGH,
    analyses: [{
      label: 'DGK-0007', expectedPending: 0,
      decisions: [{ref: 'C5', units: ['R316'], action: 'exclude', codes: ['EXC_ERRATIC'],
        other: 'Erratic terminal phase at a dose with solubility-limited absorption; λz-dependent ' +
          'parameters not reported for this animal.'}],
    }],
  },
  {
    ref: 4, label: 'Demo · 28-day TK (rat)', file: 'demo_study4_dgk0007_rat_tk28d.csv',
    frameName: '28-day TK — DGK-0007 (simulated)', sparse: false,
    groupBy: ['Dose_mg_kg', 'Sex', 'StudyDay'], passthrough: PASSTHROUGH, derived: {tauH: 24},
    analyses: [{
      label: 'DGK-0007', expectedPending: 3,
      // C6: a lab-flagged 8 h result, out of the λz regression only (PKNCA `exclude_half.life`).
      overrides: [{subject: 'T010_D28', times: [8], action: 'force-out'}],
      // C4′: the vehicle-control TK animals (PD-6), confirmatory sampling, all BLQ.
      decisions: [{ref: 'C4′',
        units: ['T001', 'T002', 'T003', 'T004', 'T005', 'T006'].flatMap((a) => [`${a}_D1`, `${a}_D28`]),
        action: 'acknowledge', codes: ['NOEXP_EXPECTED'], other: ''}],
      // C7 (delayed absorption, 100 mg/kg Day 1: Tmax 4 h vs 0.5–2 h, λz on 8–24 h) and C8 (noisy
      // terminal phase): the user reviews them.
      leaveToReview: ['T019_D1', 'T022_D1', 'T016_D28'],
    }],
  },
  {
    ref: 5, label: 'Demo · IV infusion (rat)', file: 'demo_study5_dgk0007_rat_iv_infusion.csv',
    frameName: 'Rat IV infusion — DGK-0007 (simulated)', sparse: false,
    groupBy: ['Dose_mg_kg', 'Sex'], passthrough: PASSTHROUGH, infusion: {col: 'InfusionDuration_hr', hours: 0.5},
    analyses: [{
      label: 'DGK-0007', expectedPending: 0,
      // C9 (replaced at build): the auto-fit kept the last three points of a terminal phase that is
      // log-linear from 4 h; the analyst extends λz to 4–24 h.
      overrides: [{subject: 'R510', times: [4, 6], action: 'force-in'}],
    }],
  },
];

/**
 * C11 — the composite verdicts of Study 2, measured once on the shipped engine and pinned
 * (`REGRESSION_PIN`, architect #7): the seed applies the decision for this verdict and asserts the
 * measured verdict still equals it. Keyed by analysis label.
 */
const C11_DECISION: DemoDecision = {
  ref: 'C11', units: [], action: 'accept', codes: ['ACC_VISUAL_OK'],
  other: 'Destructive design, n = 3 per time point: the composite terminal phase has 3 points (2–8 h; ' +
    '24 h below LLOQ in every animal) — t½ indicative.',
};
/** Measured 2026-11-04 on sci-comp 0.11.0: both composites Flagged, `N points (λz) ≤ 3` only
 *  (measurable composite points 0.25–8 h, 24 h BLQ in every mouse). */
export const DEMO_SPARSE_VERDICTS: Readonly<Record<string, {verdict: string; decision: DemoDecision | null}>> = {
  'DGK-0007 · F': {verdict: 'Flagged', decision: C11_DECISION},
  'DGK-0007 · M': {verdict: 'Flagged', decision: C11_DECISION},
};

/**
 * The Study 1 order of F (research §B: 0007 45 % > 0004 43 % > 0003 36 % > 0002 27 % > 0005 23 % >
 * 0001 14 % > 0006 5 %), highest first. AC-2.5 asserts it only for ADJACENT compounds whose F differ
 * by more than 2 × SE at n = 3 (R1/F-7); DGK-0008 is left out (its F rests on an extrapolated AUC∞,
 * R2/N-5a). The design order, not a measured pin — so the check is not circular (architect #7).
 */
export const DEMO_F_ORDER: ReadonlyArray<string> =
  ['DGK-0007', 'DGK-0004', 'DGK-0003', 'DGK-0002', 'DGK-0005', 'DGK-0001', 'DGK-0006'];
