/*
 * The Analyze engine: a study's long table (animal, group, day, volume) analysed as Plan pre-specifies it, so the
 * test that decides is the one the design powered: each animal's least-squares slope of log V from randomisation
 * (day 0) to the evaluation day, volumes below the smallest a caliper records raised to it (Crux stats.tgi's
 * per-animal rate, rate_floor), every treated arm against the control on those slopes with their own variances
 * (Welch's t test, Crux stats.posthoc) and the declared correction, Bonferroni by default; beside it, what one variance
 * for every group gives (Crux stats.contrast, rate ~ group, treatment vs control: Dunnett exact). TGI and
 * T/C on the evaluation day (the NCI definitions Plan declares effects in; TGI on each animal's change from its own
 * day 0, which gives its interval), the time to the endpoint volume and
 * mRECIST come from the same stats.tgi call. The digital twin of each treated animal (Crux stats.tumor_twin, the
 * line's model from the engine priors) gives what it would measure untreated, its paired effect on the evaluation
 * day, and each control left out against the others. With a body weight column, Crux stats.tumor_body_weight reads
 * whether the treatments were tolerated (each group's weight change, measured and without the tumor's mass, its
 * weight-loss grades and the animals removed at a humane endpoint); the tumor's mass from the line's density in the
 * atlas (files/tumour_density_v0.json), calibrated by the tumors' ex vivo weights when the study has them. This file
 * only translates and combines.
 */
import * as DG from 'datagrok-api/dg';

import {controlByName} from './arm-name';
import {crux, CruxResult, EmptyState, MAX_SETS, named, NoControl, rows} from './crux';
import {groupsByAnimal} from './animal-group';
import {BH_SCREENING, CORRECTION_NAME, Correction, EffectMetric, ENDPOINT, HUMANE, LAB_WEEKDAYS, StudyDeclaration,
  StudyDesign, correctionOf,
  designStudy,
  finite, isFiniteNumber as ok, median} from './engine';
import {PlanInputs, declarationOf} from './frames';
import {groupColors} from './colors';
import {gpu} from './gpu';
import {BandFactors, LineDensity, bandFactors, linePrior, modelName, tumourDensity} from './priors';
import {PLAN_TAG, SavedPlan, parsePlan, powerAt} from './plan-file';
import {registerHistoryOf} from './study-file';
import {ChecksInput, PlanChecks, Range, planChecks, widenedEnd, withPrinted} from './plan-checks';
import {mean, maxOf, minOf, minus} from './format';
import {pValue} from './ui';
import {mm3Text} from './reading-check';

/** `volume`: NaN where the row has none; `weight`: the body weight (g) of each row, null when the study has none;
 *  `mass`: the tumor's ex vivo weight (g) on the row of the day it was excised, null when the study has none;
 *  `skipped`: the table's rows left out for lacking their animal, group or day */
export interface StudyData {animal: string[], group: string[], day: number[], volume: number[],
  weight: (number | null)[] | null, mass: (number | null)[] | null, removal: string[] | null, skipped: number}

export interface AnalysisDeclaration {
  data: StudyData;
  /** the cell line: its model of control growth for the twin; '': none picked (the study's and its plan's: no
   *  twin, no atlas checks, the rest the same; never a model assumed, QA 0.1.3 A2) */
  line: string;
  /** the control (vehicle) group */
  control: string;
  /** mm3: the evaluation day is the first measurement day the median control reaches it */
  endpoint: number;
  /** the evaluation day when given (else from the endpoint) */
  day: number | null;
  alpha: number;
  correction: Correction;
  /** the plan the study ran under (null: none linked; the controls and the end are then read against the atlas) */
  plan: SavedPlan | null;
  /** humane endpoints: the weight loss (%) and the tumor volume (mm3) at which an animal is removed */
  weightLimit: number;
  humane: number;
  /** mg per mm³: the tumor's density for the tumor-corrected weight (null: the cell line's, from the atlas) */
  density: number | null;
  /** false: no digital twin (the longest call): twins, effects and checks empty, the rest the same */
  twin?: boolean;
}

/** pmxAnalyzeStudy's inputs as the form and Grokky give them: columns, an optional plan file, empty optionals. */
export interface AnalysisInputs {
  study: DG.DataFrame; animal: DG.Column; group: DG.Column; day: DG.Column; volume: DG.Column;
  /** the text of a .pmxplan file */
  plan?: string | null; line: string | null; control: string; endpoint: number; correction: string;
  evaluationDay?: number; alpha?: number; bodyWeight?: DG.Column | null; weightLimit?: number; humane?: number;
  tumorWeight?: DG.Column | null; density?: number; removalReason?: DG.Column | null; withTwin?: boolean;
}

/** A column's values, an empty cell as null. */
const valuesOf = (c: DG.Column): any[] => Array.from({length: c.length}, (_, i) => c.isNone(i) ? null : c.get(i));

/** The declaration the engine reads: the rows with their animal, group and day, the control found when not named,
 *  the plan from its file or else the study table's tag, each empty optional at its default. */
export function analysisDeclarationOf(x: AnalysisInputs): AnalysisDeclaration {
  // a model not in the atlas: its history, carried by the study or its plan, registered before the twin needs it
  registerHistoryOf(x.study);
  if (x.plan?.trim()) {
    try {
      parsePlan(x.plan);
    } catch {
      // read again below, where a damaged plan says so
    }
  }
  const [animal, day, volume] = [x.animal, x.day, x.volume].map(valuesOf);
  const group = groupsByAnimal(animal.map((a) => a === null ? null : String(a)), valuesOf(x.group));
  const keep = animal.map((a, i) => a !== null && group[i] !== null && ok(day[i]));
  const pick = <T>(xs: T[]) => xs.filter((_, i) => keep[i]);
  const groups = pick(group).map(String);
  const names = [...new Set(groups)];
  // never a treated group taken for the control: with none named so, the user picks it (no verdict until then)
  const control = x.control?.trim() || controlByName(names);
  if (!control)
    throw new NoControl('', NO_CONTROL);
  const planText = x.plan?.trim() || x.study.getTag(PLAN_TAG);
  const plan = planText ? parsePlan(planText) : null;
  // the model picked, else the plan's; none: no model assumed
  return {line: x.line || plan?.predicted.line || '', control, endpoint: x.endpoint,
    day: ok(x.evaluationDay) ? x.evaluationDay : null, plan, alpha: ok(x.alpha) ? x.alpha : 0.05,
    correction: correctionOf(x.correction),
    weightLimit: ok(x.weightLimit) ? x.weightLimit : 20, humane: ok(x.humane) ? x.humane : HUMANE,
    density: ok(x.density) && x.density > 0 ? x.density : null, twin: x.withTwin !== false,
    data: {animal: pick(animal).map(String), group: groups, day: pick(day),
      volume: pick(volume).map((v) => ok(v) ? v : NaN),
      weight: x.bodyWeight ? pick(valuesOf(x.bodyWeight)) : null,
      mass: x.tumorWeight ? pick(valuesOf(x.tumorWeight)) : null,
      removal: x.removalReason ? pick(valuesOf(x.removalReason)).map((r) => r === null ? '' : String(r)) : null,
      skipped: keep.filter((k) => !k).length}};
}

/** The smallest volume a caliper records (mm3): the design's v_lod (stats.tumor_study_design's default). */
const V_LOD = 3.4;
export const COMPARED = 'per-animal slope of log V from randomization (Welch t test)';

type Contrast = {term: string, estimate: number, lower: number, upper: number, p_value: number, level?: number,
  p_raw?: number};

/** stats.posthoc's p_adjust for each method (statsmodels' multipletests names: Crux computes BH too). */
const P_ADJUST: {[c in Correction]: string} = {bonferroni: 'bonferroni', holm: 'holm', bh: 'fdr_bh', none: 'none'};

/** Each arm's interval level, compatible with its adjusted p: its interval excludes 0 exactly when its adjusted p is
 *  below α, so p, interval and verdict agree. `raw`/`adjusted`: the arms' p-values. Bonferroni: 1 − α/k; none (or
 *  one arm): 1 − α; BH: the FCR-adjusted level of Benjamini & Yekutieli (2005), 1 − R·α/k with R arms declared
 *  (R = 0: 1 − α/k); Holm: raw p ascending, step i at 1 − α/(k − i + 1) while the steps pass; from the first that
 *  fails, every later arm at that failed step's level (at its own, a later arm's interval could exclude 0 while
 *  Holm does not declare it). */
export function armLevels(raw: number[], adjusted: number[], alpha: number, correction: Correction): number[] {
  const k = raw.length;
  if (correction === 'none' || k < 2)
    return raw.map(() => 1 - alpha);
  const declared = adjusted.filter((p) => p < alpha).length;
  if (correction !== 'holm')
    return raw.map(() => 1 - (correction === 'bh' ? Math.max(declared, 1) : 1) * alpha / k);
  const levels = new Array<number>(k);
  raw.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0])
    .forEach(([, arm], j) => levels[arm] = 1 - alpha / (k - Math.min(j, declared)));
  return levels;
}

/** Each treated group (g1, g2, ...) against the control (g0) with their own variances: Welch's t test, its p-values
 *  raw and adjusted over the comparisons (Crux stats.posthoc), and each difference's interval at its arm's level
 *  (armLevels: the one compatible with the adjusted p). A group with fewer than two slopes has no row. */
async function welchContrasts(x: {slope: Float64Array, group: string[]}, correction: Correction,
  alpha: number): Promise<Contrast[]> {
  const by = new Map<string, number[]>();
  x.group.forEach((g, i) => (by.get(g) ?? by.set(g, []).get(g)!).push(x.slope[i]));
  const of = (g: string) => {
    const v = by.get(g) ?? [];
    const m = mean(v);
    return {n: v.length, mean: m, sem2: v.reduce((t, y) => t + (y - m) ** 2, 0) / (v.length - 1) / v.length};
  };
  const arms = [...by.keys()].filter((g) => g !== 'g0' && by.get(g)!.length >= MIN_TESTED);
  if ((by.get('g0')?.length ?? 0) < MIN_TESTED || arms.length === 0)
    return [];
  const keep = x.group.map((g) => g === 'g0' || arms.includes(g));
  const posthoc = (pAdjust: string) => crux('stats.posthoc', {x: x.slope.filter((_, i) => keep[i]),
    group: x.group.filter((_, i) => keep[i])}, {test: 'ttest', equal_var: false, control: 'g0',
    p_adjust: pAdjust, alpha}).then((r) => rows(r.main as any));
  const [tested, untouched] = await Promise.all([posthoc(P_ADJUST[correction]), posthoc('none')]);
  const pOf = (t: any[], g: string) => Number(t.find((r) => r.group_a === g)?.pvalue ?? NaN);
  const [adjusted, raw] = [arms.map((g) => pOf(tested, g)), arms.map((g) => pOf(untouched, g))];
  const levels = armLevels(raw, adjusted, alpha, correction);
  const c = of('g0');
  return Promise.all(arms.map(async (g, j) => {
    const level = 1 - levels[j];
    const a = of(g);
    const [estimate, se2] = [a.mean - c.mean, a.sem2 + c.sem2];
    // Satterthwaite's degrees of freedom, as the test's
    const df = se2 * se2 / (a.sem2 * a.sem2 / (a.n - 1) + c.sem2 * c.sem2 / (c.n - 1));
    const t = se2 > 0 ? Number((await crux('stats.dist_eval', {x: Float64Array.from([1 - level / 2])},
      {dist: 't', func: 'ppf', shape1: df})).main.col['value'][0]) : 0;
    const half = t * Math.sqrt(se2);
    return {term: `${g} - g0`, estimate, lower: estimate - half, upper: estimate + half, level: levels[j],
      p_value: adjusted[j], p_raw: raw[j]};
  }));
}

export interface ArmAnalysis {
  arm: string;
  /** treated animals with a slope (the test's) */
  n: number;
  /** treated animals measured on the evaluation day (TGI's and T/C's) */
  nDay: number;
  /** percent, on the evaluation day */
  tgi: number | null;
  tc: number | null;
  /** mean slope of ln V per day, the arm's and the control's, and their difference with its (adjusted) interval at
   *  `confidence` (its method's level, armLevels): the one interval the result is read on */
  slope: number;
  slopeControl: number;
  difference: number;
  lower: number;
  upper: number;
  confidence: number;
  /** adjusted for the comparisons with the control (pRaw: before the adjustment) */
  p: number;
  pRaw: number;
  /** the sensitivity analysis: what one variance for every group gives (Dunnett's p-value; with one treated group,
   *  the two-sample t test's); null: the group has none */
  pCommon: number | null;
  /** Kaplan-Meier median time to the endpoint volume (null: not reached) */
  medianTte: number | null;
  /** how many days later than the control its median reaches the endpoint volume (T − C; null: either never does) */
  delay: number | null;
  /** the twin: animals measured on the evaluation day, the median of their paired effects (percent) and how many
   *  are smaller than their twin beyond doubt (the 90 % interval of the effect above 0) */
  twinAnimals: number;
  pairedEffect: number | null;
  responders: number;
}

/** What the test on the slopes says of a treated group: it grows slower than the control (its interval below 0),
 *  faster (above 0), or the study cannot tell it from the control. */
export type Activity = 'Slows growth' | 'Grows faster' | 'No clear effect' | typeof TOO_FEW;
/** A group the test cannot be run on (fewer than MIN_TESTED animals with a slope, in it or in the control): no p, no
 *  verdict (QA 0.1.3 round 3, chaos 4: such a group read "grows 0.0%/day faster, p < 0.001"). */
export const TOO_FEW = 'Too few animals to test';
/** The animals with a slope a group (and the control) needs to be tested: the one threshold the note, the cells and
 *  the test share (round 4, analyst A4: a group "too few for a test" was tested on 2). */
export const MIN_TESTED = 3;
/** Why a group is not tested, in a cell's tooltip and the notes: "1 animal; at least 3 needed". */
export const tooFewWhy = (n: number) => `${n} animal${n === 1 ? '' : 's'} with a slope; at least ${MIN_TESTED} needed`;
export const activityOf = (a: Pick<ArmAnalysis, 'lower' | 'upper' | 'p'>): Activity =>
  ![a.lower, a.upper, a.p].every(Number.isFinite) ? TOO_FEW :
    a.upper < 0 ? 'Slows growth' : a.lower > 0 ? 'Grows faster' : 'No clear effect';

/** A slope difference (ln V per day) as the change of the daily growth in %/day, signed: −7.1 grows 7.1% less each
 *  day than the control. */
export const pctPerDay = (d: number) => 100 * (Math.exp(d) - 1);

/** A number to one decimal with a true minus ("−7.1"; never "−0.0"). */
export const oneDecimal = (x: number) => {
  const s = x.toFixed(1);
  return /^-0\.0$/.test(s) ? '0.0' : minus(s);
};

/** The groups table's growth cell: the change of the daily growth against the control in %/day and its interval,
 *  low to high: "−7.1 (−9.6 to −4.5)"; null: no test. */
export function growthVsControl(d: number, lower: number, upper: number): string | null {
  if (![d, lower, upper].every(Number.isFinite))
    return null;
  return `${oneDecimal(pctPerDay(d))} (${oneDecimal(pctPerDay(lower))} to ${oneDecimal(pctPerDay(upper))})`;
}

/** The interval's level and method, for its header (round 3, analyst W2), and the p column's name: the one place
 *  Analyze, Compare, the CTD, the reports and the files take them from. "99% CI, Bonferroni over 5 groups" and
 *  "p (Bonferroni)"; Holm's level varies by arm: "CI by Holm step" (each arm's own on its cell's hover). Each
 *  method's interval agrees with its p (armLevels); one switch of `correction` changes both. `confidence`: the
 *  study's level (with Holm, its first step's). */
export function intervalWords(confidence: number, correction: unknown,
  k: number): {level: string, header: string, tip: string, p: string} {
  const level = levelText(confidence);
  if (correction === 'none' || k < 2) {
    return {level, header: level, p: 'p', tip: `${level}: ${k < 2 ? 'one comparison with the control' :
      'not adjusted for the comparisons (exploratory)'}.`};
  }
  const c = correctionOf(correction);
  const p = `p (${c === 'bh' ? 'BH-adjusted' : CORRECTION_NAME[c]})`;
  const agree = `the interval and the ${p} agree, and the result is read on them.`;
  if (c === 'holm') {
    const header = `CI by Holm step over ${k} groups`;
    return {level, header, p, tip: `${header} with the control: the ` +
      `groups in order of their unadjusted p, the first at 1 − α/${k}, the next at 1 − α/${k - 1}, and so on; from ` +
      `the first step not passed, the rest at that step's level. Each group's level is on its cell's hover; ${agree}`};
  }
  const header = `${level}, ${c === 'bh' ? 'FDR-adjusted (BH)' : 'Bonferroni'} over ${k} groups`;
  return {level, header, p, tip: `${header} with the control (` + (c === 'bh' ? `Benjamini and Yekutieli: each ` +
    `at 1 − R·α/${k}, R the groups declared); ${agree} Benjamini-Hochberg ${BH_SCREENING}.` :
    `each at 1 − α/${k}); ${agree}`)};
}

/** The comparisons an analysis made, as its correction counts them and the selector's hint says them: the groups
 *  tested, not every treated group (one with too few animals, or none measured past day 0, is not tested;
 *  fewanimals and day0only read 3 and 2 where 1 was made, QA 0.1.3 r4). An analysis saved without it: the treated. */
export const comparisonsMade = (analysis: DG.DataFrame): number => {
  const tested = analysis.col('tested_arms')?.get(0);
  return Number.isFinite(tested) ? tested : analysis.col('treated_arms')?.get(0) ?? 0;
};

/** An interval's level as its header and words say it ("99% CI"); never "0%": a level not computed (no group
 *  tested) reads as the unadjusted 95% (round 4, analyst A3). */
export const levelText = (confidence: number) =>
  `${+(100 * (Number.isFinite(confidence) && confidence > 0 ? confidence : 0.95)).toFixed(1)}% CI`;

/** The level the cells' intervals are at, from the comparisons actually made (k): 1 − α/k with a correction, 1 − α
 *  without (welchContrasts' own). */
export const intervalLevel = (alpha: number, correction: unknown, k: number) =>
  1 - (correction === 'none' || k < 2 ? alpha : alpha / k);

/** A slope difference `d` (ln V per day) as the change of the daily growth: 100 (1 − e^d) % slower, or
 *  100 (e^d − 1) % faster. */
const perDay = (d: number) => d <= 0 ? {side: 'slower', x: 100 * (1 - Math.exp(d))} :
  {side: 'faster', x: 100 * (Math.exp(d) - 1)};

/** The test the result is read on, in words: "grows 4.1%/day slower (95% CI 2.0–6.3%), p = 0.003", from the
 *  difference of the mean slopes of ln V per day (treated − control), its interval at `confidence` (the adjusted
 *  one) and its adjusted p. An interval across no difference names each end's side. null: no test. */
export function growthRateWords(d: number, lower: number, upper: number, p: number,
  confidence: number): string | null {
  if (activityOf({lower, upper, p}) === TOO_FEW || !Number.isFinite(d))
    return null;
  const [m, lo, hi] = [perDay(d), perDay(lower), perDay(upper)];
  const f = (x: number) => x.toFixed(1);
  // the interval in the side of the estimate, its smaller end first (round 4, analyst W2), else both sides named
  const ci = lo.side !== hi.side ? `${f(lo.x)}% slower to ${f(hi.x)}% faster` : m.side === 'slower' ?
    `${f(hi.x)} to ${f(lo.x)}% slower` : `${f(lo.x)} to ${f(hi.x)}% faster`;
  const pText = pValue(p);
  return `Tumors grew ${f(m.x)}% per day ${m.side} than in the control (${levelText(confidence)} ${ci}; ` +
    `p ${pText.startsWith('<') ? pText : `= ${pText}`})`;
}

/** What the growth-rate cell means, on hover after its sentence (the method is in the protocol). */
export const GROWTH_RATE_TIP = 'This is the test that decides the result; TGI and T/C describe the evaluation day ' +
  'only.';

/** What the view says when no group is named as the control. */
export const NO_CONTROL = 'No control group found: pick it in Control group.';

/** A control whose measurements mostly fall outside its band from the other controls. */
export interface OddControl {animal: string, outside: number, measurements: number, below: boolean}

export interface StudyAnalysis {
  declaration: AnalysisDeclaration;
  /** whether the twins were computed (false: their tables empty, asked for later) */
  twinned: boolean;
  treated: string[];
  /** the treated groups tested against the control (the comparisons made, the correction's k) and the level their
   *  intervals are at (armLevels; with Holm, the widest: its first step's) */
  tested: number;
  confidence: number;
  evaluationDay: number;
  /** the evaluation day came from the endpoint (false: given, or the last day when the controls never reach it) */
  fromEndpoint: boolean;
  nControl: number;
  /** controls measured on the evaluation day */
  nControlDay: number;
  controlMedianTte: number | null;
  arms: ArmAnalysis[];
  /** per animal: group, animal, slope (ln V per day), time to endpoint and whether reached, mRECIST */
  animals: {[k: string]: any}[];
  /** per group and day: animals measured, median volume and the observed range */
  volumes: {group: string, day: number, n: number, median: number, lower: number, upper: number}[];
  /** stats.tumor_twin: per treated animal and day the twin's band, per animal the paired effect, the control check */
  twins: {[k: string]: any}[];
  effects: {[k: string]: any}[];
  checks: {[k: string]: any}[];
  /** share of the controls' measurements inside their band from the others (about 0.9 when the line's model fits) */
  controlCoverage: number | null;
  oddControls: OddControl[];
  versusPlan: PlanComparison;
  notes: string[];
}

/** An effect against the range a study run as planned shows (95% of such studies fall inside); inconclusive: inside
 *  it, but the study cannot tell the group from the control. */
export type EffectVerdict = 'As planned' | 'Inconclusive' | 'Weaker than planned' | 'Stronger than planned';

/** The study against its plan's predictive checks (or, with none, its controls and end against the line's atlas):
 *  each observed figure against the 95% range a study run as planned shows. */
export interface PlanComparison {
  /** the plan's name; null: none linked */
  plan: string | null;
  /** a study's median control per day, its 95% range: the plan's, else the atlas's */
  controls: PlanChecks['controls'];
  /** the controls' mean slope of log V per day to the evaluation day */
  controlSlope: {observed: number | null, range: Range,
    verdict: 'As predicted' | 'Faster than predicted' | 'Slower than predicted' | null};
  /** the first measurement day the median control reached the endpoint */
  end: {observed: number | null, range: Range, notReached: number,
    /** null: nothing predicted to read it against (no model nor plan) */
    verdict: 'As predicted' | 'Earlier than predicted' | 'Later than predicted' | 'Endpoint not reached' | null};
  /** each treated group: the plan's arm it ran as, the effect assumed, the one observed on the evaluation day (no
   *  interval: the test on the slopes is the result's), the range a study run as planned shows, and the verdict */
  groups: {group: string, plannedAs: string | null, metric: EffectMetric, assumed: number | null,
    observed: number | null, expected: Range | null,
    verdict: EffectVerdict | null,
    /** the test on the slopes (adjusted p), the group's analysable animals, the power the study had for the planned
     *  effect with them, and the animals the plan needed for its power (null: no plan, or a plan without them) */
    activity: Activity, p: number, n: number, power: number | null, nNeeded: number | null,
    /** with no plan, for a group with no clear effect: the smallest TGI (%) a study with its animals detects in 80
     *  of 100 labs (null: out of reach, or not asked) */
    mde: number | null}[];
  /** the plan's power target (0–1; null: no plan) */
  powerTarget: number | null;
  notes: string[];
}

/** The end's range Analyze checks: Plan's widened one (no factor: as it was) and, unless the study was reanchored (its
 *  range is predicted from another start), with the plan's printed range inside it. */
export function endRange(end: Range, f: BandFactors | null, p: {durationQ10: number, durationQ90: number | null} | null,
  reanchored: boolean): Range {
  const w = widenedEnd(end, f);
  return reanchored || !p || w === end ? w : withPrinted(w, p.durationQ10, p.durationQ90);
}

/** The end of the study against its range: `reached` the first day the median control was at the endpoint, null when
 *  it never was (later than planned only if the study was measured beyond the range's top). */
export function endVerdict(reached: number | null, lastDay: number, r: Range): PlanComparison['end']['verdict'] {
  if (reached === null)
    return r.q975 !== null && lastDay > r.q975 ? 'Later than predicted' : 'Endpoint not reached';
  return within(reached, r) ? 'As predicted' : reached < r.q025 ? 'Earlier than predicted' : 'Later than predicted';
}

/** Whether x falls in a 95% range (an open top: no upper limit). */
const within = (x: number, r: Range) => x >= r.q025 && (r.q975 === null || x <= r.q975);

/** The controls' growth rate against its predicted range, in Analyze's words (Compare reads a control the same way). */
export const slopeSide = (x: number, r: Range) => within(x, r) ? 'As predicted' as const :
  x > (r.q975 ?? Infinity) ? 'Faster than predicted' as const : 'Slower than predicted' as const;

/** The days with a tumor measured: a study weighed on other days has rows with no volume, and those days are not
 *  ones an analysis can be on. */
const volumeDays = (data: StudyData) => data.day.filter((_, i) => Number.isFinite(data.volume[i]));

/** Plan's evaluation day, observed: the first day the median control is measured at or above the endpoint; the
 *  last measurement day when it never is. */
export function evaluationDay(d: AnalysisDeclaration): {day: number, fromEndpoint: boolean, incomplete?: number[],
  leftOut?: string[]} {
  const {data} = d;
  const on = new Map<number, Set<string>>();
  data.day.forEach((t, i) => Number.isFinite(data.volume[i]) && t >= 0 &&
    (on.get(t) ?? on.set(t, new Set()).get(t)!).add(data.group[i]));
  // a group measured on its first day only (no growth to read: "too few animals to test") waits for no day; the
  // others must all have a volume on it (QA 0.1.3 round 4 chaos N4: one such group made the evaluation day 0)
  const first = minOf(on.keys());
  const later = new Set([...on].filter(([t]) => t > first).flatMap(([, gs]) => [...gs]));
  const groups = new Set([...new Set(data.group)].filter((g) => later.has(g) || g === d.control));
  const leftOut = [...new Set(data.group)].filter((g) => !groups.has(g));
  // a day some group has no volume on (half entered) is never the evaluation day: TGI and T/C read 0% there. A group
  // already gone the measurement day before (every animal lost early) does not hold it back: the others are read on
  // the control's endpoint day, that group shown as not measured then (QA 0.1.3 round 4, analyst A2)
  const sorted = [...on.keys()].sort((a, b) => a - b);
  const complete = (t: number) => {
    const before = on.get(sorted[sorted.indexOf(t) - 1]);
    return [...groups].every((g) => on.get(t)?.has(g) || (before !== undefined && !before.has(g)));
  };
  const incomplete: number[] = [];
  const found = evaluationDayOf(d, on, complete, incomplete);
  return leftOut.length ? {...found, leftOut} : found;
}
function evaluationDayOf(d: AnalysisDeclaration, on: Map<number, Set<string>>, complete: (t: number) => boolean,
  incomplete: number[]): {day: number, fromEndpoint: boolean, incomplete?: number[]} {
  const {data} = d;
  const sorted = [...on.keys()].sort((a, b) => a - b);
  const byDay = new Map<number, number[]>();
  data.group.forEach((g, i) => {
    if (g === d.control && data.day[i] >= 0 && Number.isFinite(data.volume[i])) {
      const own = byDay.get(data.day[i]);
      if (own)
        own.push(data.volume[i]);
      else
        byDay.set(data.day[i], [data.volume[i]]);
    }
  });
  for (const day of [...byDay.keys()].sort((a, b) => a - b)) {
    if (median(byDay.get(day)!) < d.endpoint)
      continue;
    if (complete(day))
      return {day, fromEndpoint: true, incomplete};
    incomplete.push(day);
  }
  const last = [...sorted].reverse().find(complete);
  if (last === undefined)
    return {day: maxOf(volumeDays(data)), fromEndpoint: false};
  incomplete.push(...sorted.filter((t) => t > last && !incomplete.includes(t)));
  return {day: last, fromEndpoint: false, incomplete: incomplete.sort((a, b) => a - b)};
}

/** The groups with no animal measured on or after `day` while the study went on (every animal lost before it), by
 *  their last measurement day: shown as not measured on the evaluation day, never moving it. */
export function groupsGone(data: StudyData, day: number): Map<string, number> {
  const last = new Map<string, number>();
  data.group.forEach((g, i) => Number.isFinite(data.volume[i]) && data.day[i] >= 0 &&
    last.set(g, Math.max(last.get(g) ?? -Infinity, data.day[i])));
  return new Map([...last].filter(([, t]) => t < day));
}

/** Crux stats.tumor_twin's model of the line: the engine priors' study prior, noise and start spread with the twin's
 *  own fit (between-animal spread, process noise, slow animals, skew). */
export async function twinModel(line: string) {
  const lp = await linePrior(line);
  const [sp, nz, tw] = [lp.study_prior, lp.noise, lp.twin];
  if (!sp || !nz || !tw)
    throw new Error(`PMX Studio has no growth model for ${line} yet: pick another tumor model.`);
  return {mean_log_alpha: sp.mean_log_alpha, sd_study_log_alpha: sp.sd_log_alpha, mean_log_beta: sp.mean_log_beta,
    sd_study_log_beta: sp.sd_log_beta, rho_study: sp.rho, log_v_ref: lp.log_v_ref, start_sd: lp.start_sd,
    me_c0: nz.c0, me_c1: nz.c1, noise_df: nz.df, sd_session: nz.sd_session, v_lod: V_LOD,
    sd_log_alpha: tw.sd_log_alpha, sd_log_beta: tw.sd_log_beta, sd_process: tw.sd_process,
    first_factor: tw.first_factor, start_df: tw.start_df, slow_share: tw.slow_share, slow_shift: tw.slow_shift,
    skew: tw.skew};
}

type Window = {animal: string[], group: string[], day: Float64Array, volume: Float64Array};
type Twins = {twins: {[k: string]: any}[], effects: {[k: string]: any}[], checks: {[k: string]: any}[]};

/** Crux stats.tumor_twin in its three stages, each split across the worker sets and joined into one call's tables:
 *  the controls' log-likelihoods on the posterior's support; their mean curves given the posterior (lw: lp plus each
 *  control's, added in index order, as Crux adds them); then the twins and the checks (each animal on its own random
 *  stream, so the split changes nothing). Splitting the support search too (measured) saved less than its transfers
 *  cost. */
async function stagedTwin(window: Window, params: {[k: string]: any}): Promise<Twins> {
  // the first two stages split the controls, at most one a worker (each also finds the support: an idle one only
  // slows the rest, measured); the last splits every animal
  const controls = new Set(window.animal.filter((_, i) => window.group[i] === params.control)).size;
  const [few, shards] = [Math.max(1, Math.min(MAX_SETS, controls)), MAX_SETS];
  const call = (stage: string, extra: (shard: number) => object, n = few) => Promise.all(Array.from({length: n},
    (_, shard) => crux('stats.tumor_twin', window, {...params, stage, shard, shards: n, ...extra(shard)})));
  const of = (rs: CruxResult[], kind: string): Float64Array[] => rs.flatMap((r) => named(r, 'stage'))
    .filter((x) => x.kind === kind).sort((a, b) => a.index - b.index).map((x) => x.values);
  const csv = (vs: ArrayLike<number>[]) => vs.map((v) => Array.from(v).join(',')).filter((s) => s).join(',');
  const l = await call('loglik', () => ({}));
  const [[points], [lp], ll] = [of(l, 'points'), of(l, 'lp'), of(l, 'll')];
  const lw = Float64Array.from(lp, (v, i) => ll.reduce((acc, x) => acc + x[i], v));
  const shared = {points: csv([points]), lw: csv([lw])};
  const means = csv(of(await call('means', () => shared), 'mean'));
  const t = await call('twins', (s) => ({...shared, means,
    ll: params.check === false ? '' : csv(ll.filter((_, i) => i % shards === s))}), shards);
  // in the one call's order: animals as the study lists them, then days
  const order = new Map(window.animal.map((a, i): [string, number] => [a, i]).reverse());
  const sorted = (xs: {[k: string]: any}[]) =>
    xs.sort((a, b) => order.get(a.animal)! - order.get(b.animal)! || (a.day ?? 0) - (b.day ?? 0));
  return {twins: sorted(t.flatMap((r) => rows(r.main as any))), effects: sorted(t.flatMap((r) => named(r, 'effects'))),
    checks: sorted(t.flatMap((r) => named(r, 'checks')))};
}

/** The twin on the browser's WebGPU (Crux's variant webgpu, one call: the likelihoods, convolutions and draws on the
 *  device, the same posterior, the bands within Monte Carlo error); where there is no WebGPU or the device fails,
 *  the staged twin on the workers. */
let warmLate = false;

async function fastTwin(window: Window, params: {[k: string]: any}): Promise<Twins> {
  // after the warm-up: the call then finds its set free and warm, and opens no other (cold) one; at most 5 s, once a
  // page (a warm-up that hangs only leaves the twin slower). The device's question is the warm-up's first step, so it
  // is under the same limit
  if (!warmLate) {
    const late = new Promise<'late'>((r) => setTimeout(() => r('late'), 5000));
    if (await Promise.race([warmTwin().catch(() => null), late]) === 'late')
      warmLate = true;
  }
  if (await gpu()) {
    try {
      const r = await crux('stats.tumor_twin', window, {...params, device: 'webgpu'});
      return {twins: rows(r.main as any), effects: named(r, 'effects'), checks: named(r, 'checks')};
    } catch (e) {
      console.warn('PMX Studio: the twin on WebGPU failed, on the CPU workers instead', e);
    }
  }
  return stagedTwin(window, params);
}

let warming: Promise<void> | null = null;
/** Opens the twin's WebGPU device in Crux's worker and compiles its shaders, once a page (Analyze asks for it as it
 *  opens, and the twin waits for it): the first twin of a page took 2.2 s for them, then ~0.18 s (measured). A tiny
 *  study's twin, once with series of up to 8 measurements and once up to 16 (the likelihood's shader is built for the
 *  longest series, rounded up to 8). */
export const warmTwin = (): Promise<void> => warming ??= warm();

async function warm(): Promise<void> {
  if (!await gpu())
    return;
  const m = await twinModel('MC38');
  for (const n of [4, 12]) {
    const t: Window = {animal: [], group: [], day: new Float64Array(3 * n), volume: new Float64Array(3 * n)};
    ['c1', 'c2', 't1'].forEach((a, k) => Array.from({length: n}, (_, i) => {
      t.animal.push(a);
      t.group.push(a[0]);
      t.day[k * n + i] = 3 * i;
      t.volume[k * n + i] = 100 * Math.exp(0.1 * 3 * i);
    }));
    await crux('stats.tumor_twin', t, {...m, control: 'c', draws: 100, check: false, device: 'webgpu'});
  }
}

/** The last studies' twins by their inputs: a rerun that only changes the correction, α or the limits keeps them (they
 *  depend on nothing else), and going back to an evaluation day or endpoint just seen finds its twins. */
const twinOf = remembered(([window, params]: [Window, {[k: string]: any}]) => fastTwin(window, params), 8);

export async function analyzeStudy(d: AnalysisDeclaration): Promise<StudyAnalysis> {
  const groups = [...new Set(d.data.group)];
  if (!groups.includes(d.control)) {
    throw new NoControl(d.control, `There is no group "${d.control}" in the study table: pick the control among ` +
      `${groups.join(', ')}.`);
  }
  const treated = groups.filter((g) => g !== d.control);
  if (treated.length === 0) // (the user's study as it is, not a failure: said, not logged as an error)
    throw new EmptyState('The study has no treated group');
  const end: ReturnType<typeof evaluationDay> = d.day === null ? evaluationDay(d) : {day: d.day, fromEndpoint: false};
  // the tests and the twins read the volumes of that very day, and a growth rate takes two days from randomization
  const measured = [...new Set(volumeDays(d.data).filter((t) => t >= 0))].sort((a, b) => a - b);
  if (!measured.includes(end.day) || end.day <= measured[0]) {
    throw new Error(`Evaluation day ${end.day} cannot be analyzed: choose a measurement day after the first, ` +
      `one of ${measured.slice(1).join(', ')}`);
  }

  // from randomisation to the evaluation day: the window the design's slope and effects are defined on
  const keep = d.data.day.map((t) => t >= 0 && t <= end.day);
  const pick = <T>(xs: T[]) => xs.filter((_, i) => keep[i]);
  const window = {animal: pick(d.data.animal), group: pick(d.data.group),
    day: Float64Array.from(pick(d.data.day)), volume: Float64Array.from(pick(d.data.volume))};
  // the twin is the analysis's longest call and needs nothing from the others
  // (no model: no twin, which is the model's prediction)
  const twinned = d.twin !== false && !!d.line;
  const twinning = !twinned ? null :
    twinModel(d.line).then((m) => twinOf([window, {...m, control: d.control, day: end.day}]));
  twinning?.catch(() => null); // (awaited below; a failure of the calls before it is the one reported)
  const tgi = await crux('stats.tgi', window, {control: d.control, change: 'paired', rate_floor: V_LOD,
    endpoint_volume: d.endpoint, days: String(end.day)});

  // Hather's rate is the slope of log10 V: in ln V per day, the scale of the design's growth rates
  const animals = named(tgi, 'animals').map((a): {[k: string]: any} =>
    ({...a, slope: a.rate === null ? null : a.rate * Math.LN10}));
  const sloped = animals.filter((a) => a.slope !== null);
  // the groups coded g0 (the control), g1, ...: a contrast's label is then "gk - g0" whatever the groups' names
  const code = new Map([d.control, ...treated].map((g, k) => [g, `g${k}`]));
  const coded = {slope: Float64Array.from(sloped, (a) => a.slope), group: sloped.map((a) => code.get(a.group)!)};
  // beside the test, what one variance for every group gives (Dunnett is lsmeans' mvt: exact for any arm sizes, as
  // scipy and the box plot); it fails where the test has no row either
  const [contrasts, sensitivity] = await Promise.all([welchContrasts(coded, d.correction, d.alpha),
    crux('stats.contrast', coded, {formula: 'slope ~ group', specs: 'group', contrast: 'trt.vs.ctrl', ref: 'g0',
      adjust: treated.length < 2 ? 'none' : 'mvt', alpha: d.alpha})
      .then((r) => rows(r.main as any) as Contrast[], () => [] as Contrast[])]);

  // each treated animal against its twin, each control against the others
  const {twins, effects, checks} = await twinning ?? {twins: [], effects: [], checks: []};

  const table = rows(tgi.main as any);
  const delay = named(tgi, 'delay');
  // each animal's time to the endpoint over the whole study, as Time to endpoint draws it (stats.tgi's stops at the
  // evaluation day, and its file said "not reached" of an animal the curves show reaching it later)
  const all = d.data.day.map((t, i) => Number.isFinite(t) && Number.isFinite(d.data.volume[i]));
  const of = <T>(xs: T[]) => xs.filter((_, i) => all[i]);
  const survival = await crux('stats.tumor_survival', {animal: of(d.data.animal), group: of(d.data.group),
    day: Float64Array.from(of(d.data.day)), volume: Float64Array.from(of(d.data.volume))},
  {control: d.control, endpoint_volume: d.endpoint}).then((r) => named(r, 'animals'), () => []);
  const reaching = new Map(survival.map((a) => [String(a.animal), a]));
  for (const a of animals) {
    const s = reaching.get(String(a.animal));
    if (s)
      Object.assign(a, {time_to_endpoint: s.time, event: !!s.event});
  }
  // the animals each group has on the evaluation day: those its TGI and T/C are read on
  const onDay = (g: string) => new Set(window.animal.filter((_, i) => window.group[i] === g &&
    window.day[i] === end.day && Number.isFinite(window.volume[i]))).size;
  const meanSlope = (g: string) => {
    const s = sloped.filter((a) => a.group === g).map((a) => a.slope as number);
    return s.reduce((t, x) => t + x, 0) / s.length;
  };
  const arms: ArmAnalysis[] = treated.map((g, k) => {
    const e = table.find((r) => r.group === g) ?? {};
    // a group with no animal measured twice has no contrast: no test (said in the notes)
    const c = contrasts.find((r) => String(r.term).startsWith(`g${k + 1} `)) ??
      {estimate: NaN, lower: NaN, upper: NaN, p_value: NaN, level: NaN, p_raw: NaN};
    const own = effects.filter((r) => r.group === g);
    // (none measured that day: no TGI nor T/C, never a 0% read on no animal: chaos 4)
    const nDay = onDay(g);
    return {arm: g, n: sloped.filter((a) => a.group === g).length, nDay, tgi: nDay ? finite(e.tgi_change) : null,
      tc: nDay ? finite(e.tc) : null, slope: meanSlope(g), slopeControl: meanSlope(d.control),
      difference: c.estimate, lower: c.lower, upper: c.upper, confidence: c.level ?? NaN,
      p: finite(c.p_value) ?? NaN, pRaw: finite(c.p_raw) ?? NaN,
      pCommon: finite(sensitivity.find((r) => String(r.term).startsWith(`g${k + 1} `))?.p_value),
      medianTte: finite(delay.find((r) => r.group === g)?.median),
      delay: finite(delay.find((r) => r.group === g)?.t_minus_c),
      twinAnimals: own.length, pairedEffect: own.length ? 100 * median(own.map((r) => r.effect)) : null,
      responders: own.filter((r) => r.effect_lo > 0).length};
  });

  // a control that mostly falls outside its band grows unlike the others: shown, never dropped
  const byControl = new Map<string, {[k: string]: any}[]>();
  checks.forEach((r) => byControl.set(r.animal, [...(byControl.get(r.animal) ?? []), r]));
  const oddControls: OddControl[] = [...byControl].map(([animal, rs]) => {
    const out = rs.filter((r) => !r.inside);
    return {animal, outside: out.length, measurements: rs.length,
      below: out.filter((r) => r.observed < r.q05).length > out.length / 2};
  }).filter((c) => c.outside >= 2 && c.outside >= c.measurements / 2);

  // each group's median and range per measurement day (tgi's group_days holds only the evaluation day)
  const byGroupDay = new Map<string, {group: string, day: number, v: number[]}>();
  window.group.forEach((g, i) => {
    if (!Number.isFinite(window.volume[i]))
      return;
    const key = `${g}\u0000${window.day[i]}`;
    if (!byGroupDay.has(key))
      byGroupDay.set(key, {group: g, day: window.day[i], v: []});
    byGroupDay.get(key)!.v.push(window.volume[i]);
  });
  const volumes = [...byGroupDay.values()].map(({group, day, v}) =>
    ({group, day, n: v.length, median: median(v), lower: minOf(v), upper: maxOf(v)}));

  const notes: string[] = [];
  if (d.data.skipped > 0) {
    notes.push(`${d.data.skipped} row${d.data.skipped > 1 ? 's' : ''} without an animal, a group or a day left out ` +
      'of the analysis');
  }
  const leftOut = end.leftOut ?? [];
  if (d.day === null && leftOut.length) {
    const first = minOf(volumeDays(d.data).filter((t) => t >= 0));
    const many = leftOut.length > 1;
    notes.push(`${many ? 'Groups' : 'Group'} ${leftOut.join(', ')} ${many ? 'have' : 'has'} no measurements after ` +
      `day ${first}, so ${many ? 'they are' : 'it is'} left out. Effects are read on day ${end.day}.`);
  }
  if (d.day === null && end.incomplete?.length) {
    const ds = end.incomplete;
    // the groups with no volume on each such day, named (one left out above is not waited for)
    const groups = [...new Set(d.data.group)].filter((g) => !leftOut.includes(g));
    const missingOn = (t: number) => groups.filter((g) => !d.data.day.some((x, i) => x === t &&
      d.data.group[i] === g && Number.isFinite(d.data.volume[i]))).join(', ');
    const whose = ds.length > 1 ? ds.map((t) => `${missingOn(t)} on day ${t}`).join('; ') : missingOn(ds[0]);
    // only the last day can be one still being entered: an earlier one is a finished day
    const lastDay = maxOf(measured);
    notes.push(`Day${ds.length > 1 ? 's' : ''} ${ds.join(', ')} not used as the evaluation day: no volume ` +
      `${ds.length > 1 ? '' : 'on it '}for ${whose}${ds.includes(lastDay) ? ' (not entered yet?)' : ''}. ` +
      `The effects are read on day ${end.day}` + (end.fromEndpoint ? '' : ', the last day every group was measured'));
  }
  if (d.day === null && !end.fromEndpoint && !end.incomplete?.length) {
    // the median volume of the controls still measured each day, not the median day to the endpoint: animals taken
    // out at the endpoint leave that median, so it can stay under the endpoint after half the controls reached it
    const own = d.data.day.map((t, i) => d.data.group[i] === d.control && Number.isFinite(d.data.volume[i]) && t >= 0 ?
      t : NaN);
    const medians = [...new Set(own.filter(Number.isFinite))].map((t) => ({t, m: median(d.data.volume
      .filter((_, i) => own[i] === t))}));
    const top = medians.reduce((b, x) => x.m > b.m ? x : b, medians[0]);
    notes.push(`The control median never reached ${mm3Text(d.endpoint)} mm³` + (top ? ` (largest: ${mm3Text(top.m)} mm³ ` +
      `on day ${top.t})` : '') + `, so effects are read on day ${end.day}, the last measurement day`);
  }
  // a group every animal of which left before the evaluation day: one sentence, no TGI nor T/C that day
  // (one note per fact: a group the automatic day left out is said so once, above)
  const saidOut = (g: string) => d.day === null && leftOut.includes(g);
  const gone = groupsGone(d.data, end.day);
  for (const [g, t] of gone) {
    if (g !== d.control && !saidOut(g))
      notes.push(`${g} has no animals measured after day ${t}: no TGI nor T/C on day ${end.day}`);
  }
  // animals measured once: said once, by group when a whole group is (round 4, analyst U4/W6)
  const unslopedAll = animals.filter((a) => a.slope === null);
  const whole = treated.filter((g) => !saidOut(g) && !sloped.some((a) => a.group === g) &&
    unslopedAll.some((a) => a.group === g));
  for (const g of whole) {
    const days = [...new Set(d.data.day.filter((t, i) => d.data.group[i] === g && t >= 0 &&
      Number.isFinite(d.data.volume[i])))].sort((a, b) => a - b);
    notes.push(`${g} was measured only on day${days.length > 1 ? 's' : ''} ${days.join(', ') || '0'} and is left ` +
      `out of the test. Effects read on day ${end.day}`);
  }
  // (a group said once already, measured only on day 0 or left out by the automatic day, gets no other note: the
  // day0only file read "left out", "8 animals with fewer than two measurements" and "not tested, 0 animals with a
  // slope" for one group, QA 0.1.3 browser check)
  const saidOnce = (g: string) => whole.includes(g) || saidOut(g);
  const unsloped = unslopedAll.filter((a) => !saidOnce(a.group)).length;
  if (unsloped > 0)
    notes.push(`${unsloped} animal${unsloped > 1 ? 's' : ''} with fewer than two measurements left out of the test`);
  const nControl = sloped.filter((a) => a.group === d.control).length;
  const nControlDay = onDay(d.control);
  // (a group with fewer animals on the evaluation day than in the test: its Analyzed cell says "7 of 10" and why on
  // hover, no sentence here: round 4, analyst U1/W8)
  if (nControl < MIN_TESTED)
    notes.push(`${d.control}: ${tooFewWhy(nControl)}: no group is tested`);
  for (const a of arms) {
    if (a.n < MIN_TESTED && !saidOnce(a.arm))
      notes.push(`${a.arm}: not tested, ${tooFewWhy(a.n)}`);
    // the sensitivity analysis speaks when it disagrees
    if (a.pCommon !== null && Number.isFinite(a.p) && (a.p < d.alpha) !== (a.pCommon < d.alpha)) {
      notes.push(`${a.arm}: p = ${a.p.toPrecision(2)} with each group's own variance (the test); with one variance ` +
        `for every group (${treated.length < 2 ? 'the t test' : 'Dunnett'}) it would be ${a.pCommon.toPrecision(2)}`);
    }
  }
  for (const c of oddControls) {
    notes.push(`Control ${c.animal} grows ${c.below ? 'slower' : 'faster'} than the other controls predict: ` +
      `${c.outside} of its ${c.measurements} measurements fall outside their range. It stays in the analysis`);
  }
  const inside = checks.filter((r) => r.inside).length;
  const versusPlan = await compareWithPlan(d, treated, arms, nControl);
  const tested = contrasts.length;
  // the arms' level (Holm: the widest, its first step's); none tested: the method's at k
  const levels = contrasts.map((c) => c.level ?? NaN).filter(Number.isFinite);
  const confidence = levels.length ? minOf(levels) : intervalLevel(d.alpha, d.correction, tested);
  return {declaration: d, twinned, treated, tested, confidence, evaluationDay: end.day, fromEndpoint: end.fromEndpoint,
    versusPlan,
    nControl, nControlDay, controlMedianTte: finite(delay.find((r) => r.group === d.control)?.median), arms,
    animals, volumes,
    twins, effects, checks, controlCoverage: checks.length ? inside / checks.length : null, oddControls, notes};
}

/** fn's results by its input (a failure not kept, so the next call tries again), the `keep` most recently used: a
 *  rerun (another correction, a tab opened again) does not compute them again. */
export function remembered<A, R>(fn: (a: A) => Promise<R>, keep = 16): (a: A) => Promise<R> {
  const seen = new Map<string, Promise<R>>();
  return (a) => {
    const key = JSON.stringify(a);
    let r = seen.get(key);
    if (r)
      seen.delete(key); // the most recent last: the oldest goes first
    else {
      r = fn(a);
      const mine = r;
      r.catch(() => seen.get(key) === mine && seen.delete(key));
    }
    seen.set(key, r);
    if (seen.size > keep)
      seen.delete(seen.keys().next().value!);
    return r;
  };
}

/** The designs an analysis with no plan reads its study against, and the checks of a study's own start. */
const noPlanDesign = remembered<StudyDeclaration, StudyDesign>(designStudy);
const checksAt = remembered<ChecksInput, PlanChecks>(planChecks);

/** The atlas's controls for `line` (Compare: each study's control against them): as Analyze predicts them for a
 *  study without a plan, randomized at `vRand`, measured to `lastDay`. */
export const atlasControls = (line: string, vRand: number, lastDay: number, nControl: number,
  endpoint = ENDPOINT) =>
  checksAt({line, vRand, endpoint, humane: HUMANE, weekdays: LAB_WEEKDAYS, startWeekday: 0,
    endpointDay: Math.max(lastDay, 7), nControl, arms: []});

/** The study's volume at randomization: the median of every animal's day-0 volume (before any treatment, so every
 *  group's); null when none was measured. */
function startVolume(d: AnalysisDeclaration): number | null {
  const day0 = d.data.volume.filter((v, i) => d.data.day[i] === 0 && Number.isFinite(v));
  return day0.length ? median(day0) : null;
}

/** A study randomized at a volume this far from the plan's (either way) has its controls and its end predicted
 *  again from its own start, where its early controls fall outside the band in 4 to 5 of 100 studies of normal growth
 *  whatever the start. Against the plan's band (MC38, CT26, B16F10): at most 8 from 10% off, 12 from 15%, 18 from
 *  20%. Lower costs nothing but more studies shown recomputed numbers (the observed start alone varies ~8%). */
export const REANCHOR = 1.10;

/** A group's or an arm's name as they are matched by: any case, spaces and punctuation set aside ("p,o." is
 *  "p.o.", round 3 analyst 16), a decimal kept whichever its separator ("0,3" is "0.3", never "03" nor "3"), "+"
 *  kept (a combination). */
export const armKey = (s: string) => String(s ?? '').normalize('NFKC').toLowerCase()
  .replace(/(\d)[.,](?=\d)/g, '$1ǀ').replace(/[^\p{L}\p{N}+]+/gu, '');

/** Each treated group with the plan's arm it ran as: the same name (any case), else, for the groups left, the arms
 *  left in their order when they are as many (said in a note). */
export function pairGroups(treated: string[], arms: SavedPlan['predicted']['arms'], notes: string[]) {
  const key = armKey;
  const byName = new Map(arms.map((a) => [key(a.name), a]));
  const pairs = new Map(treated.flatMap((g) => byName.has(key(g)) ? [[g, byName.get(key(g))!] as const] : []));
  const used = new Set(pairs.values());
  const [groupsLeft, armsLeft] = [treated.filter((g) => !pairs.has(g)), arms.filter((a) => !used.has(a))];
  if (groupsLeft.length > 0 && groupsLeft.length === armsLeft.length) {
    groupsLeft.forEach((g, i) => pairs.set(g, armsLeft[i]));
    notes.push('Matched to the plan by their order: ' +
      groupsLeft.map((g, i) => `${g} as ${armsLeft[i].name}`).join(', '));
  } else {
    for (const g of groupsLeft)
      notes.push(`${g} is not in the plan`);
  }
  return pairs;
}

/** The observed effect inside the plan's 95% range, or beyond it on the strong or the weak side (a larger TGI, a
 *  smaller T/C is stronger). Inside it while the test on the slopes cannot tell the group from the control, the
 *  study cannot tell the planned effect from none: inconclusive, not as planned (one test decides both). */
function verdictOf(metric: EffectMetric, observed: number | null, activity: Activity,
  expected: Range): EffectVerdict | null {
  if (observed === null || activity === TOO_FEW)
    return null;
  if (within(observed, expected))
    return activity === 'Slows growth' ? 'As planned' : 'Inconclusive';
  return (metric === 'TGI' ? observed > (expected.q975 ?? Infinity) : observed < expected.q025) ?
    'Stronger than planned' : 'Weaker than planned';
}

/** The study against its plan's predictive checks: the controls' growth, the end and each treated group's effect,
 *  each against the 95% range a study run as planned shows. With no plan, the controls and the end against the
 *  line's atlas: the checks of Plan's design for this line and endpoint, randomized at the study's median on day 0,
 *  with this study's controls (ponytail: on Mon / Wed / Fri; the study's own weekdays if the calendar matters). */
async function compareWithPlan(d: AnalysisDeclaration, treated: string[], arms: ArmAnalysis[],
  controls: number): Promise<PlanComparison> {
  const notes: string[] = [];
  let p = d.plan?.predicted;
  const mde = new Map<string, number | null>();
  const start = startVolume(d);
  if (!p && !d.line) {
    const at = evaluationDay(d);
    return modelless(arms, at.fromEndpoint ? at.day : null);
  }
  if (!p) {
    const vRand = start ?? 100;
    const design = (n: number | null) => noPlanDesign({line: d.line, vRand, endpoint: d.endpoint, power: 0.8,
      alpha: d.alpha, correction: d.correction, n, treatedArms: arms.length,
      arms: [{name: 'x', therapy: 'small_molecule', metric: 'TGI', effect: 60}],
      preLoss: 0, postLoss: 0, lossRule: 'confident', confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0,
      humane: HUMANE, schedules: false, lastDay: null});
    // a group the test cannot tell from the control: the smallest effect a study with its animals (and the study's
    // comparisons) detects, which only it shows (a search of seconds: not for the groups with an effect); all at once
    const unclear = arms.filter((x) => activityOf(x) === 'No clear effect');
    const [s, ...searched] = await Promise.all([design(null), ...unclear.map((a) => design(Math.max(3, a.n)))]);
    unclear.forEach((a, i) => mde.set(a.arm, searched[i].arms[0].mde));
    const checks = await checksAt({line: d.line, vRand, endpoint: d.endpoint, humane: HUMANE,
      weekdays: LAB_WEEKDAYS, startWeekday: 0, endpointDay: s.endpointDay, nControl: controls, arms: []});
    p = {line: d.line, endpoint: d.endpoint, vRand, studyEnd: s.endpointDay, endpointDay: s.endpointDay,
      durationQ10: s.durationQ10, durationQ90: s.durationQ90, nControl: controls, arms: [], checks};
  } else {
    if (p.line !== d.line)
      notes.push(`The plan is for ${modelName(p.line)}; this analysis uses the ${modelName(d.line)} model`);
    if (p.endpoint !== d.endpoint)
      notes.push(`The plan's endpoint is ${p.endpoint} mm³; this analysis uses ${d.endpoint} mm³`);
  }
  // a study randomized larger or smaller than planned grows apart from the plan's band from day 0: its controls and
  // its end predicted again from its own start, by the plan's model, line and schedule (the .pmxplan unchanged; each
  // arm's range stays the plan's: it moves little with the start, and the plan keeps no kill rate to compute it again)
  let c = p.checks;
  // (a median control that never reaches the endpoint has no day)
  const dayOf = (d: number | null) => Number.isFinite(d) ? `about day ${Math.round(d!)}` : 'not within the study';
  let reanchored = false;
  if (d.plan && start !== null && Math.abs(Math.log(start / p.vRand)) > Math.log(REANCHOR)) {
    const pd = declarationOf(d.plan.inputs as unknown as PlanInputs);
    const own = await checksAt({line: p.line, vRand: start, endpoint: p.endpoint, humane: pd.humane,
      weekdays: pd.weekdays, startWeekday: pd.startWeekday, endpointDay: p.endpointDay,
      nControl: p.nControl ?? controls, arms: []});
    notes.push(`Randomized at ${Math.round(start)} mm³ (median), the plan at ${p.vRand} mm³: the controls and the ` +
      `end are predicted from ${Math.round(start)} mm³, the median control reaching ${p.endpoint} mm³ ` +
      `${dayOf(own.end.median)} (the plan: ${dayOf(c.end.median)})`);
    c = {...c, end: own.end, notReached: own.notReached, controlSlope: own.controlSlope, controls: own.controls};
    reanchored = true;
  }
  // the end as widened in the lines the validation measured, and no narrower than the range Plan printed (its 8 in 10)
  c = {...c, end: endRange(c.end, await bandFactors(d.line), d.plan ? p : null, reanchored)};

  // the controls' growth: their mean slope of log V to the evaluation day (the test's own quantity)
  const slope = arms.length ? arms[0].slopeControl : null;
  const observedSlope = slope !== null && Number.isFinite(slope) ? slope : null;
  const controlSlope: PlanComparison['controlSlope'] = {observed: observedSlope, range: c.controlSlope,
    verdict: observedSlope === null ? null : slopeSide(observedSlope, c.controlSlope)};

  // the first measurement day the median control reached the endpoint, whatever evaluation day was chosen; a study
  // that never did is later than planned only if it was measured beyond the range's top
  const at = evaluationDay({...d, endpoint: p.endpoint});
  const reached = at.fromEndpoint ? at.day : null;
  const end: PlanComparison['end'] = {observed: reached, range: c.end, notReached: c.notReached,
    verdict: endVerdict(reached, maxOf(volumeDays(d.data)), c.end)};

  const pairs = d.plan ? pairGroups(treated, p.arms, notes) : new Map<string, SavedPlan['predicted']['arms'][0]>();
  const groups = arms.map((a) => {
    const planned = pairs.get(a.arm);
    const metric: EffectMetric = planned?.metric ?? 'TGI';
    const observed = metric === 'TGI' ? a.tgi : a.tc;
    const own = planned ? c.arms.find((x) => x.name === planned.name) : undefined;
    const expected = own ? (metric === 'TGI' ? own.tgi : own.tc) : null;
    const activity = activityOf(a);
    return {group: a.arm, plannedAs: planned?.name ?? null, metric, assumed: planned?.effect ?? null, observed,
      expected, verdict: expected ? verdictOf(metric, observed, activity, expected) : null,
      activity, p: a.p, n: a.n, power: planned?.power?.length ? powerAt(planned.power, a.n) : null,
      nNeeded: planned?.nNeeded ?? null, mde: mde.get(a.arm) ?? null};
  });
  // the plan's power, as its form gave it: a percent, or a share
  const target = Number(d.plan?.inputs.power);
  const powerTarget = Number.isFinite(target) && target > 0 ? (target > 1 ? target / 100 : target) : null;
  return {plan: d.plan?.name ?? null, controls: c.controls, controlSlope, end, groups, powerTarget, notes};
}

/** The study with no model nor plan: nothing predicted to read it against (no range, no verdict), each group's
 *  observed effect and test as they are; `reached`: the first day the median control reached the endpoint (null:
 *  never), as the header's evaluation day says it (round 3, analyst 4: "not reached" beside "evaluated on day 16"). */
export function modelless(arms: ArmAnalysis[], reached: number | null = null): PlanComparison {
  const none = {q025: NaN, median: NaN, q975: null};
  const slope = arms.length ? arms[0].slopeControl : null;
  return {plan: null, controls: [], controlSlope: {observed: slope !== null && Number.isFinite(slope) ? slope : null,
    range: none, verdict: null}, end: {observed: reached, range: none, notReached: NaN,
    verdict: reached === null ? 'Endpoint not reached' : null},
  groups: arms.map((a) => ({group: a.arm, plannedAs: null, metric: 'TGI' as EffectMetric, assumed: null,
    observed: a.tgi, expected: null, verdict: null, activity: activityOf(a),
    p: a.p, n: a.n, power: null, nNeeded: null, mde: null})), powerTarget: null, notes: []};
}

/** Days for a tumor growing at `slope` (ln V per day) to double; null when it does not grow. */
export const doublingOf = (slope: number | null | undefined) => slope && slope > 0 ? Math.LN2 / slope : null;

export interface AnalysisFrames {analysis: DG.DataFrame, comparisons: DG.DataFrame, animals: DG.DataFrame,
  volumes: DG.DataFrame, twins: DG.DataFrame, effects: DG.DataFrame, checks: DG.DataFrame, versusPlan: DG.DataFrame,
  predicted: DG.DataFrame, tolerability: DG.DataFrame, weights: DG.DataFrame}

/** The two readings of the body weight (stats.tumor_body_weight's basis), as the Tolerability tab names them. */
export const BASES: {[basis: string]: string} = {measured: 'Measured', tumor_corrected: 'Tumor-corrected'};

/** The NCI limits for a treatment in mice (Plowman et al. 1997): a group's mean weight loss at its nadir, as weighed,
 *  and the share of its animals dead or removed for their weight, both in percent (Crux toxic_nadir, toxic_deaths). */
export const NCI = {nadir: 20, deaths: 10};

/** The tolerability table's tag: the density the tumor-corrected weight starts from, the line's (with the atlas
 *  papers and range behind it) or the one declared. */
export const DENSITY_TAG = 'pmx-density';
export interface DensityUsed extends LineDensity {line: string, declared: boolean}

/** The study's body weight through Crux stats.tumor_body_weight: per group and reading its weight-loss grades, mean
 *  nadir, recoveries and removals (with the animals removed, in words), and per group, reading and day the mean change
 *  from day 0 and its range. The tumor's mass is its volume times the line's density (the atlas's), or the declared
 *  one; with the tumors' ex vivo weights Crux calibrates it to the study and its groups, and the row a tumor was
 *  excised takes off its own mass. Empty tables when the study has no body weight. */
export async function tolerabilityFrames(d: AnalysisDeclaration): Promise<{tolerability: DG.DataFrame,
  weights: DG.DataFrame, animals: {[k: string]: any}[]}> {
  const w = d.data.weight;
  if (!w?.some((x) => typeof x === 'number' && x > 0))
    return {tolerability: DG.DataFrame.create(0), weights: DG.DataFrame.create(0), animals: []};
  const num = (xs: (number | null)[]) => Float64Array.from(xs, (x) => x ?? NaN);
  const atlas = await tumourDensity();
  const line = atlas.lines[d.line] ?? atlas.all;
  const used: DensityUsed = {...line, line: atlas.lines[d.line] ? d.line : 'all lines',
    density: d.density ?? line.density, declared: d.density !== null};
  const s = atlas.spread_log;
  const bw = await crux('stats.tumor_body_weight', {animal: d.data.animal, group: d.data.group,
    day: num(d.data.day), weight: num(w), volume: num(d.data.volume),
    ...(d.data.mass ? {tumor_mass: num(d.data.mass)} : {}), ...(d.data.removal ? {removal: d.data.removal} : {})},
  {weight_limit: d.weightLimit, humane: d.humane, density: used.density, lab_sd: s.between_labs,
    group_sd: s.between_groups, within_sd: Math.hypot(s.between_animals, s.caliper), toxic_nadir: NCI.nadir,
    toxic_deaths: NCI.deaths});
  const why: {[r: string]: string} = {death: 'dead', weight: 'weight', volume: 'tumor', other: 'other reason'};
  const animals = named(bw, 'animals');
  const removed = animals.filter((a) => a.removal !== 'none');
  const tolerability = DG.DataFrame.fromObjects(rows(bw.main as any).map((g) => ({...g, basis: BASES[g.basis],
    removed_animals: removed.filter((a) => a.group === g.group).map((a) =>
      `${a.animal} (day ${Math.round(a.last_day)}, ${why[a.removal]})`).join(', ')})))!;
  tolerability.setTag(DENSITY_TAG, JSON.stringify(used));
  const weights = DG.DataFrame.fromObjects(named(bw, 'days').map((x) => ({...x, basis: BASES[x.basis]})))!;
  weights.col('day')!.meta.units = 'd';
  for (const q of ['mean', 'min', 'max']) {
    weights.col(q)!.meta.units = '%';
    weights.col(q)!.meta.format = '0.0';
  }
  // the grades' edges and the humane endpoint, as the chart's reference lines
  for (const [y, title] of [[-5, '−5%'], [-10, '−10%'], [-d.weightLimit, `−${d.weightLimit}%, humane endpoint`]] as
    [number, string][])
    weights.meta.formulaLines.addLine({title, formula: `\${mean} = ${y}`, color: '#9aa5b1'});
  return {tolerability, weights, animals};
}

/** The predicted control series and the observed one, as the band chart draws them. */
export const PREDICTED = 'Predicted controls';
export const OBSERVED = 'Observed controls';

type ControlBand = PlanComparison['controls'][number];
/** The predicted controls to `day`: the plan's measurement days up to it, and `day` itself between two of them (its
 *  volumes geometric between theirs: the plan measured days 0, 2, 4, 7 and the study is read on day 6). */
export function controlsTo(controls: ControlBand[], day: number): ControlBand[] {
  const kept = controls.filter((b) => b.day <= day);
  const [a, b] = [kept[kept.length - 1], controls.find((x) => x.day > day)];
  if (!a || !b || a.day === day)
    return kept;
  const t = (day - a.day) / (b.day - a.day);
  const at = (x: number, y: number) => x > 0 && y > 0 ? Math.exp(Math.log(x) + t * (Math.log(y) - Math.log(x))) :
    x + t * (y - x);
  return [...kept, {day, q025: at(a.q025, b.q025), median: at(a.median, b.median), q975: at(a.q975, b.q975)}];
}

/** The columns of `cols` the table has as volumes: in mm³, to one decimal in every viewer's tooltip and cell (Crux's
 *  floats carry a dozen). */
export function inMm3(df: DG.DataFrame, cols: readonly string[]): void {
  for (const q of cols) {
    const c = df.col(q);
    if (c) {
      c.meta.units = 'mm³';
      c.meta.format = '0.0';
    }
  }
}

/** The analysis (one row), the comparisons (one row per treated arm), the animals, the observed volumes, and the
 *  twin's bands, paired effects and control check. */
export async function analysisFrames(d: AnalysisDeclaration): Promise<AnalysisFrames> {
  const [s, {tolerability, weights}] = await Promise.all([analyzeStudy(d), tolerabilityFrames(d)]);
  const vp = s.versusPlan;
  const analysis = DG.DataFrame.fromObjects([{
    line: d.line, twinned: s.twinned, control: d.control, n_control: s.nControl, n_control_day: s.nControlDay,
    treated_arms: s.treated.length, tested_arms: s.tested, confidence: s.confidence,
    evaluation_day: s.evaluationDay, from_endpoint: s.fromEndpoint, endpoint: d.endpoint,
    control_median_tte: s.controlMedianTte, alpha: d.alpha, correction: d.correction, compared: COMPARED,
    v_lod: V_LOD, control_coverage: s.controlCoverage, odd_controls: s.oddControls.length,
    notes: s.notes.join(' | '),
    plan: vp.plan, control_slope: vp.controlSlope.observed, control_slope_q025: vp.controlSlope.range.q025,
    control_slope_median: vp.controlSlope.range.median, control_slope_q975: vp.controlSlope.range.q975,
    control_verdict: vp.controlSlope.verdict, observed_end: vp.end.observed, end_q025: vp.end.range.q025,
    end_median: vp.end.range.median, end_q975: vp.end.range.q975, end_not_reached: vp.end.notReached,
    end_verdict: vp.end.verdict, power_target: vp.powerTarget,
    plan_notes: vp.notes.join(' | '),
  }])!;
  const comparisons = DG.DataFrame.fromObjects(s.arms.map((a) => ({
    arm: a.arm, n: a.n, n_day: a.nDay, tgi: a.tgi, tc: a.tc, slope_per_day: a.slope,
    slope_control_per_day: a.slopeControl, slope_difference: a.difference,
    slope_difference_lower: a.lower, slope_difference_upper: a.upper,
    slope_difference_confidence: a.confidence, correction: d.correction, p_value: a.p, p_raw: a.pRaw,
    significant: a.p < d.alpha,
    p_common_variance: a.pCommon,
    activity: activityOf(a),
    median_time_to_endpoint: a.medianTte, growth_delay_days: a.delay, twin_animals: a.twinAnimals,
    paired_effect: a.pairedEffect,
    responders: a.responders,
  })))!;
  const animals = DG.DataFrame.fromObjects(s.animals.map((a) => ({group: a.group, animal: a.animal,
    slope_per_day: a.slope, doubling_time: doublingOf(a.slope), time_to_endpoint: a.time_to_endpoint,
    endpoint_reached: a.event, mrecist: a.mrecist,
    // mRECIST's best response (% change of the volume from baseline, its minimum from day 10) and best average:
    // a waterfall's axis
    best_response: a.best_response, best_avg_response: a.best_avg_response})))!;
  const volumes = DG.DataFrame.fromObjects(s.volumes)!;
  // (no rows, as without the twin: an empty table, which fromObjects does not build)
  const table = (xs: object[]) => xs.length ? DG.DataFrame.fromObjects(xs)! : DG.DataFrame.create(0);
  const twins = table(s.twins);
  const effects = table(s.effects.map((e) => ({group: e.group, animal: e.animal,
    effect: 100 * e.effect, effect_lower: 100 * e.effect_lo, effect_upper: 100 * e.effect_hi,
    responder: e.effect_lo > 0})));
  const checks = table(s.checks);
  // each treated group against the plan (with none: its effect, nothing assumed)
  const versusPlan = DG.DataFrame.fromObjects(vp.groups.map((g) => ({group: g.group, planned_as: g.plannedAs,
    metric: g.metric, assumed: g.assumed, observed: g.observed,
    expected_q025: g.expected?.q025 ?? null, expected_median: g.expected?.median ?? null,
    expected_q975: g.expected?.q975 ?? null, verdict: g.verdict, activity: g.activity, p_value: g.p, n: g.n,
    power: g.power, n_needed: g.nNeeded, mde: g.mde})))!;
  // a study's median control as predicted (its 95% range) to the evaluation day, then the observed median control
  const predicted = DG.DataFrame.fromObjects([
    ...controlsTo(vp.controls, s.evaluationDay).map((b) => ({series: PREDICTED, day: b.day,
      lower: b.q025, upper: b.q975, median: b.median})),
    ...s.volumes.filter((v) => v.group === d.control).sort((a, b) => a.day - b.day).map((v) => ({series: OBSERVED,
      day: v.day, lower: v.median, upper: v.median, median: v.median}))])!;

  // the automatic colours: the view gives them the study's own (colors.ts follow)
  const colors = groupColors(s.treated, d.control);
  for (const [df, col] of [[comparisons, 'arm'], [animals, 'group'], [volumes, 'group'], [twins, 'group'],
    [effects, 'group'], [tolerability, 'group'], [weights, 'group']] as const)
    df.col(col)?.meta.colors.setCategorical(colors);
  predicted.col('series')!.meta.colors.setCategorical({[PREDICTED]: '#9aa5b1', [OBSERVED]: colors['Control']});
  for (const df of [volumes, twins, checks, predicted]) {
    const c = df.col('day');
    if (c)
      c.meta.units = 'd';
  }
  for (const [df, cols] of [[volumes, ['median', 'lower', 'upper']], [twins, ['observed', 'twin_q05', 'twin_median',
    'twin_q95']], [checks, ['observed', 'q05', 'median', 'q95']], [predicted, ['lower', 'upper', 'median']]] as const)
    inMm3(df, cols);
  for (const q of ['effect', 'effect_lower', 'effect_upper']) {
    const c = effects.col(q);
    if (c)
      c.meta.units = '%';
  }
  volumes.meta.formulaLines.addLine({title: `Endpoint ${d.endpoint} mm³`, formula: `\${median} = ${d.endpoint}`});
  volumes.meta.formulaLines.addLine({title: `Evaluation, day ${s.evaluationDay}`,
    formula: `\${day} = ${s.evaluationDay}`, style: 'solid', color: '#4d5261'});
  if (s.effects.length)
    effects.meta.formulaLines.addLine({title: 'No effect', formula: '${effect} = 0', color: '#4d5261'});
  const tag = `${s.treated.length + 1} groups, day ${s.evaluationDay}`;
  for (const [df, name] of [[analysis, 'Analysis'], [comparisons, 'Comparisons'], [animals, 'Animals'],
    [volumes, 'Volumes'], [twins, 'Twins'], [effects, 'Paired effects'], [checks, 'Control check'],
    [versusPlan, 'Against the plan'], [predicted, 'Controls against the prediction'],
    [tolerability, 'Tolerability'], [weights, 'Body weight']] as const)
    df.name = `${name} · ${tag}`;
  return {analysis, comparisons, animals, volumes, twins, effects, checks, versusPlan, predicted, tolerability,
    weights};
}
