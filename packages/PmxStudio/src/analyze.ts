/*
 * The Analyze engine: a study's long table (animal, group, day, volume) analysed as Plan pre-specifies it, so the
 * test that decides is the one the design powered: each animal's least-squares slope of log V from randomisation
 * (day 0) to the evaluation day, volumes below the smallest a caliper records raised to it (Crux stats.tgi's
 * per-animal rate, rate_floor), every treated arm against the control by the two-sample t test on those slopes
 * with the declared correction (Crux stats.contrast, rate ~ group, treatment vs control: Dunnett exact). TGI and
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

import {crux, named, rows} from './crux';
import {Correction, EffectMetric, HUMANE, LAB_WEEKDAYS, SEED, StudyDeclaration, StudyDesign, correctionOf,
  designStudy, finite, isFiniteNumber as ok, lineModel} from './engine';
import {armColors} from './frames';
import {LineDensity, linePrior, tumourDensity} from './priors';
import {PLAN_TAG, SavedPlan, parsePlan, powerAt, savedPlan} from './plan-file';
import {PlanChecks, Range, planChecks} from './plan-checks';

/** `volume`: NaN where the row has none; `weight`: the body weight (g) of each row, null when the study has none;
 *  `mass`: the tumor's ex vivo weight (g) on the row of the day it was excised, null when the study has none;
 *  `skipped`: the table's rows left out for lacking their animal, group or day */
export interface StudyData {animal: string[], group: string[], day: number[], volume: number[],
  weight: (number | null)[] | null, mass: (number | null)[] | null, removal: string[] | null, skipped: number}

export interface AnalysisDeclaration {
  data: StudyData;
  /** the cell line: its model of control growth for the twin */
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
}

/** pmxAnalyzeStudy's inputs as the form and Grokky give them: columns, an optional plan file, empty optionals. */
export interface AnalysisInputs {
  study: DG.DataFrame; animal: DG.Column; group: DG.Column; day: DG.Column; volume: DG.Column;
  /** the text of a .pmxplan file */
  plan?: string | null; line: string; control: string; endpoint: number; correction: string;
  evaluationDay?: number; alpha?: number; bodyWeight?: DG.Column | null; weightLimit?: number; humane?: number;
  tumorWeight?: DG.Column | null; density?: number; removalReason?: DG.Column | null;
}

/** A column's values, an empty cell as null. */
const valuesOf = (c: DG.Column): any[] => Array.from({length: c.length}, (_, i) => c.isNone(i) ? null : c.get(i));

/** The declaration the engine reads: the rows with their animal, group and day, the control found when not named,
 *  the plan from its file or else the study table's tag, each empty optional at its default. */
export function analysisDeclarationOf(x: AnalysisInputs): AnalysisDeclaration {
  const [animal, group, day, volume] = [x.animal, x.group, x.day, x.volume].map(valuesOf);
  const keep = animal.map((a, i) => a !== null && group[i] !== null && ok(day[i]));
  const pick = <T>(xs: T[]) => xs.filter((_, i) => keep[i]);
  const groups = pick(group).map(String);
  const names = [...new Set(groups)];
  const control = x.control?.trim() || names.find((g) => /vehicle|control|pbs|saline|untreated/i.test(g)) ||
    names[0];
  const planText = x.plan?.trim() || x.study.getTag(PLAN_TAG);
  return {line: x.line || 'MC38', control, endpoint: x.endpoint, day: ok(x.evaluationDay) ? x.evaluationDay : null,
    plan: planText ? parsePlan(planText) : null, alpha: ok(x.alpha) ? x.alpha : 0.05,
    correction: correctionOf(x.correction),
    weightLimit: ok(x.weightLimit) ? x.weightLimit : 20, humane: ok(x.humane) ? x.humane : HUMANE,
    density: ok(x.density) && x.density > 0 ? x.density : null,
    data: {animal: pick(animal).map(String), group: groups, day: pick(day),
      volume: pick(volume).map((v) => ok(v) ? v : NaN),
      weight: x.bodyWeight ? pick(valuesOf(x.bodyWeight)) : null,
      mass: x.tumorWeight ? pick(valuesOf(x.tumorWeight)) : null,
      removal: x.removalReason ? pick(valuesOf(x.removalReason)).map((r) => r === null ? '' : String(r)) : null,
      skipped: keep.filter((k) => !k).length}};
}

/** The smallest volume a caliper records (mm3): the design's v_lod (stats.tumor_study_design's default). */
const V_LOD = 3.4;
const ADJUST: {[c in Correction]: string} = {dunnett: 'dunnettx', bonferroni: 'bonferroni', none: 'none'};
export const COMPARED = 'per-animal slope of log V from randomization (two-sample t test)';

export interface ArmAnalysis {
  arm: string;
  /** treated animals with a slope */
  n: number;
  /** percent, on the evaluation day */
  tgi: number | null;
  tgiLower: number | null;
  tgiUpper: number | null;
  tc: number | null;
  tcLower: number | null;
  tcUpper: number | null;
  /** mean slope of ln V per day, the arm's and the control's, and their difference with its (adjusted) interval */
  slope: number;
  slopeControl: number;
  difference: number;
  lower: number;
  upper: number;
  /** adjusted for the comparisons with the control */
  p: number;
  /** Kaplan-Meier median time to the endpoint volume (null: not reached) */
  medianTte: number | null;
  /** the twin: animals measured on the evaluation day, the median of their paired effects (percent) and how many
   *  are smaller than their twin beyond doubt (the 90 % interval of the effect above 0) */
  twinAnimals: number;
  pairedEffect: number | null;
  responders: number;
}

/** What the test on the slopes says of a treated group: it grows slower than the control (its interval below 0),
 *  faster (above 0), or the study cannot tell it from the control. */
export type Activity = 'Slows growth' | 'Grows faster' | 'No clear effect';
export const activityOf = (a: ArmAnalysis): Activity =>
  a.upper < 0 ? 'Slows growth' : a.lower > 0 ? 'Grows faster' : 'No clear effect';

/** A control whose measurements mostly fall outside its band from the other controls. */
export interface OddControl {animal: string, outside: number, measurements: number, below: boolean}

export interface StudyAnalysis {
  declaration: AnalysisDeclaration;
  treated: string[];
  evaluationDay: number;
  /** the evaluation day came from the endpoint (false: given, or the last day when the controls never reach it) */
  fromEndpoint: boolean;
  nControl: number;
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
    verdict: 'As predicted' | 'Earlier than predicted' | 'Later than predicted' | 'Endpoint not reached'};
  /** each treated group: the plan's arm it ran as, the effect assumed, the one observed with its 95% interval, the
   *  range a study run as planned shows, and the verdict */
  groups: {group: string, plannedAs: string | null, metric: EffectMetric, assumed: number | null,
    observed: number | null, lower: number | null, upper: number | null, expected: Range | null,
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

/** Whether x falls in a 95% range (an open top: no upper limit). */
const within = (x: number, r: Range) => x >= r.q025 && (r.q975 === null || x <= r.q975);

/** The median (the mean of the middle two for an even count). */
export const median = (xs: number[]) => {
  const s = [...xs].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};

/** Plan's evaluation day, observed: the first day the median control is measured at or above the endpoint; the
 *  last measurement day when it never is. */
function evaluationDay(d: AnalysisDeclaration): {day: number, fromEndpoint: boolean} {
  const {data} = d;
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
    if (median(byDay.get(day)!) >= d.endpoint)
      return {day, fromEndpoint: true};
  }
  return {day: Math.max(...data.day), fromEndpoint: false};
}

/** Crux stats.tumor_twin's model of the line: the engine priors' study prior, noise and start spread with the twin's
 *  own fit (between-animal spread, process noise, slow animals, skew). */
async function twinModel(line: string) {
  const lp = await linePrior(line);
  const [sp, nz, tw] = [lp.study_prior, lp.noise, lp.twin];
  if (!sp || !nz || !tw)
    throw new Error(`The engine priors have no digital twin model for ${line}`);
  return {mean_log_alpha: sp.mean_log_alpha, sd_study_log_alpha: sp.sd_log_alpha, mean_log_beta: sp.mean_log_beta,
    sd_study_log_beta: sp.sd_log_beta, rho_study: sp.rho, log_v_ref: lp.log_v_ref, start_sd: lp.start_sd,
    me_c0: nz.c0, me_c1: nz.c1, noise_df: nz.df, sd_session: nz.sd_session, v_lod: V_LOD,
    sd_log_alpha: tw.sd_log_alpha, sd_log_beta: tw.sd_log_beta, sd_process: tw.sd_process,
    first_factor: tw.first_factor, start_df: tw.start_df, slow_share: tw.slow_share, slow_shift: tw.slow_shift,
    skew: tw.skew};
}

export async function analyzeStudy(d: AnalysisDeclaration): Promise<StudyAnalysis> {
  const groups = [...new Set(d.data.group)];
  if (!groups.includes(d.control))
    throw new Error(`No group "${d.control}" in the study: the control must be one of ${groups.join(', ')}`);
  const treated = groups.filter((g) => g !== d.control);
  if (treated.length === 0)
    throw new Error('The study has no treated group');
  const end = d.day === null ? evaluationDay(d) : {day: d.day, fromEndpoint: false};

  // from randomisation to the evaluation day: the window the design's slope and effects are defined on
  const keep = d.data.day.map((t) => t >= 0 && t <= end.day);
  const pick = <T>(xs: T[]) => xs.filter((_, i) => keep[i]);
  const window = {animal: pick(d.data.animal), group: pick(d.data.group),
    day: Float64Array.from(pick(d.data.day)), volume: Float64Array.from(pick(d.data.volume))};
  const tgi = await crux('stats.tgi', window, {control: d.control, change: 'paired', rate_floor: V_LOD,
    endpoint_volume: d.endpoint, days: String(end.day)});

  // Hather's rate is the slope of log10 V: in ln V per day, the scale of the design's growth rates
  const animals = named(tgi, 'animals').map((a): {[k: string]: any} =>
    ({...a, slope: a.rate === null ? null : a.rate * Math.LN10}));
  const sloped = animals.filter((a) => a.slope !== null);
  // the groups coded g0 (the control), g1, ...: a contrast's label is then "gk - g0" whatever the groups' names
  const code = new Map([d.control, ...treated].map((g, k) => [g, `g${k}`]));
  const contrasts = rows((await crux('stats.contrast', {slope: Float64Array.from(sloped, (a) => a.slope),
    group: sloped.map((a) => code.get(a.group)!)}, {formula: 'slope ~ group', specs: 'group',
    contrast: 'trt.vs.ctrl', ref: 'g0', adjust: ADJUST[d.correction], alpha: d.alpha})).main as any);

  // each treated animal against its twin, each control against the others
  const tw = await crux('stats.tumor_twin', window, {...await twinModel(d.line), control: d.control, day: end.day});
  const twins = rows(tw.main as any);
  const effects = named(tw, 'effects');
  const checks = named(tw, 'checks');

  const table = rows(tgi.main as any);
  const delay = named(tgi, 'delay');
  const meanSlope = (g: string) => {
    const s = sloped.filter((a) => a.group === g).map((a) => a.slope as number);
    return s.reduce((t, x) => t + x, 0) / s.length;
  };
  const arms: ArmAnalysis[] = treated.map((g, k) => {
    const e = table.find((r) => r.group === g) ?? {};
    // a group with no animal measured twice has no contrast: no test (said in the notes)
    const c = contrasts.find((r) => String(r.term).startsWith(`g${k + 1} `)) ??
      {estimate: NaN, lower: NaN, upper: NaN, p_value: NaN};
    const own = effects.filter((r) => r.group === g);
    return {arm: g, n: sloped.filter((a) => a.group === g).length, tgi: finite(e.tgi_change),
      tgiLower: finite(e.tgi_change_lower), tgiUpper: finite(e.tgi_change_upper), tc: finite(e.tc),
      tcLower: finite(e.tc_lower), tcUpper: finite(e.tc_upper), slope: meanSlope(g), slopeControl: meanSlope(d.control),
      difference: c.estimate, lower: c.lower, upper: c.upper, p: c.p_value,
      medianTte: finite(delay.find((r) => r.group === g)?.median),
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
    ({group, day, n: v.length, median: median(v), lower: Math.min(...v), upper: Math.max(...v)}));

  const notes: string[] = [];
  if (d.data.skipped > 0) {
    notes.push(`${d.data.skipped} row${d.data.skipped > 1 ? 's' : ''} without an animal, a group or a day left out ` +
      'of the analysis');
  }
  if (d.day === null && !end.fromEndpoint) {
    notes.push(`The median control stays below ${d.endpoint} mm³: the effects are on the last measurement day, ` +
      `day ${end.day}`);
  }
  const unsloped = animals.length - sloped.length;
  if (unsloped > 0)
    notes.push(`${unsloped} animal${unsloped > 1 ? 's' : ''} with fewer than two measurements left out of the test`);
  for (const a of arms) {
    if (a.n < 3)
      notes.push(`${a.arm}: ${a.n} animal${a.n === 1 ? '' : 's'} with a slope, too few for a test`);
    if (a.twinAnimals < a.n) {
      notes.push(`${a.arm}: ${a.n - a.twinAnimals} animal${a.n - a.twinAnimals > 1 ? 's' : ''} not measured on day ` +
        `${end.day}, so not compared with ${a.n - a.twinAnimals > 1 ? 'their twins' : 'its twin'}`);
    }
  }
  for (const c of oddControls) {
    notes.push(`Control ${c.animal} grows ${c.below ? 'slower' : 'faster'} than the other controls predict: ` +
      `${c.outside} of its ${c.measurements} measurements fall outside their range. It stays in the analysis`);
  }
  const inside = checks.filter((r) => r.inside).length;
  const versusPlan = await compareWithPlan(d, treated, arms, sloped.filter((a) => a.group === d.control).length);
  return {declaration: d, treated, evaluationDay: end.day, fromEndpoint: end.fromEndpoint, versusPlan,
    nControl: sloped.filter((a) => a.group === d.control).length,
    controlMedianTte: finite(delay.find((r) => r.group === d.control)?.median), arms, animals, volumes,
    twins, effects, checks, controlCoverage: checks.length ? inside / checks.length : null, oddControls, notes};
}

/** The designs an analysis with no plan reads its study against: the same line, endpoint, volume, comparisons and
 *  animals give the same design, so a rerun (another correction, a tab opened again) does not search again. */
const noPlanDesigns = new Map<string, Promise<StudyDesign>>();

function noPlanDesign(d: StudyDeclaration): Promise<StudyDesign> {
  const key = JSON.stringify(d);
  let design = noPlanDesigns.get(key);
  if (!design) {
    design = designStudy(d);
    noPlanDesigns.set(key, design);
    design.catch(() => noPlanDesigns.delete(key));
  }
  return design;
}

/** Each treated group with the plan's arm it ran as: the same name (any case), else, for the groups left, the arms
 *  left in their order when they are as many (said in a note). */
function pairGroups(treated: string[], arms: SavedPlan['predicted']['arms'], notes: string[]) {
  const key = (s: string) => s.trim().toLowerCase();
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
  if (observed === null)
    return null;
  if (within(observed, expected))
    return activity === 'Slows growth' ? 'As planned' : 'Inconclusive';
  return (metric === 'TGI' ? observed > (expected.q975 ?? Infinity) : observed < expected.q025) ?
    'Stronger than planned' : 'Weaker than planned';
}

/** The study against its plan's predictive checks: the controls' growth, the end and each treated group's effect,
 *  each against the 95% range a study run as planned shows. With no plan, the controls and the end against the
 *  line's atlas: the checks of Plan's design for this line and endpoint, randomized at the controls' median on day
 *  0, with this study's controls (ponytail: on Mon / Wed / Fri; the study's own weekdays if the calendar matters). */
async function compareWithPlan(d: AnalysisDeclaration, treated: string[], arms: ArmAnalysis[],
  controls: number): Promise<PlanComparison> {
  const notes: string[] = [];
  let p = d.plan?.predicted;
  const mde = new Map<string, number | null>();
  if (!p) {
    const day0 = d.data.volume.filter((v, i) => d.data.group[i] === d.control && d.data.day[i] === 0 &&
      Number.isFinite(v));
    const vRand = day0.length ? median(day0) : 100;
    const design = (n: number | null) => noPlanDesign({line: d.line, vRand, endpoint: d.endpoint, power: 0.8,
      alpha: d.alpha, correction: d.correction, n, treatedArms: arms.length,
      arms: [{name: 'x', therapy: 'small_molecule', metric: 'TGI', effect: 60}],
      preLoss: 0, postLoss: 0, lossRule: 'confident', confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0,
      humane: HUMANE, schedules: false, lastDay: null});
    const s = await design(null);
    // a group the test cannot tell from the control: the smallest effect a study with its animals (and the study's
    // comparisons) detects, which only it shows (a search of seconds: not for the groups with an effect)
    for (const a of arms.filter((x) => activityOf(x) === 'No clear effect'))
      mde.set(a.arm, (await design(Math.max(3, a.n))).arms[0].mde);
    const checks = await planChecks({line: d.line, vRand, endpoint: d.endpoint, humane: HUMANE,
      weekdays: LAB_WEEKDAYS, startWeekday: 0, endpointDay: s.endpointDay, nControl: controls, arms: []});
    p = {line: d.line, endpoint: d.endpoint, vRand, studyEnd: s.endpointDay, endpointDay: s.endpointDay,
      durationQ10: s.durationQ10, durationQ90: s.durationQ90, nControl: controls, arms: [], checks};
  } else {
    if (p.line !== d.line)
      notes.push(`The plan is for ${p.line}; this analysis uses the ${d.line} model`);
    if (p.endpoint !== d.endpoint)
      notes.push(`The plan's endpoint is ${p.endpoint} mm³; this analysis uses ${d.endpoint} mm³`);
  }
  // the band starts from the plan's volume at randomization: a study randomized much larger or smaller grows apart
  // from it from day 0 (said, not corrected)
  const atStart = d.data.volume.filter((v, i) => d.data.group[i] === d.control && d.data.day[i] === 0 &&
    Number.isFinite(v));
  if (d.plan && atStart.length) {
    const start = median(atStart);
    if (start / p.vRand > 1.25 || start / p.vRand < 0.8) {
      notes.push(`The controls were randomized at ${Math.round(start)} mm³ (median), the plan at ${p.vRand} mm³: ` +
        'their growth is predicted from another volume');
    }
  }
  const c = p.checks;

  // the controls' growth: their mean slope of log V to the evaluation day (the test's own quantity)
  const slope = arms.length ? arms[0].slopeControl : null;
  const observedSlope = slope !== null && Number.isFinite(slope) ? slope : null;
  const controlSlope: PlanComparison['controlSlope'] = {observed: observedSlope, range: c.controlSlope,
    verdict: observedSlope === null ? null : within(observedSlope, c.controlSlope) ? 'As predicted' :
      observedSlope > (c.controlSlope.q975 ?? Infinity) ? 'Faster than predicted' : 'Slower than predicted'};

  // the first measurement day the median control reached the endpoint, whatever evaluation day was chosen; a study
  // that never did is later than planned only if it was measured beyond the range's top
  const reached = evaluationDay({...d, endpoint: p.endpoint});
  const lastDay = Math.max(...d.data.day);
  const end: PlanComparison['end'] = {observed: reached.fromEndpoint ? reached.day : null, range: c.end,
    notReached: c.notReached, verdict: !reached.fromEndpoint ?
      (c.end.q975 !== null && lastDay > c.end.q975 ? 'Later than predicted' : 'Endpoint not reached') :
      within(reached.day, c.end) ? 'As predicted' :
        reached.day < c.end.q025 ? 'Earlier than predicted' : 'Later than predicted'};

  const pairs = d.plan ? pairGroups(treated, p.arms, notes) : new Map<string, SavedPlan['predicted']['arms'][0]>();
  const groups = arms.map((a) => {
    const planned = pairs.get(a.arm);
    const metric: EffectMetric = planned?.metric ?? 'TGI';
    const [observed, lower, upper] = metric === 'TGI' ? [a.tgi, a.tgiLower, a.tgiUpper] :
      [a.tc, a.tcLower, a.tcUpper];
    const own = planned ? c.arms.find((x) => x.name === planned.name) : undefined;
    const expected = own ? (metric === 'TGI' ? own.tgi : own.tc) : null;
    return {group: a.arm, plannedAs: planned?.name ?? null, metric, assumed: planned?.effect ?? null, observed,
      lower, upper, expected, verdict: expected ? verdictOf(metric, observed, activityOf(a), expected) : null,
      activity: activityOf(a), p: a.p, n: a.n, power: planned?.power?.length ? powerAt(planned.power, a.n) : null,
      nNeeded: planned?.nNeeded ?? null, mde: mde.get(a.arm) ?? null};
  });
  // the plan's power, as its form gave it: a percent, or a share
  const target = Number(d.plan?.inputs.power);
  const powerTarget = Number.isFinite(target) && target > 0 ? (target > 1 ? target / 100 : target) : null;
  return {plan: d.plan?.name ?? null, controls: c.controls, controlSlope, end, groups, powerTarget, notes};
}

/** The analysis (one row), the comparisons (one row per treated arm), the animals, the observed volumes, and the
 *  twin's bands, paired effects and control check. */
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
async function tolerabilityFrames(d: AnalysisDeclaration): Promise<{tolerability: DG.DataFrame,
  weights: DG.DataFrame}> {
  const w = d.data.weight;
  if (!w?.some((x) => typeof x === 'number' && x > 0))
    return {tolerability: DG.DataFrame.create(0), weights: DG.DataFrame.create(0)};
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
  const removed = named(bw, 'animals').filter((a) => a.removal !== 'none');
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
  return {tolerability, weights};
}

/** The predicted control series and the observed one, as the band chart draws them. */
export const PREDICTED = 'Predicted controls';
export const OBSERVED = 'Observed controls';

export async function analysisFrames(d: AnalysisDeclaration): Promise<AnalysisFrames> {
  const s = await analyzeStudy(d);
  const {tolerability, weights} = await tolerabilityFrames(d);
  const vp = s.versusPlan;
  const analysis = DG.DataFrame.fromObjects([{
    line: d.line, control: d.control, n_control: s.nControl, treated_arms: s.treated.length,
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
    arm: a.arm, n: a.n, tgi: a.tgi, tgi_lower: a.tgiLower, tgi_upper: a.tgiUpper, tc: a.tc, tc_lower: a.tcLower,
    tc_upper: a.tcUpper,
    slope_per_day: a.slope, slope_control_per_day: a.slopeControl, slope_difference: a.difference,
    slope_difference_lower: a.lower, slope_difference_upper: a.upper, p_value: a.p, significant: a.p < d.alpha,
    activity: activityOf(a),
    median_time_to_endpoint: a.medianTte, twin_animals: a.twinAnimals, paired_effect: a.pairedEffect,
    responders: a.responders,
  })))!;
  const animals = DG.DataFrame.fromObjects(s.animals.map((a) => ({group: a.group, animal: a.animal,
    slope_per_day: a.slope, time_to_endpoint: a.time_to_endpoint, endpoint_reached: a.event, mrecist: a.mrecist})))!;
  const volumes = DG.DataFrame.fromObjects(s.volumes)!;
  const twins = DG.DataFrame.fromObjects(s.twins)!;
  const effects = DG.DataFrame.fromObjects(s.effects.map((e) => ({group: e.group, animal: e.animal,
    effect: 100 * e.effect, effect_lower: 100 * e.effect_lo, effect_upper: 100 * e.effect_hi,
    responder: e.effect_lo > 0})))!;
  const checks = DG.DataFrame.fromObjects(s.checks)!;
  // each treated group against the plan (with none: its effect, nothing assumed)
  const versusPlan = DG.DataFrame.fromObjects(vp.groups.map((g) => ({group: g.group, planned_as: g.plannedAs,
    metric: g.metric, assumed: g.assumed, observed: g.observed, lower: g.lower, upper: g.upper,
    expected_q025: g.expected?.q025 ?? null, expected_median: g.expected?.median ?? null,
    expected_q975: g.expected?.q975 ?? null, verdict: g.verdict, activity: g.activity, p_value: g.p, n: g.n,
    power: g.power, n_needed: g.nNeeded, mde: g.mde})))!;
  // a study's median control as predicted (its 95% range) to the evaluation day, then the observed median control
  const predicted = DG.DataFrame.fromObjects([
    ...vp.controls.filter((b) => b.day <= s.evaluationDay).map((b) => ({series: PREDICTED, day: b.day,
      lower: b.q025, upper: b.q975, median: b.median})),
    ...s.volumes.filter((v) => v.group === d.control).sort((a, b) => a.day - b.day).map((v) => ({series: OBSERVED,
      day: v.day, lower: v.median, upper: v.median, median: v.median}))])!;

  // the control in Plan's slate, the treated arms in its palette, by declaration order
  const colors = armColors(s.treated);
  colors[d.control] = colors['Control'];
  for (const [df, col] of [[comparisons, 'arm'], [animals, 'group'], [volumes, 'group'], [twins, 'group'],
    [effects, 'group'], [tolerability, 'group'], [weights, 'group']] as const)
    df.col(col)?.meta.colors.setCategorical(colors);
  predicted.col('series')!.meta.colors.setCategorical({[PREDICTED]: '#9aa5b1', [OBSERVED]: colors['Control']});
  for (const df of [volumes, twins, checks, predicted])
    df.col('day')!.meta.units = 'd';
  for (const [df, cols] of [[volumes, ['median', 'lower', 'upper']], [twins, ['observed', 'twin_q05', 'twin_median',
    'twin_q95']], [checks, ['observed', 'q05', 'median', 'q95']], [predicted, ['lower', 'upper', 'median']]] as const) {
    for (const q of cols) {
      df.col(q)!.meta.units = 'mm³';
      df.col(q)!.meta.format = '0.0';
    }
  }
  for (const q of ['effect', 'effect_lower', 'effect_upper'])
    effects.col(q)!.meta.units = '%';
  volumes.meta.formulaLines.addLine({title: `Endpoint ${d.endpoint} mm³`, formula: `\${median} = ${d.endpoint}`});
  volumes.meta.formulaLines.addLine({title: `Evaluation, day ${s.evaluationDay}`,
    formula: `\${day} = ${s.evaluationDay}`, style: 'solid', color: '#4d5261'});
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

/** The example's body weight, an illustrative shape (Crux stats.tumor_body_weight_simulate; the atlas has no body
 *  weights): per group its median loss at the nadir (%), the stronger arm losing more. */
const EXAMPLE_NADIR = [0, 4, 10];

/** An example study run as Plan designs it: MC38, a vehicle and two small-molecule arms (TGI 40% and 70%), each
 *  with the animals Plan asks for its power, measured 3 times a week up to a week past the design's end of study.
 *  One lab of the line's atlas, drawn with the design's seed (whatever it shows: nothing is picked by its outcome),
 *  as Crux stats.tumor_study_simulate measures it: one call per arm with the same seed gives the same lab and the
 *  same controls, and each treated group takes animals of its own (no shared mice). Each animal is weighed with
 *  every measurement: an illustrative body weight (EXAMPLE_NADIR), removed at the weight endpoint like a lab would;
 *  its tumor is weighed when excised, after its last measurement (tumor_weight, on that row). */
export async function exampleStudy(): Promise<DG.DataFrame> {
  const line = 'MC38';
  const arms = [40, 70].map((effect) => ({name: `Cpd TGI ${effect}%`, therapy: 'small_molecule' as const,
    metric: 'TGI' as const, effect}));
  const design = await designStudy({line, vRand: 100, endpoint: 1500, power: 0.8, alpha: 0.05, correction: 'dunnett',
    n: null, arms, preLoss: 0, postLoss: 0, lossRule: 'confident', confidence: 0.9, weekdays: LAB_WEEKDAYS,
    startWeekday: 0, humane: HUMANE, schedules: false, lastDay: null});
  const sizes = design.arms.map((a) => a.n);
  if (design.nControl === null || sizes.some((x) => x === null))
    throw new Error('The example needs more than 40 animals per arm');
  const nControl = design.nControl;
  const n = sizes as number[];
  // arm k's treated animals: the k-th run of the draws, so that no two groups share a mouse
  const offsets = n.map((_, k) => n.slice(0, k).reduce((t, x) => t + x, 0));
  const {prior, params} = await lineModel(line);
  // a week past the plan's end: the study is read on the first day its median control reaches the endpoint
  const shared = {...params, v_rand: 100, n: Math.max(nControl, offsets[n.length - 1] + n[n.length - 1]),
    last_day: design.studyEnd + 7, per_week: 3, study: 'prior', reps: 1, seed: SEED};
  const out: {animal: string, group: string, day: number, volume: number, body_weight: number,
    tumor_weight: number | null}[] = [];
  const density = (await tumourDensity()).lines[line].density;
  let g = 0;
  for (const [k, a] of design.arms.entries()) {
    const sim = rows((await crux('stats.tumor_study_simulate', prior, {...shared, kill_rate: a.killRate,
      onset: a.archetype.onsetDays, kappa_sd: a.archetype.kappaSd, responders: a.archetype.responders})).main as any);
    const take = (arm: string, from: number, count: number) =>
      sim.filter((r) => r.arm === arm && r.animal >= from && r.animal < from + count);
    const groups: [string, {[k: string]: any}[]][] = [[a.declaration.name, take('treated', offsets[k], n[k])]];
    if (k === 0)
      groups.unshift(['Vehicle', take('control', 0, nControl)]);
    // each group weighed with its measurements (its own seed), each tumor weighed once excised (the animal's last
    // measurement, at the line's density); an animal removed at the weight endpoint is no longer measured
    for (const [group, rs] of groups) {
      const animal = rs.map((r) => `${group} #${r.animal}`);
      const bw = rows((await crux('stats.tumor_body_weight_simulate', {animal, day: Float64Array.from(rs, (r) => r.day),
        volume: Float64Array.from(rs, (r) => r.volume)}, {nadir: EXAMPLE_NADIR[g], density, seed: SEED + g}))
        .main as any);
      const weighed = new Map(bw.map((r) => [`${r.animal}|${r.day}`, r]));
      rs.forEach((r, i) => {
        const w = weighed.get(`${animal[i]}|${r.day}`);
        if (w !== undefined) {
          out.push({animal: animal[i], group, day: r.day, volume: r.volume, body_weight: w.weight,
            tumor_weight: Number.isFinite(w.tumor_mass) ? w.tumor_mass : null});
        }
      });
      g++;
    }
  }
  const df = DG.DataFrame.fromObjects(out)!;
  df.name = `Example ${line} study`;
  // its plan in a tag: the plan Plan would save for it, with its predictive checks and each arm's power by n
  const checks = await planChecks({line, vRand: 100, endpoint: 1500, humane: HUMANE, weekdays: LAB_WEEKDAYS,
    startWeekday: 0, endpointDay: design.endpointDay, nControl, arms: design.arms.map((a, k) => ({
      name: a.declaration.name, therapy: a.declaration.therapy, killRate: a.killRate, n: n[k]}))});
  const plan = savedPlan(`Example ${line} plan`, {line, endpoint: 1500, arms, power: 80, correction: 'dunnett',
    n: null, vRand: 100, measurementDays: 'Mon, Wed, Fri', startWeekday: 'Monday', humane: HUMANE, alpha: 0.05,
    preLoss: 0, postLoss: 0}, {line, endpoint: 1500, vRand: 100, studyEnd: design.studyEnd,
    endpointDay: design.endpointDay, durationQ10: design.durationQ10, durationQ90: design.durationQ90, nControl,
    arms: arms.map((a, k) => ({...a, n: n[k], nNeeded: design.arms[k].nAssurance,
      power: design.arms[k].power.flatMap((r) => r.assurance === null ? [] : [{n: r.n, assurance: r.assurance}])})),
    checks});
  df.setTag(PLAN_TAG, JSON.stringify(plan));
  df.col('day')!.meta.units = 'd';
  df.col('volume')!.meta.units = 'mm³';
  df.col('volume')!.meta.format = '0.0';
  df.col('body_weight')!.meta.units = 'g';
  df.col('body_weight')!.meta.format = '0.0';
  df.col('tumor_weight')!.meta.units = 'g';
  df.col('tumor_weight')!.meta.format = '0.00';
  return df;
}
