/* A tumor model the atlas has never seen, planned from its own history: an implantation pilot, the provider's pilot,
   or earlier studies of the same model. The history is weighed against PMX Studio's PDX reference
   (files/pdx_history_v1.json, priors.ts pdxHistoryAnchor): what the user is told, each study's growth rate, the
   weighing, and the growth prior the engine takes. Reference implementation with its self-test:
   preclinical_bundle/analysis/out_pdx/d1_r7_rerun/history_anchor.py (update). */
import type {HistoryAnchor, LinePrior} from './priors';
import {mean} from './format';

export type Answer = 'yes' | 'no' | 'unknown';

/** What the user says about where the history comes from (the dialog's questions). */
export interface HistorySource {
  /** run in the lab that will run the new study */
  sameLab: Answer;
  /** the same mouse strain and implant site as the new study */
  sameHost: Answer;
  /** the passage the history was measured at: P0–P1, P2 or later */
  passage: 'early' | 'later' | 'unknown';
  /** the tumors came from frozen material */
  frozen: boolean;
  /** group means without animals (a provider's sheet); `animalsPerGroup` null: not stated */
  meansOnly: boolean;
  animalsPerGroup: number | null;
}

/** The history counts as another lab's (its lab's shift adds to the ranges) unless the user says it is this lab's,
 *  with the same strain and implant site (D694: a strain or site not known counts as another lab's). Validated
 *  (D685): a same-lab history, by an independent local contrast (D681) and the real engine's contrast (D685). Not
 *  validated: the lab's shift tau_L 0.15 of another lab's history, which no data at hand can check; it keeps its own
 *  note. */
export const otherLab = (s: HistorySource) => s.sameLab !== 'yes' || s.sameHost !== 'yes';

/** Days the history's sustained growth was validated for (D681/D682): longer plans get a note. */
export const HORIZON_DAYS = 28;

/** Every note Plan, the protocol and the report show for a plan made from a model's history, in that order: where
 *  the growth comes from (`weight`: the history's share against the reference, 0–1), then why the ranges are wider,
 *  then what may make the model grow otherwise. `nAssumed`: the animals per group assumed for means without them. */
export function historyNotes(s: HistorySource, studies: number, weight: number, nAssumed: number,
  plannedDays?: number): string[] {
  const notes = [`Growth from this model’s own history (${studies} ${studies === 1 ? 'study' : 'studies'}), ` +
    `weighted ${Math.round(100 * weight)}% against PMX Studio’s PDX reference.`];
  if (s.sameLab === 'no')
    notes.push('The history comes from another lab, so the ranges are wider.');
  else if (s.sameLab === 'unknown')
    notes.push('It is not known whether the history comes from this lab, so the ranges are wider.');
  if (s.sameHost === 'no')
    notes.push('The history used another mouse strain or implant site, so the ranges are wider.');
  else if (s.sameHost === 'unknown') {
    notes.push('It is not known whether the history used the same mouse strain and implant site, so the ranges are ' +
      'wider.');
  }
  if (s.passage === 'early')
    notes.push('The history is from an early passage (P0–P1): the tumor usually grows faster later.');
  if (s.frozen)
    notes.push('The tumors came from frozen material: their first growth may be slower.');
  if (s.meansOnly && s.animalsPerGroup === null)
    notes.push(`The history gives group means without the number of animals: ${nAssumed} per group assumed.`);
  if ((plannedDays ?? 0) > HORIZON_DAYS) {
    notes.push('The model’s growth is taken from its history and assumed to keep its pace. That has been checked ' +
      `for studies of up to about ${HORIZON_DAYS / 7} weeks; in a longer study the controls may end up smaller ` +
      'than shown.');
  }
  return notes;
}

/** A measurement of the history's control group: an animal's, or (group means) the experiment's mean curve. */
export interface HistoryRow {animal: string, day: number, volume: number}

/** One study of the history as the weighing takes it: ln of its control growth rate (1/day) and its animals. */
export interface HistoryStudy {lnRate: number, n: number}

const GAP = 21; // days: a curve is cut at its first gap longer than this (the atlas's recipe, explore_tau_t100.py)
const MIN_POINTS = 4;
const MIN_RATE = 1e-3; // per day: slower is no growth (no take), left out as in the atlas

/** The exponential growth rate of one curve: the least-squares slope of ln V on day over its positive points, cut at
 *  the first gap > GAP days; null with fewer than MIN_POINTS points or a rate <= MIN_RATE. */
export function curveRate(points: {day: number, volume: number}[]): number | null {
  const p = points.filter((x) => Number.isFinite(x.day) && x.volume > 0).sort((a, b) => a.day - b.day);
  const cut = p.findIndex((x, i) => i > 0 && x.day - p[i - 1].day > GAP);
  const q = cut < 0 ? p : p.slice(0, cut);
  if (q.length < MIN_POINTS)
    return null;
  const [t, y] = [q.map((x) => x.day), q.map((x) => Math.log(x.volume))];
  const [mt, my] = [mean(t), mean(y)];
  const sxx = t.reduce((s, x) => s + (x - mt) ** 2, 0);
  const rate = t.reduce((s, x, i) => s + (x - mt) * (y[i] - my), 0) / sxx;
  return sxx > 0 && rate > MIN_RATE ? rate : null;
}

/** A study's control group as the weighing takes it: the mean over its animals of ln of each one's rate (the atlas's
 *  per-animal recipe; ponytail: without the atlas's fit-quality filter, RMSE <= 0.30 of the best growth model). Group
 *  means (`meansOnly`): each "animal" is one experiment's mean curve, its n `animalsPerGroup` (null: `nAssumed`). Null
 *  when no curve has a rate. */
export function historyStudy(rows: HistoryRow[], meansOnly: boolean, animalsPerGroup: number | null,
  nAssumed: number): HistoryStudy | null {
  const byAnimal = new Map<string, HistoryRow[]>();
  for (const r of rows)
    byAnimal.set(r.animal, [...byAnimal.get(r.animal) ?? [], r]);
  const ln = [...byAnimal.values()].map(curveRate).filter((r): r is number => r !== null).map(Math.log);
  return ln.length ? {lnRate: mean(ln), n: meansOnly ? animalsPerGroup ?? nAssumed : ln.length} : null;
}

/** The history weighed against the PDX reference, in ln rate (the scale d1 (a) validated): the model's own mean
 *  u ~ N(mu, tau_model²); each study ~ N(u, tau_study² + sigma²/n). Another lab's history shares one lab shift
 *  (tau_lab², counted once) and the new study has its own. Returns the next study's ln rate (mean, sd: between
 *  studies, without the animals, which the engine draws itself) and the history's weight. ponytail: closed form, not
 *  Crux stats.map_prior (its fixed tau covers the same lab; the shared lab shift is not expressible there). */
export function weighHistory(a: HistoryAnchor, studies: HistoryStudy[], otherLab: boolean):
  {mean: number, sd: number, weight: number} {
  const lab = otherLab ? a.tau_lab ** 2 : 0;
  if (!studies.length)
    return {mean: a.mu, sd: Math.sqrt(a.tau_model ** 2 + a.tau_study ** 2), weight: 0};
  const prec = studies.map((s) => 1 / (a.tau_study ** 2 + a.sigma_animal ** 2 / s.n));
  const total = prec.reduce((t, x) => t + x, 0);
  const ybar = studies.reduce((t, s, i) => t + prec[i] * s.lnRate, 0) / total;
  const v = 1 / total + lab;
  const weight = a.tau_model ** 2 / (a.tau_model ** 2 + v);
  const post = 1 / (1 / a.tau_model ** 2 + 1 / v);
  return {mean: a.mu + weight * (ybar - a.mu), sd: Math.sqrt(post + a.tau_study ** 2 + lab), weight};
}

/** Each numeric leaf the mean over `xs` (objects recursed; anything else the first's). */
function meanOf<T>(xs: T[]): T {
  const first = xs[0] as any;
  if (typeof first === 'number')
    return mean(xs as number[]) as T;
  if (!first || typeof first !== 'object' || Array.isArray(first))
    return first;
  return Object.fromEntries(Object.keys(first).map((k) => [k, meanOf(xs.map((x: any) => x?.[k]))])) as T;
}

// (Converting the weighed rate to the engine's log alpha, alpha = rate + beta c, was REJECTED by an independent local
// contrast, D681, as was the delta method before it, D674: that route is closed and its code, never called, is in
// git history, rateToLogAlpha. historyPrior uses the exponential.)

/** log beta of a history's growth: a deceleration so small (about 8e-7 per day) that the engine's reduced Gompertz
 *  is the exponential the history measured, ln V(t) = ln V0 + alpha t (it departs by under 1e-4 in ln V over four
 *  weeks); Crux needs beta > 0, and its closed form stays accurate down to beta ~ 1e-12. */
export const HISTORY_LOG_BETA = -14;
/** The study sd of that log beta: a small positive value, since Crux tumor_twin refuses 0 (D683); Plan's growth stays
 *  as it was. Whether the twin fits well with beta held this low is not checked: 2e's Docker test before Monitor and
 *  Analyze take a model with a history. */
export const HISTORY_SD_LOG_BETA = 0.01;

/** The growth prior of a model planned from its history (D682): the history's sustained exponential growth, as
 *  validated by an independent local contrast (D681). A new study's log alpha is the weighed ln rate itself
 *  (alpha = rate, beta ~ 0: log_v_ref no longer matters), its between-animal sd the reference's sigma_animal; the
 *  measurement error and start spread are the atlas PDX lines' mean, without their shift per session (the validated
 *  predictive had none; it cancels between arms). The non-progressing share stays the lines': it sizes the animals
 *  to implant, outside what the validation measured. Holds for studies of up to about four weeks (historyNotes says
 *  so past 28 days). Plan, Monitor and Analyze find it by the model's key (history-models.ts). */
export function historyPrior(lines: LinePrior[], a: HistoryAnchor, w: {mean: number, sd: number},
  animals: number): LinePrior {
  if (!lines.length || lines.some((l) => !l.study_prior || !l.noise))
    throw new Error('PMX Studio has no PDX growth reference to weigh this history against.');
  const base = {...lines[0], ...meanOf(lines.map((l) => ({log_v_ref: l.log_v_ref, sd_within: l.sd_within,
    noise: l.noise, twin: l.twin, take: l.take, start_sd: l.start_sd, sd_residual_log_v: l.sd_residual_log_v})))};
  const [la, sd] = [w.mean, w.sd];
  const twin = base.twin ? {...base.twin} : undefined;
  delete twin?.conformal; // the atlas lines' validations, not this model's
  delete twin?.next_ranges;
  return {...base, papers: 0, animals, log_beta: HISTORY_LOG_BETA, twin, calibration: undefined,
    sd_within: {...base.sd_within, log_alpha: a.sigma_animal, log_beta: 0}, facet_shifts_log_rate: [],
    noise: {...base.noise!, sd_session: 0},
    study_prior: {mean_log_alpha: la, sd_log_alpha: sd, mean_log_beta: HISTORY_LOG_BETA,
      sd_log_beta: HISTORY_SD_LOG_BETA, rho: 0},
    log_alpha: {...base.log_alpha, mixture: [{kind: 'history', w: 1, a: la, b: sd}], map_mean: la, map_sd: sd,
      map_q025: la - 1.959964 * sd, map_q975: la + 1.959964 * sd}};
}
