/*
 * The Plan engine: Crux's stats.tumor_study_design (crate tgi, [decision:tumor-study-design]), called once per
 * treated arm. This file only translates and combines: the cell line's prior (files/engine_priors_v3.json: the labs'
 * joint distribution of a new study's log alpha and log beta, and the heavy-tailed observation error) and the arm's
 * therapy (archetypes.ts) become the
 * function's table and parameters; every arm is designed with the same seed, so the studies drawn, the study end
 * and the controls are the same in every call, and the study's animals per arm are the largest arm's.
 */
import * as DG from 'datagrok-api/dg';

import {crux, CruxResult, MAX_SETS, named, rows} from './crux';
import {MeasurementError, PriorFit, enginePriors, linePrior} from './priors';
import {ARCHETYPES, Archetype, Therapy} from './archetypes';

export type EffectMetric = 'TGI' | 'TC';
export type Correction = 'holm' | 'bonferroni' | 'none';

/** A form's correction choice, Holm when it is none of them. */
export const correctionOf = (x: unknown): Correction =>
  x === 'bonferroni' || x === 'none' ? x : 'holm';

/** The test of Plan, Monitor and Analyze: each treated group against the control with their own variances (Welch's t
 *  test on the per-animal slopes). With one variance for every group, an arm less spread than another lost 13 to 38
 *  points of the power its plan gave it (preclinical bundle, analysis/test_unequal_spread.py). */
export const TEST = 'welch';

/** Whether x is a finite number; `finite` gives it or null (Crux's nulls, a form's empty optionals, and a column's
 *  empty cell where the platform gives its sentinel and not null: INT_NULL, FLOAT_NULL). */
export const isFiniteNumber = (x: unknown): x is number => typeof x === 'number' && Number.isFinite(x);
export const finite = (x: unknown): number | null =>
  isFiniteNumber(x) && x !== DG.INT_NULL && x !== DG.FLOAT_NULL ? x : null;
/** The median (the mean of the middle two for an even count). */
export const median = (xs: number[]) => {
  const s = [...xs].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};

/** Monday, Wednesday and Friday (0 Monday … 6 Sunday): the measurement weekdays when none are given. */
export const LAB_WEEKDAYS = [0, 2, 4];

/** A treated arm as the user declares it. */
export interface ArmDeclaration {
  name: string;
  therapy: Therapy;
  metric: EffectMetric;
  /** percent: TGI or T/C on the endpoint day (the median control reaches the endpoint) */
  effect: number;
}

export interface StudyDeclaration {
  line: string;
  /** mm3 at randomisation (start of treatment) */
  vRand: number;
  /** mm3: by default the study ends when the median control reaches it (single animals may pass it) */
  endpoint: number;
  /** mm3: the welfare limit; an animal measured above it is removed, whatever the study's day */
  humane: number;
  /** the study's last day, days from randomisation, when chosen (null: the endpoint day) */
  lastDay: number | null;
  /** the lab's measurement weekdays, 0 Monday … 6 Sunday (at least one) */
  weekdays: number[];
  /** day of the week of randomization (0 Monday … 6 Sunday): the measurement weekdays stay the lab's */
  startWeekday: number;
  power: number;
  /** two-sided significance level of the family of comparisons with the control */
  alpha: number;
  correction: Correction;
  /** analysable animals per arm when fixed by the user (its power and minimum detectable effect) */
  n: number | null;
  arms: ArmDeclaration[];
  /** treated arms the correction counts, when not the arms designed (one arm standing for a study's several) */
  treatedArms?: number;
  /** share of implanted animals not randomised: tumours that do not take or fall outside the randomisation window */
  preLoss: number;
  /** share of randomised animals not analysable: ulceration, deaths, technical failures (at least the atlas's
   *  non-progressing share) */
  postLoss: number;
  /** probability that the animals planned leave a group its analysable n: each group's own, and the implanted
   *  pool's, not that of all of them at once (three groups at 0.9 keep all their animals about 0.7 of the time) */
  confidence: number;
  /** 'confident': enough with that probability (binomial); 'expected': the average loss only, n / (1 - loss) */
  lossRule: 'confident' | 'expected';
  /** whether to scan the measurement schedules (the costliest alternative; off: none, the rest the same) */
  schedules: boolean;
  /** false: no search of the minimum detectable effect at a fixed n (the costliest step; none, the rest the same) */
  mde?: boolean;
  /** false: the groups are not sized under the line's other fits (a caller that reads the kill rates alone) */
  fits?: boolean;
}

/** The pools the trials draw from: 100 animals per arm and study (Crux's defaults: 1000, and 300 by last day and
 *  schedule). A random subset of a pool is a sample of the study, so the power is unbiased at any size that fills n;
 *  measured over 20 seeds and four cases (Crux spike tumor-design-shards, 2026-11-04): the assurance's sd at the
 *  reference n 0.0090-0.0112 at 1000, 0.0095-0.0113 at 100, n as often the reference's, the design about half the
 *  time. More studies do not help (the noise floor is elsewhere); 60 starts to show. The kill rate is calibrated over
 *  every study drawn (400 labs, Crux's default 200, the same 8000 animals): its sd over seeds 24-35 % lower in all
 *  four cases, at no cost. */
const POOLS = {animals: 100, mde_animals: 100, cal_studies: 400};

/** One arm's design. */
export interface ArmDesign {
  declaration: ArmDeclaration;
  archetype: Archetype;
  killRate: number;
  /** the effect the kill rate gives (percent, as declared) */
  effect: number;
  /** analysable animals this arm needs for the target assurance (null: not up to 40); on a line with several fits
   *  of its growth spreads, for the assurance averaged over them */
  nAssurance: number | null;
  /** what each fit that describes the atlas's controls as well needs alone: the prior's own first, alone unless the
   *  line has several (the engine priors' calibration.identification) */
  fits: (number | null)[];
  /** the arm's analysable animals (the user's n when given, else nAssurance, at least MIN_GROUP) and the animals to
   *  randomise for them */
  n: number | null;
  randomise: number | null;
  /** at the arm's n */
  assurance: number | null;
  /** minimum detectable effect at the user's n (percent; null without one, or out of reach) */
  mde: number | null;
  /** null where Crux gives no power */
  power: {n: number, implant: number, assurance: number | null}[];
}

export interface StudyDesign {
  declaration: StudyDeclaration;
  /** the study's last day: the one chosen, else the first measurement day from the endpoint day */
  studyEnd: number;
  /** the median over the atlas's labs of the day their controls reach the endpoint: the day the declared effects
   *  refer to */
  endpointDay: number;
  durationQ10: number;
  /** null: the controls of some plausible labs level off below the endpoint */
  durationQ90: number | null;
  pNotReached: number;
  /** the control's analysable animals: the largest arm's (null: an arm needs more than 40), and its randomised ones */
  nControl: number | null;
  randomiseControl: number | null;
  /** every group's analysable and randomised animals, summed (null: some group has no answer) */
  nTotal: number | null;
  randomiseTotal: number | null;
  /** animals to implant so that the randomised ones are there: one pool for every group */
  implant: number | null;
  /** the chance that the implanted animals give the randomised ones and every group then keeps its analysable
   *  animals, all at once (null: not sized) */
  confidenceStudy: number | null;
  /** the post-randomisation loss used: the declared one, or the atlas's non-progressing share when larger */
  postLoss: number;
  /** per-comparison level of the largest arm */
  alphaArm: number;
  calendar: number[];
  arms: ArmDesign[];
  /** day, arm ('Control' or the arm's name), volume quantiles (mm3) */
  curves: {day: number, arm: string, q05: number, q25: number, median: number, q75: number, q95: number}[];
  measurement: MeasurementError;
  nonprogressing: number;
  atlas: {papers: number, animals: number, ess: number};
  /** per possible last day: each arm's power (arms with an n) and the share of the controls past the endpoint and
   *  over the humane volume */
  duration: {day: number, assurance: {[arm: string]: number | null}, pastEndpoint: number, pastHumane: number}[];
  /** per measurements a week: the measurement days, each arm's power (arms with an n) and the analysable animals the
   *  study needs with that schedule (null: some arm needs more than 40) */
  schedules: {perWeek: number, measurements: number, assurance: {[arm: string]: number | null},
    nTotal: number | null}[];
  /** warnings about the groups' sizes, and about the study's end (the controls' growth) */
  notes: {groups: string[], growth: string[]};
}

export const SEED = 20261021;

/** The lab's measurement weekdays (0 Monday … 6 Sunday) as Crux's `weekdays` parameter: "mon,wed,fri". */
export const weekdaysParam = (days: number[]) =>
  days.map((i) => ['mon', 'tue', 'wed', 'thu', 'fri', 'sat', 'sun'][i]).join(',');

/** mm3: where a study starts until the lab sets its own limit (study-file.ts labHumane): the humane limit under the
 *  guideline's, the endpoint below it. Against 1,500 and 2,000, a default study needs no more animals and ends 2 to 9
 *  days sooner (preclinical bundle, analysis/humane_default_impact.py). */
export const HUMANE = 1500;
export const ENDPOINT = 1000;
/** mm3: Workman et al. 2010 (Br J Cancer 102:1555) advise a mean diameter of 1.5 cm at most in therapy studies; as a
 *  volume by length x width^2 / 2. */
export const GUIDE_LIMIT = 1690;
/** The mean diameter (cm) of a tumor of `volume` mm3, by length x width^2 / 2 with both equal. */
export const diameterOf = (volume: number) => Math.cbrt(2 * volume) / 10;

/** Fewest analyzable animals in any group, whatever the power needs: with fewer, losing one animal can leave a
 *  comparison without a usable variance (a design rule; usual xenograft groups are 8-10). */
export const MIN_GROUP = 5;

/** The fewest animals that leave at least `n` with probability `confidence` when each stays with probability `keep`:
 *  P(Binom(m, keep) >= n) >= confidence, i.e. n + the confidence quantile of the failures before the n-th success,
 *  NB(n, keep) (Crux stats.dist_eval: scipy's nbinom). */
async function enough(n: number | null, keep: number, confidence: number,
  rule: StudyDeclaration['lossRule'] = 'confident'): Promise<number | null> {
  if (n === null)
    return null;
  if (keep >= 1)
    return n;
  // the classic correction: the expected number left is n (a group can end up short)
  if (rule === 'expected')
    return Math.ceil(n / keep - 1e-9);
  const r = await crux('stats.dist_eval', {x: Float64Array.from([confidence])},
    {dist: 'nbinom', func: 'ppf', shape1: n, shape2: keep});
  return n + Number(r.main.col['value'][0]);
}

/** A line's growth model as Crux stats.tumor_study_design / tumor_study_simulate take it: the table of a new study's
 *  log alpha (one component, the labs' bivariate normal) and the model parameters. `fit`: another fit of the line's
 *  growth spreads in place of the prior's own (a line the atlas's controls do not pin down has several). */
export async function lineModel(line: string, fit?: PriorFit) {
  const [lp, atlasMe] = [await linePrior(line), (await enginePriors()).measurement_error];
  const nz = lp.noise;
  if (!lp.study_prior || !nz)
    throw new Error(`The engine priors have no study prior for ${line}`);
  const [sp, within] = [{...lp.study_prior, ...fit?.study_prior}, fit?.sd_within ?? lp.sd_within];
  const prior = {w: Float64Array.of(1), a: Float64Array.of(sp.mean_log_alpha), b: Float64Array.of(sp.sd_log_alpha)};
  const me: MeasurementError = {...atlasMe, model: `${atlasMe.model}, Student t (df ${nz.df.toFixed(1)})`,
    c0: nz.c0, c1: nz.c1};
  const params = {
    log_beta: sp.mean_log_beta, log_v_ref: lp.log_v_ref, sd_log_alpha: within.log_alpha,
    sd_log_beta: within.log_beta, sd_study_log_beta: sp.sd_log_beta,
    rho_study: sp.rho, start_sd: lp.start_sd, nonprogressing: lp.take.nonprogressing_map,
    me_c0: me.c0, me_c1: me.c1, noise_df: nz.df, sd_process: 0,
  };
  return {lp, me, prior, params};
}

export interface Scenario {tgi: number, design: StudyDesign}

/** The study sized for each assumed TGI (percent), every declared arm at that TGI: what the design needs if the effect
 *  is that large, which says nothing about whether a treatment reaches it. */
export function designScenarios(d: StudyDeclaration, tgis: number[]): Promise<Scenario[]> {
  return Promise.all(tgis.map(async (tgi) => ({tgi, design: await designStudy({...d, n: null, schedules: true,
    mde: false, arms: d.arms.map((a) => ({...a, metric: 'TGI' as const, effect: tgi}))})})));
}

/** Calls one design is split into, each over its block of the studies drawn on a worker set of its own (crux.ts opens
 *  as many): wasm runs one thread per worker, so the shards are the design's parallelism. */
const SHARDS = MAX_SETS;

type Row = {[k: string]: any};
interface ArmTables {plan: Row, power: Row[], duration: Row[], schedules: Row[], arms: Row[]}

/** One arm's stats.tumor_study_design: one call, or `shards` calls (Crux's shard / shards) whose partial sums give
 *  the assurance curves the single call would; the rest of every shard's design is the single call's, the conditional
 *  power shard 0's. The n chosen, its level and its power by last day and schedule follow from the summed curves. */
function designArm(prior: object, params: Row, shards: number, power: number,
  nFixed: number | null): Promise<ArmTables> {
  // the seed is fixed: the same inputs give the same tables, so a change recomputes only the arms it touches
  const key = JSON.stringify([prior, params, shards, power, nFixed],
    (_, v) => ArrayBuffer.isView(v) ? Array.from(v as any) : v);
  let r = armCache.get(key);
  if (r)
    armCache.delete(key); // re-inserted below: most recently used last
  else {
    const mine = r = designArmNow(prior, params, shards, power, nFixed);
    r.catch(() => armCache.get(key) === mine && armCache.delete(key));
    // ponytail: oldest-used out at 512 arms (three arms on a line with eleven fits, in a plan and its eight
    // scenarios, are 3 x 11 x 9 = 297); some MB of tables at most
    if (armCache.size >= 512)
      armCache.delete(armCache.keys().next().value!);
  }
  armCache.set(key, r);
  return r;
}

/** designArm's results by their inputs (read only by designStudy: never mutated once resolved). */
const armCache = new Map<string, Promise<ArmTables>>();

async function designArmNow(prior: object, params: Row, shards: number, power: number,
  nFixed: number | null): Promise<ArmTables> {
  const calls: CruxResult[] = await Promise.all(Array.from({length: shards}, (_, shard) =>
    crux('stats.tumor_study_design', prior, shards > 1 ? {...params, shard, shards} : params)));
  const r = calls[0];
  const t: ArmTables = {plan: rows(r.main as any)[0], power: named(r, 'power'), duration: named(r, 'duration'),
    schedules: named(r, 'schedules'), arms: named(r, 'arms')};
  if (shards === 1)
    return t;
  // per (part, key, n): the power summed over every shard's studies, and how many
  const sums = new Map<string, {sum: number, count: number}>();
  for (const c of calls) {
    for (const p of named(c, 'partials')) {
      const k = `${p.part} ${p.key} ${p.n}`;
      const s = sums.get(k) ?? {sum: 0, count: 0};
      sums.set(k, {sum: s.sum + p.sum, count: s.count + p.count});
    }
  }
  const mean = (part: string, key: number | null, n: number) => {
    const s = sums.get(`${part} ${key} ${n}`);
    return s && s.count > 0 ? s.sum / s.count : null;
  };
  const reaches = (part: string, key: number | null) =>
    t.power.find((x) => (mean(part, key, x.n) ?? -1) >= power)?.n ?? null;
  for (const x of t.power)
    x.assurance = mean('assurance', null, x.n);
  const nAssurance = reaches('assurance', null);
  const chosen = nFixed ?? nAssurance ?? t.power[t.power.length - 1].n;
  // (a given n outside the power table, 2 or over 40, keeps the α its own call gave)
  const alphaArm = t.power.find((x) => x.n === chosen)?.alpha ?? t.plan.alpha_arm;
  t.plan = {...t.plan, n_assurance: nAssurance, alpha_arm: alphaArm,
    ...(nFixed ? {assurance_at_n: mean('assurance', null, nFixed)} : {})};
  for (const x of t.duration)
    x.assurance = mean('duration', x.day, chosen);
  for (const x of t.schedules)
    Object.assign(x, {assurance: mean('schedule', x.per_week, chosen), n_assurance: reaches('schedule', x.per_week)});
  return t;
}

/** One arm's stats.tumor_study_design parameters in study `d`, over the line's model parameters `model`. */
function armCallParams(d: StudyDeclaration, arm: ArmDeclaration, model: Row): Row {
  const a = ARCHETYPES[arm.therapy];
  return {...model, v_rand: d.vRand,
    endpoint: d.endpoint, power: d.power, alpha: d.alpha, treated_arms: d.treatedArms ?? d.arms.length,
    correction: d.correction, test: TEST,
    weekdays: weekdaysParam(d.weekdays),
    // (the smallest detectable effect is searched only at a given n: without one the flag would only keep Plan's
    // designs and Monitor's or Analyze's of the same study apart in the cache)
    start_weekday: d.startWeekday, humane: d.humane, schedules: d.schedules, mde: !!d.n && d.mde !== false,
    seed: SEED,
    ...POOLS, ...(d.n ? {n: d.n} : {}), ...(d.lastDay ? {last_day: d.lastDay} : {}),
    metric: arm.metric.toLowerCase(), effect: arm.effect, onset: a.onsetDays, kappa_sd: a.kappaSd,
    responders: a.responders};
}

/** The share of studies whose controls may level off below the endpoint before Plan warns of it. */
export const ENDPOINT_RISK = 0.02;

/** The highest round endpoint below the declared one (50 mm³ steps, from 90 % of it down to 40 %) at which the
 *  controls of at most `limit` of the studies level off before reaching it; null: none. The share is the control's,
 *  whatever the arm: one call per candidate, the first arm's at a given n (no n search, no shards), in parallel. */
export async function safeEndpoint(d: StudyDeclaration, limit = ENDPOINT_RISK): Promise<number | null> {
  const step = 50;
  const xs = [...new Set([0.9, 0.8, 0.7, 0.6, 0.5, 0.4].map((f) => Math.floor(d.endpoint * f / step) * step))]
    .filter((x) => x > d.vRand && x < d.endpoint);
  if (!xs.length || !d.arms.length)
    return null;
  const {prior, params} = await lineModel(d.line);
  const ps = await Promise.all(xs.map((x) => designArm(prior,
    armCallParams({...d, endpoint: x, n: MIN_GROUP, schedules: false, mde: false}, d.arms[0], params), 1, d.power,
    MIN_GROUP).then((t) => Number(t.plan['p_not_reached']))));
  return xs.find((_, i) => ps[i] <= limit) ?? null;
}

export async function designStudy(d: StudyDeclaration): Promise<StudyDesign> {
  if (d.arms.length === 0)
    throw new Error('Add at least one treatment arm');
  const {lp, me, prior, params} = await lineModel(d.line);
  const armParams = (arm: ArmDeclaration, model: Row) => armCallParams(d, arm, model);
  const armError = (arm: ArmDeclaration) => (e: any) => {
    throw new Error(`${arm.name}: ${e?.message ?? e}`);
  };
  // the minimum detectable effect searches over every study at once: such a design is one call; else the worker sets
  // shared out over the arms (each shard repeats its arm's calibration: more shards than sets only queue)
  const shards = d.n && d.mde !== false ? 1 : Math.max(1, Math.floor(SHARDS / d.arms.length));
  // a line whose growth spreads the published controls do not pin down: each arm is also designed under every other
  // fit that describes them as well, for its power by n alone (one call each: no schedules, no n, no detectable
  // effect, so a fixed n reads the same calls). The fits count only with the prior's own marked among them: without
  // it, which one the prior carries is not known. Every fit, not some: three of Panc02's eleven gave one animal more
  // than all of them in half the designs tried. Only their average is used, so the other fits share one design's
  // studies: fit j simulates the j-th of their blocks (Crux's shard), disjoint, at one design's cost in all (checked
  // in Crux: the blocks are disjoint and a lone block's power table is its studies' mean). The average is per fit, not
  // per study: a block one study short when the fits do not divide the studies weighs as much as the others, which
  // only moves its own fit's noise, the same for every fit
  const fitted = lp.calibration?.identification?.solutions ?? [];
  const others = d.fits === false || !fitted.some((s) => s.best) ? [] : fitted.filter((s) => !s.best);
  const powerOnly = (p: Row, shard: number) => {
    const q: Row = {...p, schedules: false, mde: false, ...(others.length > 1 ? {shard, shards: others.length} : {})};
    delete q.n;
    return q;
  };
  const [results, underOthers]: [ArmTables[], ArmTables[][]] = await Promise.all([
    Promise.all(d.arms.map((arm) => designArm(prior, armParams(arm, params), shards, d.power, d.n)
      .catch(armError(arm)))),
    Promise.all(others.map(async (fit, j) => {
      const m = await lineModel(d.line, fit);
      return Promise.all(d.arms.map((arm) =>
        designArm(m.prior, powerOnly(armParams(arm, m.params), j), 1, d.power, null).catch(armError(arm))));
    }))]);

  // every arm with its own n (the treatments differ in what they need); the control as large as the largest arm,
  // so each comparison has at least the control its own design assumed (equal groups): conservative
  const plans = results.map((r) => r.plan);
  const arms: ArmDesign[] = d.arms.map((decl, k) => {
    const p = plans[k];
    const under = underOthers.map((f) => f[k]);
    const curves = under.map((t) => new Map(t.power.map((x) => [x.n as number, finite(x.assurance)])));
    // the power at each n, averaged over the fits as it already is over the atlas's labs: which fit is the true one
    // is one more thing not known (the mean, not the most any needs: that would grow with the fits tried)
    const power = results[k].power.map((x) => {
      const each = [finite(x.assurance), ...curves.map((c) => c.get(x.n) ?? null)];
      return {n: x.n as number, implant: x.implant as number, assurance: each.some((a) => a === null) ? null :
        (each as number[]).reduce((t, a) => t + a, 0) / each.length};
    });
    const fits = [finite(p['n_assurance']), ...under.map((t) => finite(t.plan['n_assurance']))];
    const nAssurance = under.length ? power.find((r) => (r.assurance ?? -1) >= d.power)?.n ?? null : fits[0];
    const n = d.n ?? (nAssurance === null ? null : Math.max(nAssurance, MIN_GROUP));
    const row = n === null ? undefined : power.find((r) => r.n === n);
    return {declaration: decl, archetype: ARCHETYPES[decl.therapy], killRate: p['kill_rate'], effect: p['effect'],
      nAssurance, fits, n, randomise: null,
      // (a given n the averaged table has no power for, as one over 40: the prior's own fit's, as its α)
      assurance: d.n ? (under.length ? row?.assurance : null) ?? finite(p['assurance_at_n']) : row?.assurance ?? null,
      mde: d.n ? finite(p['mde']) : null, power};
  });
  const sizes = arms.map((a) => a.n);
  const nControl = sizes.some((x) => x === null) ? null : Math.max(...(sizes as number[]));
  const largest = nControl === null ? 0 : sizes.indexOf(nControl);
  // two stages of loss: randomised animals that end up not analysable (per group), then implanted animals that never
  // reach randomisation (one pool for every group); each covered with the declared confidence
  const postLoss = Math.max(d.postLoss, lp.take.nonprogressing_map);
  // every group's at once (each a worker round trip)
  const [randomiseControl, ...randomise] = await Promise.all([nControl, ...sizes].map((n) =>
    enough(n, 1 - postLoss, d.confidence, d.lossRule)));
  arms.forEach((a, k) => a.randomise = randomise[k]);
  const sum = (xs: (number | null)[]) =>
    xs.some((x) => x === null) ? null : (xs as number[]).reduce((t, x) => t + x, 0);
  const nTotal = sum([nControl, ...sizes]);
  const randomiseTotal = sum([randomiseControl, ...arms.map((a) => a.randomise)]);
  const implant = await enough(randomiseTotal, 1 - d.preLoss, d.confidence, d.lossRule);
  // each count covers its own loss; what they give together: P(Binom(m, keep) >= n) of the pool and of every group
  const kept = async (m: number | null, n: number | null, keep: number) => m === null || n === null ? null :
    keep >= 1 ? 1 : Number((await crux('stats.dist_eval', {x: Float64Array.from([n - 1])},
      {dist: 'binom', func: 'sf', shape1: m, shape2: keep})).main.col['value'][0]);
  const each = await Promise.all([kept(implant, randomiseTotal, 1 - d.preLoss),
    kept(randomiseControl, nControl, 1 - postLoss), ...arms.map((a) => kept(a.randomise, a.n, 1 - postLoss))]);
  const confidenceStudy = each.some((p) => p === null) ? null : (each as number[]).reduce((t, p) => t * p, 1);
  const first = plans[0];

  // the control from the first call (the same in every call), then each arm's treated rows under its name
  type Curve = StudyDesign['curves'][number];
  const curves: Curve[] = results[0].arms.filter((r) => r.arm === 'control')
    .map((r) => ({...r, arm: 'Control'}) as Curve);
  results.forEach((r, k) => curves.push(...r.arms.filter((x) => x.arm === 'treated')
    .map((x) => ({...x, arm: d.arms[k].name}) as Curve)));

  // by last day: each arm's power at Crux's n for it (the user's n, else its n_assurance, which the minimum group can
  // raise and a line's other fits move: this table and the schedules' are the prior's own fit's) and the controls'
  // state, from the first call (its shares move a little between calls: simulation)
  const byDay = results.map((r) => r.duration);
  // each arm's value in row i of its own table, for the arms designed (an arm with no n has none)
  const perArm = (tables: {[k: string]: any}[][], i: number) =>
    Object.fromEntries(d.arms.flatMap((a, k) => arms[k].n === null ? [] : [[a.name, finite(tables[k][i].assurance)]]));
  const duration = byDay[0].map((row, i) => ({day: row.day, pastEndpoint: row.past_endpoint,
    pastHumane: row.past_humane, assurance: perArm(byDay, i)}));

  // by measurements a week: each arm's power as above, and the study's animals by the same rules as the design's
  const bySchedule = results.map((r) => r.schedules);
  // ponytail: what the other fits move an arm by at the study's own schedule is carried to every schedule; exact
  // would be every fit designed with its schedules (the costliest part, once per fit)
  const moved = arms.map((a) => d.n || a.n === null || a.fits[0] === null ? 0 :
    a.n - Math.max(a.fits[0], MIN_GROUP));
  const schedules = bySchedule[0].map((row, i) => {
    const each = bySchedule.map((t, k) => {
      const needed = finite(t[i].n_assurance);
      return d.n ?? (needed === null ? null : Math.max(Math.max(needed, MIN_GROUP) + moved[k], MIN_GROUP));
    });
    return {perWeek: row.per_week, measurements: row.measurements,
      assurance: perArm(bySchedule, i),
      nTotal: each.some((x) => x === null) ? null : Math.max(...(each as number[])) + sum(each)!};
  });

  const notes = {groups: [] as string[], growth: [] as string[]};
  const end = duration.find((r) => r.day === first['evaluation_day']);
  if (d.lastDay && end && end.pastHumane > 0.1) {
    notes.growth.push(`By day ${end.day}, ${Math.round(end.pastHumane * 100)}% of the controls are over the humane ` +
      `limit (${d.humane} mm³) and leave the study before its end`);
  }
  if (d.humane > GUIDE_LIMIT) {
    notes.growth.push(`The humane limit (${d.humane} mm³, a mean diameter of about ` +
      `${diameterOf(d.humane).toFixed(1)} cm) is above the 1.5 cm that Workman et al. 2010 advise for therapy ` +
      `studies (about ${GUIDE_LIMIT} mm³)`);
  }
  if (d.endpoint >= d.humane) {
    notes.growth.push(`The endpoint (${d.endpoint} mm³) is not below the humane limit (${d.humane} mm³): ` +
      'controls are ' +
      'removed before the median reaches it. Lower the endpoint');
  }
  const pNot = first['p_not_reached'];
  if (pNot > ENDPOINT_RISK) {
    notes.growth.push(`In ${Math.round(pNot * 100)}% of studies like this one the controls level off below ` +
      `${d.endpoint} mm³: a lower endpoint is safer`);
  }
  for (const a of arms) {
    if (!d.n && a.nAssurance !== null && a.nAssurance < MIN_GROUP) {
      notes.groups.push(`${a.declaration.name}: ${a.nAssurance} animals already give ` +
        `${Math.round(d.power * 100)}% power; ` +
        `${MIN_GROUP} are set as the minimum per group, so that losing one animal cannot compromise the comparison`);
    }
    if (a.n === null) {
      notes.groups.push(`${a.declaration.name}: even 40 animals per group do not reach ${Math.round(d.power * 100)}% ` +
        'power; this effect is too small to detect reliably');
    }
    // the line's other fits, said only where they disagree on the group: fits that agree need no words, and one animal
    // apart is agreement (each other fit's own n comes from a block of the studies: one animal is its noise)
    if (!d.n && a.fits.length > 1) {
      const each = a.fits.map((x) => x === null ? Infinity : Math.max(x, MIN_GROUP));
      const [least, most] = [Math.min(...each), Math.max(...each)];
      if (most - least > 1) {
        notes.groups.push(`${a.declaration.name}: the published ${d.line} controls do not pin down how much the ` +
          `animals of one study differ. Fits that describe them equally well need about ${least} to ` +
          `${Number.isFinite(most) ? most : 'more than 40'} animals per group each; the power here is their average`);
      }
    }
  }
  for (const a of arms) {
    if (d.n && a.assurance !== null && a.assurance < d.power) {
      notes.groups.push(`${a.declaration.name}: ${d.n} animals per arm give ${Math.round(a.assurance * 100)}% power, ` +
        `below the ${Math.round(d.power * 100)}% declared`);
    }
  }
  return {
    declaration: d, studyEnd: first['evaluation_day'], endpointDay: first['endpoint_day'], duration, schedules,
    durationQ10: first['duration_q10'],
    durationQ90: finite(first['duration_q90']), pNotReached: pNot, nControl, randomiseControl, nTotal, randomiseTotal,
    implant, postLoss, confidenceStudy,
    alphaArm: plans[largest]['alpha_arm'], calendar: String(first['calendar']).split(',').map(Number), arms, curves,
    measurement: me, nonprogressing: lp.take.nonprogressing_map,
    atlas: {papers: lp.papers, animals: lp.animals, ess: lp.log_alpha.ess}, notes,
  };
}
