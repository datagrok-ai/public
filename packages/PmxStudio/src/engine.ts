/*
 * The Plan engine: Crux's stats.tumor_study_design (crate tgi, [decision:tumor-study-design]), called once per
 * treated arm. This file only translates and combines: the cell line's prior (files/engine_priors_v2.json: the labs'
 * joint distribution of a new study's log alpha and log beta, and the heavy-tailed observation error; an older file's
 * MAP prior of log alpha, Crux stats.map_prior) and the arm's therapy (archetypes.ts) become the
 * function's table and parameters; every arm is designed with the same seed, so the studies drawn, the study end
 * and the controls are the same in every call, and the study's animals per arm are the largest arm's.
 */
import {crux, CruxResult, named, rows} from './crux';
import {MeasurementError, enginePriors, linePrior, mapPrior} from './priors';
import {ARCHETYPES, Archetype, Therapy} from './archetypes';

export type EffectMetric = 'TGI' | 'TC';
export type Correction = 'dunnett' | 'bonferroni' | 'none';

/** A form's correction choice, Dunnett when it is none of them. */
export const correctionOf = (x: unknown): Correction =>
  x === 'bonferroni' || x === 'none' ? x : 'dunnett';

/** Whether x is a finite number; `finite` gives it or null (Crux's nulls, a form's empty optionals). */
export const isFiniteNumber = (x: unknown): x is number => typeof x === 'number' && Number.isFinite(x);
export const finite = (x: unknown): number | null => isFiniteNumber(x) ? x : null;

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
  /** probability that the animals planned leave every group its analysable n */
  confidence: number;
  /** 'confident': enough with that probability (binomial); 'expected': the average loss only, n / (1 - loss) */
  lossRule: 'confident' | 'expected';
  /** whether to scan the measurement schedules (the costliest alternative; off: none, the rest the same) */
  schedules: boolean;
}

/** One arm's design. */
export interface ArmDesign {
  declaration: ArmDeclaration;
  archetype: Archetype;
  killRate: number;
  /** the effect the kill rate gives (percent, as declared) */
  effect: number;
  /** analysable animals this arm needs for the target assurance (null: not up to 40) */
  nAssurance: number | null;
  nConditional: number | null;
  /** the arm's analysable animals (the user's n when given, else nAssurance) and the animals to randomise for them */
  n: number | null;
  randomise: number | null;
  /** at the arm's n */
  assurance: number | null;
  conditional: number | null;
  /** minimum detectable effect at the user's n (percent; null without one, or out of reach) */
  mde: number | null;
  /** null where Crux gives no power */
  power: {n: number, implant: number, assurance: number | null, conditional: number | null}[];
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

/** mm3: the humane limit when none is declared (Crux's default) */
export const HUMANE = 2000;

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
 *  log alpha (one component with v2's labs' bivariate normal; v0 / v1: the MAP prior) and the model parameters. */
export async function lineModel(line: string) {
  const [lp, atlasMe] = [await linePrior(line), (await enginePriors()).measurement_error];
  // v2: a new study from the labs' bivariate normal and the heavy-tailed error; v0 / v1: the MAP prior of log alpha
  const sp = lp.study_prior;
  const nz = lp.noise;
  let prior: {w: ArrayLike<number>, a: ArrayLike<number>, b: ArrayLike<number>};
  if (sp)
    prior = {w: Float64Array.of(1), a: Float64Array.of(sp.mean_log_alpha), b: Float64Array.of(sp.sd_log_alpha)};
  else {
    const map = await mapPrior(line);
    prior = {w: map.main.col['w'], a: map.main.col['a'], b: map.main.col['b']};
  }
  const me: MeasurementError = nz ? {...atlasMe, model: `${atlasMe.model}, Student t (df ${nz.df.toFixed(1)})`,
    c0: nz.c0, c1: nz.c1} : atlasMe;
  const params = {
    log_beta: sp?.mean_log_beta ?? lp.log_beta, log_v_ref: lp.log_v_ref, sd_log_alpha: lp.sd_within.log_alpha,
    sd_log_beta: lp.sd_within.log_beta, sd_study_log_beta: sp?.sd_log_beta ?? lp.sd_study_log_beta,
    rho_study: sp?.rho ?? 0, start_sd: lp.start_sd, nonprogressing: lp.take.nonprogressing_map,
    me_c0: me.c0, me_c1: me.c1, noise_df: nz?.df ?? 0, sd_process: nz ? 0 : lp.sd_process ?? 0,
  };
  return {lp, me, prior, params};
}

export async function designStudy(d: StudyDeclaration): Promise<StudyDesign> {
  if (d.arms.length === 0)
    throw new Error('Add at least one treatment arm');
  const {lp, me, prior, params} = await lineModel(d.line);
  const shared = {
    ...params,
    v_rand: d.vRand,
    endpoint: d.endpoint, power: d.power, alpha: d.alpha, treated_arms: d.treatedArms ?? d.arms.length,
    correction: d.correction,
    weekdays: weekdaysParam(d.weekdays),
    start_weekday: d.startWeekday, humane: d.humane, schedules: d.schedules, seed: SEED,
    ...(d.n ? {n: d.n} : {}), ...(d.lastDay ? {last_day: d.lastDay} : {}),
  };
  const results: CruxResult[] = [];
  for (const arm of d.arms) {
    const a = ARCHETYPES[arm.therapy];
    try {
      results.push(await crux('stats.tumor_study_design', prior, {...shared, metric: arm.metric.toLowerCase(),
        effect: arm.effect, onset: a.onsetDays, kappa_sd: a.kappaSd, responders: a.responders}));
    } catch (e: any) {
      throw new Error(`${arm.name}: ${e?.message ?? e}`);
    }
  }

  // every arm with its own n (the treatments differ in what they need); the control as large as the largest arm,
  // so each comparison has at least the control its own design assumed (equal groups): conservative
  const plans = results.map((r) => rows(r.main as any)[0]);
  const powers = results.map((r) => named(r, 'power').map((x) => ({n: x.n, implant: x.implant,
    assurance: finite(x.assurance), conditional: finite(x.conditional)})));
  const arms: ArmDesign[] = d.arms.map((decl, k) => {
    const p = plans[k];
    const nAssurance = finite(p['n_assurance']);
    const n = d.n ?? (nAssurance === null ? null : Math.max(nAssurance, MIN_GROUP));
    const row = n === null ? undefined : powers[k].find((r) => r.n === n);
    return {declaration: decl, archetype: ARCHETYPES[decl.therapy], killRate: p['kill_rate'], effect: p['effect'],
      nAssurance, nConditional: finite(p['n_conditional']), n, randomise: null,
      assurance: d.n ? finite(p['assurance_at_n']) : row?.assurance ?? null,
      conditional: d.n ? finite(p['conditional_at_n']) : row?.conditional ?? null,
      mde: d.n ? finite(p['mde']) : null, power: powers[k]};
  });
  const sizes = arms.map((a) => a.n);
  const nControl = sizes.some((x) => x === null) ? null : Math.max(...(sizes as number[]));
  const largest = nControl === null ? 0 : sizes.indexOf(nControl);
  // two stages of loss: randomised animals that end up not analysable (per group), then implanted animals that never
  // reach randomisation (one pool for every group); each covered with the declared confidence
  const postLoss = Math.max(d.postLoss, lp.take.nonprogressing_map);
  for (const a of arms)
    a.randomise = await enough(a.n, 1 - postLoss, d.confidence, d.lossRule);
  const randomiseControl = await enough(nControl, 1 - postLoss, d.confidence, d.lossRule);
  const sum = (xs: (number | null)[]) =>
    xs.some((x) => x === null) ? null : (xs as number[]).reduce((t, x) => t + x, 0);
  const nTotal = sum([nControl, ...sizes]);
  const randomiseTotal = sum([randomiseControl, ...arms.map((a) => a.randomise)]);
  const implant = await enough(randomiseTotal, 1 - d.preLoss, d.confidence, d.lossRule);
  const first = plans[0];

  // the control from the first call (the same in every call), then each arm's treated rows under its name
  type Curve = StudyDesign['curves'][number];
  const curves: Curve[] = named(results[0], 'arms').filter((r) => r.arm === 'control')
    .map((r) => ({...r, arm: 'Control'}) as Curve);
  results.forEach((r, k) => curves.push(...named(r, 'arms').filter((x) => x.arm === 'treated')
    .map((x) => ({...x, arm: d.arms[k].name}) as Curve)));

  // by last day: each arm's power at Crux's n for it (the user's n, else its n_assurance, which the minimum group can
  // raise) and the controls' state (the same in every call)
  const byDay = results.map((r) => named(r, 'duration'));
  // each arm's value in row i of its own table, for the arms designed (an arm with no n has none)
  const perArm = (tables: {[k: string]: any}[][], i: number) =>
    Object.fromEntries(d.arms.flatMap((a, k) => arms[k].n === null ? [] : [[a.name, finite(tables[k][i].assurance)]]));
  const duration = byDay[0].map((row, i) => ({day: row.day, pastEndpoint: row.past_endpoint,
    pastHumane: row.past_humane, assurance: perArm(byDay, i)}));

  // by measurements a week: each arm's power as above, and the study's animals by the same rules as the design's
  const bySchedule = results.map((r) => named(r, 'schedules'));
  const schedules = bySchedule[0].map((row, i) => {
    const each = bySchedule.map((t) => {
      const needed = finite(t[i].n_assurance);
      return d.n ?? (needed === null ? null : Math.max(needed, MIN_GROUP));
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
  if (d.endpoint >= d.humane) {
    notes.growth.push(`The endpoint (${d.endpoint} mm³) is not below the humane limit (${d.humane} mm³): ` +
      'controls are ' +
      'removed before the median reaches it. Lower the endpoint');
  }
  const pNot = first['p_not_reached'];
  if (pNot > 0.02) {
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
    implant, postLoss,
    alphaArm: plans[largest]['alpha_arm'], calendar: String(first['calendar']).split(',').map(Number), arms, curves,
    measurement: me, nonprogressing: lp.take.nonprogressing_map,
    atlas: {papers: lp.papers, animals: lp.animals, ess: lp.log_alpha.ess}, notes,
  };
}
