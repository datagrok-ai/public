/* Monitor's engine side: a running study's measurements and its plan as Crux's stats.tumor_monitor reads them, its
   results as the tables pmxMonitorStudy returns, and who needs attention; stats.tumor_allocate for the groups. The
   forecasts and the checks are Crux's; here only the translation and the reasons a person reads. */
import * as DG from 'datagrok-api/dg';

import {AnalysisInputs, analysisDeclarationOf, pairGroups, tolerabilityFrames, twinModel} from './analyze';
import {ARCHETYPES} from './archetypes';
import {crux, CruxResult, named, rows} from './crux';
import {HUMANE, LAB_WEEKDAYS, TEST, designStudy, finite} from './engine';
import {gpu} from './gpu';
import {SavedPlan} from './plan-file';
import {oneDecimal} from './study-file';
import {WEEKDAY_NAMES, weekdaysOf} from './frames';
import {LinePrior, linePrior} from './priors';
import {metaOf, weekdayOf} from './study-meta';

type Study = {animal: string[], group: string[], day: Float64Array, volume: Float64Array};

/** stats.tumor_monitor on the twin's WebGPU device where there is one (its study posterior is the twin's grid work,
 *  most of the call), else on the CPU workers. */
async function fastMonitor(study: Study, params: {[k: string]: any}): Promise<CruxResult> {
  if (await gpu()) {
    try {
      return await crux('stats.tumor_monitor', study, {...params, device: 'webgpu'});
    } catch (e) {
      console.warn('PMX Studio: the Monitor on WebGPU failed, on the CPU workers instead', e);
    }
  }
  return crux('stats.tumor_monitor', study, params);
}

/** The last forecasts by their inputs, as twinOf keeps the twins: a value set back (a limit tried and undone) costs
 *  nothing. */
const monitorCache = new Map<string, Promise<CruxResult>>();
function monitorOf(study: Study, params: {[k: string]: any}): Promise<CruxResult> {
  const key = JSON.stringify([study.animal, study.group, Array.from(study.day), Array.from(study.volume), params]);
  let result = monitorCache.get(key);
  if (result)
    monitorCache.delete(key); // re-inserted below: most recently used last
  else {
    const mine = result = fastMonitor(study, params);
    result.catch(() => monitorCache.get(key) === mine && monitorCache.delete(key));
    if (monitorCache.size >= 8) // ponytail: 8 forecasts, the least recently used out
      monitorCache.delete(monitorCache.keys().next().value!);
  }
  monitorCache.set(key, result);
  return result;
}

/** Chance of passing the limit by the next measurement from which an animal is listed. */
export const LIMIT_ALERT = 0.2;
/** How often a session of `n` animals asks to measure a correct value again, in words. `rate`: the share of correct
 *  values outside the "measure again" range (checkChance). */
export function promptsByChance(n: number, rate: number): string {
  const each = Math.max(n, 1) * rate;
  // (from two thirds of a prompt a session, "one every 1 sessions" is one a session)
  return Math.round(1 / each) < 2 ? `about ${Math.max(1, Math.round(each))} a session` :
    `about one every ${Math.round(1 / each)} sessions`;
}

/** What "measure again" costs by chance: how often a correct value fell outside the range shown before it was
 *  measured, on the atlas's published controls at Crux's own levels, and the share of those above the range (the
 *  engine priors' twin.next_ranges, an average per study; null: not measured).
 *  The line's fitted levels (twin.conformal) are not passed to Crux: they were fitted on the predictive of an animal's
 *  last value, built after its session, and Crux takes one pair of levels for that flag and for the ranges shown
 *  before measuring, which they would leave too narrow (a correct value outside 2.1% of the time for the 1% asked,
 *  on the same controls). */
export function checkChance(ranges: NonNullable<LinePrior['twin']>['next_ranges']):
  {chance: number, above: number} | null {
  const v = ranges?.observed['99'];
  return v && v.miss_lo + v.miss_hi > 0 ? {chance: v.miss_lo + v.miss_hi, above: v.miss_hi / (v.miss_lo + v.miss_hi)} :
    null;
}
/** The share of the weight-loss limit from which an animal is listed (15% for the usual 20%). Our criterion, not a
 *  published one: a margin to act before the limit. */
export const WEIGHT_WATCH = 0.75;

export interface MonitorInputs extends AnalysisInputs {
  /** the planned reviews: study days, comma-separated */
  looks?: string | null;
}

export interface MonitorFrames {
  /** one row: the study today */
  status: DG.DataFrame;
  /** per animal: its next measurement's range and why it needs attention */
  animals: DG.DataFrame;
  /** per group and day, measured and to come */
  days: DG.DataFrame;
  /** per treated arm: the chance of success, on a planned review */
  decisions: DG.DataFrame;
  /** the body weight by group, as Tolerability reads it; empty without a body weight column */
  tolerability: DG.DataFrame;
  weights: DG.DataFrame;
}

/** The planned reviews from their text: the study days, ascending; anything else is an error. */
export function parseLooks(text: string | null | undefined): number[] {
  const days = (text ?? '').split(/[,;\s]+/).filter((x) => x).map(Number);
  if (days.some((d) => !Number.isFinite(d))) {
    throw new Error(`Planned reviews: "${text}" is not a list of study days. Write the days separated by commas, ` +
      'such as 14, 21.');
  }
  return [...new Set(days)].sort((a, b) => a - b);
}

/** The days to come: the plan's measurement weekdays; without a plan, the study's own usual gap (ponytail: a
 *  2-3-2 week reads as every 2 days). */
export function futureDays(days: number[], plan?: SavedPlan | null, start?: string | null): number[] {
  // with a plan, its measurement weekdays from day 0's (the study's date, else the plan's start weekday)
  if (plan && days.length) {
    const on = new Set(weekdaysOf(plan.inputs.measurementDays as string | undefined));
    const day0 = start ? weekdayOf(start) :
      Math.max(0, WEEKDAY_NAMES.indexOf(String(plan.inputs.startWeekday ?? 'Monday')));
    const from = Math.max(...days);
    return Array.from({length: 60}, (_, i) => from + i + 1).filter((d) => on.has((day0 + d) % 7));
  }
  const seen = [...new Set(days)].sort((a, b) => a - b);
  const gaps = seen.slice(1).map((d, i) => d - seen[i]).sort((a, b) => a - b);
  const gap = gaps.length ? Math.max(gaps[Math.floor(gaps.length / 2)], 1) : 3;
  const last = seen[seen.length - 1];
  return Array.from({length: Math.ceil(60 / gap)}, (_, i) => last + gap * (i + 1));
}

/** Each treated group as an arm of the plan: its kill rate (the design's for the plan's effect) and its therapy's
 *  archetype. A group not in the plan, or a study with no plan: no assumed effect (a kill rate of 0), and no chance of
 *  success is given for it. */
async function armsOf(plan: SavedPlan | null, treated: string[], notes: string[]) {
  const pairs = plan ? pairGroups(treated, plan.predicted.arms, notes) :
    new Map<string, SavedPlan['predicted']['arms'][0]>();
  const kill = new Map<string, number>();
  if (plan && plan.predicted.arms.every((a) => typeof a.killRate === 'number')) {
    // as planned: the rate the plan gave each arm, kept in its file
    plan.predicted.arms.forEach((a) => kill.set(a.name, a.killRate!));
  } else if (plan) {
    // a plan saved before its rates were kept: designed again as the plan declared it (its endpoint and limit, not
    // the ones in Monitor's form: the declared effect is at the plan's endpoint), which the engine remembers
    const arms = plan.predicted.arms.map((a) => ({name: a.name, therapy: a.therapy, metric: a.metric,
      effect: a.effect}));
    const d = await designStudy({line: plan.predicted.line, vRand: plan.predicted.vRand,
      endpoint: plan.predicted.endpoint, power: 0.8, alpha: 0.05, correction: 'holm', n: null, arms, preLoss: 0,
      postLoss: 0, lossRule: 'expected', confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0,
      humane: finite(plan.inputs.humane) ?? HUMANE, schedules: false, lastDay: null, fits: false});
    d.arms.forEach((a) => kill.set(a.declaration.name, a.killRate));
  }
  return treated.map((g) => {
    const p = pairs.get(g);
    const arch = p ? ARCHETYPES[p.therapy] : ARCHETYPES.small_molecule;
    return {name: g, planned: !!p, kill_rate: p ? kill.get(p.name) ?? 0 : 0, onset: arch.onsetDays,
      kappa_sd: arch.kappaSd, responders: arch.responders};
  });
}

/** Why an animal needs attention today, most urgent first; none: not listed. */
export function reasonsOf(a: {p_limit_next?: number | null, last_flag?: number | null}, latest: number, limit: number,
  lastDay: number, lastVolume: number, previous: number | null = null): string[] {
  const out: string[] = [];
  // an animal gone since an earlier session was already seen missing: only the session it first misses lists it
  if (lastDay < latest && lastDay !== previous)
    return out;
  if (lastVolume >= limit)
    out.push(`Over the ${limit} mm³ limit`);
  else if (lastDay < latest)
    out.push(`Not measured on day ${latest}`);
  else if ((a.p_limit_next ?? 0) >= LIMIT_ALERT)
    out.push(`May pass the ${limit} mm³ limit by the next measurement`);
  // Crux's one decision, the same range the sheet shows: below it (-1), inside (0), above (+1); none on a first value
  if (lastDay === latest && a.last_flag)
    out.push('Last value unexpected: measure again');
  return out;
}

/** Why an animal's body weight needs attention today. `change`: its % change from baseline at its last weighing, as
 *  weighed (the tumor's mass, an estimate, stays in Welfare). `latest`, `previous`: the study's last two weighing days;
 *  an animal not weighed on the latest is listed once, and only while it is still measured. */
export function weightReasons(w: {change: number | null, lastDay: number}, latest: number, previous: number | null,
  limit: number, present: boolean): string[] {
  if (w.lastDay < latest)
    return present && w.lastDay === previous ? [`Not weighed on day ${latest}`] : [];
  if (w.change === null)
    return [];
  const loss = Math.round(-w.change);
  if (-w.change >= limit)
    return [`Over the ${limit}% weight-loss limit (lost ${loss}%)`];
  return -w.change >= limit * WEIGHT_WATCH ? [`Losing weight: ${loss}%, limit ${limit}%`] : [];
}

const emptyDecisions = () => DG.DataFrame.fromColumns([DG.Column.string('arm', 0), DG.Column.bool('planned', 0),
  DG.Column.float('success', 0), DG.Column.float('tgi_lower', 0), DG.Column.float('tgi_median', 0),
  DG.Column.float('tgi_upper', 0), DG.Column.float('n_end', 0)]);

export async function monitorFrames(x: MonitorInputs): Promise<MonitorFrames> {
  const d = analysisDeclarationOf(x);
  const looks = parseLooks(x.looks);
  const keep = d.data.volume.map((v) => Number.isFinite(v));
  const pick = <T>(xs: T[]) => xs.filter((_, i) => keep[i]);
  const [animal, group, day, volume] = [pick(d.data.animal), pick(d.data.group), pick(d.data.day), pick(d.data.volume)];
  if (day.length === 0) {
    // rows held back for a missing animal, group or day are not "no measurements": say what they lack
    throw new Error(d.data.skipped > 0 ? 'No row of the study table has its animal, group, day and volume yet: ' +
      'fill them in, and the forecasts follow.' : 'The study has no measurements yet: paste or type the volumes ' +
      'into the study table, and the forecasts follow.');
  }
  const names = [...new Set(group)];
  if (!names.includes(d.control)) {
    throw new Error(`There is no group "${d.control}" in the study table: pick the control among ` +
      `${names.join(', ')}.`);
  }
  const notes: string[] = [];
  // the plan's arms, the line's twin and the body weights do not wait for one another
  const [arms, model, {tolerability, weights, animals: weighed}, byChance] = await Promise.all([
    armsOf(d.plan, names.filter((g) => g !== d.control), notes), twinModel(d.line),
    tolerabilityFrames(d), linePrior(d.line).then((lp) => checkChance(lp.twin?.next_ranges))]);
  const latest = Math.max(...day);
  const future = futureDays(day, d.plan, metaOf(x.study)?.start);
  const r = await monitorOf({animal, group, day: Float64Array.from(day),
    volume: Float64Array.from(volume)}, {...model, control: d.control,
    ...arms.length ? {arms: JSON.stringify(arms)} : {}, future: future.join(','), endpoint: d.endpoint,
    humane: d.humane, alpha: d.alpha, correction: d.correction, test: TEST,
    ...d.day === null ? {} : {last_day: d.day}});

  // each animal's last measurement, which the forecast does not repeat
  const last = new Map<string, {day: number, volume: number}>();
  animal.forEach((a, i) => {
    const l = last.get(a);
    if (!l || day[i] > l.day)
      last.set(a, {day: day[i], volume: volume[i]});
  });
  const previous = [...new Set(day)].sort((p, q) => q - p)[1] ?? null;
  // an id with its group in it shows the group on every list, hidden or not
  if ([...new Set(animal.map((a, i) => a.toLowerCase().includes(group[i].toLowerCase())))].includes(true))
    notes.push('Some animal ids contain their group\'s name, so the lists below show the group');
  // the body weight is optional: without it, none of this says anything
  const baseline = new Map(weighed.map((w) => [String(w.animal), w.baseline as number | null]));
  const lastWeight = new Map<string, {day: number, weight: number}>();
  const weighDays = new Set<number>();
  d.data.weight?.forEach((w, i) => {
    const [a, wd] = [d.data.animal[i], d.data.day[i]];
    if (typeof w === 'number' && w > 0) {
      weighDays.add(wd);
      if (!(lastWeight.get(a)?.day! >= wd))
        lastWeight.set(a, {day: wd, weight: w});
    }
  });
  const [latestWeighed, previousWeighed] = [...weighDays].sort((p, q) => q - p);
  const weightOf = (a: string, lastMeasured: number) => {
    const [w, b] = [lastWeight.get(a), baseline.get(a)];
    if (!w)
      return {reasons: [] as string[], fields: weighDays.size ? {weight_change: null} : {}};
    const change = b ? 100 * (w.weight / b - 1) : null;
    return {reasons: weightReasons({change, lastDay: w.day}, latestWeighed, previousWeighed ?? null, d.weightLimit,
      lastMeasured >= latestWeighed), fields: {weight_change: change}};
  };
  const animals = rows(r.main as any).map((a) => {
    const l = last.get(String(a.animal))!;
    const w = weightOf(String(a.animal), l.day);
    return {animal: String(a.animal), group: a.group, last_day: l.day, last_volume: l.volume, next_lower: a.band_lower,
      next_median: a.next_median, next_upper: a.band_upper, check_lower: a.check_lower, check_upper: a.check_upper,
      p_limit: a.p_limit_next, limit_day: a.limit_day, last_flag: a.last_flag,
      ...w.fields,
      attention: [...reasonsOf(a, latest, d.humane, l.day, l.volume, previous), ...w.reasons].join('; ')};
  });
  const study = (r.named as any).study.value;
  // the chance of success only on a planned review: between them, an empty table (MONITOR_DESIGN.md §3)
  // a review planned for a day that is not measured (day 10, a Thursday, of a Monday-Wednesday-Friday study) is
  // held on the first measurement from that day on
  const atLook = looks.some((l) => l <= latest && (previous === null || l > previous));
  const decisions = named(r, 'decisions').map((a, i) => ({arm: arms[i].name, planned: arms[i].planned,
    success: a.success, tgi_lower: a.tgi_lo, tgi_median: a.tgi_median, tgi_upper: a.tgi_hi, n_end: a.n_end}));
  const days = DG.DataFrame.fromObjects(named(r, 'days').map((g) => ({group: g.group, day: g.day, lower: g.q05,
    median: g.median, upper: g.q95, animals: g.animals})))!;
  days.col('day')!.meta.units = 'd';
  // volumes to one decimal in every viewer's tooltip and cell (Crux's floats carry a dozen)
  const animalsDf = DG.DataFrame.fromObjects(animals)!;
  for (const [df, cols] of [[days, ['lower', 'median', 'upper']], [animalsDf, ['last_volume', 'next_lower',
    'next_median', 'next_upper', 'check_lower', 'check_upper']]] as const) {
    for (const q of cols) {
      const c = df.col(q);
      if (c) {
        c.meta.units = 'mm³';
        c.meta.format = '0.0';
      }
    }
  }
  return {
    status: DG.DataFrame.fromObjects([{line: d.line, control: d.control, day: latest, next_day: future[0],
      end_median: study.end_median, end_lower: study.end_q10, end_upper: study.end_q90,
      planned_end: d.plan?.predicted.studyEnd ?? null, plan: d.plan?.name ?? null,
      attention: animals.filter((a) => a.attention).length, animals: animals.length,
      next_look: looks.find((l) => l > latest) ?? null, at_look: atLook, looks: looks.join(', '), humane: d.humane,
      check_chance: byChance?.chance ?? null, check_above: byChance?.above ?? null, notes: notes.join(' | ')}])!,
    animals: animalsDf,
    days,
    decisions: atLook && decisions.length ? DG.DataFrame.fromObjects(decisions)! : emptyDecisions(),
    tolerability,
    weights,
  };
}

/** The values just entered for the next measurement day that fall outside what the last run expected ("measure
 *  again"). Only that day: the band is for it, and another day's value has no band yet. */
export function outsideBand(rows: {animal: string, day: number, volume: number}[], nextDay: number,
  bands: Map<string, {check_lower: number | null, check_upper: number | null}>):
  {animal: string, volume: number, lower: number, upper: number}[] {
  const out = [];
  for (const r of rows) {
    const b = bands.get(r.animal);
    if (r.day === nextDay && Number.isFinite(r.volume) && b && b.check_lower !== null && b.check_upper !== null &&
      (r.volume < b.check_lower || r.volume > b.check_upper))
      out.push({animal: r.animal, volume: r.volume, lower: b.check_lower, upper: b.check_upper});
  }
  return out;
}

/** Balanced groups for these animals: who goes where, and each group's balance. */
export async function allocateGroups(animal: string[], volume: number[], groups: {name: string, n: number}[],
  seed: number): Promise<{allocation: DG.DataFrame, balance: DG.DataFrame}> {
  const r = await crux('stats.tumor_allocate', {animal, volume: Float64Array.from(volume)},
    {groups: JSON.stringify(groups), seed});
  return {allocation: oneDecimal(DG.DataFrame.fromObjects(rows(r.main as any))!),
    balance: oneDecimal(DG.DataFrame.fromObjects(named(r, 'balance'))!)};
}
