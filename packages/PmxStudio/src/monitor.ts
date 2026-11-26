/* Monitor's engine side: a running study's measurements and its plan as Crux's stats.tumor_monitor reads them, its
   results as the tables pmxMonitorStudy returns, and who needs attention; stats.tumor_allocate for the groups. The
   forecasts and the checks are Crux's; here only the translation and the reasons a person reads. */
import * as DG from 'datagrok-api/dg';

import {AnalysisInputs, analysisDeclarationOf, pairGroups, tolerabilityFrames, twinModel} from './analyze';
import {ARCHETYPES} from './archetypes';
import {crux, CruxResult, EmptyState, named, NoControl, rows} from './crux';
import {HUMANE, LAB_WEEKDAYS, TEST, designStudy, finite} from './engine';
import {gpu} from './gpu';
import {SavedPlan} from './plan-file';
import {FINISHED_TAG, oneDecimal} from './study-file';
import {WEEKDAY_NAMES, weekdaysOf} from './frames';
import {LinePrior, bandFactors, linePrior} from './priors';
import {monitorEnd} from './plan-checks';
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
  // (a call that leaves the field out may hand it over as no text at all)
  text = typeof text === 'string' ? text : typeof text === 'number' ? String(text) : '';
  const days = text.split(/[,;\s]+/).filter((x) => x).map(Number);
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

/** Why an animal needs attention today, most urgent first; none: not listed. `removal`: why Welfare reads the animal
 *  as gone (Crux stats.tumor_body_weight: a removal column, else one not measured on its group's last day is removed),
 *  so a missed day says it is counted so there too; none without body weights, where Welfare is empty; PENDING: a
 *  running study's session may still be under way, so the animal is only not measured yet (notYetMeasured). */
export function reasonsOf(a: {p_limit_next?: number | null, last_flag?: number | null}, latest: number, limit: number,
  lastDay: number, lastVolume: number, previous: number | null = null, removal: string | null = null): string[] {
  const out: string[] = [];
  // an animal gone since an earlier session was already seen missing: only the session it first misses lists it
  if (lastDay < latest && lastDay !== previous)
    return out;
  if (lastVolume >= limit)
    out.push(`Over the ${limit} mm³ limit`);
  else if (lastDay < latest && removal === PENDING)
    out.push(`Not measured yet on day ${latest}`);
  else if (lastDay < latest)
    out.push(`Not measured on day ${latest}` + (removal && removal !== 'none' ? ': counted as removed early' : ''));
  else if ((a.p_limit_next ?? 0) >= LIMIT_ALERT)
    out.push(`May pass the ${limit} mm³ limit by the next measurement`);
  // Crux's one decision, the same range the sheet shows: below it (-1), inside (0), above (+1); none on a first value
  if (lastDay === latest && a.last_flag)
    out.push('Last value unexpected: measure again');
  return out;
}

/** reasonsOf's `removal` for an animal a running study has not measured yet on its latest day. */
export const PENDING = 'pending';

/** A running study's animals Crux reads as gone only because today's session has not reached them: Crux's rule (one
 *  not weighed on its group's last day is removed) cannot tell a session under way from a removal. Those it counts
 *  gone for no reason of its own ('other': not over a limit, none given in the table) and last weighed on the session
 *  before (`previous`) go back to the study in Welfare's counts and list (`tolerability`, changed in place); the
 *  animals returned. Crux itself is not changed: a finished study keeps its rule. */
export function notYetMeasured(animals: {[k: string]: any}[], tolerability: DG.DataFrame, previous: number | null,
  given: Set<string>): Set<string> {
  const back = animals.filter((a) => a.removal === 'other' && a.last_day === previous && !given.has(String(a.animal)));
  const out = new Set(back.map((a) => String(a.animal)));
  const [count, list] = [tolerability.col('removed_other'), tolerability.col('removed_animals')];
  for (let i = 0; count && i < tolerability.rowCount; i++) {
    const own = back.filter((a) => a.group === tolerability.get('group', i)).map((a) => String(a.animal));
    if (!own.length)
      continue;
    count.set(i, count.get(i) - own.length);
    // (removed_animals: "M3 (day 7, other reason), M5 (day 4, dead)": split only before an ID, a reason has ", ")
    list?.set(i, String(list.get(i) ?? '').split(/, (?=[^,()]+ \()/)
      .filter((s) => !own.some((a) => s.startsWith(`${a} (`))).join(', '));
  }
  return out;
}

/** Why an animal's body weight needs attention today. `change`: its % change from baseline at its last weighing, as
 *  weighed (the tumor's mass, an estimate, stays in Welfare). `latest`, `previous`: the study's last two weighing days;
 *  an animal not weighed on the latest is listed once, and only while it is still measured; one over the limit at its
 *  last weighing says so on that session too, gone or not (it is why Welfare counts it removed for its weight). */
export function weightReasons(w: {change: number | null, lastDay: number}, latest: number, previous: number | null,
  limit: number, present: boolean): string[] {
  // whole percents, but near the limit one decimal, rounded (19.57% lost reads 19.6%, not "20%, limit 20%" nor a
  // truncated 19.5%); only a loss a hair under the limit, which would round onto it (19.97%), is cut down to 19.9%
  const lost = -(w.change ?? 0);
  const tenth = Math.round(lost * 10) / 10;
  const loss = Math.round(lost) >= limit && lost < limit ?
    (tenth >= limit ? Math.floor(lost * 10) / 10 : tenth).toFixed(1) : Math.round(lost);
  const over = w.change !== null && -w.change >= limit ? [`Over the ${limit}% weight-loss limit (lost ${loss}%)`] : [];
  if (w.lastDay < latest)
    return w.lastDay !== previous ? [] : over.length ? over : present ? [`Not weighed on day ${latest}`] : [];
  if (w.change === null)
    return [];
  if (over.length)
    return over;
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
    throw new EmptyState(d.data.skipped > 0 ? 'No row of the study table has its animal, group, day and volume yet: ' +
      'fill them in, and the forecasts follow.' : 'The study has no measurements yet: paste or type the volumes ' +
      'into the study table, and the forecasts follow.');
  }
  const names = [...new Set(group)];
  if (!names.includes(d.control)) {
    throw new NoControl(d.control, `There is no group "${d.control}" in the study table: pick the control among ` +
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
  const removal = new Map(weighed.map((w) => [String(w.animal), String(w.removal)]));
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
  // a running study (not saved from Analyze): its latest session may be under way, so an animal it has not reached
  // yet is not removed, here nor in Welfare
  const running = x.study?.getTag(FINISHED_TAG) !== 'true';
  if (running) {
    const given = new Set(d.data.removal ? d.data.animal.filter((_, i) => d.data.removal![i].trim()) : []);
    for (const a of notYetMeasured(weighed, tolerability, previousWeighed ?? null, given))
      removal.set(a, PENDING);
  }
  const removalOf = (a: string) => running && (removal.get(a) ?? 'none') === 'none' ? PENDING :
    removal.get(a) ?? null;
  const weightOf = (a: string, lastMeasured: number) => {
    const [w, b] = [lastWeight.get(a), baseline.get(a)];
    if (!w)
      return {reasons: [] as string[], fields: weighDays.size ? {weight_change: null, last_weight: null} : {}};
    const change = b ? 100 * (w.weight / b - 1) : null;
    return {reasons: weightReasons({change, lastDay: w.day}, latestWeighed, previousWeighed ?? null, d.weightLimit,
      lastMeasured >= latestWeighed), fields: {weight_change: change, last_weight: w.weight}};
  };
  const animals = rows(r.main as any).map((a) => {
    const l = last.get(String(a.animal))!;
    const w = weightOf(String(a.animal), l.day);
    return {animal: String(a.animal), group: a.group, last_day: l.day, last_volume: l.volume, next_lower: a.band_lower,
      next_median: a.next_median, next_upper: a.band_upper, check_lower: a.check_lower, check_upper: a.check_upper,
      p_limit: a.p_limit_next, limit_day: a.limit_day, last_flag: a.last_flag,
      ...w.fields,
      // one sentence: each reason after the first goes on in lower case ("…limit; last value unexpected")
      attention: [...reasonsOf(a, latest, d.humane, l.day, l.volume, previous, removalOf(String(a.animal))),
        ...w.reasons].map((x, i) => i ? x[0].toLowerCase() + x.slice(1) : x).join('; ')};
  });
  const study = (r.named as any).study.value;
  const end = monitorEnd({median: study.end_median, q10: study.end_q10, q90: study.end_q90 ?? null},
    await bandFactors(d.line), d.plan?.predicted);
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
      // (the likely range never starts before the day the study is on)
      end_median: study.end_median, end_lower: Number.isFinite(end.q10) ?
        Math.min(Math.max(end.q10, latest), end.q90 ?? Infinity) : end.q10,
      end_upper: end.q90,
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

/** Each animal's latest volume, and its body weight: the latest valid (finite, positive) one on or before that day,
 *  from any row, with or without a volume (a duplicate day does not lose it). */
export function latestReadings(rows: {animal: string, day: number, volume: number, weight: number | null}[]):
  {animal: string[], volume: number[], weight: (number | null)[]} {
  const last = new Map<string, {day: number, v: number}>();
  for (const r of rows) {
    if (Number.isFinite(r.volume) && !(last.get(r.animal)?.day! >= r.day))
      last.set(r.animal, {day: r.day, v: r.volume});
  }
  const w = new Map<string, {day: number, w: number}>();
  for (const r of rows) {
    const at = last.get(r.animal);
    if (at && r.weight !== null && Number.isFinite(r.weight) && r.weight > 0 && r.day <= at.day &&
      !(w.get(r.animal)?.day! >= r.day))
      w.set(r.animal, {day: r.day, w: r.weight});
  }
  const animal = [...last.keys()];
  return {animal, volume: animal.map((a) => last.get(a)!.v), weight: animal.map((a) => w.get(a)?.w ?? null)};
}

/** How many animals with a usable volume (positive, finite) have no usable body weight (finite, positive): the rule
 *  the engine applies before it balances on weight at all (mirrors Crux's `x > 0.0` input filter in fns.rs). */
export function weightMissing(volume: number[], weight: (number | null)[]): number {
  return volume.filter((v, i) => Number.isFinite(v) && v > 0 && !(Number.isFinite(weight[i]) && weight[i]! > 0)).length;
}

/** What the balance area says about body weight; null when the table has no weight column. Complete weights the
 *  engine does not use yet get said too. */
export function balanceNote(hasWeightColumn: boolean, missing: number, used: boolean): string | null {
  if (used)
    return 'Groups balanced on tumor volume first, then on body weight.';
  if (!hasWeightColumn)
    return null;
  return missing > 0 ?
    `Body weight missing for ${missing} animal${missing === 1 ? '' : 's'}: balanced on tumor volume only.` :
    'Balanced on tumor volume only: body weight is not used yet.';
}

/** Balanced groups for these animals: who goes where, and each group's balance. With `weight` (the animals' body
 *  weights, in the animals' order), the engine balances on weight too when every animal has one; `weightUsed` says
 *  whether it did (an engine without it, or a missing weight, gives volume only). */
export async function allocateGroups(animal: string[], volume: number[], groups: {name: string, n: number}[],
  seed: number, weight?: (number | null)[]):
  Promise<{allocation: DG.DataFrame, balance: DG.DataFrame, weightUsed: boolean}> {
  const params = {groups: JSON.stringify(groups), seed};
  const vol = Float64Array.from(volume);
  // weights the engine would ignore are not sent
  const send = weight && weightMissing(volume, weight) === 0;
  const r = await crux('stats.tumor_allocate',
    send ? {animal, volume: vol, weight: Float64Array.from(weight!.map((w) => w ?? NaN))} : {animal, volume: vol},
    params);
  const balance = oneDecimal(DG.DataFrame.fromObjects(named(r, 'balance'))!);
  const mean = balance.col('weight_mean');
  const weightUsed = !!mean && Array.from({length: balance.rowCount}, (_, i) => !mean.isNone(i)).some((x) => x);
  // empty columns would show: only a balance that used weight keeps them
  if (!weightUsed) {
    for (const c of ['weight_mean', 'weight_sd']) {
      if (balance.col(c))
        balance.columns.remove(c);
    }
  }
  return {allocation: oneDecimal(DG.DataFrame.fromObjects(rows(r.main as any))!), balance, weightUsed};
}

/** Animal IDs in natural order (T3-2 before T3-10), a run of digits with a leading zero read as text: T03-01 and
 *  T3-01 are two series, not one interleaved (numeric collation alone reads 03 as 3). */
export function naturalCompare(a: string, b: string): number {
  const [x, y] = [a.match(/\d+|\D+/g) ?? [], b.match(/\d+|\D+/g) ?? []];
  for (let i = 0; i < Math.min(x.length, y.length); i++) {
    const [p, q] = [x[i], y[i]];
    const numeric = /^[1-9]\d*$/.test(p) && /^[1-9]\d*$/.test(q);
    const c = numeric ? Number(p) - Number(q) : p.localeCompare(q);
    if (c)
      return c;
  }
  return x.length - y.length;
}
