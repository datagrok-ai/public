/* Monitor's engine side: a running study's measurements and its plan as Crux's stats.tumor_monitor reads them, its
   results as the tables pmxMonitorStudy returns, and who needs attention; stats.tumor_allocate for the groups. The
   forecasts and the checks are Crux's; here only the translation and the reasons a person reads. */
import {maxOf} from './format';
import * as DG from 'datagrok-api/dg';

import {AnalysisInputs, analysisDeclarationOf, inMm3, pairGroups, remembered, tolerabilityFrames, twinModel}
  from './analyze';
import {ARCHETYPES} from './archetypes';
import {crux, CruxResult, EmptyState, named, NoControl, rows} from './crux';
import {HUMANE, LAB_WEEKDAYS, TEST, cruxCorrection, designStudy, finite} from './engine';
import {gpu} from './gpu';
import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';
import {FINISHED_TAG, FLAGGED, MARK_REASON, groupsShown, oneDecimal, ticked} from './study-file';
import {withoutExcluded} from './study-table';
import {checkReading, likelyTypo, mm3Text} from './reading-check';
export {TYPO_FOLD, checkReading, likelyTypo, mm3Text} from './reading-check';
export type {ReadingCheck} from './reading-check';
import {WEEKDAY_NAMES, weekdaysOf} from './frames';
import {LinePrior, bandFactors, linePrior} from './priors';
import {monitorEnd} from './plan-checks';
import {lossText, metaOf, weekdayOf} from './study-meta';

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
const monitorOf = remembered(([study, params]: [Study, {[k: string]: any}]) => fastMonitor(study, params), 8);

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

/** Day 0's weekday (Monday 0): the study's date's, else the plan's start weekday (Monday with neither). */
export const day0Weekday = (start: string | null | undefined, plan: SavedPlan | null | undefined): number =>
  start ? weekdayOf(start) : Math.max(0, WEEKDAY_NAMES.indexOf(String(plan?.inputs.startWeekday ?? 'Monday')));

/** The days to come: the plan's measurement weekdays; without a plan, the study's own usual gap (ponytail: a
 *  2-3-2 week reads as every 2 days). */
export function futureDays(days: number[], plan?: SavedPlan | null, start?: string | null): number[] {
  // with a plan, its measurement weekdays from day 0's
  if (plan && days.length) {
    const on = new Set(weekdaysOf(plan.inputs.measurementDays as string | undefined));
    const day0 = day0Weekday(start, plan);
    const from = maxOf(days);
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
      endpoint: plan.predicted.endpoint, power: 0.8, alpha: 0.05, correction: 'bonferroni', n: null, arms, preLoss: 0,
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

/** Welfare's reason for an animal it counts as removed early, as Crux gives it ('' when it does not: none given, or
 *  not measured yet on a day under way). */
export const earlyRemoval = (removal: string | null): string =>
  !removal || removal === 'none' || removal === PENDING ? '' : removal;

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
  const loss = lossText(-(w.change ?? 0), limit);
  const over = w.change !== null && -w.change >= limit ? [`Over the ${limit}% weight-loss limit (lost ${loss}%)`] : [];
  if (w.lastDay < latest)
    return w.lastDay !== previous ? [] : over.length ? over : present ? [`Not weighed on day ${latest}`] : [];
  if (w.change === null)
    return [];
  if (over.length)
    return over;
  return -w.change >= limit * WEIGHT_WATCH ? [`Lost ${loss}% of body weight (limit ${limit}%)`] : [];
}

const emptyDecisions = () => DG.DataFrame.fromColumns([DG.Column.string('arm', 0), DG.Column.bool('planned', 0),
  DG.Column.float('success', 0), DG.Column.float('tgi_lower', 0), DG.Column.float('tgi_median', 0),
  DG.Column.float('tgi_upper', 0), DG.Column.float('n_end', 0)]);

// (the checks of one value: reading-check.ts, kept there out of an import cycle with study-table.ts)

/** The control as the plan has it: the one group of the table that is none of the plan's arms (undefined: no plan, or
 *  not exactly one such group). */
export function planControl(x: AnalysisInputs): string | undefined {
  const text = x.plan?.trim() || x.study.getTag(PLAN_TAG);
  let arms: Set<string>;
  try {
    arms = new Set(text ? parsePlan(text).predicted.arms.map((a) => a.name) : []);
  } catch {
    return undefined;
  }
  if (!arms.size)
    return undefined;
  const others = [...new Set(Array.from({length: x.group.length}, (_, i) => x.group.isNone(i) ? '' :
    String(x.group.get(i)).trim()))].filter((g) => g && !arms.has(g));
  return others.length === 1 ? others[0] : undefined;
}
/** A reading's key among likelyTypos'. */
export const typoKey = (animal: string, day: number) => `${animal}\u0000${day}`;

export type Typo = {animal: string, day: number, volume: number, t: string, trusted: {day: number, volume: number},
  last: boolean};

/** Each animal's likely typos, by typoKey: a volume about tenfold off the animal's last trusted one (the next one is
 *  read against the trusted value, not the typo); `last`: its latest reading. The one rule Monitor's welfare and
 *  Compare's ranges leave them out by. */
export function likelyTypos(rows: {animal: string, day: number, volume: number}[]): Map<string, Typo> {
  const series = new Map<string, {day: number, volume: number}[]>();
  for (const r of rows)
    (series.get(r.animal) ?? series.set(r.animal, []).get(r.animal)!).push({day: r.day, volume: r.volume});
  const out = new Map<string, Typo>();
  series.forEach((s, animal) => {
    s.sort((p, q) => p.day - q.day);
    let trusted: {day: number, volume: number} | null = null;
    for (const m of s) {
      const t = trusted && likelyTypo(trusted.volume, m.volume);
      if (!t)
        trusted = m;
      else
        out.set(typoKey(animal, m.day), {animal, ...m, t, trusted: trusted!, last: m === s[s.length - 1]});
    }
  });
  return out;
}

/** Said until the study has a tumor model: never one assumed, nor the last study's. */
export const LINE_NOT_SET = 'Tumor model: not set. Pick the study\'s model in Tumor model, and the forecasts follow.';

/** The study table's own removal column, by its name (a reason, a fate), never PMX Studio's mark reason. */
export function removalColumn(df: DG.DataFrame): DG.Column | null {
  return df.columns.toList().find((c) => c.type === DG.TYPE.STRING && c.name !== MARK_REASON &&
    /reason|removal|fate|cause/i.test(c.name)) ?? null;
}

/** The animals the table says left the study, each with the first day a removal reason is given on and the reason. */
export function removedAnimals(df: DG.DataFrame, animal: DG.Column, day: DG.Column):
  Map<string, {day: number, why: string}> {
  const out = new Map<string, {day: number, why: string}>();
  const c = removalColumn(df);
  for (let i = 0; c && i < df.rowCount; i++) {
    const why = c.isNone(i) ? '' : String(c.get(i)).trim();
    if (!why || animal.isNone(i) || day.isNone(i))
      continue;
    const a = String(animal.get(i));
    if (!(out.get(a)?.day! <= day.get(i)))
      out.set(a, {day: day.get(i), why});
  }
  return out;
}

/** Two readings of one animal on one day (two volumes, or two body weights), in words; null when there are none. The
 *  engine cannot read them, and only the user can tell which one stands. */
export function twoReadings(rows: {animal: string, day: number, volume: number | null, weight: number | null}[]):
  string | null {
  const seen = new Map<string, {volume: boolean, weight: boolean}>();
  for (const r of rows) {
    const key = `${r.animal}\u0000${r.day}`;
    const [v, w] = [Number.isFinite(r.volume) && r.volume !== null, Number.isFinite(r.weight) && r.weight !== null];
    const was = seen.get(key);
    if (was && ((v && was.volume) || (w && was.weight))) {
      return `${r.animal} has two ${v && was.volume ? 'volumes' : 'body weights'} on day ${r.day}: keep one row per ` +
        'animal and day in the study table (correct or delete the other).';
    }
    seen.set(key, {volume: v || !!was?.volume, weight: w || !!was?.weight});
  }
  return null;
}

export async function monitorFrames(all: MonitorInputs): Promise<MonitorFrames> {
  // the rows a person excluded never reach the engine; an animal excluded whole is still read for its welfare (below):
  // the exclusion is about the analysis, not the animal
  const kept = withoutExcluded(all);
  // nothing measured yet (a new study, one just followed): that, not a control to pick
  if (!Array.from({length: kept.volume.length}, (_, i) => i).some((i) => !kept.volume.isNone(i))) {
    throw new EmptyState('The study has no measurements yet: paste or type the volumes into the study table, ' +
      'and the forecasts follow.');
  }
  // every forecast is the model's: none until the study has one
  if (!kept.line?.trim())
    throw new EmptyState(LINE_NOT_SET);
  const x = kept.control?.trim() ? kept : {...kept, control: planControl(kept) ?? kept.control};
  const d = analysisDeclarationOf(x);
  const twice = twoReadings(d.data.animal.map((a, i) => ({animal: a, day: d.data.day[i], volume: d.data.volume[i],
    weight: d.data.weight?.[i] ?? null})));
  if (twice)
    throw new Error(twice);
  const looks = parseLooks(x.looks);
  // a negative volume is no measurement: said as impossible, kept from the forecasts
  const keep = d.data.volume.map((v) => Number.isFinite(v) && v >= 0);
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
  const key = typoKey;
  const impossible = new Map<string, {day: number, volume: number}>();
  d.data.volume.forEach((v, i) => {
    const a = d.data.animal[i];
    if (v < 0 && !(impossible.get(a)?.day! > d.data.day[i]))
      impossible.set(a, {day: d.data.day[i], volume: v});
  });
  // each animal's volumes by day; a value about tenfold off the last trusted one is a likely typo (likelyTypos)
  const series = new Map<string, {day: number, volume: number}[]>();
  animal.forEach((a, i) => (series.get(a) ?? series.set(a, []).get(a)!).push({day: day[i], volume: volume[i]}));
  series.forEach((s) => s.sort((p, q) => p.day - q.day));
  const found = likelyTypos(animal.map((a, i) => ({animal: a, day: day[i], volume: volume[i]})));
  const typos = new Set(found.keys());
  const typoNow = new Map<string, string>();
  for (const m of found.values()) {
    if (m.last) {
      typoNow.set(m.animal, `Likely a typo: ${mm3Text(m.volume)} mm³ on day ${m.day} is ${m.t} day ` +
        `${m.trusted.day}'s ${mm3Text(m.trusted.volume)} mm³: check it`);
    }
  }
  // the tumor-corrected weight takes the tumor's mass off: not from a volume flagged or a likely typo (a 4,201 mm³
  // slip is 4 g of "tumor" off a 20 g mouse)
  const flagged = x.study.col(FLAGGED);
  for (let i = 0; flagged && i < x.study.rowCount; i++) {
    if (ticked(flagged, i) && !x.animal.isNone(i) && !x.day.isNone(i))
      typos.add(key(String(x.animal.get(i)), x.day.get(i)));
  }
  const forWelfare = {...d, data: {...d.data, volume: d.data.volume.map((v, i) =>
    v < 0 || typos.has(key(d.data.animal[i], d.data.day[i])) ? NaN : v)}};
  const notes: string[] = [];
  // the plan's arms, the line's twin and the body weights do not wait for one another
  const [arms, model, {tolerability, weights, animals: weighed}, byChance] = await Promise.all([
    armsOf(d.plan, names.filter((g) => g !== d.control), notes), twinModel(d.line),
    tolerabilityFrames(forWelfare), linePrior(d.line).then((lp) => checkChance(lp.twin?.next_ranges))]);
  const latest = maxOf(day);
  const future = futureDays(day, d.plan, metaOf(x.study)?.start);
  const sent = {animal, group, day: Float64Array.from(day), volume: Float64Array.from(volume)};
  const params = {...model, control: d.control, future: future.join(','), endpoint: d.endpoint, humane: d.humane,
    alpha: d.alpha, correction: cruxCorrection(d.correction), test: TEST, ...d.day === null ? {} : {last_day: d.day}};
  // while the groups are hidden, each animal's ranges are the study's as one group: a group's own would tell the
  // groups apart on the sheet (a failure leaves no ranges, never a group's)
  const blinded = !groupsShown(all.study);
  const POOLED = 'all';
  const [r, pooled] = await Promise.all([
    monitorOf([sent, {...params, ...arms.length ? {arms: JSON.stringify(arms)} : {}}]),
    blinded ? monitorOf([{...sent, group: animal.map(() => POOLED)}, {...params, control: POOLED}])
      .catch((e) => (console.warn('PMX Studio: the pooled ranges failed', e), null)) : null]);
  const pooledOf = new Map(pooled ? rows(pooled.main as any).map((p) => [String(p.animal), p]) : []);
  const own = (a: {[k: string]: any}): {[k: string]: any} => blinded ? pooledOf.get(String(a.animal)) ?? {} : a;

  // each animal's last measurement, which the forecast does not repeat
  const lastOf = (a: string) => {
    const s = series.get(a)!;
    return s[s.length - 1];
  };
  const [, previous = null, before = null] = [...new Set(day)].sort((p, q) => q - p);
  // an id with its group in it shows the group on every list, hidden or not
  if (animal.some((a, i) => a.toLowerCase().includes(group[i].toLowerCase())))
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
  const [latestWeighed, previousWeighed, beforeWeighed] = [...weighDays].sort((p, q) => q - p);
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
  // a day under way: groups measured on the session before and not one animal of them yet on the latest. Their
  // animals are read as of their own last session (their alerts stay), and counted once as still to measure
  const groupsOn = (dd: number | null) => new Set(group.filter((_, i) => day[i] === dd));
  const measuredNow = groupsOn(latest);
  const waitingGroups = running && previous !== null ? [...groupsOn(previous)].filter((g) => !measuredNow.has(g)) :
    [];
  // an animal gone before the session before, in a group the day under way has not reached: read as of that
  // session, as its group is, so the session it first missed still lists it (welfare alerts and all) and a few
  // values of the new day in other groups do not drop it
  const readAsBefore = (g: string, lastDay: number) => previous !== null && lastDay < previous &&
    waitingGroups.includes(g);
  const weightOf = (a: string, lastMeasured: number, waiting: boolean, asBefore: boolean) => {
    const [w, b] = [lastWeight.get(a), baseline.get(a)];
    if (!w)
      return {reasons: [] as string[], fields: weighDays.size ? {weight_change: null, last_weight: null} : {}};
    const change = b ? 100 * (w.weight / b - 1) : null;
    const at = {change, lastDay: w.day};
    const reasons = asBefore ? previousWeighed === undefined ? [] :
      weightReasons(at, previousWeighed, beforeWeighed ?? null, d.weightLimit, lastMeasured >= previousWeighed) :
      waiting ? weightReasons(at, previousWeighed ?? latestWeighed, null, d.weightLimit, true) :
        weightReasons(at, latestWeighed, previousWeighed ?? null, d.weightLimit, lastMeasured >= latestWeighed);
    return {reasons, fields: {weight_change: change, last_weight: w.weight}};
  };
  // one sentence: each reason after the first goes on in lower case ("…limit; last value unexpected")
  const sentence = (xs: string[]) => xs.map((s, i) => i ? s[0].toLowerCase() + s.slice(1) : s).join('; ');
  // an animal the table says has left the study (a removal reason on its rows) needs no attention any more
  const removedOn = removedAnimals(all.study, all.animal, all.day);
  const off = (a: string, lastDay: number) => (removedOn.get(a)?.day ?? Infinity) <= lastDay;
  const animals: {[k: string]: any}[] = rows(r.main as any).map((a) => {
    const name = String(a.animal);
    const l = lastOf(name);
    if (off(name, l.day)) {
      return {animal: name, group: a.group, last_day: l.day, last_volume: l.volume, next_lower: null,
        next_median: null, next_upper: null, check_lower: null, check_upper: null, p_limit: null, limit_day: null,
        last_flag: null, waiting: false, ...weightOf(name, l.day, false, false).fields, removed_early: '',
        attention: ''};
    }
    const o = own(a);
    const waiting = l.day === previous && waitingGroups.includes(String(a.group));
    const asBefore = readAsBefore(String(a.group), l.day);
    const w = weightOf(name, l.day, waiting, asBefore);
    const bad = impossible.get(name);
    const typo = typoNow.get(name);
    const reasons = waiting ? reasonsOf(o, l.day, d.humane, l.day, l.volume) : asBefore ?
      reasonsOf(o, previous!, d.humane, l.day, l.volume, before, removalOf(name)) :
      reasonsOf(o, latest, d.humane, l.day, l.volume, previous, removalOf(name));
    return {animal: name, group: a.group, last_day: l.day, last_volume: l.volume, next_lower: o.band_lower ?? null,
      next_median: o.next_median ?? null, next_upper: o.band_upper ?? null, check_lower: o.check_lower ?? null,
      check_upper: o.check_upper ?? null, p_limit: o.p_limit_next ?? null, limit_day: o.limit_day ?? null,
      last_flag: o.last_flag ?? null, waiting, ...w.fields,
      // Welfare counts it removed (Crux: not measured on its group's last day): the day's entry lists it greyed,
      // not counted (QA 0.1.3 round 3: technician F22)
      removed_early: earlyRemoval(removalOf(name)),
      attention: sentence([...bad && bad.day >= l.day ? [`Impossible volume on day ${bad.day} ` +
        `(${mm3Text(bad.volume)} mm³): correct it`] : [], ...typo ? [typo] : [],
      // (a likely typo already says to check the value; the limits wait until it is confirmed or corrected)
      ...typo ? reasons.filter((s) => !/^(Last value unexpected|Over the|May pass)/.test(s)) : reasons, ...w.reasons])};
  });
  // an animal excluded whole from the analysis: its humane limits are still read, from its own rows
  if (kept !== all) {
    const full = analysisDeclarationOf(all).data;
    const rowsOfAnimal = new Map<string, number[]>();
    full.animal.forEach((a, i) => (rowsOfAnimal.get(a) ?? rowsOfAnimal.set(a, []).get(a)!).push(i));
    for (const [a, is] of rowsOfAnimal) {
      if (series.has(a))
        continue;
      is.sort((p, q) => full.day[p] - full.day[q]);
      const vols = is.filter((i) => Number.isFinite(full.volume[i]) && full.volume[i] >= 0);
      const lv = vols[vols.length - 1];
      // gone before the session before: as reasonsOf, not listed (the one before that while its group waits)
      const asBefore = lv !== undefined && readAsBefore(full.group[lv], full.day[lv]);
      if (lv === undefined || (previous !== null && full.day[lv] < (asBefore ? before ?? -Infinity : previous)) ||
        off(a, full.day[lv]))
        continue;
      const ws = full.weight ? is.filter((i) => typeof full.weight![i] === 'number' && full.weight![i]! > 0) : [];
      // ponytail: its first weighing as its baseline (Crux's own baseline is read only for the animals it is given)
      const [b, lw] = [ws.length ? full.weight![ws[0]]! : null, ws[ws.length - 1]];
      const change = b !== null && lw !== undefined ? 100 * (full.weight![lw]! / b - 1) : null;
      const [wl, wp] = asBefore ? [previousWeighed, beforeWeighed] : [latestWeighed, previousWeighed];
      const wr = lw === undefined || wl === undefined ? [] : weightReasons({change,
        lastDay: full.day[lw]}, wl, wp ?? null, d.weightLimit, full.day[lv] >= wl);
      animals.push({animal: a, group: full.group[lv], last_day: full.day[lv], last_volume: full.volume[lv],
        next_lower: null, next_median: null, next_upper: null, check_lower: null, check_upper: null, p_limit: null,
        limit_day: null, last_flag: null, waiting: false, removed_early: '',
        ...weighDays.size ? {weight_change: change, last_weight: lw === undefined ? null : full.weight![lw]} : {},
        attention: sentence([...full.volume[lv] >= d.humane ? [`Over the ${d.humane} mm³ limit`] : [], ...wr])});
    }
  }
  const study = (r.named as any).study.value;
  const end = monitorEnd({median: study.end_median, q10: study.end_q10, q90: study.end_q90 ?? null},
    await bandFactors(d.line), d.plan?.predicted);
  // the chance of success only on a planned review: between them, an empty table (MONITOR_DESIGN.md §3)
  // a review planned for a day that is not measured (day 10, a Thursday, of a Monday-Wednesday-Friday study) is
  // held on the first measurement from that day on; a review on a day under way waits for its last group
  const lookDue = looks.some((l) => l <= latest && (previous === null || l > previous));
  const atLook = lookDue && waitingGroups.length === 0;
  const decisions = named(r, 'decisions').map((a, i) => ({arm: arms[i].name, planned: arms[i].planned,
    success: a.success, tgi_lower: a.tgi_lo, tgi_median: a.tgi_median, tgi_upper: a.tgi_hi, n_end: a.n_end}));
  const days = DG.DataFrame.fromObjects(named(r, 'days').map((g) => ({group: g.group, day: g.day, lower: g.q05,
    median: g.median, upper: g.q95, animals: g.animals})))!;
  days.col('day')!.meta.units = 'd';
  // volumes to one decimal in every viewer's tooltip and cell (Crux's floats carry a dozen)
  const animalsDf = DG.DataFrame.fromObjects(animals)!;
  inMm3(days, ['lower', 'median', 'upper']);
  inMm3(animalsDf, ['last_volume', 'next_lower', 'next_median', 'next_upper', 'check_lower', 'check_upper']);
  return {
    // (the days to come, for the next measurement from today on: nextMeasurement)
    status: DG.DataFrame.fromObjects([{line: d.line, control: d.control, day: latest, next_day: future[0],
      future_days: future.slice(0, 30).join(','),
      // (the likely range never starts before the day the study is on)
      end_median: study.end_median, end_lower: Number.isFinite(end.q10) ?
        Math.min(Math.max(end.q10, latest), end.q90 ?? Infinity) : end.q10,
      end_upper: end.q90,
      planned_end: d.plan?.predicted.studyEnd ?? null, plan: d.plan?.name ?? null,
      attention: animals.filter((a) => a.attention).length, animals: animals.length,
      // every animal off the study: the study has ended (its header says so, not an expected end)
      ended: animals.length > 0 && animals.every((a) => off(String(a.animal), a.last_day)),
      // the day under way: its animals still to measure, their groups, and the last day every group has
      waiting: animals.filter((a) => a.waiting).length, waiting_groups: waitingGroups.join(', '),
      complete_day: waitingGroups.length ? previous : latest, look_waits: lookDue && !atLook,
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
 *  again"). Only that day: the band is for it, and another day's value has no band yet; a negative value is
 *  impossible, not one to measure again (said apart). */
export function outsideBand(rows: {animal: string, day: number, volume: number}[], nextDay: number,
  bands: Map<string, {check_lower: number | null, check_upper: number | null}>):
  {animal: string, volume: number, lower: number, upper: number}[] {
  const out = [];
  for (const r of rows) {
    const b = bands.get(r.animal);
    if (r.day === nextDay && b && checkReading(r.volume, {lower: b.check_lower, upper: b.check_upper})?.kind ===
      'again')
      out.push({animal: r.animal, volume: r.volume, lower: b.check_lower!, upper: b.check_upper!});
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

/** Why Propose groups cannot place the plan's animals, in words: how many it needs, how many have a tumor volume,
 *  and which of the table's animals have none (the engine's own "60 needed, 22 available" is no one's to read). */
export function shortOfAnimals(needed: number, withVolume: string[], inTable: string[]): string {
  const have = new Set(withVolume);
  const without = [...new Set(inTable)].filter((a) => !have.has(a)).sort(naturalCompare);
  const n = withVolume.length;
  return `The plan needs ${needed} animals, and ${n} ${n === 1 ? 'has' : 'have'} a tumor volume` +
    (without.length ? ` (no volume yet: ${without.slice(0, 8).join(', ')}${without.length > 8 ? '…' : ''})` : '') +
    ': enter the missing animals\' volumes, or size the groups to the animals you have in Plan.';
}

/** The next measurement as Monitor says it: the first of the days to come (`future`, after the last day measured)
 *  from today's study day on, never one already past, and the last planned day before today left unmeasured
 *  (`missed`, null: none). Without today's day (no dated day 0, or a date the data contradicts): the first to come.
 *  (QA 0.1.3 round 3: technician N12, designer W1, chaos 7) */
export function nextMeasurement(future: number[], today: number | null):
  {next: number | undefined, missed: number | null} {
  if (today === null)
    return {next: future[0], missed: null};
  return {next: future.find((d) => d >= today), missed: future.filter((d) => d < today).pop() ?? null};
}

/** What follows an animal's reason as what to do about it ("…: measure again", "…: check it", "…: correct it"). */
const TO_DO = /^(.*): (measure again|check it|correct it|counted as removed early)$/;

/** An animal's reasons to need attention ("; "-joined) as two columns: why, and what to do now (each reason's own
 *  closing words, "Measure again"; none for one that is only a fact). */
export function whyAndStatus(attention: string): {why: string, status: string} {
  const parts = attention.split('; ').filter((s) => s).map((s) => TO_DO.exec(s) ?? [s, s, '']);
  const up = (s: string) => s ? s[0].toUpperCase() + s.slice(1) : s;
  return {why: parts.map((p) => p[1]).join('; '),
    status: up([...new Set(parts.map((p) => p[2]).filter((s) => s))].join('; '))};
}

/** The animal IDs typed or pasted for day 0 (one a line, or separated by tabs, commas or semicolons): each once, in
 *  the order given, those typed twice named (whatever their case: M6 and m6 are one cage card), and a line counting
 *  them against the animals the plan randomizes, saying what to do when short (QA 0.1.3 round 3: technician N5;
 *  round 4 designer D4). */
export function animalIds(text: string, planned: number): {ids: string[], twice: string[], line: string} {
  const all = text.split(/[\n\t,;]+/).map((s) => s.trim()).filter((s) => s);
  const key = (a: string) => a.toLowerCase();
  const first = new Map<string, string>();
  for (const a of all)
    first.set(key(a), first.get(key(a)) ?? a);
  const ids = [...first.values()];
  const spellings = (a: string) => all.filter((x) => key(x) === key(a));
  const twice = ids.filter((a) => spellings(a).length > 1);
  const k = ids.length;
  const count = k < planned ? `${k} of ${planned} – add ${planned - k} or resize the groups in Plan` :
    k > planned ? `${k} of ${planned} – ${k - planned} more than the plan randomizes` : `${k} of ${planned}`;
  const line = [count, ...twice.map((a) => {
    const as = spellings(a);
    const ways = [...new Set(as)];
    return ways.length > 1 ? `${ways.join(' and ')} are the same ID` :
      `${a} appears ${as.length === 2 ? 'twice' : `${as.length} times`}`;
  })].join(' · ');
  return {ids, twice, line};
}

/** What entering day 0 asks when its planned date is not today (null: no date, or today): "Enter day 0 values
 *  (planned for Mon 30 Nov; today is Sat 28 Nov)" (QA 0.1.3 round 4 designer D5: saved under the planned date
 *  without a word). */
export function dayZeroAsk(start: string | null | undefined, now = new Date()): string | null {
  if (!start || !/^\d{4}-\d{2}-\d{2}$/.test(start))
    return null;
  const today = new Date(Date.UTC(now.getFullYear(), now.getMonth(), now.getDate()));
  const planned = new Date(`${start}T00:00:00Z`);
  if (planned.getTime() === today.getTime())
    return null;
  const short = (d: Date) => d.toLocaleDateString('en-GB', {weekday: 'short', day: 'numeric', month: 'short',
    timeZone: 'UTC'}).replace(',', '');
  return `Enter day 0 values (planned for ${short(planned)}; today is ${short(today)})`;
}
