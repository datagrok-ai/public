/* Monitor's engine side: a running study's measurements and its plan as Crux's stats.tumor_monitor reads them, its
   results as the tables pmxMonitorStudy returns, and who needs attention; stats.tumor_allocate for the groups. The
   forecasts and the checks are Crux's; here only the translation and the reasons a person reads. */
import * as DG from 'datagrok-api/dg';

import {AnalysisInputs, analysisDeclarationOf, pairGroups, twinModel} from './analyze';
import {ARCHETYPES} from './archetypes';
import {crux, named, rows} from './crux';
import {LAB_WEEKDAYS, designStudy} from './engine';
import {SavedPlan} from './plan-file';

/** Chance of passing the limit by the next measurement from which an animal is listed. */
export const LIMIT_ALERT = 0.2;
/** A last value this far out of its forecast (two-sided) is listed as unexpected. */
export const UNEXPECTED = 0.005;

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
}

/** The planned reviews from their text: the study days, ascending; anything else is an error. */
export function parseLooks(text: string | null | undefined): number[] {
  const days = (text ?? '').split(/[,;\s]+/).filter((x) => x).map(Number);
  if (days.some((d) => !Number.isFinite(d)))
    throw new Error(`Planned reviews: "${text}" is not a list of study days`);
  return [...new Set(days)].sort((a, b) => a - b);
}

/** The days to come. ponytail: the study's own usual gap, not the plan's weekdays (a 2-3-2 week reads as every 2
 *  days); the calendar when a study's dates are recorded. */
export function futureDays(days: number[]): number[] {
  const seen = [...new Set(days)].sort((a, b) => a - b);
  const gaps = seen.slice(1).map((d, i) => d - seen[i]).sort((a, b) => a - b);
  const gap = gaps.length ? Math.max(gaps[Math.floor(gaps.length / 2)], 1) : 3;
  const last = seen[seen.length - 1];
  return Array.from({length: Math.ceil(60 / gap)}, (_, i) => last + gap * (i + 1));
}

/** Each treated group as an arm of the plan: its kill rate (the design's for the plan's effect) and its therapy's
 *  archetype. A group not in the plan, or a study with no plan: no assumed effect (a kill rate of 0), and no chance of
 *  success is given for it. */
async function armsOf(plan: SavedPlan | null, treated: string[], endpoint: number, humane: number, notes: string[]) {
  const pairs = plan ? pairGroups(treated, plan.predicted.arms, notes) :
    new Map<string, SavedPlan['predicted']['arms'][0]>();
  const kill = new Map<string, number>();
  if (plan) {
    const arms = plan.predicted.arms.map((a) => ({name: a.name, therapy: a.therapy, metric: a.metric, effect: a.effect}));
    const d = await designStudy({line: plan.predicted.line, vRand: plan.predicted.vRand, endpoint, power: 0.8,
      alpha: 0.05, correction: 'dunnett', n: null, arms, preLoss: 0, postLoss: 0, lossRule: 'confident',
      confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0, humane, schedules: false, lastDay: null});
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
export function reasonsOf(a: {p_limit_next?: number | null, last_pit?: number | null}, latest: number, limit: number,
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
  if (lastDay === latest && a.last_pit !== null && a.last_pit !== undefined &&
    (a.last_pit < UNEXPECTED || a.last_pit > 1 - UNEXPECTED))
    out.push('Last value unexpected: measure again');
  return out;
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
  if (day.length === 0)
    throw new Error('The study has no measurements yet');
  const names = [...new Set(group)];
  if (!names.includes(d.control))
    throw new Error(`No group "${d.control}" in the study: the control must be one of ${names.join(', ')}`);
  const notes: string[] = [];
  const arms = await armsOf(d.plan, names.filter((g) => g !== d.control), d.endpoint, d.humane, notes);
  const latest = Math.max(...day);
  const future = futureDays(day);
  const r = await crux('stats.tumor_monitor', {animal, group, day: Float64Array.from(day),
    volume: Float64Array.from(volume)}, {...await twinModel(d.line), control: d.control,
    ...arms.length ? {arms: JSON.stringify(arms)} : {}, future: future.join(','), endpoint: d.endpoint,
    humane: d.humane, alpha: d.alpha, correction: d.correction, ...d.day === null ? {} : {last_day: d.day}});

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
  const animals = rows(r.main as any).map((a) => {
    const l = last.get(String(a.animal))!;
    return {animal: String(a.animal), group: a.group, last_day: l.day, last_volume: l.volume, next_lower: a.next_q05,
      next_median: a.next_median, next_upper: a.next_q95, check_lower: a.next_q005, check_upper: a.next_q995,
      p_limit: a.p_limit_next, limit_day: a.limit_day, last_pit: a.last_pit,
      attention: reasonsOf(a, latest, d.humane, l.day, l.volume, previous).join('; ')};
  });
  const study = (r.named as any).study.value;
  // the chance of success only on a planned review: between them, an empty table (MONITOR_DESIGN.md §3)
  const atLook = looks.includes(latest);
  const decisions = named(r, 'decisions').map((a, i) => ({arm: arms[i].name, planned: arms[i].planned,
    success: a.success, tgi_lower: a.tgi_lo, tgi_median: a.tgi_median, tgi_upper: a.tgi_hi, n_end: a.n_end}));
  const days = DG.DataFrame.fromObjects(named(r, 'days').map((g) => ({group: g.group, day: g.day, lower: g.q05,
    median: g.median, upper: g.q95, animals: g.animals})))!;
  days.col('day')!.meta.units = 'd';
  return {
    status: DG.DataFrame.fromObjects([{line: d.line, control: d.control, day: latest, next_day: future[0],
      end_median: study.end_median, end_lower: study.end_q10, end_upper: study.end_q90,
      planned_end: d.plan?.predicted.studyEnd ?? null, plan: d.plan?.name ?? null,
      attention: animals.filter((a) => a.attention).length, animals: animals.length,
      next_look: looks.find((l) => l > latest) ?? null, at_look: atLook, looks: looks.join(', '), humane: d.humane,
      notes: notes.join(' | ')}])!,
    animals: DG.DataFrame.fromObjects(animals)!,
    days,
    decisions: atLook && decisions.length ? DG.DataFrame.fromObjects(decisions)! : emptyDecisions(),
  };
}

/** Balanced groups for these animals: who goes where, and each group's balance. */
export async function allocateGroups(animal: string[], volume: number[], groups: {name: string, n: number}[],
  seed: number): Promise<{allocation: DG.DataFrame, balance: DG.DataFrame}> {
  const r = await crux('stats.tumor_allocate', {animal, volume: Float64Array.from(volume)},
    {groups: JSON.stringify(groups), seed});
  return {allocation: DG.DataFrame.fromObjects(rows(r.main as any))!,
    balance: DG.DataFrame.fromObjects(named(r, 'balance'))!};
}
