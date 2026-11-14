/* The study's strict Kaplan–Meier: each animal's days to the endpoint volume as measured (Crux stats.tumor_survival,
   about 2 ms), an animal that never reached it censored at its last measurement; the curves drawn as steps and each
   treated group's log-rank test against the control. Projected (on request, labelled as such): each animal still
   under the endpoint carried on at the log-linear rate of its last measurements to a horizon (projectedRows, closed
   form), the same Crux KM over it; the log-rank tests stay the measured study's. */
import * as DG from 'datagrok-api/dg';

import {crux, named, rows} from './crux';
import {numbersOf, type StudyInputs} from './explore';

/** The curves as a step line (two points per step), the medians and the log-rank tests. */
export interface Survival {
  steps: DG.DataFrame;
  medians: {group: string, median: number | null, n: number, events: number}[];
  logrank: {group: string, statistic: number, p_value: number}[];
  /** the last day measured; projected: the horizon followed to */
  lastDay: number;
  horizon?: number;
  /** the groups, the control first, and their colours */
  groups: string[];
  colors: {[g: string]: string};
}

/** A treated group against the control with the animals removed for another reason as a competing event. */
export interface Competing {group: string, shr: number, lower: number, upper: number, p: number, competing: number}

/** Fine-Gray (Crux survival.fine_gray, as cmprsk's crr): the endpoint is the event, an animal removed before it for
 *  another reason (`removal`: the column of the reason, any text) competes with it, the rest are censored at their
 *  last measurement; each treated group against the control, its subdistribution hazard ratio (under 1: reaches
 *  the endpoint less). null: no removal column, or no animal removed for another reason. */
export async function competingOf(s: StudyInputs, removal: string | null, control: string, endpoint: number):
  Promise<Competing[] | null> {
  const why = removal ? s.study.col(removal) : null;
  if (!why)
    return null;
  const [an, gr] = [s.animal, s.group].map((c) => s.study.col(c)!);
  const [dy, vol] = [s.day, s.volume].map((c) => numbersOf(s.study.col(c)!));
  type Animal = {group: string, last: number, reached: number | null, removed: boolean};
  const animals = new Map<string, Animal>();
  for (let i = 0; i < s.study.rowCount; i++) {
    const d = dy[i];
    if (!Number.isFinite(d))
      continue;
    const a = String(an.get(i));
    const x = animals.get(a) ?? {group: String(gr.get(i)), last: -Infinity, reached: null, removed: false};
    x.last = Math.max(x.last, d);
    if (vol[i] >= endpoint)
      x.reached = x.reached === null ? d : Math.min(x.reached, d);
    // (a status column says of most animals that they finished: that is no removal)
    const r = why.get(i);
    if (r !== null && r !== undefined && String(r).trim() !== '' && !FINISHED.test(String(r)))
      x.removed = true;
    animals.set(a, x);
  }
  const all = [...animals.values()];
  const status = (x: Animal) => x.reached !== null ? '1' : x.removed ? '2' : '0';
  if (!all.some((x) => status(x) === '2'))
    return null;
  // one fit per treated group, all at once
  const fits = await Promise.all([...new Set(all.map((x) => x.group))].filter((g) => g !== control).map(async (g) => {
    const sub = all.filter((x) => x.group === control || x.group === g);
    try {
      const r = await crux('survival.fine_gray', {duration: Float64Array.from(sub, (x) => x.reached ?? x.last),
        status: sub.map(status), treated: Float64Array.from(sub, (x) => x.group === g ? 1 : 0)},
      {failcode: '1', cencode: '0'});
      const k = rows(r.main as any)[0];
      return {group: g, shr: k.subdistribution_hazard_ratio, lower: k.shr_ci_lower, upper: k.shr_ci_upper, p: k.p,
        competing: sub.filter((x) => status(x) === '2').length} as Competing;
    } catch {
      return null; // a pair the fit cannot take (no event in it): left out
    }
  }));
  return fits.filter((x): x is Competing => x !== null);
}

/** A measurement row, as tumor_survival takes them. */
export interface VolumeRow {animal: string, group: string, day: number, volume: number}

/** A removal column's entry for an animal that was not removed: it finished the study (what Tolerability reads
 *  as not removed too). */
const FINISHED = /\b(none|completed?|terminal|scheduled)\b|end of study/i;

/** How many of an animal's last measurements give its rate. */
const RATE_FROM = 3;

/** The rows that carry each animal still under `endpoint` on to `horizon` at the log-linear rate of its last
 *  measurements (least squares of ln V on the last RATE_FROM): a row on the day it would reach the endpoint, if it
 *  does by the horizon; else one on the horizon, under it (a tumor that does not grow stays as it is). An animal with
 *  fewer than two measurements, or past the endpoint already, adds nothing. */
export function projectedRows(rows: VolumeRow[], endpoint: number, horizon: number): VolumeRow[] {
  const by = new Map<string, VolumeRow[]>();
  for (const r of rows) {
    if (Number.isFinite(r.day) && r.volume > 0) {
      const own = by.get(r.animal);
      own ? own.push(r) : by.set(r.animal, [r]);
    }
  }
  const out: VolumeRow[] = [];
  // an animal that left before the study's last measurement (found dead, removed) is not carried on: it was not
  // there to be followed
  const studyLast = Math.max(...[...by.values()].flat().map((r) => r.day));
  for (const own of by.values()) {
    own.sort((a, b) => a.day - b.day);
    const last = own[own.length - 1];
    if (own.some((r) => r.volume >= endpoint) || own.length < 2 || horizon <= last.day || last.day < studyLast)
      continue;
    const fit = own.slice(-RATE_FROM);
    const mx = fit.reduce((s, r) => s + r.day, 0) / fit.length;
    const my = fit.reduce((s, r) => s + Math.log(r.volume), 0) / fit.length;
    const sxx = fit.reduce((s, r) => s + (r.day - mx) ** 2, 0);
    const rate = sxx > 0 ? fit.reduce((s, r) => s + (r.day - mx) * (Math.log(r.volume) - my), 0) / sxx : 0;
    const reach = rate > 0 ? last.day + Math.log(endpoint / last.volume) / rate : Infinity;
    out.push(reach <= horizon ? {...last, day: reach, volume: endpoint} :
      {...last, day: horizon, volume: last.volume * Math.exp(Math.max(rate, 0) * (horizon - last.day))});
  }
  return out;
}

/** How far before a step its drop starts (days): invisible, but its own x. */
const DROP = 1e-3;

export const KM = {day: 'Day after randomization', share: 'Animals below the endpoint (%)', group: 'Group'};

/** The study's KM on `endpoint` (mm³), `control` first, the groups in `colors`; `project`: carried on to its
 *  horizon (null: the last day measured plus two weeks). */
export async function survivalOf(s: StudyInputs, control: string, endpoint: number,
  colors: {[g: string]: string}, project?: {horizon: number | null}): Promise<Survival> {
  const [an, gr] = [s.animal, s.group].map((c) => s.study.col(c)!);
  const [dy, vol] = [s.day, s.volume].map((c) => numbersOf(s.study.col(c)!));
  const keep = Array.from({length: s.study.rowCount}, (_, i) => i)
    .filter((i) => Number.isFinite(dy[i]) && Number.isFinite(vol[i]));
  const measured: VolumeRow[] = keep.map((i) => ({animal: String(an.get(i)), group: String(gr.get(i)),
    day: dy[i], volume: vol[i]}));
  const lastDay = measured.reduce((m, x) => Math.max(m, x.day), 0);
  const horizon = project ? project.horizon ?? lastDay + 14 : undefined;
  const all = horizon === undefined ? measured : [...measured, ...projectedRows(measured, endpoint, horizon)];
  const r = await crux('stats.tumor_survival', {animal: all.map((x) => x.animal), group: all.map((x) => x.group),
    day: Float64Array.from(all, (x) => x.day), volume: Float64Array.from(all, (x) => x.volume)},
  {control, endpoint_volume: endpoint});
  // each group from 100% on day 0, each step flat to its time, then down: the line chart averages the points of one
  // x, so the drop's top sits a hair (DROP) before it
  const out: {[c: string]: any}[] = [];
  const last = new Map<string, number>();
  for (const k of rows(r.main as any)) {
    if (!last.has(k.group))
      out.push({[KM.group]: k.group, [KM.day]: 0, [KM.share]: 100});
    const before = last.get(k.group) ?? 100;
    if (k.time > 0)
      out.push({[KM.group]: k.group, [KM.day]: k.time - DROP, [KM.share]: before});
    out.push({[KM.group]: k.group, [KM.day]: k.time, [KM.share]: 100 * k.survival});
    last.set(k.group, 100 * k.survival);
  }
  const steps = DG.DataFrame.fromObjects(out)!;
  steps.col(KM.group)!.meta.colors.setCategorical(colors);
  const groups = [control, ...Object.keys(colors).filter((g) => g !== control)];
  steps.col(KM.group)!.setCategoryOrder(groups);
  steps.col(KM.share)!.meta.format = '0';
  // what lies beyond it is projected
  if (horizon !== undefined)
    steps.meta.formulaLines.addLine({title: `Last measured, day ${lastDay}`, formula: `\${${KM.day}} = ${lastDay}`});
  return {steps, medians: named(r, 'medians') as Survival['medians'],
    logrank: named(r, 'logrank') as Survival['logrank'], lastDay, horizon, groups, colors};
}
