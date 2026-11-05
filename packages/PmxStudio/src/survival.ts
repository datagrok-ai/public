/* The study's strict Kaplan–Meier: each animal's days to the endpoint volume as measured (Crux stats.tumor_survival,
   about 2 ms), an animal that never reached it censored at its last measurement; the curves drawn as steps and each
   treated group's log-rank test against the control. Projected (on request, labelled as such): each animal still
   under the endpoint carried on at the log-linear rate of its last measurements to a horizon (projectedRows, closed
   form), the same Crux KM over it; the log-rank tests stay the measured study's. */
import * as DG from 'datagrok-api/dg';

import {crux, named, rows} from './crux';
import type {StudyInputs} from './explore';

/** The curves as a step line (two points per step), the medians and the log-rank tests. */
export interface Survival {
  steps: DG.DataFrame;
  medians: {group: string, median: number | null, n: number, events: number}[];
  logrank: {group: string, statistic: number, p_value: number}[];
  /** the last day measured; projected: the horizon followed to */
  lastDay: number;
  horizon?: number;
}

/** A measurement row, as tumor_survival takes them. */
export interface VolumeRow {animal: string, group: string, day: number, volume: number}

/** How many of an animal's last measurements give its rate. */
const RATE_FROM = 3;

/** The rows that carry each animal still under `endpoint` on to `horizon` at the log-linear rate of its last
 *  measurements (least squares of ln V on the last RATE_FROM): a row on the day it would reach the endpoint, if it
 *  does by the horizon; else one on the horizon, under it (a tumor that does not grow stays as it is). An animal with
 *  fewer than two measurements, or past the endpoint already, adds nothing. */
export function projectedRows(rows: VolumeRow[], endpoint: number, horizon: number): VolumeRow[] {
  const by = new Map<string, VolumeRow[]>();
  for (const r of rows) {
    if (Number.isFinite(r.day) && r.volume > 0)
      by.set(r.animal, [...by.get(r.animal) ?? [], r]);
  }
  const out: VolumeRow[] = [];
  for (const own of by.values()) {
    own.sort((a, b) => a.day - b.day);
    const last = own[own.length - 1];
    if (own.some((r) => r.volume >= endpoint) || own.length < 2 || horizon <= last.day)
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
  const [an, gr, dy, vol] = [s.animal, s.group, s.day, s.volume].map((c) => s.study.col(c)!);
  const keep = Array.from({length: s.study.rowCount}, (_, i) => i)
    .filter((i) => Number.isFinite(dy.get(i)) && Number.isFinite(vol.get(i)));
  const measured: VolumeRow[] = keep.map((i) => ({animal: String(an.get(i)), group: String(gr.get(i)),
    day: dy.get(i), volume: vol.get(i)}));
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
  steps.col(KM.group)!.setCategoryOrder([control, ...Object.keys(colors).filter((g) => g !== control)]);
  steps.col(KM.share)!.meta.format = '0';
  // what lies beyond it is projected
  if (horizon !== undefined)
    steps.meta.formulaLines.addLine({title: `Last measured, day ${lastDay}`, formula: `\${${KM.day}} = ${lastDay}`});
  return {steps, medians: named(r, 'medians') as Survival['medians'],
    logrank: named(r, 'logrank') as Survival['logrank'], lastDay, horizon};
}
