/* The study's strict Kaplan–Meier: each animal's days to the endpoint volume as measured (Crux stats.tumor_survival,
   about 2 ms), an animal that never reached it censored at its last measurement; the curves drawn as steps and each
   treated group's log-rank test against the control. Nothing projected: what is drawn is what was measured. */
import * as DG from 'datagrok-api/dg';

import {crux, named, rows} from './crux';
import type {StudyInputs} from './explore';

/** The curves as a step line (two points per step), the medians and the log-rank tests. */
export interface Survival {
  steps: DG.DataFrame;
  medians: {group: string, median: number | null, n: number, events: number}[];
  logrank: {group: string, statistic: number, p_value: number}[];
}

/** How far before a step its drop starts (days): invisible, but its own x. */
const DROP = 1e-3;

export const KM = {day: 'Day after randomization', share: 'Animals below the endpoint (%)', group: 'Group'};

/** The study's KM on `endpoint` (mm³), `control` first, the groups in `colors`. */
export async function survivalOf(s: StudyInputs, control: string, endpoint: number,
  colors: {[g: string]: string}): Promise<Survival> {
  const [an, gr, dy, vol] = [s.animal, s.group, s.day, s.volume].map((c) => s.study.col(c)!);
  const keep = Array.from({length: s.study.rowCount}, (_, i) => i)
    .filter((i) => Number.isFinite(dy.get(i)) && Number.isFinite(vol.get(i)));
  const r = await crux('stats.tumor_survival', {animal: keep.map((i) => String(an.get(i))),
    group: keep.map((i) => String(gr.get(i))), day: Float64Array.from(keep, (i) => dy.get(i)),
    volume: Float64Array.from(keep, (i) => vol.get(i))}, {control, endpoint_volume: endpoint});
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
  return {steps, medians: named(r, 'medians') as Survival['medians'],
    logrank: named(r, 'logrank') as Survival['logrank']};
}
