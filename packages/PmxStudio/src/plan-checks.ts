/* The plan's predictive checks: what a study run as planned should observe (Crux stats.tumor_plan_checks,
   [decision:tumor-plan-checks]), from the line's model and each arm's kill rate with the animals each group really
   has: the day its median control reaches the endpoint, its controls' growth, each treated arm's TGI and T/C there,
   and its median control day by day, as 95% ranges over the studies. Computed when a plan is saved (and for the
   example study, or an analysis with no plan), not on every replan. This file only translates. */
import {crux, CruxResult, named, rows} from './crux';
import {lineModel, weekdaysParam} from './engine';
import {ARCHETYPES, Therapy} from './archetypes';
import {linePrior} from './priors';

/** 2.5 / 50 / 97.5% quantiles over the studies (q975 null: more than 2.5% never got there). */
export interface Range {q025: number, median: number, q975: number | null}

export interface PlanChecks {
  /** the evaluation day: the first measurement day the median control is at or above the endpoint */
  end: Range;
  /** share of the studies whose controls had not reached the endpoint by the last day simulated */
  notReached: number;
  /** the controls' mean slope of log volume per day to the evaluation day */
  controlSlope: Range;
  /** a study's median control volume (mm³) on each measurement day */
  controls: {day: number, q025: number, median: number, q975: number}[];
  /** each treated arm's TGI and T/C (%) on the evaluation day */
  arms: {name: string, tgi: Range, tc: Range}[];
}

export interface ChecksInput {
  line: string; vRand: number; endpoint: number; humane: number;
  /** the lab's measurement weekdays (0 Monday … 6 Sunday) and the weekday of randomization */
  weekdays: number[]; startWeekday: number;
  /** the design's endpoint day: a study is simulated to twice it */
  endpointDay: number;
  nControl: number;
  /** none: only the controls' and the end's checks */
  arms: {name: string, therapy: Therapy, killRate: number, n: number}[];
}

const range = (r: {[k: string]: any}, k: string): Range => ({q025: r[`${k}_q025`], median: r[`${k}_median`],
  q975: r[`${k}_q975`] ?? null});

export async function planChecks(x: ChecksInput): Promise<PlanChecks> {
  const [{prior, params}, lp] = [await lineModel(x.line), await linePrior(x.line)];
  const shared = {...params, v_rand: x.vRand, endpoint: x.endpoint, humane: x.humane,
    weekdays: weekdaysParam(x.weekdays), start_weekday: x.startWeekday, sd_session: lp.noise?.sd_session ?? 0,
    n_control: Math.max(2, x.nControl), last_day: Math.min(365, Math.max(2 * x.endpointDay, x.endpointDay + 14))};
  // with no treated arm, one with no effect: the controls and the end are the same in every call (one seed)
  const calls = x.arms.length ? x.arms : [{name: '', therapy: 'small_molecule' as Therapy, killRate: 0, n: 2}];
  const results: CruxResult[] = [];
  for (const arm of calls) {
    const a = ARCHETYPES[arm.therapy];
    results.push(await crux('stats.tumor_plan_checks', prior, {...shared, kill_rate: arm.killRate,
      n_treated: Math.max(2, arm.n), onset: a.onsetDays, kappa_sd: a.kappaSd, responders: a.responders}));
  }
  const first = rows(results[0].main as any)[0];
  return {
    end: range(first, 'end'), notReached: first.not_reached, controlSlope: range(first, 'control_slope'),
    controls: named(results[0], 'controls').map((d) => ({day: d.day, q025: d.q025, median: d.median, q975: d.q975})),
    arms: x.arms.map((arm, k) => {
      const r = rows(results[k].main as any)[0];
      return {name: arm.name, tgi: range(r, 'tgi'), tc: range(r, 'tc')};
    }),
  };
}
