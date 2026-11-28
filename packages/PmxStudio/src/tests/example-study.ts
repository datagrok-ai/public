/* A study for the package's tests only, not in the app: MC38 (or any line) run as Plan designs it, with its plan in
 * the pmx.plan tag. Moved from analyze.ts as it was (same numbers, same seed) when the app's examples were removed. */
import * as DG from 'datagrok-api/dg';

import {crux, rows} from '../crux';
import {ENDPOINT, EffectMetric, HUMANE, LAB_WEEKDAYS, SEED, designStudy, lineModel} from '../engine';
import type {Therapy} from '../archetypes';
import {tumourDensity} from '../priors';
import {PLAN_TAG, savedPlan} from '../plan-file';
import {planChecks} from '../plan-checks';

/** The example's body weight, an illustrative shape (Crux stats.tumor_body_weight_simulate; the atlas has no body
 *  weights): per group its median loss at the nadir (%), the stronger arm losing more. */
const EXAMPLE_NADIR = [0, 4, 10];

/** An example's study: its tumor model and treated arms (the control is the vehicle), and the atlas lab drawn
 *  (`seed`, else the package's); `name` names its table and plan. */
export interface ExampleSpec {
  line: string;
  arms: {name: string, therapy: Therapy, metric: EffectMetric, effect: number}[];
  seed?: number;
  name?: string;
}

/** MC38, a vehicle and Cpd-A at 10 and 30 mg/kg (small molecule, TGI 40% and 70%): the example Monitor and Analyze
 *  offer, and the Overview's demo study. */
export const MC38_EXAMPLE: ExampleSpec = {line: 'MC38', arms: [[10, 40], [30, 70]].map(([dose, effect]) =>
  ({name: `Cpd-A ${dose} mg/kg`, therapy: 'small_molecule' as const, metric: 'TGI' as const, effect}))};

/** Plan's form for an example: what Plan reopens it with, and what its study's plan carries. */
export const exampleInputs = (spec: ExampleSpec) => ({line: spec.line, endpoint: ENDPOINT, arms: spec.arms, power: 80,
  correction: 'holm', n: null, vRand: 100, measurementDays: 'Mon, Wed, Fri', startWeekday: 'Monday',
  humane: HUMANE, alpha: 0.05, preLoss: 0, postLoss: 0});

/** An example's design, as Plan makes it from exampleInputs (its end, its animals per arm). */
export const exampleDesign = (spec: ExampleSpec) => designStudy({line: spec.line, vRand: 100, endpoint: ENDPOINT,
  power: 0.8, alpha: 0.05, correction: 'holm', n: null, arms: spec.arms, preLoss: 0, postLoss: 0,
  lossRule: 'confident', confidence: 0.9, weekdays: LAB_WEEKDAYS, startWeekday: 0, humane: HUMANE, schedules: false,
  lastDay: null});

/** An example study run as Plan designs it: MC38, a vehicle and two small-molecule arms (TGI 40% and 70%), each
 *  with the animals Plan asks for its power, measured 3 times a week up to a week past the design's end of study.
 *  One lab of the line's atlas, drawn with the design's seed (whatever it shows: nothing is picked by its outcome),
 *  as Crux stats.tumor_study_simulate measures it: one call per arm with the same seed gives the same lab and the
 *  same controls, and each treated group takes animals of its own (no shared mice). Each animal is weighed with
 *  every measurement: an illustrative body weight (EXAMPLE_NADIR), removed at the weight endpoint like a lab would;
 *  its tumor is weighed when excised, after its last measurement (tumor_weight, on that row). */
export async function exampleStudy(spec: ExampleSpec = MC38_EXAMPLE): Promise<DG.DataFrame> {
  const {line, arms} = spec;
  const seed = spec.seed ?? SEED;
  const design = await exampleDesign(spec);
  const sizes = design.arms.map((a) => a.n);
  if (design.nControl === null || sizes.some((x) => x === null))
    throw new Error('The example needs more than 40 animals per arm');
  const nControl = design.nControl;
  const n = sizes as number[];
  // arm k's treated animals: the k-th run of the draws, so that no two groups share a mouse
  const offsets = n.map((_, k) => n.slice(0, k).reduce((t, x) => t + x, 0));
  const {prior, params} = await lineModel(line);
  // a week past the plan's end: the study is read on the first day its median control reaches the endpoint
  const shared = {...params, v_rand: 100, n: Math.max(nControl, offsets[n.length - 1] + n[n.length - 1]),
    last_day: design.studyEnd + 7, per_week: 3, study: 'prior', reps: 1, seed};
  const out: {animal: string, group: string, day: number, volume: number, body_weight: number,
    tumor_weight: number | null}[] = [];
  // the line's density, or every line's where the atlas has none for it (as tolerabilityFrames)
  const atlas = await tumourDensity();
  const density = (atlas.lines[line] ?? atlas.all).density;
  // every arm simulated at once (each call is seeded: the order they finish in changes nothing)
  const sims = await Promise.all(design.arms.map((a) => crux('stats.tumor_study_simulate', prior, {...shared,
    kill_rate: a.killRate, onset: a.archetype.onsetDays, kappa_sd: a.archetype.kappaSd,
    responders: a.archetype.responders}).then((r) => rows(r.main as any))));
  const take = (k: number, arm: string, from: number, count: number) =>
    sims[k].filter((r) => r.arm === arm && r.animal >= from && r.animal < from + count);
  const groups: [string, {[k: string]: any}[]][] = [['Vehicle', take(0, 'control', 0, nControl)],
    ...design.arms.map((a, k) => [a.declaration.name, take(k, 'treated', offsets[k], n[k])] as [string, any[]])];
  // each group weighed with its measurements (its own seed, g in the groups' order), each tumor weighed once excised
  // (the animal's last measurement, at the line's density); an animal removed at the weight endpoint is no longer
  // measured
  const weighings = await Promise.all(groups.map(([group, rs], g) => crux('stats.tumor_body_weight_simulate',
    {animal: rs.map((r) => `${group} #${r.animal}`), day: Float64Array.from(rs, (r) => r.day),
      volume: Float64Array.from(rs, (r) => r.volume)},
    {nadir: EXAMPLE_NADIR[Math.min(g, EXAMPLE_NADIR.length - 1)], density, seed: seed + g})
    .then((r) => rows(r.main as any))));
  groups.forEach(([group, rs], g) => {
    const weighed = new Map(weighings[g].map((r) => [`${r.animal}|${r.day}`, r]));
    rs.forEach((r) => {
      const animal = `${group} #${r.animal}`;
      const w = weighed.get(`${animal}|${r.day}`);
      if (w !== undefined) {
        out.push({animal, group, day: r.day, volume: r.volume, body_weight: w.weight,
          tumor_weight: Number.isFinite(w.tumor_mass) ? w.tumor_mass : null});
      }
    });
  });
  // lab ids that say nothing of the group (as a blinded study's), in an order drawn from the seed: the same animal
  // has the same id at every look at the study
  let draw = seed;
  const ids = [...new Set(out.map((r) => r.animal))]
    .map((a) => [a, draw = (Math.imul(draw, 1103515245) + 12345) >>> 0] as const).sort((a, b) => a[1] - b[1]);
  const coded = new Map(ids.map(([a], i) => [a, `M${String(i + 1).padStart(3, '0')}`]));
  out.forEach((r) => r.animal = coded.get(r.animal)!);
  const df = DG.DataFrame.fromObjects(out)!;
  df.name = spec.name ?? `Example ${line} study`;
  // its plan in a tag: the plan Plan would save for it, with its predictive checks and each arm's power by n
  const checks = await planChecks({line, vRand: 100, endpoint: ENDPOINT, humane: HUMANE, weekdays: LAB_WEEKDAYS,
    startWeekday: 0, endpointDay: design.endpointDay, nControl, arms: design.arms.map((a, k) => ({
      name: a.declaration.name, therapy: a.declaration.therapy, killRate: a.killRate, n: n[k]}))});
  const plan = savedPlan(spec.name ? `${spec.name} plan` : `Example ${line} plan`, exampleInputs(spec), {line,
    endpoint: ENDPOINT, vRand: 100, studyEnd: design.studyEnd,
    endpointDay: design.endpointDay, durationQ10: design.durationQ10, durationQ90: design.durationQ90, nControl,
    arms: arms.map((a, k) => ({...a, n: n[k], nNeeded: design.arms[k].nAssurance, killRate: design.arms[k].killRate,
      power: design.arms[k].power.flatMap((r) => r.assurance === null ? [] : [{n: r.n, assurance: r.assurance}])})),
    checks});
  df.setTag(PLAN_TAG, JSON.stringify(plan));
  df.col('day')!.meta.units = 'd';
  df.col('volume')!.meta.units = 'mm³';
  df.col('volume')!.meta.format = '0.0';
  df.col('body_weight')!.meta.units = 'g';
  df.col('body_weight')!.meta.format = '0.0';
  df.col('tumor_weight')!.meta.units = 'g';
  df.col('tumor_weight')!.meta.format = '0.00';
  return df;
}
