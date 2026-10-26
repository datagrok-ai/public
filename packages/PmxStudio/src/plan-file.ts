/* A saved plan: a .pmxplan file (JSON) with Plan's declaration, to reopen it, and what Plan predicted, frozen, so
   that Analyze compares a study with the plan as it was made, not with one recomputed on a later atlas. It is saved
   to the user's Home (shared from Files like any file) or downloaded, and opens in Plan from Files or a drop
   (package.ts's file handler); a study table can carry its plan in a tag (the example study does). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Arm, armsOf} from './arms-input';
import {PlanFrames, PlanInputs, declarationOf} from './frames';
import {PlanChecks, planChecks} from './plan-checks';
import {finite} from './engine';

export const PLAN_EXT = 'pmxplan';
/** The tag of a study table that carries its plan (the plan's JSON). */
export const PLAN_TAG = 'pmx.plan';
const FORMAT = 'pmx-plan';

export interface SavedPlan {
  format: typeof FORMAT;
  version: 1;
  name: string;
  /** ISO date and the login of who saved it */
  created: string;
  author: string;
  /** pmxPlanStudy's inputs as its form gives them, the arms as a list: Plan reopens with them */
  inputs: {[param: string]: unknown, arms: Arm[]};
  /** what Plan predicted, as it was shown */
  predicted: {
    line: string; endpoint: number; vRand: number;
    /** the study's last day, the day the median control reaches the endpoint, and 8 in 10 studies' range of it */
    studyEnd: number; endpointDay: number; durationQ10: number; durationQ90: number | null;
    nControl: number | null;
    /** each arm with its animals, the animals it needed for the plan's power, and its power (assurance) by animals
     *  per arm (a plan saved before these were kept has neither) */
    arms: (Arm & {n: number | null, nNeeded?: number | null, power?: PowerAt[]})[];
    /** what a study run as planned should observe, as 95% ranges: the end, the controls' growth, each arm's effect */
    checks: PlanChecks;
  };
}

/** An arm's power with n analysable animals per group: the share of the atlas's labs whose study detects it. */
export interface PowerAt {n: number, assurance: number}

/** The power at n animals, read off an arm's curve (between two of its n, the line between them; beyond its ends,
 *  the end). */
export function powerAt(curve: PowerAt[], n: number): number | null {
  const c = [...curve].sort((a, b) => a.n - b.n);
  if (c.length === 0)
    return null;
  if (n <= c[0].n)
    return c[0].assurance;
  const k = c.findIndex((x) => x.n >= n);
  if (k < 0)
    return c[c.length - 1].assurance;
  const [a, b] = [c[k - 1], c[k]];
  return b.n === n ? b.assurance : a.assurance + (b.assurance - a.assurance) * (n - a.n) / (b.n - a.n);
}

/** A plan from its parts: the inputs (the arms as a table or a list) and what was predicted. */
export function savedPlan(name: string, inputs: {[param: string]: unknown},
  predicted: SavedPlan['predicted']): SavedPlan {
  const arms = Array.isArray(inputs.arms) ? inputs.arms as Arm[] : armsOf(inputs.arms as DG.DataFrame);
  return {format: FORMAT, version: 1, name, created: new Date().toISOString(),
    author: grok.shell.user?.login ?? '', inputs: {...inputs, arms}, predicted};
}

/** The plan Plan shows: the form's inputs (`call`) and its outputs, with its predictive checks for the animals
 *  each group has (an arm that needs more than 40: 40). */
export async function planOf(name: string, call: DG.FuncCall, f: PlanFrames): Promise<SavedPlan> {
  const inputs = Object.fromEntries(call.func.inputs.filter((p) => p.name !== 'compareSchedules')
    .map((p) => [p.name, call.getParamValue(p.name)]));
  const v = (c: string) => f.plan.get(c, 0);
  // (a column with no value reads as an empty string)
  const n = (c: string) => finite(v(c));
  const a = f.armResults;
  const pc = f.powerCurve;
  // the points with a power (Crux gives none where it cannot compute it)
  const curve = (arm: string) => Array.from({length: pc.rowCount}, (_, i) => i)
    .filter((i) => pc.get('arm', i) === arm && finite(pc.get('assurance', i)) !== null)
    .map((i) => ({n: pc.get('n', i), assurance: pc.get('assurance', i)}));
  const arms = armsOf(inputs.arms as DG.DataFrame).map((arm, i) => ({...arm, n: finite(a.get('n', i)),
    nNeeded: finite(a.get('n_needed', i)), power: curve(a.get('arm', i))}));
  const d = declarationOf(inputs as unknown as PlanInputs);
  const checks = await planChecks({line: d.line, vRand: d.vRand, endpoint: d.endpoint, humane: d.humane,
    weekdays: d.weekdays, startWeekday: d.startWeekday, endpointDay: v('endpoint_day'), nControl: n('n_control') ?? 40,
    arms: arms.map((arm, i) => ({name: arm.name, therapy: arm.therapy, killRate: a.get('kill_rate_per_day', i),
      n: arm.n ?? 40}))});
  return savedPlan(name, inputs, {line: v('line'), endpoint: v('endpoint'), vRand: v('v_rand'),
    studyEnd: v('study_end'), endpointDay: v('endpoint_day'), durationQ10: v('duration_q10'),
    durationQ90: n('duration_q90'), nControl: n('n_control'), arms, checks});
}

/** A plan read from a file someone may have edited or sent: checked before use. */
export function parsePlan(text: string): SavedPlan {
  let p: any;
  try {
    p = JSON.parse(text);
  } catch {
    throw new Error('Not a PMX Studio plan: the file is not JSON');
  }
  if (p?.format !== FORMAT)
    throw new Error('Not a PMX Studio plan');
  if (p.version !== 1)
    throw new Error(`A plan of version ${p.version}, which this PMX Studio cannot read: update the package`);
  const fin = (x: unknown) => typeof x === 'number' && Number.isFinite(x);
  const range = (r: any) => fin(r?.q025) && fin(r.median) && (r.q975 === null || fin(r.q975));
  const pr = p.predicted;
  const ch = pr?.checks;
  const ok = pr && typeof pr.line === 'string' && fin(pr.endpoint) && fin(pr.studyEnd) && Array.isArray(pr.arms) &&
    pr.arms.every((a: any) => typeof a?.name === 'string' && ['TGI', 'TC'].includes(a.metric) && fin(a.effect)) &&
    ch && range(ch.end) && range(ch.controlSlope) && Array.isArray(ch.controls) &&
    ch.controls.every((b: any) => fin(b?.day) && fin(b.median)) && Array.isArray(ch.arms) &&
    ch.arms.every((a: any) => typeof a?.name === 'string' && range(a.tgi) && range(a.tc)) &&
    p.inputs && typeof p.inputs === 'object';
  if (!ok)
    throw new Error('The plan file is damaged: save it again from Plan');
  return p as SavedPlan;
}

/** The folder of the user's plans, in their Home. */
export const plansFolder = () => `${grok.shell.user.project.name}:Home/PMX Studio/plans/`;

const fileName = (name: string) => `${name.replace(/[\\/:*?"<>|]+/g, ' ').trim() || 'Plan'}.${PLAN_EXT}`;
const textOf = (p: SavedPlan) => JSON.stringify(p, null, 2);

/** Where a plan of this name is saved, in the user's plans folder. */
export const planPath = (name: string) => plansFolder() + fileName(name);

/** Whether a plan of this name is already saved. */
export const planSaved = (name: string) => grok.dapi.files.exists(planPath(name));

/** Saves the plan to the user's plans folder (replacing one of the same name); returns its path. */
export async function savePlan(p: SavedPlan): Promise<string> {
  const path = planPath(p.name);
  await grok.dapi.files.writeAsText(path, textOf(p));
  return path;
}

/** Gives the plan as a file to save on the computer, to send to someone. */
export const downloadPlan = (p: SavedPlan) => DG.Utils.download(fileName(p.name), textOf(p), 'application/json');

/** The plans in the user's plans folder (none yet: an empty list). */
export async function myPlans(): Promise<DG.FileInfo[]> {
  return await grok.dapi.files.exists(plansFolder()) ?
    await grok.dapi.files.list(plansFolder(), false, PLAN_EXT) : [];
}
