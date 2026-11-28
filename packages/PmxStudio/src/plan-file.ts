/* A saved plan: a .pmxplan file (JSON) with Plan's declaration, to reopen it, and what Plan predicted, frozen, so
   that Analyze compares a study with the plan as it was made, not with one recomputed on a later atlas. It is saved
   to the user's Home (shared from Files like any file) or downloaded, and opens in Plan from Files or a drop
   (package.ts's file handler); a study table can carry its plan in a tag (the example study does). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Arm, armsOf} from './arms-input';
import {noBidi, parseArmName} from './arm-name';
import {PlanFrames, PlanInputs, declarationOf} from './frames';
import {PlanChecks, planChecks} from './plan-checks';
import {MAX_ARMS, StudyDeclaration, finite, isFiniteNumber as fin, tgiProblem, tooManyArms} from './engine';
import {StudyMeta, onlyGroups, parseMeta} from './study-meta';
import {_package} from './package';
import {priorsText} from './priors';
import {cruxManifestText} from './crux';
import {canonical, sameJson} from './format';
import {SavedHistory, historyOf, registerHistory, savedHistory} from './history-models';

export const PLAN_EXT = 'pmxplan';
/** The tag of a study table that carries its plan (the plan's JSON). */
export const PLAN_TAG = 'pmx.plan';
const FORMAT = 'pmx-plan';
export {FORMAT as PLAN_FORMAT};

/** What a file was made with: the sha256 of the atlas file (hex) and the app's version. Both only say where the
 *  file's numbers come from; they change no calculation. */
export interface Provenance {atlas: string, app: string}

/** The atlas file's hash: its text with LF line ends (as git stores it, so a CRLF checkout gives the same), sha256. */
export const atlasHashOf = (text: string) => sha256Hex(text.replace(/\r\n/g, '\n'));

/** The sha256 of a text's UTF-8 bytes, in hex. */
export async function sha256Hex(text: string): Promise<string> {
  const d = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
  return Array.from(new Uint8Array(d), (b) => b.toString(16).padStart(2, '0')).join('');
}

let atlasHash: Promise<string> | null = null;

/** The current atlas and app, the hash computed once; a failed read is not kept. The atlas is hashed with LF line
 *  ends, as git stores it, so a package built on Windows (CRLF checkout) gives the same hash. */
export async function provenance(): Promise<Provenance> {
  atlasHash ??= priorsText().then(atlasHashOf).catch((e) => {
    atlasHash = null;
    throw e;
  });
  return {atlas: await atlasHash, app: _package.version};
}

/** A file's provenance, as far as it is well formed (a file without it, or with it damaged, has none). */
export const provenanceOf = (x: any): Provenance | undefined =>
  typeof x?.atlas === 'string' && typeof x.app === 'string' ? {atlas: x.atlas, app: x.app} : undefined;

/** What a file made with another atlas says: Plan recomputes what it opens, Monitor keeps the plan's predictions. */
export const PLAN_ATLAS_NOTE = 'Saved with a different version of the atlas; the numbers shown are computed with ' +
  'the current one.';
export const FROZEN_ATLAS_NOTE = 'Made with a different version of the atlas: the plan\'s effects and end are kept ' +
  'as saved; the forecasts use the current atlas.';

/** `note` for a file made with another atlas than the current one; null when it is the same, says nothing, or the
 *  atlas cannot be read. */
export async function atlasNote(p: Provenance | undefined, note: string): Promise<string | null> {
  const now = await provenance().catch(() => undefined);
  return p && now && p.atlas !== now.atlas ? note : null;
}

export interface SavedPlan {
  format: typeof FORMAT;
  version: 1;
  name: string;
  /** ISO date and the login of who saved it */
  created: string;
  author: string;
  /** the study's ID, owner, day 0 and what each group receives (older plans: none) */
  study?: StudyMeta;
  /** the atlas and app it was made with (older plans: none) */
  provenance?: Provenance;
  /** the history of a model not in the atlas (its line is the key `new:pdx:<id>`; atlas models: none) */
  history?: SavedHistory;
  /** pmxPlanStudy's inputs as its form gives them, the arms as a list: Plan reopens with them */
  inputs: {[param: string]: unknown, arms: Arm[]};
  /** what Plan predicted, as it was shown */
  predicted: {
    line: string; endpoint: number; vRand: number;
    /** the study's last day, the day the median control reaches the endpoint, and 8 in 10 studies' range of it */
    studyEnd: number; endpointDay: number; durationQ10: number; durationQ90: number | null;
    nControl: number | null;
    /** the control's animals to randomize, the analyzable ones and those lost after randomization (older plans:
     *  none, read as the analyzable) */
    randomiseControl?: number | null;
    /** each arm with its animals (`randomise`: to randomize, as randomiseControl), the animals it needed for the
     *  plan's power, its power (assurance) by animals per arm, and the kill rate that gives its declared effect: Monitor's chance of success rests on that rate, so
     *  it is kept as planned (a plan saved before these were kept has none of them) */
    arms: (Arm & {n: number | null, randomise?: number | null, nNeeded?: number | null, power?: PowerAt[],
      killRate?: number | null})[];
    /** what a study run as planned should observe, as 95% ranges: the end, the controls' growth, each arm's effect */
    checks: PlanChecks;
    /** the plan as Plan showed it, to show it again at once: its tables (d42, base64, by PlanFrames key), the inputs
     *  they were computed from (the file's own, as JSON) and the declaration Crux was given. Plan recomputes when
     *  the inputs, the atlas or the app differ (a plan saved before these were kept has none of them: computed as it
     *  opens) */
    shown?: SavedShown;
  };
}

export interface SavedShown {frames: {[key: string]: string}, asked: string, declared: StudyDeclaration,
  /** what computed them (engineFingerprint): reopened elsewhere, they are computed again */
  engine: string}

/** Bumped when this package's maths of a plan (frames.ts, engine.ts) change what its tables hold: tables kept by
 *  an earlier one are computed again. */
export const SHOWN_VERSION = 1;

let engineHash: Promise<string> | null = null;

/** What a plan's tables are computed with, as one hash: SHOWN_VERSION, the app's version, the atlas, the bundled
 *  Crux modules (their manifest carries the wasm's hash), the end's band calibration and the PDX history; computed
 *  once, a failed read not kept. */
export function engineFingerprint(): Promise<string> {
  const lf = (t: string) => t.replace(/\r\n/g, '\n');
  engineHash ??= Promise.all([provenance(), cruxManifestText(), _package.files.readAsText('band_calibration_v1.json'),
    _package.files.readAsText('pdx_history_v1.json')]).then(([p, ...texts]) =>
    sha256Hex([SHOWN_VERSION, p.app, p.atlas, ...texts.map(lf)].join('\n'))).catch((e) => {
    engineHash = null;
    throw e;
  });
  return engineHash;
}

/** A saved plan's text for a study table's tag: without the tables the .pmxplan keeps (the study reads its plan's
 *  predictions, not its tables). */
export function planTag(p: SavedPlan | string): string {
  const x = typeof p === 'string' ? (() => {
    try {
      return JSON.parse(p);
    } catch {
      return null;
    }
  })() : p;
  if (!x?.predicted?.shown)
    return typeof p === 'string' ? p : JSON.stringify(p);
  const {shown: _, ...predicted} = x.predicted;
  return JSON.stringify({...x, predicted});
}

const FRAME_KEYS: (keyof PlanFrames)[] = ['plan', 'armResults', 'powerCurve', 'curves', 'duration', 'schedules',
  'calendar', 'scenarios'];

/** A plan's tables as a file keeps them. */
export function shownOf(f: PlanFrames, asked: string, declared: StudyDeclaration, engine: string): SavedShown {
  const b64 = (b: Uint8Array) => {
    let s = '';
    for (let i = 0; i < b.length; i += 0x8000)
      s += String.fromCharCode(...b.subarray(i, i + 0x8000));
    return btoa(s);
  };
  return {frames: Object.fromEntries(FRAME_KEYS.map((k) => [k, b64(f[k].toByteArray())])), asked, declared, engine};
}

/** The tables a file kept, read back; null when any is missing, does not read or lacks what Plan reads first (Plan
 *  then computes them). */
export function framesOf(s: SavedShown | undefined): PlanFrames | null {
  if (!s?.frames || typeof s.asked !== 'string' || !s.declared || typeof s.engine !== 'string' ||
    !FRAME_KEYS.every((k) => typeof s.frames[k] === 'string'))
    return null;
  try {
    const f = Object.fromEntries(FRAME_KEYS.map((k) => [k, DG.DataFrame.fromByteArray(
      Uint8Array.from(atob(s.frames[k]), (c) => c.charCodeAt(0)))])) as unknown as PlanFrames;
    const has = (df: DG.DataFrame, cols: string[]) => df.rowCount > 0 && cols.every((c) => !!df.col(c));
    return has(f.plan, ['line', 'study_end', 'endpoint', 'n_control']) && has(f.armResults, ['arm']) ? f : null;
  } catch {
    return null;
  }
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
  predicted: SavedPlan['predicted'], study?: StudyMeta): SavedPlan {
  const arms = Array.isArray(inputs.arms) ? inputs.arms as Arm[] : armsOf(inputs.arms as DG.DataFrame);
  return {format: FORMAT, version: 1, name, created: new Date().toISOString(),
    author: grok.shell.user?.login ?? '', ...study ? {study} : {}, inputs: {...inputs, arms}, predicted};
}

/** The plan Plan shows: the form's inputs (`call`) and its outputs, with its predictive checks for the animals
 *  each group has (an arm that needs more than 40: 40). */
export async function planOf(name: string, call: DG.FuncCall, f: PlanFrames, study?: StudyMeta):
  Promise<SavedPlan> {
  const inputs = Object.fromEntries(call.func.inputs
    .filter((p) => p.name !== 'compareSchedules' && p.name !== 'detectableEffect')
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
  // (each arm's kind as planned: one read from its name follows the model's host, frames.ts planFrames)
  const arms = armsOf(inputs.arms as DG.DataFrame).map((arm, i) => ({...arm,
    therapy: (a.col('therapy_key') ? a.get('therapy_key', i) : arm.therapy) as Arm['therapy'],
    n: finite(a.get('n', i)), randomise: a.col('randomise') ? finite(a.get('randomise', i)) : null,
    nNeeded: finite(a.get('n_needed', i)), power: curve(a.get('arm', i)),
    killRate: finite(a.get('kill_rate_per_day', i))}));
  const d = declarationOf(inputs as unknown as PlanInputs);
  const checks = await planChecks({line: d.line, vRand: d.vRand, endpoint: d.endpoint, humane: d.humane,
    weekdays: d.weekdays, startWeekday: d.startWeekday, endpointDay: v('endpoint_day'), nControl: n('n_control') ?? 40,
    arms: arms.map((arm, i) => ({name: arm.name, therapy: arm.therapy, killRate: a.get('kill_rate_per_day', i),
      n: arm.n ?? 40}))});
  const history = savedHistory(d.line);
  return {...savedPlan(name, inputs, {line: v('line'), endpoint: v('endpoint'), vRand: v('v_rand'),
    studyEnd: v('study_end'), endpointDay: v('endpoint_day'), durationQ10: v('duration_q10'),
    durationQ90: n('duration_q90'), nControl: n('n_control'),
    randomiseControl: f.plan.col('randomise_control') ? n('randomise_control') : null, arms, checks}, study),
  ...await provenance().then((x) => ({provenance: x}), () => ({})), ...history ? {history} : {}};
}

/** A plan as one text whatever its file's spacing, its keys' order or its study's details since: its name and inputs,
 *  each object's keys sorted (a plan read back from its row's json comes with its keys reordered) (StudyLink.plan). */
export const planKey = (p: {name: string, inputs: unknown} | null | undefined): string | null =>
  p ? JSON.stringify(canonical([p.name, p.inputs])) : null;
/** Whether the inputs a plan's kept tables were computed from (`asked`, their JSON text) are `inputs`. */
export function sameAsked(asked: string, inputs: unknown): boolean {
  try {
    return sameJson(JSON.parse(asked), inputs);
  } catch {
    return false;
  }
}

/** A study's plan text as Plan opens it from that study: with the study's own ID, not the one its plan came with (QA
 *  0.1.3 round 4 chaos N3: a study saved from the demo opened Plan as "… · PMX-DEMO-CT26"). Text not read: as is. */
export function planOfStudy<T extends string | null | undefined>(text: T, id: string | null | undefined): T | string {
  if (!text || !id)
    return text;
  try {
    const p = JSON.parse(text);
    if (p?.study && p.study.id !== id)
      return JSON.stringify({...p, study: {...p.study, id}});
  } catch {} // (a damaged plan: Plan says so)
  return text;
}

/** The plan key of a study table's plan (null: none, or not readable). */
export function planKeyOf(df: DG.DataFrame | null | undefined): string | null {
  try {
    const text = df?.getTag(PLAN_TAG);
    return text ? planKey(parsePlan(text)) : null;
  } catch {
    return null;
  }
}

/** The groups a plan randomizes, the control as Vehicle: each with its animals to randomize (the analyzable and those
 *  lost after randomization, as Plan's "to randomize" counts them), else, a plan saved before these were kept, its
 *  analyzable ones. The day-0 IDs and Propose groups both count by it (QA 0.1.3 round 4 designer D4: 60 asked, 78
 *  planned). */
export function randomizedGroups(plan: SavedPlan): {name: string, n: number}[] {
  const p = plan.predicted;
  const many = (r: number | null | undefined, n: number | null | undefined) => fin(r) && r > 0 ? r : n ?? 0;
  return [{name: 'Vehicle', n: many(p.randomiseControl, p.nControl)},
    ...p.arms.map((a) => ({name: a.name, n: many(a.randomise, a.n)}))];
}

/** A plan read from a file someone may have edited or sent: checked before use. */
export function parsePlan(text: string): SavedPlan {
  let p: any;
  try {
    // (no text-direction controls in any of its names: they reverse them on screen and in the documents)
    p = JSON.parse(noBidi(text, true));
  } catch {
    throw new Error('Not a PMX Studio plan: the file cannot be read (cut short, or another kind of file). ' +
      'Open the .pmxplan file saved from Plan.');
  }
  if (p?.format !== FORMAT)
    throw new Error('Not a PMX Studio plan. Open the .pmxplan file saved from Plan.');
  if (p.version !== 1) {
    throw new Error('This plan was saved by another version of PMX Studio, which this one cannot read. Ask your ' +
      'administrator to update the PMX Studio package.');
  }
  const range = (r: any) => fin(r?.q025) && fin(r.median) && (r.q975 === null || fin(r.q975));
  // the end: none of it when the median control does not reach the endpoint in the study (a slow model, a high
  // endpoint: Panc02 at 2000 mm³), as the engine says it
  const open = (x: unknown) => x === null || fin(x);
  const endOk = (r: any) => !!r && open(r.q025) && open(r.median) && open(r.q975);
  const pr = p.predicted;
  const ch = pr?.checks;
  const ok = pr && typeof pr.line === 'string' && fin(pr.endpoint) && fin(pr.studyEnd) && Array.isArray(pr.arms) &&
    pr.arms.every((a: any) => typeof a?.name === 'string' && ['TGI', 'TC'].includes(a.metric) && fin(a.effect)) &&
    ch && endOk(ch.end) && range(ch.controlSlope) && Array.isArray(ch.controls) &&
    ch.controls.every((b: any) => fin(b?.day) && fin(b.median)) && Array.isArray(ch.arms) &&
    ch.arms.every((a: any) => typeof a?.name === 'string' && range(a.tgi) && range(a.tc)) &&
    p.inputs && typeof p.inputs === 'object';
  if (!ok) {
    throw new Error('The plan file is damaged: some of what Plan predicted is missing. Open another copy of it, or ' +
      'make the plan again in Plan.');
  }
  // values Plan could not have saved (a file edited by hand): refused with what is wrong, never shown as "NaN%" or
  // cut to fit (QA 0.1.3 final designer #2, #11)
  const bad = badValues(p);
  if (bad)
    throw new Error(`This plan cannot be opened: ${bad}. Correct the file, or make the plan again in Plan.`);
  // (a plan saved before a removed arm left its details: its doses dropped here, for every reader)
  const study = parseMeta(p.study);
  if (study)
    p.study = onlyGroups(study, ['Control', ...pr.arms.map((a: any) => a.name)]);
  else
    delete p.study;
  const prov = provenanceOf(p.provenance);
  if (prov)
    p.provenance = prov;
  else
    delete p.provenance;
  // a model not in the atlas: its history registered, so every stage finds it by its key (the key the history gives
  // now, should the file have been edited)
  const history = historyOf(p.history);
  if (history) {
    p.history = history;
    const key = registerHistory(history);
    if (p.inputs.line === pr.line)
      p.inputs.line = key;
    pr.line = key;
  } else
    delete p.history;
  return p as SavedPlan;
}

/** What in a plan's file Plan could not have saved, in words ('' when nothing): more arms than a study holds, an
 *  endpoint or volume at randomization not above 0, an arm's TGI not above 0 and below complete regression (tgiMax)
 *  or its T/C not above 0 to 100%. */
export function badValues(p: any): string {
  const arms: any[] = Array.isArray(p?.inputs?.arms) ? p.inputs.arms : [];
  const most = Math.max(arms.length, Array.isArray(p?.predicted?.arms) ? p.predicted.arms.length : 0);
  if (most > MAX_ARMS) {
    const extra = arms.slice(MAX_ARMS).map((a) => `"${a?.name ?? ''}"`).join(', ');
    return tooManyArms(most).replace(/\.$/, '') + (extra ? ` (past the ${MAX_ARMS}th: ${extra})` : '');
  }
  const positive = (x: unknown, what: string) => x != null && x !== '' && !(Number(x) > 0) ?
    `the ${what} is ${x}, not above 0` : '';
  // (Plan predicts for the model of its inputs: two models are a file edited by hand, its preview and Plan disagreeing)
  if (typeof p?.inputs?.line === 'string' && p.inputs.line !== p?.predicted?.line)
    return `its model is ${p.inputs.line}, but its predictions are for ${p.predicted?.line}`;
  const said = positive(p?.inputs?.endpoint, 'endpoint volume') ||
    positive(p?.predicted?.endpoint, 'endpoint volume') ||
    positive(p?.inputs?.vRand, 'volume at randomization');
  if (said)
    return said;
  // a TGI above 0 and below complete regression (tgiMax: past 100% the tumors shrink), a T/C above 0 to 100%
  const [endpoint, v0] = [Number(p?.inputs?.endpoint ?? p?.predicted?.endpoint ?? 1000),
    Number(p?.inputs?.vRand ?? p?.predicted?.vRand ?? 100)];
  const why = (a: any) => a?.metric === 'TC' ? Number(a?.effect) > 0 && Number(a?.effect) <= 100 ? null :
    'it must be above 0 and at most 100%' : tgiProblem(Number(a?.effect), endpoint, v0)
    ?.replace('Expected TGI: ', 'it must be ') ?? null;
  const effect = [...arms, ...p?.predicted?.arms ?? []].find((a) => why(a) !== null);
  if (effect)
    return `arm "${effect.name}": expected ${effect.metric === 'TC' ? 'T/C' : 'TGI'} is ${effect.effect}%, ` +
      why(effect);
  // a dose below 0 (or no number), in an arm's name or a group's details: Plan's own check, its words (QA 0.1.3
  // round 4 technician R4-3: dropped silently, "dose from the group details")
  const named = arms.flatMap((a) => parseArmName(String(a?.name ?? '')).filter((x) => x.negative))[0];
  if (named)
    return `${named.drug}: a dose cannot be below 0`;
  const dosed = (Array.isArray(p?.study?.groups) ? p.study.groups : []).find((g: any) => g?.dose != null &&
    !(typeof g.dose === 'number' && Number.isFinite(g.dose) && g.dose >= 0));
  return dosed ? `group "${dosed.group}": a dose cannot be below 0 (it is ${dosed.dose})` : '';
}

/** The folder of the user's plans, in their Home. */
export const plansFolder = () => `${grok.shell.user.project.name}:Home/PMX Studio/plans/`;
/** A path in the user's Home as the app shows it: "your Home/…", never the account's login (a shared screen, a
 *  recorded demo) */
export const shownPath = (path: string) => {
  const own = `${grok.shell.user.project.name}:Home/`;
  // another project's Home (shared with the user) keeps its name
  return path.startsWith(own) ? 'your Home/' + path.slice(own.length) : path;
};

/** A name as a file name: without the characters a path cannot hold, `fallback` when nothing is left. */
export const fileNameOf = (name: string, fallback: string, ext: string) =>
  `${name.replace(/[\\/:*?"<>|]+/g, ' ').trim() || fallback}.${ext}`;
const fileName = (name: string) => fileNameOf(name, 'Plan', PLAN_EXT);
const textOf = (p: SavedPlan) => JSON.stringify(p, null, 2);

/** Where a plan of this name is saved, in the user's plans folder. */
export const planPath = (name: string) => plansFolder() + fileName(name);

/** Whether a plan of this name is already saved. */
export const planSaved = (name: string) => grok.dapi.files.exists(planPath(name));

/** Saves the plan to the user's plans folder (replacing one of the same name); returns its path. */
export async function savePlan(p: SavedPlan): Promise<string> {
  const path = planPath(p.name);
  await grok.dapi.files.writeAsText(path, textOf(p)); // (its index row, now older than the file, is read again)
  return path;
}

/** Gives the plan as a file to save on the computer, to send to someone. */
export const downloadPlan = (p: SavedPlan) => DG.Utils.download(fileName(p.name), textOf(p), 'application/json');

/** The plans in the user's plans folder (none yet: an empty list). */
export async function myPlans(): Promise<DG.FileInfo[]> {
  // (no folder yet lists as nothing; files.exists says false of any path ending in a slash)
  return grok.dapi.files.list(plansFolder(), false, PLAN_EXT);
}
