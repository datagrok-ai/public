/* An arm's effect from the user's own earlier study of the same drug (the preclinical bundle's
   analysis/out_pk/plan_pk/ESPEC.md, v1): the effect that study measured at its dose, scaled to the new dose with the
   same schedule, route, tumor model and reading day, as a point and a range. Plan sizes the arm on the low end. Where
   the scaling was not validated Plan says why and leaves the effect to the user: what Plan can tell from the names,
   the model and the days it checks here, without calling Crux; the rest Crux's function checks. */
import {ArmPart, controlByName, parseArmName} from './arm-name';
import {crux, rows} from './crux';
import {finite, isFiniteNumber as num, median, SEED} from './engine';
import type {DrugPk, SavedExposure} from './exposure';
import {StudyMeta, armParts} from './study-meta';

/** An arm linked to an earlier study: the study and its group, the doses, and what Crux made of them (percent TGI),
 *  or why Plan stays silent. */
export interface FromStudy {
  study: string;
  arm: string;
  dosePrev: number;
  doseNew: number;
  /** the doses' unit, as both names write it */
  unit: string;
  /** null with `why`, or while it is computed */
  point: number | null;
  /** the 10th and 90th percentiles */
  range: [number, number] | null;
  why?: string;
  /** percent: the effect the earlier study measured, as Crux read it */
  anchorTgi?: number;
  /** "AUC at 10 mg/kg: 1,230 h*ng/mL", when NCA gives it */
  exposure?: string;
}

/** A link as a file or a table cell may give it (null: none, or damaged). */
export function fromStudyOf(x: unknown): FromStudy | null {
  const f = x as Partial<FromStudy> | null;
  if (!f || typeof f !== 'object' || typeof f.study !== 'string' || typeof f.arm !== 'string' ||
    !num(f.dosePrev) || !num(f.doseNew))
    return null;
  const range = Array.isArray(f.range) && f.range.length === 2 && f.range.every(num) ?
    [f.range[0], f.range[1]] as [number, number] : null;
  return {study: f.study, arm: f.arm, dosePrev: f.dosePrev!, doseNew: f.doseNew!,
    unit: typeof f.unit === 'string' ? f.unit : '',
    point: num(f.point) ? f.point! : null, range, ...typeof f.why === 'string' && f.why ? {why: f.why} : {},
    ...num(f.anchorTgi) ? {anchorTgi: f.anchorTgi} : {},
    ...typeof f.exposure === 'string' && f.exposure ? {exposure: f.exposure} : {}};
}

/** An earlier study as Plan reads it: its name, tumor model, group details and measurements. */
export interface PreviousStudy {
  name: string;
  line: string | null;
  meta: StudyMeta | null;
  rows: {animal: string, group: string | null, day: number | null, volume: number | null}[];
  /** the exposure Analyze kept with it (the patients' and NCA's per group; older studies: none) */
  exposure?: SavedExposure;
}

/** A group of an earlier study that received the arm's drug alone, with its control. */
export interface Anchor {
  study: PreviousStudy;
  arm: string;
  control: string;
  /** what the group received, read from its details or else its name */
  part: ArmPart;
  /** the last day with at least 3 animals measured in the group and in the control (null: none) */
  readDay: number | null;
  /** every such day, in order: the days the effect can be read on (the plan's own day must be one) */
  days: number[];
  /** percent: the analysis's TGI when the study kept it (Analyze, as Compare shows it), else 1 - median treated /
   *  median control on readDay */
  tgi: number | null;
  /** whether `tgi` is the analysis's, read on `tgiDay` */
  analyzed: boolean;
  tgiDay: number | null;
}

/** Fewest animals in the group and in its control on the day the effect is read (as Crux's function). */
const MIN_ANIMALS = 3;

/** The control of a study's groups: the one its details say receives a vehicle (a compound that names no drug), else
 *  the one named so, as Analyze finds it; undefined: none (Plan then stays silent, it does not guess). */
export function controlOf(groups: string[], meta: StudyMeta | null): string | undefined {
  const declared = meta?.groups.find((g) => groups.includes(g.group) && g.compound.trim() &&
    parseArmName(g.compound).length === 0)?.group;
  return declared ?? controlByName(groups);
}

/** The groups of a study's rows, in their order. */
const groupsIn = (study: PreviousStudy) =>
  [...new Set(study.rows.map((r) => r.group).filter((g): g is string => !!g))];

/** The day of the first dose: day 0, the app's convention (study-meta.ts: day 0 is randomization and first dose). */
export const FIRST_DOSE_DAY = 0;

/** The day the effect is read: the last from the first dose with enough distinct animals measured in both groups (as
 *  Crux's function reads it). */
export function readDayOf(study: PreviousStudy, arm: string, control: string): number | null {
  const days = readDaysOf(study, arm, control);
  return days.length ? days[days.length - 1] : null;
}

/** The days from the first dose with enough distinct animals measured in both groups, in order: the effect can be
 *  read on any of them (Crux reads the last of the rows it is given, so Plan gives it the rows up to the day). */
export function readDaysOf(study: PreviousStudy, arm: string, control: string): number[] {
  const animals = new Map<number, {t: Set<string>, c: Set<string>}>();
  for (const r of study.rows) {
    if (r.day === null || r.day < FIRST_DOSE_DAY || finite(r.volume) === null ||
      (r.group !== arm && r.group !== control))
      continue;
    const n = animals.get(r.day) ?? {t: new Set<string>(), c: new Set<string>()};
    (r.group === arm ? n.t : n.c).add(r.animal);
    animals.set(r.day, n);
  }
  return [...animals].filter(([, n]) => n.t.size >= MIN_ANIMALS && n.c.size >= MIN_ANIMALS).map(([d]) => d)
    .sort((a, b) => a - b);
}

/** Of a study's groups that can anchor `armName`, the one whose dose is nearest the arm's on a log scale (the dose
 *  ratio is what Crux scales by, within a quarter to four times), whatever their order in the file; null: none. */
export function nearestAnchor(anchors: Anchor[], armName: string): Anchor | null {
  const dose = parseArmName(armName)[0]?.dose?.value;
  const far = (a: Anchor) => dose && dose > 0 ? Math.abs(Math.log(a.part.dose!.value / dose)) : 0;
  return anchors.reduce<Anchor | null>((best, a) => best === null || far(a) < far(best) ? a : best, null);
}

/** The groups of `study` that received the drug of `armName` alone, with a dose (none: the arm is no single drug). */
export function anchorsIn(study: PreviousStudy, armName: string): Anchor[] {
  const mine = parseArmName(armName);
  if (mine.length !== 1)
    return [];
  const groups = groupsIn(study);
  const control = controlOf(groups, study.meta);
  if (control === undefined)
    return [];
  return groups.filter((g) => g !== control).flatMap((arm): Anchor[] => {
    // what the group received: its details when they name a compound (empty cells from its name), else its name
    const parts = armParts(arm, study.meta);
    if (parts.length !== 1 || parts[0].key !== mine[0].key || !parts[0].dose)
      return [];
    const days = readDaysOf(study, arm, control);
    const readDay = days.length ? days[days.length - 1] : null;
    const at = (g: string) => study.rows.filter((r) => r.group === g && r.day === readDay)
      .map((r) => finite(r.volume)).filter((v): v is number => v !== null);
    // the analysis's own TGI first, so Plan says what Analyze said of that study (the medians' ratio is another number)
    const kept = study.exposure?.effect;
    const own = kept?.tgi?.[arm];
    if (num(own))
      return [{study, arm, control, part: parts[0], readDay, days, tgi: own, analyzed: true,
        tgiDay: kept!.day ?? null}];
    const tgi = readDay === null ? null : 100 * (1 - median(at(arm)) / median(at(control)));
    return [{study, arm, control, part: parts[0], readDay, days, tgi: Number.isFinite(tgi) ? tgi : null,
      analyzed: false, tgiDay: readDay}];
  });
}

/** A group of `study` that received the drug of `armName` alone but whose details and name say no dose (null: none):
 *  Plan cannot scale from it until the study says the dose, so it says so rather than drop the study. */
export function undosedIn(study: PreviousStudy, armName: string): string | null {
  const mine = parseArmName(armName);
  const groups = groupsIn(study);
  const control = controlOf(groups, study.meta);
  if (mine.length !== 1 || control === undefined)
    return null;
  return groups.find((g) => {
    const p = g === control ? [] : armParts(g, study.meta);
    return p.length === 1 && p[0].key === mine[0].key && !p[0].dose;
  }) ?? null;
}

/** The note under an arm whose drug's earlier study says no dose. */
export const undosedText = (drug: string, study: string, group: string) =>
  `You have a study of ${drug} ("${study}", group ${group}) that does not say its dose: add it in that study's ` +
  'details to estimate this arm from it.';

/** Why Plan stays silent before asking Crux, in plain words (null: nothing it can see). The rest (the dose ratio, the
 *  doses with PK, the schedule and route, the anchor's effect and animals) is the function's. */
export function silence(armName: string, anchor: Anchor, planLine: string | null,
  evaluationDay: number | null, hostOf: (line: string | null) => string | null = () => null): string | null {
  const parts = parseArmName(armName);
  if (parts.length !== 1)
    return 'a combination: your study measured the drug alone';
  const dose = parts[0].dose;
  if (!dose)
    return 'write the new dose in the arm\'s name, e.g. "20 mg/kg"';
  if (dose.unit.toLowerCase() !== anchor.part.dose!.unit.toLowerCase())
    return `the dose is not in the unit of your study (${anchor.part.dose!.unit})`;
  // the study's days counted from its day 0, the first dose (FIRST_DOSE_DAY)
  if (!anchor.study.meta?.start)
    return 'your study does not say its day 0 (randomization and first dose)';
  if (anchor.study.rows.some((r) => r.day !== null && r.day < FIRST_DOSE_DAY))
    return 'your study has days before day 0; Plan reads day 0 as the first dose';
  if (!anchor.study.line)
    return 'your study does not say its tumor model';
  if (hostOf(planLine) !== hostOf(anchor.study.line))
    return 'your study\'s mice are not the kind of host this plan uses';
  if (planLine !== anchor.study.line)
    return `your study was run in ${anchor.study.line}, not ${planLine ?? 'this model'}`;
  if (anchor.readDay === null)
    return `fewer than ${MIN_ANIMALS} animals measured in its groups on any one day`;
  if (evaluationDay === null)
    return 'the day this plan reads its effect is not known yet';
  // the same day for the earlier study and the plan, as in every pair the scaling was validated on: the plan's day
  // must be one the study measured (any of them, not only its last: Plan gives Crux the rows up to that day)
  const day = Math.round(evaluationDay);
  if (anchor.days.includes(day))
    return null;
  const before = anchor.days.filter((d) => d < day).pop();
  const after = anchor.days.find((d) => d > day);
  if (after === undefined)
    return `your study's last day with ${MIN_ANIMALS} animals in each group is day ${before}, this plan reads ` +
      `its effect on day ${day}`;
  return `your study did not measure day ${day}, the day this plan reads its effect` +
    (before === undefined ? ` (its first day is ${after})` : ` (it measured days ${before} and ${after})`);
}

/** Whether the new arm keeps the anchor's schedule or route: a name that says none keeps them. */
const same = (mine: string | null, theirs: string | null) =>
  mine === null || mine.replace(/\s+/g, '').toLowerCase() === (theirs ?? '').replace(/\s+/g, '').toLowerCase();

/** The earlier study's dose as the texts say it: "10 mg/kg". */
export const doseText = (x: Pick<FromStudy, 'dosePrev' | 'unit'>) => `${x.dosePrev} ${x.unit}`.trim();

/** What Plan assumed of a name that says no schedule or route, for the effect's tooltip: "Assumes the same schedule
 *  and route as your 10 mg/kg study." (empty: the name says both) */
export function assumedText(armName: string, x: Pick<FromStudy, 'dosePrev' | 'unit'>): string {
  const p = parseArmName(armName)[0];
  const missing = [!p?.schedule && 'schedule', !p?.route && 'route'].filter((w) => w);
  return missing.length ? `Assumes the same ${missing.join(' and ')} as your ${doseText(x)} study.` : '';
}

/** The doses' unit as the names and NCA write it, compared. */
export const unitKey = (u: string) => u.replace(/[µμ]/g, 'u').toLowerCase();

/** Said with every range: Crux scales it with its 2 structures, without the drug's exposure over time (the curve from
 *  NCA failed its sealed check, plan_pk/v1/SELLO_NCA_1CMT.md). */
export const SIMPLIFIED = 'The range comes from a simplified model and may be somewhat narrow.';

/** What to do when Plan has no PK of the drug to scale within: Plan reads the concentration tables open in the
 *  workspace, and a saved NCA analysis opened in the NCA app is one (NCA offers no function to find saved ones). */
export const noPkText = (drug: string) =>
  `open your NCA analysis of ${drug} (or its plasma concentrations by dose) in Datagrok to estimate this arm from it`;

/** The one call to Crux: stats.tumor_effect_from_study (plan_pk/v1/CRUX_V1.md) on the anchor and its control, within
 *  the doses of the drug's PK (`pk`, as NCA read it: exposure.ts drugPkInWorkspace). No PK in that unit: silent,
 *  saying what to add, without calling Crux (the doses validated are those with PK). Proportionality: NCA's AUC per
 *  dose (exposure.ts proportionalOf). Crux's module without the function yet, or any failure, is a silence. */
export async function effectFromStudy(anchor: Anchor, armName: string, pk: DrugPk | null,
  day: number | null = anchor.readDay): Promise<Pick<FromStudy, 'point' | 'range' | 'why' | 'anchorTgi' | 'exposure'>> {
  const mine = parseArmName(armName)[0];
  if (!pk || unitKey(pk.unit) !== unitKey(anchor.part.dose!.unit))
    return {point: null, range: null, why: noPkText(anchor.part.drug)};
  // the exposure at the earlier study's dose, for the tooltip, when NCA gives it
  const auc = pk.auc.get(anchor.part.dose!.value);
  const exposure = auc ? {exposure: `AUC at ${anchor.part.dose!.value} ${anchor.part.dose!.unit}: ${auc}`} : {};
  // stats.tgi's long table; the function reads the anchor and its control only; an animal is its group's (xlsx
  // studies number their animals 1..n in every group); up to the plan's day, the last Crux reads (silence checked
  // the study measured it)
  const r = anchor.study.rows.filter((x) => (x.group === anchor.arm || x.group === anchor.control) &&
    x.day !== null && (day === null || x.day <= day) && finite(x.volume) !== null);
  try {
    const out = await crux('stats.tumor_effect_from_study', {animal: r.map((x) => `${x.group}|${x.animal}`),
      group: r.map((x) => x.group!), day: Float64Array.from(r, (x) => x.day!),
      volume: Float64Array.from(r, (x) => x.volume!)}, {anchor: anchor.arm, control: anchor.control,
      first_dose_day: FIRST_DOSE_DAY, dose_prev: anchor.part.dose!.value, dose_new: mine.dose!.value,
      pk_doses: pk.doses.join(','), pk_proportional: pk.proportional,
      same_schedule: same(mine.schedule, anchor.part.schedule), same_route: same(mine.route, anchor.part.route),
      seed: SEED});
    const row = rows(out.main)[0] ?? {};
    // (percent, as Crux gives them)
    const [point, p10, p90] = [finite(row.tgi_point), finite(row.tgi_p10), finite(row.tgi_p90)];
    const anchorTgi = finite(row.tgi_anchor) ?? undefined;
    if (row.reason || point === null || p10 === null || p90 === null)
      return {point: null, range: null, why: plainReason(String(row.reason ?? '')), anchorTgi, ...exposure};
    // the day Crux read must be the one Plan checked against this plan's (silence)
    if (day !== null && finite(row.read_day) !== day)
      return {point: null, range: null, why: `your study's effect could not be read on day ${day}`};
    return {point, range: [p10, p90], anchorTgi, ...exposure};
  } catch (e) {
    if (/unknown function/.test(String((e as Error)?.message ?? e)))
      return {point: null, range: null, why: 'this version of PMX Studio cannot scale from an earlier study yet'};
    console.error('PMX Studio: estimating from an earlier study failed', e);
    return {point: null, range: null, why: 'the estimate failed; try again, or set the effect by hand'};
  }
}

/** Crux's reasons (CRUX_V1.md), in the user's words. */
const REASONS: {[code: string]: string} = {
  schedule_differs: 'the schedule is not your study\'s',
  route_differs: 'the route is not your study\'s',
  ratio_out_of_range: 'the new dose is more than 4 times, or less than a quarter of, your study\'s',
  pk_not_proportional: 'exposure does not rise in step with the dose',
  dose_outside_pk_window: 'the new dose is outside the doses your PK study measured (from half the lowest to twice ' +
    'the highest)',
  too_few_animals: `fewer than ${MIN_ANIMALS} animals measured in its groups on any one day`,
  anchor_effect_out_of_range: 'the effect in your study was too small or too large to scale',
  control_fit_failed: 'the growth of your study\'s controls could not be described',
};

/** Crux's reason code in the user's words; one it does not know, the general one. */
export const plainReason = (code: string) => REASONS[code.trim()] ?? 'no basis to scale from your study';

const round = Math.round;

/** A range and its point (percent TGI) as `metric` says them, [low, high, point]: in T/C, 100 - TGI, the range turned
 *  round. */
export const rangeAs = ([lo, hi]: [number, number], point: number, metric: 'TGI' | 'TC'): [number, number, number] =>
  metric === 'TC' ? [100 - hi, 100 - lo, 100 - point] : [lo, hi, point];

/** What the effect field shows for a linked arm: "≈ 65%, between 55 and 80%, from your 10 mg/kg study", in T/C when
 *  the study's effects are (T/C = 100 - TGI, the range turned round); or why Plan is silent. */
export function effectText(x: FromStudy, metric: 'TGI' | 'TC'): string {
  if (!x.why && (!x.range || x.point === null))
    return `Estimating from your ${doseText(x)} study…`;
  if (!x.range || x.point === null)
    return `No basis to scale from your ${doseText(x)} study: ${x.why}`;
  const [lo, hi, point] = rangeAs(x.range, x.point, metric);
  return `≈ ${round(point)}%, between ${round(lo)} and ${round(hi)}%, from your ${doseText(x)} study`;
}

/** The effect Plan sizes a linked arm on: the low end of its range (percent, in the study's metric, within the
 *  field's 5-95); null when the arm is not linked or Plan is silent. */
export function sizingEffect(x: FromStudy | null | undefined, metric: 'TGI' | 'TC'): number | null {
  if (!x?.range)
    return null;
  // ponytail: a low end under 5% is sized as 5%, the field's floor; such an arm is barely detectable anyway
  const tgi = Math.min(95, Math.max(5, x.range[0]));
  return metric === 'TC' ? 100 - tgi : tgi;
}

/** The assumed TGIs a plan with linked arms adds to its scenarios: each one's low end and point. */
export const linkedScenarios = (arms: {fromStudy?: FromStudy | null}[]) =>
  arms.flatMap((a) => a.fromStudy?.range && a.fromStudy.point !== null ?
    [a.fromStudy.range[0], a.fromStudy.point].map(round).filter((t) => t > 0 && t < 100) : []);

/** The note under an arm: "You have a study of LD4172 at 10 mg/kg (TGI 45%, day 21)." TGI only when it is the
 *  analysis's; a study never analyzed says its medians' ratio, which is not Analyze's TGI. */
export function offerText(a: Anchor, planLine: string | null = null): string {
  const effect = a.tgi === null ? '' : a.analyzed ? `TGI ${round(a.tgi)}%` :
    `median volume ${round(100 - a.tgi)}% of the control's`;
  const said = [effect, a.tgiDay === null ? '' : `day ${a.tgiDay}`].filter((s) => s).join(', ');
  return `You have a study of ${a.part.drug} at ${a.part.dose!.value} ${a.part.dose!.unit}` +
    // the study's model when it is not the plan's: an effect seen in another model is not this one's (nor, unsaid,
    // one whose model the study never recorded)
    `${a.study.line && a.study.line !== planLine ? ` in ${a.study.line}` : !a.study.line && planLine ?
      ' in a tumor model the study does not record' : ''}${said ? ` (${said})` : ''}.`;
}
