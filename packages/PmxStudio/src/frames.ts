/* The engine's results as Datagrok tables: the glue between Crux's study design (engine.ts) and the package
 * functions. */
import * as DG from 'datagrok-api/dg';

import {rows} from './crux';
import {enginePriors, mapPrior} from './priors';
import {HUMANE, LAB_WEEKDAYS, StudyDeclaration, correctionOf, designScenarios, designStudy,
  isFiniteNumber as ok} from './engine';
import {armsOf} from './arms-input';

/** pmxPlanStudy's startWeekday choices, Monday first (the engine's 0 … 6). */
export const WEEKDAY_NAMES = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'];

/** The measurement weekdays of a text such as "Mon, Wed, Fri" (any case, full names too): their indices, Monday
 *  first; none recognised: Mon / Wed / Fri. */
export function weekdaysOf(text: string | undefined): number[] {
  const words = String(text ?? '').toLowerCase().split(/[,;\s]+/).filter((w) => w.length >= 3);
  const days = WEEKDAY_NAMES.map((_, i) => i)
    .filter((i) => words.some((w) => WEEKDAY_NAMES[i].toLowerCase().startsWith(w.slice(0, 3))));
  return days.length > 0 ? days : [...LAB_WEEKDAYS];
}

/** The weekdays as the measurementDays parameter writes them: "Mon, Wed, Fri". */
export const weekdaysText = (days: number[]) => days.map((i) => WEEKDAY_NAMES[i].slice(0, 3)).join(', ');

/** Most assumed TGIs in one run: each is a whole design. */
const MAX_SCENARIOS = 8;

/** The assumed TGIs of a text such as "30, 50, 70", in order; empty: none (no scenarios are computed). */
export function tgisOf(text: string | null | undefined): number[] {
  const tgis = [...new Set(String(text ?? '').split(/[,;\s%]+/).filter((w) => w !== '').map(Number))]
    .sort((a, b) => a - b);
  if (tgis.some((x) => !(x > 0 && x < 100)))
    throw new Error(`Assumed TGIs: write percentages between 0 and 100, such as 30, 50, 70 (got "${text}")`);
  if (tgis.length > MAX_SCENARIOS)
    throw new Error(`Assumed TGIs: at most ${MAX_SCENARIOS} at a time, each is a whole design`);
  return tgis;
}

/** pmxPlanStudy's inputs as the form and Grokky give them: percentages, choice labels, empty optionals. */
export interface PlanInputs {
  line: string; endpoint: number; arms: DG.DataFrame | null; power: number; correction: string;
  measurementDays?: string; n?: number; lastDay?: number; alpha?: number; preLoss?: number; postLoss?: number;
  lossRule?: string;
  confidence?: number; vRand?: number; startWeekday?: string; humane?: number; compareSchedules?: boolean;
  detectableEffect?: boolean;
}

/** The declaration the engine reads: each empty optional at its default, percentages as shares, counts rounded and
 *  kept in range. */
export function declarationOf(x: PlanInputs): StudyDeclaration {
  const share = (v: number | undefined, d: number) => (ok(v) ? v : d) / 100;
  return {
    line: x.line, endpoint: x.endpoint, arms: armsOf(x.arms), power: x.power / 100,
    correction: correctionOf(x.correction),
    weekdays: weekdaysOf(x.measurementDays),
    startWeekday: Math.max(0, WEEKDAY_NAMES.indexOf(x.startWeekday ?? 'Monday')),
    n: ok(x.n) && x.n >= 2 ? Math.round(x.n) : null,
    lastDay: ok(x.lastDay) && x.lastDay >= 1 ? Math.round(x.lastDay) : null,
    alpha: ok(x.alpha) ? x.alpha : 0.05, vRand: ok(x.vRand) ? x.vRand : 100,
    humane: ok(x.humane) && x.humane > 0 ? x.humane : HUMANE,
    preLoss: share(x.preLoss, 20), postLoss: share(x.postLoss, 10), confidence: share(x.confidence, 90),
    lossRule: x.lossRule === 'Average loss only' ? 'expected' : 'confident',
    schedules: x.compareSchedules !== false, mde: x.detectableEffect !== false,
  };
}

export async function priorsFrame(): Promise<DG.DataFrame> {
  const p = await enginePriors();
  const df = DG.DataFrame.fromObjects(Object.entries(p.lines).map(([line, l]) => ({
    line,
    alpha_per_day: Math.exp(l.log_alpha.map_mean),
    alpha_q025: Math.exp(l.log_alpha.map_q025),
    alpha_q975: Math.exp(l.log_alpha.map_q975),
    beta_per_day: Math.exp(l.log_beta),
    tau_log_alpha: l.log_alpha.tau,
    ess: l.log_alpha.ess,
    nonprogressing: l.take.nonprogressing_map,
    papers: l.papers,
    animals: l.animals,
  })))!;
  df.name = `PMX Studio priors (${p.version})`;
  return df;
}

export async function mapPriorFrame(line: string): Promise<DG.DataFrame> {
  // (Crux gives each column decoded, not its name: read as rows, as every other result is)
  const df = DG.DataFrame.fromObjects(rows((await mapPrior(line)).main))!;
  df.name = `MAP prior log alpha ${line}`;
  return df;
}

/** The control in a neutral slate, each treated arm in the platform's categorical palette, in declaration order:
 *  set as the arm columns' colour coding, so every viewer of the tables shows an arm in one colour. */
export function armColors(arms: string[]): {[arm: string]: string} {
  return Object.fromEntries([['Control', '#5b6b7f'],
    ...arms.map((a, k) => [a, DG.Color.toHtml(DG.Color.getCategoricalColor(k))])]);
}

/** A study's groups in the colours every stage gives them: the control (named `control`) in the control's grey, the
 *  treated groups in order, as Plan colours its arms, whether or not `groups` lists the control among them. */
export function groupColors(groups: string[], control: string): {[group: string]: string} {
  const colors = armColors(groups.filter((g) => g !== control));
  colors[control] = colors['Control'];
  return colors;
}

/** The study (one row), its arms (one row each), the power curves (per arm and n), the predicted volumes (per
 *  day and arm: the control, then each treated arm), the power by last day (per arm and day) and by measurement
 *  schedule (per arm and measurements a week), the study calendar (per day from randomization) and the assumed-TGI
 *  scenarios (per TGI and measurements a week). */
export interface PlanFrames {plan: DG.DataFrame, armResults: DG.DataFrame, powerCurve: DG.DataFrame,
  curves: DG.DataFrame,
  duration: DG.DataFrame, schedules: DG.DataFrame, calendar: DG.DataFrame, scenarios: DG.DataFrame}

export async function planFrames(d: StudyDeclaration, tgis: number[] = []): Promise<PlanFrames> {
  const [s, sweep] = await Promise.all([designStudy(d), designScenarios(d, tgis)]);
  // (none asked for: an empty table, which fromObjects does not build)
  const scenarios = DG.DataFrame.fromObjects(sweep.flatMap(({tgi, design: x}) => x.schedules.map((c) => ({tgi,
    per_week: c.perWeek, n_total: c.nTotal, implant_total: x.implant, randomise_total: x.randomiseTotal,
    study_end: x.studyEnd})))) ?? DG.DataFrame.create();
  const me = s.measurement;
  const tag = `${d.line} · ${d.arms.length} arm${d.arms.length > 1 ? 's' : ''}`;
  // the first day some controls (5%) may be over the humane limit, before the end
  const watch = s.duration.find((r) => r.pastHumane >= 0.05 && r.day < s.studyEnd)?.day ?? null;
  const plan = DG.DataFrame.fromObjects([{
    line: d.line, atlas_papers: s.atlas.papers, atlas_animals: s.atlas.animals, prior_ess: s.atlas.ess,
    study_end: s.studyEnd, endpoint_day: s.endpointDay, last_day_chosen: d.lastDay !== null,
    duration_q10: s.durationQ10, duration_q90: s.durationQ90,
    p_endpoint_not_reached: s.pNotReached, endpoint: d.endpoint, v_rand: d.vRand,
    n_control: s.nControl, randomise_control: s.randomiseControl, n_total: s.nTotal, randomise_total: s.randomiseTotal,
    implant_total: s.implant, pre_loss: d.preLoss, post_loss: s.postLoss, confidence: d.confidence,
    confidence_study: s.confidenceStudy,
    loss_rule: d.lossRule,
    n_fixed: d.n !== null,
    nonprogressing: s.nonprogressing, power_target: d.power, alpha: d.alpha, correction: d.correction,
    treated_arms: d.arms.length, alpha_per_comparison: s.alphaArm,
    measurement_days: s.calendar.join(', '), measurements: s.calendar.length, per_week: d.weekdays.length,
    measurement_weekdays: weekdaysText(d.weekdays),
    start_weekday: d.startWeekday,
    humane_watch_day: watch, humane_limit: d.humane,
    compared: 'per-animal slope of log V over the calendar (Welch t test)',
    measurement_error: `${me.model}: c0 ${me.c0.toFixed(3)}, c1 ${me.c1.toFixed(2)} (atlas, ${me.triples} triples)`,
    notes_groups: s.notes.groups.join(' | '), notes_growth: s.notes.growth.join(' | '),
  }])!;
  const arms = DG.DataFrame.fromObjects(s.arms.map((a) => ({
    arm: a.declaration.name, therapy: a.archetype.label, metric: a.declaration.metric,
    effect: a.declaration.effect, effect_realised: a.effect, kill_rate_per_day: a.killRate,
    // (not Crux's power at its "median study": that lab is median in time, not in effect, so its power is not a
    // typical lab's; decision tumor-plan-checks, Revisit)
    n_needed: a.nAssurance, n: a.n, randomise: a.randomise, assurance: a.assurance, mde: a.mde,
    // what each fit of the line's spreads needs alone, the prior's own first (several only on a line the atlas's
    // controls do not pin down, where n_needed is for the power averaged over them)
    n_fits: a.fits.map((x) => x ?? '>40').join(', '),
    responders: a.archetype.responders, onset_day: a.archetype.onsetDays,
    translation: a.archetype.translation, archetype_basis: a.archetype.basis,
  })))!;
  const powerCurve = DG.DataFrame.fromObjects(s.arms.flatMap((a) => a.power.map((r) =>
    ({arm: a.declaration.name, n: r.n, implant: r.implant, assurance: r.assurance}))))!;
  const curves = DG.DataFrame.fromObjects(s.curves)!;
  // the controls' shares repeated on each arm's rows: one tidy table, coloured by arm like the others (no arm reaches
  // the target power with any n: an empty table, which fromObjects does not build)
  const duration = DG.DataFrame.fromObjects(s.duration.flatMap((r) => Object.entries(r.assurance).map(([arm, p]) =>
    ({arm, day: r.day, assurance: p, past_endpoint: r.pastEndpoint, past_humane: r.pastHumane})))) ??
    DG.DataFrame.create();
  // (not asked for: an empty table, which fromObjects does not build)
  const schedules = DG.DataFrame.fromObjects(s.schedules.flatMap((r) => Object.entries(r.assurance).map(([arm, p]) =>
    ({arm, per_week: r.perWeek, measurements: r.measurements, assurance: p, n_total: r.nTotal})))) ??
    DG.DataFrame.create();
  // the study day by day from randomization (its date is known only when the tumors reach the window): what is done,
  // and from when some controls may pass the humane limit
  const events = new Map<number, string[]>();
  const add = (day: number, what: string) => events.set(day, [...(events.get(day) ?? []), what]);
  for (const day of s.calendar) {
    add(day, day === 0 ? 'Randomization and first dose; tumor volume and body weight (baseline)' :
      'Tumor volume and body weight');
  }
  if (watch !== null)
    add(watch, `From today some controls may pass the humane limit (${d.humane} mm³): check them more often`);
  const [q10, q90] = [Math.round(s.durationQ10), s.durationQ90 === null ? null : Math.round(s.durationQ90)];
  const range = `8 in 10 studies: ${q90 === null ? `from day ${q10}` : `day ${q10} to ${q90}`}`;
  add(s.studyEnd, d.lastDay !== null ? 'End of study, as chosen' : s.studyEnd === s.endpointDay ?
    `End of study, when the median control reaches ${d.endpoint} mm³ (${range})` :
    `End of study: the first measurement day after the median control reaches ${d.endpoint} mm³ on day ` +
    `${s.endpointDay} (${range})`);
  const calendar = DG.DataFrame.fromObjects([...events].sort(([a], [b]) => a - b)
    .map(([day, what]) => ({day, activity: what.join('. ')})))!;
  calendar.col('day')!.meta.units = 'd';
  const colors = armColors(d.arms.map((a) => a.name));
  curves.col('arm')!.meta.colors.setCategorical(colors);
  powerCurve.col('arm')!.meta.colors.setCategorical(colors);
  duration.col('arm')?.meta.colors.setCategorical(colors);
  schedules.col('arm')?.meta.colors.setCategorical(colors);
  if (duration.col('day'))
    duration.col('day')!.meta.units = 'd';
  // units and the reference lines travel with the table: any viewer of it shows them
  curves.col('day')!.meta.units = 'd';
  for (const q of ['q05', 'q25', 'median', 'q75', 'q95']) {
    curves.col(q)!.meta.units = 'mm³';
    curves.col(q)!.meta.format = '0.0';
  }
  curves.meta.formulaLines.addLine({title: `Endpoint ${d.endpoint} mm³`, formula: `\${median} = ${d.endpoint}`});
  curves.meta.formulaLines.addLine({title: `End of study, day ${s.studyEnd}`, formula: `\${day} = ${s.studyEnd}`,
    style: 'solid', color: '#4d5261'});
  for (const [df, name] of [[plan, 'Plan'], [arms, 'Arms'], [powerCurve, 'Power'], [curves, 'Volumes'],
    [duration, 'Power by end of study'], [schedules, 'Power by measurement schedule'],
    [calendar, 'Study calendar'], [scenarios, 'Assumed TGI scenarios']] as const)
    df.name = `${name} ${tag}`;
  return {plan, armResults: arms, powerCurve, curves, duration, schedules, calendar, scenarios};
}
