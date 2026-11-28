/* The study against its plan, for the Analyze view: the Plan field of its form (a .pmxplan file from Datagrok's
   files, the computer or a drop, or one of the user's saved plans) and the "Plan vs observed" tab: a sentence that
   sums it up, each treated group's effect against the assumed one, the controls against their predicted growth, the
   end against the predicted day, and "Contribute to the atlas" (coming soon). The comparison itself is
   pmxAnalyzeStudy's (analyze.ts); this file only shows it. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {badge, computed, div, divH, iconButton, span, Tooltip} from '@datagrok-libraries/u2';
import type {BadgeVariant, ReadonlySignal, Signal} from '@datagrok-libraries/u2';
import {FileInput, FuncCallForm} from '@datagrok-libraries/u2/dg';

import {AnalysisFrames, PREDICTED, TOO_FEW, armKey} from './analyze';
import {modelName} from './priors';
import {controlByName} from './arm-name';
import {weekdaysOf, weekdaysText} from './frames';
import {parseLooks} from './monitor';
import {LOOKS_TAG, lineOf, roleColumns, studyNameOf} from './study-file';
import {GroupDetails, StudyMeta, emptyMeta, metaOf, setMetaOn} from './study-meta';
import {BandChart} from './band-chart';
import {CORRECTION_NAME, correctionOf, finite} from './engine';
import {PLAN_TAG, SavedPlan, parsePlan, planKey, planTag, plansFolder, shownPath} from './plan-file';
import {domainStorageOn, savedPlans} from './store';
import {axisAt, block, guide, infoNote, metricWord, modal, noteLine, notes, notesBox, okNamed, pct,
  percentRange, plainError, pValue, showCentred, signed, sized, swatch, volumeAxis, whenSettled} from './ui';
import {esc} from './sheet';

/** A treated group's effect against the plan's on a shared percent axis from `lo` to 100: the range a study run as
 *  planned shows, shaded in the group's colour; the assumed effect a dark tick; the observed one a dot (no interval:
 *  one day's description, the result is the test on the growth rates). */
function planBar(observed: number, assumed: number, expected: [number, number],
  lo: number, color: string, metric: string): HTMLElement {
  const track = div([], 'pmx-plan-bar');
  Object.assign(track.style, {position: 'relative', height: '18px'});
  const band = div([]);
  Object.assign(band.style, {position: 'absolute', left: axisAt(expected[0], lo),
    width: `calc(${axisAt(expected[1], lo)} - ${axisAt(expected[0], lo)})`, top: '3px', height: '12px',
    borderRadius: '2px', background: `color-mix(in srgb, ${color} 22%, transparent)`});
  const point = div([]);
  Object.assign(point.style, {position: 'absolute', left: `calc(${axisAt(observed, lo)} - 5px)`, top: '4px',
    width: '10px', height: '10px', borderRadius: '50%', background: color});
  const tick = div([]);
  Object.assign(tick.style, {position: 'absolute', left: `calc(${axisAt(assumed, lo)} - 1px)`, top: '0',
    width: '2px', height: '18px', background: 'var(--dg-text-color, #1f2a37)'});
  track.append(guide(0, lo), guide(50, lo), guide(100, lo), band, point, tick);
  // one sentence for the bar, the assumption first: the range alone read as the effect's interval (round 3, W3)
  const words = `Plan assumed ${metric} ${signed(assumed)}; 95 of 100 labs running this plan would see ` +
    `${percentRange(expected[0], expected[1])}; observed ${signed(observed)}`;
  Tooltip.bind(track, words);
  track.setAttribute('role', 'img');
  track.setAttribute('aria-label', words);
  return track;
}

/** The percent axis over the plan bars, in the header row: each value centred on the row's line like the headers
 *  beside it, its tick below at the guide. */
function percentAxis(lo: number): HTMLElement {
  const axis = div([], 'pmx-plan-axis');
  Object.assign(axis.style, {position: 'relative', height: '26px', fontSize: '11px', fontWeight: '400',
    color: 'var(--dg-text-color-muted, #68727f)'});
  for (const [x, shift] of [[lo, '0'], [0, '-50%'], [50, '-50%'], [100, '-100%']] as [number, string][]) {
    const label = span(`${x}%`);
    Object.assign(label.style, {position: 'absolute', top: '50%', left: axisAt(x, lo),
      transform: `translate(${shift}, -50%)`, lineHeight: '14px', whiteSpace: 'nowrap'});
    const tick = div([]);
    Object.assign(tick.style, {position: 'absolute', bottom: '0', left: axisAt(x, lo), width: '1px', height: '5px',
      background: 'var(--dg-text-color-muted, #9aa3ad)'});
    axis.append(label, tick);
  }
  return axis;
}

/** A verdict as a badge: held in green, off in amber, better than assumed in the accent, inconclusive in grey. */
function verdictBadge(text: string | null): HTMLElement {
  if (!text)
    return span('—');
  const variant: BadgeVariant = /^As /.test(text) ? 'success' : /Stronger/.test(text) ? 'accent' :
    text === 'Inconclusive' ? 'default' : 'warning';
  const el = badge(text, {variant});
  el.classList.add('pmx-badge');
  if (text === 'Inconclusive')
    Tooltip.bind(el, 'In the planned range, but the study cannot tell this group from the control');
  return el;
}

/** The test on the slopes as a badge: slower than the control in green, faster in amber, no clear effect in grey. */
export function activityBadge(activity: string, p: unknown): HTMLElement {
  const variant: BadgeVariant = activity === 'Slows growth' ? 'success' : activity === 'Grows faster' ? 'warning' :
    'default';
  const el = badge(activity, {variant});
  el.classList.add('pmx-badge');
  Tooltip.bind(el, activity === TOO_FEW ? 'Fewer than two animals with two measurements or more, in the group or ' +
    'its control: no test, no p' : `Tumor growth against the control, p ${pValue(p)}`);
  return el;
}

/** Why a group with no clear effect says little about the compound: the power the study had for the planned effect
 *  with its animals (and the animals the plan needed), or, with no plan, the smallest effect it detects. `g` reads
 *  the group's row of versusPlan, `a` the analysis row. Null for a group the test tells from the control. */
export function noEffectNote(g: (c: string) => any, a: (c: string) => any): string | null {
  if (g('activity') !== 'No clear effect')
    return null;
  const [n, power, needed, target, mde] = [finite(g('n')), finite(g('power')), finite(g('n_needed')),
    finite(a('power_target')), finite(g('mde'))];
  // the correction named beside its p (round 4, analyst W3)
  const [corr, k] = [a('correction'), finite(a('tested_arms')) ?? 0];
  const adjusted = k >= 2 && corr !== 'none' ? `, ${CORRECTION_NAME[correctionOf(corr)]}-adjusted for ${k} ` +
    'groups' : '';
  const head = `${g('group')}: no clear effect (p ${pValue(g('p_value'))}${adjusted}).`;
  if (power !== null) {
    return `${head} With ${n} animals, a study like this one detects the planned effect in ${Math.round(power * 100)}` +
      ' of 100 labs' + (needed !== null && n !== null && needed > n ?
      `; the plan needed ${needed}${target === null ? '' : ` for ${Math.round(target * 100)} of 100`}` : '') +
      '. This does not show the compound is inactive.';
  }
  if (a('plan'))
    return `${head} The study cannot tell this group from the control.`;
  return `${head} With ${n} animals, this study detects ` + (mde === null ? 'no effect reliably' :
    `a TGI of ${Math.round(mde)}% or more in 80 of 100 labs`) + ': a smaller effect can go unseen.';
}

const PLACEHOLDER = 'Optional: a plan saved in Plan';

/** What a study already holds of what a plan picked for it brings (planOnStudy). */
export interface StudyHolds {
  /** its details (metaOf), its planned reviews, its tumor model (lineOf) */
  meta: StudyMeta | null; looks: string; line: string | null;
  /** its reviews are fixed (Monitor: saved with them, or measured on): the plan's are not asked about */
  looksFixed?: boolean;
  /** the plan its table carries (null: none) */
  plan: SavedPlan | null;
  /** the groups in its data */
  groups: string[];
}

/** The study with what a plan brings: what the study lacks filled from the plan (the planned reviews, Run by, day 0,
 *  each group's dose, route and schedule by its arm, matched by name in any case, and its colour; the plan carried
 *  as the table's, for its measurement days); where the study holds another value, its own, or with `use` the
 *  plan's, each such difference said ("Day 0: the plan says …, the study has …"). `missing`: the plan's arms with
 *  no group in the data; `extra`: the data's treated groups with no arm in the plan. Nothing is invented for
 *  either. */
export function planOnStudy(p: SavedPlan, s: StudyHolds, use = false) {
  const differences: string[] = [];
  const pick = (what: string, plan: string, study: string, same = plan === study, shown = (x: string) => x): string => {
    if (!plan || !study || same)
      return study || plan;
    differences.push(`${what}: the plan says ${shown(plan)}, the study has ${shown(study)}`);
    return use ? plan : study;
  };
  // a model by the name the user knows it by, never its internal id ("new:pdx:6cr3vk9bjm": QA 0.1.3 round 4 R4-2)
  const modelShown = (x: string) => x === p.predicted.line && p.history?.name ? p.history.name : modelName(x);
  const sameLooks = (a: string, b: string) => {
    try {
      return parseLooks(a).join() === parseLooks(b).join();
    } catch {
      return a === b;
    }
  };
  const planLooks = s.looksFixed ? s.looks.trim() : typeof p.inputs.looks === 'string' ? p.inputs.looks.trim() : '';
  const looks = pick('Planned reviews', planLooks, s.looks.trim(), sameLooks(planLooks, s.looks.trim()));
  const line = pick('Tumor model', p.predicted.line ?? '', s.line ?? '', undefined, modelShown) || null;
  // the measurement days are the plan's inputs: the plan goes on the table unless it carries another
  const other = !!s.plan && planKey(s.plan) !== planKey(p);
  if (other) {
    const days = (q: SavedPlan) => weekdaysText(weekdaysOf(q.inputs.measurementDays as string | undefined));
    pick('Measurement days', days(p), days(s.plan!));
  }
  let meta = s.meta;
  const pm = p.study;
  const byKey = <T>(x: {[k: string]: T} | undefined, g: string) =>
    Object.entries(x ?? {}).find(([k]) => armKey(k) === armKey(g))?.[1];
  if (pm) {
    const m = s.meta ?? emptyMeta();
    const runBy = pick('Run by', pm.runBy.trim(), m.runBy.trim());
    const start = pick('Day 0', pm.start ?? '', m.start ?? '') || null;
    let groups = m.groups;
    const colors = {...m.colors};
    const slots = {...m.slots ?? {}};
    const rowsText = (rows: GroupDetails[]) => rows.map((r) => [r.compound, r.dose === null ? '' :
      `${r.dose} ${r.unit}`, r.route, r.schedule].map((x) => x.trim()).filter((x) => x).join(' ')).filter((x) => x)
      .join(' + ');
    for (const g of s.groups) {
      const planRows = pm.groups.filter((r) => armKey(r.group) === armKey(g)).map((r) => ({...r, group: g}));
      const own = groups.filter((r) => armKey(r.group) === armKey(g));
      if (planRows.length && pick(g, rowsText(planRows), rowsText(own)) !== rowsText(own))
        groups = [...groups.filter((r) => !own.includes(r)), ...planRows];
      const color = byKey(pm.colors, g);
      if (color && !colors[g])
        colors[g] = color;
      const slot = byKey(pm.slots, g);
      if (slot !== undefined && !(g in slots) && !Object.values(slots).includes(slot))
        slots[g] = slot;
    }
    meta = {...m, runBy, start, groups, colors, ...Object.keys(slots).length ? {slots} : {}};
  }
  // every arm the plan names, its form's and the ones it predicted (a combination built from two: chaos round 3 #10,
  // 3 of 8 arms listed), each once
  const arms = [...(p.inputs.arms ?? []).map((a) => a.name), ...(p.predicted?.arms ?? []).map((a) => a.name)]
    .filter((a, i, all) => a && all.findIndex((b) => armKey(b) === armKey(a)) === i);
  const control = controlByName(s.groups);
  const planned = [...arms, ...(pm?.groups ?? []).map((r) => r.group)].map(armKey);
  const missing = arms.filter((a) => !s.groups.some((g) => armKey(g) === armKey(a)));
  const extra = s.groups.filter((g) => g !== control && !planned.includes(armKey(g)));
  // the arms themselves, and how many, when they do not match: the data is the study's whichever is taken, so they
  // alone ask nothing (said under the field), but a question asked lists them too (it listed the model alone, QA
  // 0.1.3 round 4 R4-2)
  const treated = s.groups.filter((g) => g !== control);
  if (differences.length && s.groups.length && (missing.length || extra.length)) {
    differences.push(`Treated groups: the plan has ${arms.length} (${arms.join(', ')}), the study has ` +
      `${treated.length} (${treated.join(', ') || 'none'})`);
  }
  return {meta, looks, line, differences,
    /** the plan goes on the table as its own: it carries none, or another and the plan's values are taken */
    carry: !s.plan || (use && other), missing, extra};
}

/** What a plan and its study do not share, under the Plan field: the plan's arms with no data, the data's groups the
 *  plan lacks, each with the plan's arm it nearly matches (a typo, a dose written otherwise); '': they match. */
export function mismatchWords(missing: string[], extra: string[]): string {
  const near = (g: string) => {
    const m = nearestName(g, missing);
    return m ? `${g} (close to “${m}”?)` : g;
  };
  return [missing.length ? `Plan arms with no data: ${missing.join(', ')}` : '',
    extra.length ? `Groups not in the plan: ${extra.map(near).join(', ')}` : ''].filter((x) => x).join('. ');
}

/** Of `names`, the one `name` nearly matches once case and punctuation are set aside (a few letters apart: at most
 *  2, or a fifth of its length); null: none. */
export function nearestName(name: string, names: string[]): string | null {
  const a = armKey(name);
  let best: string | null = null;
  let bestD = Math.max(2, Math.floor(a.length / 5)) + 1;
  for (const n of names) {
    const d = editDistance(a, armKey(n));
    if (d < bestD) {
      best = n;
      bestD = d;
    }
  }
  return best;
}

function editDistance(a: string, b: string): number {
  let prev = Array.from({length: b.length + 1}, (_, j) => j);
  for (let i = 1; i <= a.length; i++) {
    const row = [i];
    for (let j = 1; j <= b.length; j++)
      row[j] = Math.min(prev[j] + 1, row[j - 1] + 1, prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
    prev = row;
  }
  return prev[b.length];
}

/** The differences between a plan and its study, asked about once: true to take the plan's, false to keep the
 *  study's (also on Cancel). */
function askDifferences(differences: string[]): Promise<boolean> {
  return new Promise((resolve) => {
    const d = modal('plan-differences', 'The plan and the study differ');
    if (!d)
      return resolve(false);
    d.add(ui.divText('This plan and the study disagree on:'));
    d.add(ui.divV(differences.map((x) => ui.divText(`• ${x}`))));
    // (Cancel, Esc and × alike keep the study's; ponytail: a dialog closed by the app leaving the page settles
    // nothing, the view goes with it)
    d.onOK(() => resolve(true)).onCancel(() => resolve(false));
    okNamed(d, 'Use the plan\'s');
    const cancel = d.getButton('CANCEL');
    if (cancel)
      cancel.textContent = 'Keep the study\'s';
    showCentred(d);
  });
}

/** The Plan field, under the form's Study fields; a plan picked writes its pre-specified analysis into the form and
 *  its text into the call (`state`: the view's stale mark, error line and last analysis). The view owns it. */
export function planInput(form: FuncCallForm, call: DG.FuncCall, state: {stale: Signal<boolean>,
  error: Signal<string | null>, result: ReadonlySignal<{analysis: DG.DataFrame} | undefined>,
  /** the plan file being read: a run waits for it, or it would analyze without the plan just picked */
  reading: Promise<unknown>,
  /** the plan picked or taken off: the view runs again (the call's plan is no field of the form, so the form does
   *  not say it changed, and the results would stay dimmed for good) */
  changed: () => void}): FileInput {
  const {stale, error, result} = state;
  const setField = (name: string, value: unknown) => {
    const field = form.getInput(name);
    if (field)
      field.value = value;
    call.setParamValue(name, value);
  };
  // the Plan field also offers the user's saved plans; a plan picked writes its pre-specified analysis into the form
  // and what it says of the study onto the study's table (planOnStudy: what the study lacks; where it holds another
  // value, asked once). A change of the study: the view marks it unsaved. Taking the plan off leaves it. The note
  // (returned) names the arms and groups that do not match
  const applyPlan = async (p: SavedPlan): Promise<string> => {
    for (const [name, value] of [['endpoint', p.predicted.endpoint], ['correction', p.inputs.correction],
      ['alpha', p.inputs.alpha]] as [string, unknown][]) {
      if (value !== undefined && value !== null)
        setField(name, value);
    }
    const df = call.getParamValue('study');
    if (!(df instanceof DG.DataFrame)) {
      if (p.predicted.line)
        setField('line', p.predicted.line);
      return '';
    }
    const groupParam = call.getParamValue('group');
    let groupCol: DG.Column | null = groupParam instanceof DG.Column ? groupParam : null;
    try {
      groupCol ??= roleColumns(df).group;
    } catch {
      // no group column yet: no group matched
    }
    // (the reviews Monitor fixed stay as they are)
    const looksField = form.getInput('looks');
    const looks = looksField ? String(looksField.value ?? '') : df.getTag(LOOKS_TAG) ?? '';
    let carried: SavedPlan | null = null;
    try {
      carried = df.getTag(PLAN_TAG) ? parsePlan(df.getTag(PLAN_TAG)!) : null;
    } catch {
      // a damaged plan on the table: as none
    }
    const holds: StudyHolds = {meta: metaOf(df), looks, looksFixed: !!looksField?.disabled, line: lineOf(df),
      plan: carried, groups: groupCol ? groupCol.categories.filter((g) => g !== '') : []};
    let r = planOnStudy(p, holds);
    if (r.differences.length && await askDifferences(r.differences))
      r = planOnStudy(p, holds, true);
    // (the study keeps its own name: a plan carried names a study after it, studyNameOf; QA 0.1.3 round 4 R4-2)
    if (r.carry)
      df.setTag(PLAN_TAG, planTag({...p, name: studyNameOf(df) || p.name}));
    if (r.meta && JSON.stringify(r.meta) !== JSON.stringify(holds.meta))
      setMetaOn(df, r.meta);
    if (r.looks !== looks.trim()) {
      df.setTag(LOOKS_TAG, r.looks);
      if (looksField)
        setField('looks', r.looks);
    }
    if (r.line)
      setField('line', r.line);
    return mismatchWords(r.missing, r.extra);
  };
  // what the plan picked did not match, under the field
  const note = div([], 'pmx-plan-note');
  Object.assign(note.style, {fontSize: '12px', color: 'var(--dg-text-color-secondary, #6b7785)'});
  // the Plan field: u2's FileInput (a path in Datagrok's files, a file from the computer, or a drop) under the Study
  // fields; the file's text goes into the call's hidden `plan` param (ponytail:
  // placed by hand, u2 funcForm has no file field yet)
  const planField = form.build(() => new FileInput({label: 'Plan', nullable: true,
    placeholder: PLACEHOLDER, tooltipText: 'The plan saved in Plan (.pmxplan), to compare the ' +
      'study with. Empty: the plan the study table carries, if any'}));
  planField.addOptions(form.build(() => iconButton('list', () => void savedPlans().then((plans) => {
    if (plans.length === 0) {
      grok.shell.info(domainStorageOn() ? 'No saved plans yet: save one in Plan' :
        `No saved plans yet: save one in Plan (it goes to ${shownPath(plansFolder())})`);
      return;
    }
    DG.Menu.popup().items(plans.map((p) => p.name), (name: string) => void plans.find((p) => p.name === name)?.read()
      .then((f) => planField.value = f,
        (e) => grok.shell.error(plainError(e, `"${esc(name)}" could not be read.`)))).show();
  }, (e) => grok.shell.error(plainError(e, 'Your saved plans could not be listed. Try again in a moment.'))),
  {tooltip: 'Choose one of your saved plans'})));
  // ponytail: u2's FileInput keeps its last file when the box is emptied (it reports the empty path as missing), so
  // the plan is taken off with a button of its own
  planField.addOptions(form.build(() => iconButton('times', () => planField.value = null,
    {tooltip: 'No plan: compare with the plan the study table carries, if any, else the atlas'})));
  // with the box empty, the plan in use says where it comes from: the study table's own, or none
  const box = planField.root.querySelector('input');
  form.effect(() => {
    // (no result yet, as a study just followed from Plan: the plan its table carries)
    void error.value; // (a run that ends empty says so here: the table is read by then)
    const used = result.value?.analysis.col('plan') ? result.value.analysis.get('plan', 0) :
      (call.getParamValue('study') as DG.DataFrame | null)?.getTag(PLAN_TAG) ? 'the plan it was followed from' : null;
    if (box)
      box.placeholder = used && !planField.value ? `The study's own: ${used}` : PLACEHOLDER;
  });
  // the newest file picked: a slower read of an older one is dropped
  let latest = 0;
  whenSettled(form, () => {
    form.getInput('volume')?.root.after(planField.root);
    planField.root.after(note);
    form.effect(() => {
      const file = planField.value;
      // (the effect's first pass is the view opening with no plan picked: nothing changed)
      const opening = latest === 0;
      const mine = ++latest;
      stale.value = result.peek() !== undefined;
      note.textContent = '';
      if (!file) {
        call.setParamValue('plan', null);
        if (!opening)
          state.changed();
        return;
      }
      state.reading = file.readAsString().then(async (text) => {
        if (mine !== latest)
          return;
        const said = await applyPlan(parsePlan(text));
        if (mine !== latest)
          return;
        note.textContent = said;
        call.setParamValue('plan', text);
        state.changed();
      }).catch((e) => error.value = plainError(e, 'The plan file could not be read. Pick it again, or pick ' +
        'another.'));
    });
  });

  return planField;
}

/** The Plan vs observed tab's content, filled from each analysis, and what the view must own. */
export function planTab(form: FuncCallForm, result: ReadonlySignal<AnalysisFrames | undefined>,
  colorsOf: (r: AnalysisFrames) => {[group: string]: string}): {content: HTMLElement[], own: {dispose(): void}[]} {
  // the study against its plan: a sentence that sums it up, each group's effect against the assumed one, the
  // controls against their predicted growth, the end against the predicted day
  const planSummary = span('', 'pmx-plan-summary');
  Object.assign(planSummary.style, {fontWeight: '500', lineHeight: '1.45'});
  const noPlan = infoNote('');
  const planEffects = div([], 'pmx-plan-effects');
  const planNotes = notesBox();
  const controlLine = divH([], 'pmx-plan-controls');
  const endLine = divH([], 'pmx-plan-end');
  const atlasLine = divH([], 'pmx-plan-atlas');
  for (const line of [controlLine, endLine, atlasLine]) {
    Object.assign(line.style, {alignItems: 'baseline', gap: 'var(--dg-space-s)', lineHeight: '1.45',
      margin: '0 0 var(--dg-space-s)', flexWrap: 'wrap'});
  }
  const effectsBlock = block('Effect of each treatment', 'Whether each treated group grows slower than the control, ' +
    'then its effect against the range a study run as planned shows', planEffects);
  // what the comparison rests on (groups matched by order, another volume at randomization, another line or
  // endpoint): under the summary, as it holds for every block
  planNotes.style.margin = '0';
  // the rows' columns: one grid the rows share (subgrid), so every cell sits in its column and the bars under the
  // axis; the text columns as wide as their text while there is room, else wrapped to their widest word, the bar
  // taking the rest
  const planTable = div([], 'pmx-plan-table');
  Object.assign(planTable.style, {display: 'grid', columnGap: 'var(--dg-space-m)',
    gridTemplateColumns: 'auto auto auto auto minmax(160px, 1fr) auto'});
  const planRow = (cells: (HTMLElement | string)[], cls: string) => {
    const row = div(cells.map((c) => typeof c === 'string' ? span(c) : c), cls);
    // every cell at its own width but the bar and its axis, which span their column; a badge wraps like the text
    for (const [k, cell] of [...row.children].entries()) {
      const el = cell as HTMLElement;
      if (k !== 4)
        el.style.justifySelf = 'start';
      if (el.classList.contains('pmx-badge')) {
        Object.assign(el.style, {whiteSpace: 'normal', height: 'auto', minHeight: '16px', lineHeight: '1.2',
          padding: '2px var(--dg-space-s)'});
      }
    }
    Object.assign(row.style, {display: 'grid', gridTemplateColumns: 'subgrid', gridColumn: '1 / -1',
      alignItems: 'center', padding: '4px var(--dg-space-s)'});
    return row;
  };
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.analysis.get(c, 0);
    // a number of the analysis row, or null (a column with no value reads as an empty string)
    const n = (c: string) => finite(v(c));
    // the verdict, then its sentence beside it
    const said = (text: string) => {
      const el = span(text);
      el.style.flex = '1 1 320px';
      return el;
    };
    const planned = v('plan') !== null && v('plan') !== undefined && v('plan') !== '';
    const source = planned ? 'the plan' : `the ${v('line')} atlas`;
    // "a study run as planned", "a study run as the MC38 atlas predicts"
    const asRun = planned ? 'planned' : `${source} predicts`;
    noPlan.hidden = planned;
    (noPlan.lastElementChild as HTMLElement).innerHTML = '<b>No plan linked.</b> To compare each group\'s effect ' +
      'with the one assumed, save the plan in Plan (<i>File › Save</i>) and pick it in the Plan field, or drop its ' +
      `.pmxplan file there. Below, each group against the control, and the controls and the end of the study ` +
      (v('line') ? `against the ${v('line')} atlas.` : 'against the atlas once the tumor model is picked.');

    // each treated group: the test against the control first, then the effect assumed and the one observed, on one
    // axis, and the verdict against the plan
    const g = r.versusPlan;
    const colors = colorsOf(r);
    const LO = -50;
    const head = planRow(['Group', 'Result', 'Assumed', 'Observed', planned ? percentAxis(LO) : '',
      planned ? 'Against the plan' : ''], 'pmx-plan-head');
    Object.assign(head.style, {fontWeight: '500', borderBottom: '1px solid var(--dg-border-color, #e1e6ec)'});
    const rows: HTMLElement[] = [head];
    let held = 0;
    let unclear = 0;
    let compared = 0;
    let slows = 0;
    const why: HTMLElement[] = [];
    for (let i = 0; i < g.rowCount; i++) {
      const x = (c: string) => g.get(c, i);
      const nx = (c: string) => finite(x(c));
      const metric = metricWord(x('metric'));
      const [obs, assumed] = [nx('observed'), nx('assumed')];
      // the range a study run as planned shows (an open top: drawn to the axis's end)
      const [eLo, eHi] = [nx('expected_q025'), nx('expected_q975') ?? 100];
      if (x('verdict')) {
        compared++;
        held += x('verdict') === 'As planned' ? 1 : 0;
        unclear += x('verdict') === 'Inconclusive' ? 1 : 0;
      }
      const name = swatch(x('group'), colors[x('group')]);
      if (x('planned_as') && x('planned_as') !== x('group'))
        Tooltip.bind(name, `Planned as ${x('planned_as')}`);
      slows += x('activity') === 'Slows growth' ? 1 : 0;
      const note = noEffectNote(x, v);
      if (note)
        why.push(noteLine(note));
      if (x('activity') === 'Grows faster')
        why.push(noteLine(`${x('group')}: tumors grow faster than in the control (p ${pValue(x('p_value'))}).`));
      rows.push(planRow([name, activityBadge(x('activity'), x('p_value')),
        assumed === null ? (planned ? 'Not in the plan' : '—') : `${metric} ${Math.round(assumed)}%`,
        obs === null ? '—' : `${metric} ${Math.round(obs)}%`,
        obs === null || assumed === null || eLo === null ? span('') :
          planBar(obs, assumed, [eLo, eHi], LO, colors[x('group')], metric),
        planned ? verdictBadge(x('verdict')) : span('')],
      'pmx-plan-row'));
    }
    const legend = span('Each bar: the shaded range is what a study run as planned shows (95 of 100 labs), the dark ' +
      'tick the effect the plan assumed, the dot the observed effect. The effect is as planned ' +
      'when the dot falls in the shaded range and the group grows slower than the control, inconclusive when it ' +
      'falls there but the study cannot tell the group from the control.');
    Object.assign(legend.style, {fontSize: '12px', color: 'var(--dg-text-color-muted, #68727f)',
      margin: 'var(--dg-space-m) 0 0', lineHeight: '1.45'});
    legend.hidden = !planned;
    // what a group with no clear effect does and does not say, under the table
    const whyBox = div(why);
    whyBox.style.margin = why.length ? 'var(--dg-space-m) 0 0' : '0';
    planTable.replaceChildren(...rows);
    planEffects.replaceChildren(planTable, whyBox, legend);
    planNotes.replaceChildren(...notes(v('plan_notes')));

    // the controls' growth a day, against a study run as planned (the slope of log volume as a daily percent)
    const daily = (s: number | null) => s === null ? '—' : pct(Math.exp(s) - 1);
    const slope = n('control_slope');
    const [sLo, sHi] = [n('control_slope_q025'), n('control_slope_q975')];
    controlLine.replaceChildren(verdictBadge(v('control_verdict') || null), said(slope === null ? 'No control ' +
      `grew long enough to compare with ${source}.` : `The controls grew ${daily(slope)} a day; a study run as ` +
      `${asRun} grows ${daily(sLo)} to ${sHi === null ? 'more' : daily(sHi)} a day.`));

    // the day the median control reached the endpoint, against a study run as planned
    // (the range only: its median falls on a later measurement day than Plan's end of study, the median lab's day,
    // whenever that day is the last measurement before a gap, and would read as a contradiction)
    const [eLo, eHi, obs] = [n('end_q025'), n('end_q975'), n('observed_end')];
    const range = eLo === null ? '' : eHi === null ? `from day ${Math.round(eLo)}` :
      `between day ${Math.round(eLo)} and day ${Math.round(eHi)}`;
    // an open top: more than 2.5% of such studies never get there, said with their share
    const never = n('end_not_reached');
    const expected = `a study run as ${asRun} reaches it ${range}` + (eLo !== null && eHi === null && never !== null &&
      never > 0 ? `, and ${pct(never)} of such studies never do` : '');
    endLine.replaceChildren(verdictBadge(v('end_verdict') || null), said(obs === null ? 'The median control never ' +
      `reached ${v('endpoint')} mm³ in the study; ${expected}.` :
      `The median control reached ${v('endpoint')} mm³ on day ${obs}; ${expected}.`));
    // no plan and no model: nothing predicted to read them against (QA 0.1.3 A2: never a model assumed)
    if (!planned && !v('line')) {
      controlLine.replaceChildren(said('Pick the tumor model to see the controls against its atlas.'));
      endLine.replaceChildren(said('Pick the tumor model to see the end of the study against its atlas.'));
    }

    // what this study's controls could add to the line's atlas (coming soon): more so when they grew unlike it
    const controls = n('n_control');
    const unlike = v('control_verdict') && !/^As /.test(v('control_verdict'));
    const soon = badge('Coming soon', {variant: 'default'});
    Tooltip.bind(soon, `Coming soon: once checked, this study's controls join your private ${v('line')} atlas ` +
      `and inform the next ${v('line')} plans`);
    atlasLine.replaceChildren(soon, said(`This study's ${controls ?? ''} controls could join your private ` +
      `${v('line')} atlas, so the next plans for ${v('line')} start from them too.` + (unlike ? ' They grew unlike ' +
      'the atlas predicts, so they would add a lab the atlas does not cover yet.' : '')));

    // the sum of it, in one sentence
    const lower = (t: unknown) => String(t ?? '').toLowerCase();
    // never reached: the largest median control and the day the effects are read on instead (round 3, W5)
    const vol = r.volumes;
    const largest = Array.from({length: vol.rowCount}, (_, i) => i)
      .filter((i) => vol.get('group', i) === v('control') && Number.isFinite(vol.get('median', i)))
      .reduce((b: number | null, i) => b === null || vol.get('median', i) > vol.get('median', b) ? i : b, null);
    const neverReached = `the control median never reached ${v('endpoint')} mm³` + (largest === null ? '' :
      ` (largest ${Math.round(vol.get('median', largest))} mm³ on day ${vol.get('day', largest)})`) +
      `, so effects are read on day ${v('evaluation_day')}`;
    planSummary.textContent = [
      `${slows} of ${g.rowCount} treated group${g.rowCount > 1 ? 's' : ''} slow${g.rowCount > 1 ? '' : 's'} growth`,
      v('control_verdict') ? `Controls grew ${lower(v('control_verdict'))}` : null,
      planned && compared ? `${held} of ${compared} group${compared > 1 ? 's' : ''} as planned` +
        (unclear ? `, ${unclear} inconclusive` : '') : null,
      obs === null ? neverReached : `the endpoint was reached on day ${obs}` +
        (range ? ` (expected ${range})` : ''),
    ].filter((x) => x).join(' · ') + (planned ? ` (plan: ${v('plan')})` : '');
  });
  // the observed median control (solid, its measurements marked) over the plan's prediction (dashed, its 95% range)
  const predictedAxis = volumeAxis();
  const predictedGrowth = new BandChart(computed(() => result.value?.predicted), {x: 'day', lower: 'lower',
    upper: 'upper', line: 'median', split: 'series', xTitle: 'Day after randomization', yTitle: 'Tumor volume',
    log: predictedAxis.log, markers: true, dashed: [PREDICTED]});

  return {own: [predictedAxis.input, predictedGrowth], content: [planSummary, planNotes, noPlan, effectsBlock,
    block('Control growth against the prediction', 'The median control (solid), against a study\'s median ' +
      'control run as planned (dashed, its range in 95 of 100 labs)', controlLine, predictedAxis.row,
    sized(predictedGrowth, 340)),
    block('End of the study', 'The day the median control reached the endpoint, against a study run as planned',
      endLine)]};
  // "Contribute to the atlas" (atlasLine) is not shown until it does something: a "Coming soon" placeholder read as
  // unfinished to testers. Put the block back here when it is built:
  // block('Contribute to the atlas', 'This study\'s controls could strengthen ...', atlasLine)
}
