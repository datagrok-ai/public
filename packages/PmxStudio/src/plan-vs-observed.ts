/* The study against its plan, for the Analyze view: the Plan field of its form (a .pmxplan file from Datagrok's
   files, the computer or a drop, or one of the user's saved plans) and the "Plan vs observed" tab: a sentence that
   sums it up, each treated group's effect against the assumed one, the controls against their predicted growth, the
   end against the predicted day, and "Contribute to the atlas" (coming soon). The comparison itself is
   pmxAnalyzeStudy's (analyze.ts); this file only shows it. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {badge, computed, div, divH, iconButton, span, Tooltip} from '@datagrok-libraries/u2';
import type {BadgeVariant, ReadonlySignal, Signal} from '@datagrok-libraries/u2';
import {FileInput, FuncCallForm} from '@datagrok-libraries/u2/src/dg/index.js';

import {AnalysisFrames} from './analyze';
import {finite} from './engine';
import {BandChart} from './band-chart';
import {SavedPlan, myPlans, parsePlan, plansFolder} from './plan-file';
import {axisAt, block, guide, infoNote, noteLine, notes, notesBox, pct, pValue, sized, swatch,
  volumeAxis} from './ui';

/** A treated group's effect against the plan's on a shared percent axis from `lo` to 100: the range a study run as
 *  planned shows, shaded in the group's colour; the assumed effect a dark tick; the observed one a dot with its 95%
 *  interval. */
function planBar(observed: number, lower: number, upper: number, assumed: number, expected: [number, number],
  lo: number, color: string, metric: string): HTMLElement {
  const track = div([], 'pmx-plan-bar');
  Object.assign(track.style, {position: 'relative', height: '18px'});
  const band = div([]);
  Object.assign(band.style, {position: 'absolute', left: axisAt(expected[0], lo),
    width: `calc(${axisAt(expected[1], lo)} - ${axisAt(expected[0], lo)})`, top: '3px', height: '12px',
    borderRadius: '2px', background: `color-mix(in srgb, ${color} 22%, transparent)`});
  const range = div([]);
  Object.assign(range.style, {position: 'absolute', left: axisAt(lower, lo),
    width: `calc(${axisAt(upper, lo)} - ${axisAt(lower, lo)})`, top: '8px', height: '2px', background: color});
  const point = div([]);
  Object.assign(point.style, {position: 'absolute', left: `calc(${axisAt(observed, lo)} - 5px)`, top: '4px',
    width: '10px', height: '10px', borderRadius: '50%', background: color});
  const tick = div([]);
  Object.assign(tick.style, {position: 'absolute', left: `calc(${axisAt(assumed, lo)} - 1px)`, top: '0',
    width: '2px', height: '18px', background: 'var(--dg-text-color, #1f2a37)'});
  track.append(guide(0, lo), guide(50, lo), guide(100, lo), band, range, point, tick);
  // each mark's tooltip under the mark itself (a tooltip opens under its element's left edge)
  const observedText = `Observed ${metric} ${Math.round(observed)}% ` +
    `(95% CI ${Math.round(lower)}–${Math.round(upper)}%)`;
  Tooltip.bind(point, observedText);
  Tooltip.bind(range, observedText);
  Tooltip.bind(tick, `The plan assumed ${metric} ${Math.round(assumed)}%`);
  Tooltip.bind(band, `A study run as this plan shows ${metric} ${Math.round(expected[0])}% to ` +
    `${Math.round(expected[1])}% in 95 of 100 labs`);
  track.setAttribute('role', 'img');
  track.setAttribute('aria-label', `Observed ${pct(observed / 100)} (95% CI ${pct(lower / 100)} to ` +
    `${pct(upper / 100)}), assumed ${pct(assumed / 100)}, a study run as planned ${pct(expected[0] / 100)} to ` +
    `${pct(expected[1] / 100)}`);
  return track;
}

/** The observed effect and its interval: on narrow rows the interval goes to the next line whole. */
export function observedCell(effect: string, interval: string): HTMLElement {
  const ci = span(interval);
  ci.style.whiteSpace = 'nowrap';
  const el = span(`${effect} `);
  el.append(ci);
  return el;
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
  if (text === 'Inconclusive')
    Tooltip.bind(el, 'In the planned range, but the study cannot tell this group from the control');
  return el;
}

/** The test on the slopes as a badge: slower than the control in green, faster in amber, no clear effect in grey. */
export function activityBadge(activity: string, p: unknown): HTMLElement {
  const variant: BadgeVariant = activity === 'Slows growth' ? 'success' : activity === 'Grows faster' ? 'warning' :
    'default';
  const el = badge(activity, {variant});
  Tooltip.bind(el, `Tumor growth against the control, p ${pValue(p)}`);
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
  const head = `${g('group')}: no clear effect (p ${pValue(g('p_value'))}).`;
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

/** The Plan field, under the form's Study fields; a plan picked writes its pre-specified analysis into the form and
 *  its text into the call (`state`: the view's stale mark, error line and last analysis). The view owns it. */
export function planInput(form: FuncCallForm, call: DG.FuncCall, state: {stale: Signal<boolean>,
  error: Signal<string | null>, result: ReadonlySignal<AnalysisFrames | undefined>}): FileInput {
  const {stale, error, result} = state;
  // the Plan field also offers the user's saved plans; a plan picked writes its pre-specified analysis into the form
  const applyPlan = (p: SavedPlan) => {
    for (const [name, value] of [['line', p.predicted.line], ['endpoint', p.predicted.endpoint],
      ['correction', p.inputs.correction], ['alpha', p.inputs.alpha]] as [string, unknown][]) {
      if (value === undefined || value === null)
        continue;
      const field = form.getInput(name);
      if (field)
        field.value.value = value;
      call.setParamValue(name, value);
    }
  };
  // the Plan field: u2's FileInput (a path in Datagrok's files, a file from the computer, or a drop) under the Study
  // fields; the file's text goes into the call's hidden `plan` param, so a saved run carries its plan (ponytail:
  // placed by hand, u2 funcForm has no file field yet)
  const planField = form.runInScope(() => new FileInput({label: 'Plan', nullable: true,
    placeholder: PLACEHOLDER, tooltipText: 'The plan saved in Plan (.pmxplan), to compare the ' +
      'study with. Empty: the plan the study table carries, if any'}));
  planField.addOptions(form.runInScope(() => iconButton('list', () => void myPlans().then((files) => {
    if (files.length === 0) {
      grok.shell.info(`No saved plans yet: save one in Plan (it goes to ${plansFolder()})`);
      return;
    }
    DG.Menu.popup().items(files.map((f) => f.name), (name: string) =>
      planField.value.value = files.find((f) => f.name === name) ?? null).show();
  }), {tooltip: 'Choose one of your saved plans'})));
  // ponytail: u2's FileInput keeps its last file when the box is emptied (it reports the empty path as missing), so
  // the plan is taken off with a button of its own
  planField.addOptions(form.runInScope(() => iconButton('times', () => planField.value.value = null,
    {tooltip: 'No plan: compare with the plan the study table carries, if any, else the atlas'})));
  // with the box empty, the plan in use says where it comes from: the study table's own, or none
  const box = planField.root.querySelector('input');
  form.effect(() => {
    const used = result.value?.analysis.get('plan', 0);
    if (box)
      box.placeholder = used && !planField.value.value ? `The study's own: ${used}` : PLACEHOLDER;
  });
  // the newest file picked: a slower read of an older one is dropped
  let latest = 0;
  void form.settled.then(() => {
    form.getInput('volume')?.root.after(planField.root);
    form.effect(() => {
      const file = planField.value.value;
      const mine = ++latest;
      stale.value = result.peek() !== undefined;
      if (!file) {
        call.setParamValue('plan', null);
        return;
      }
      void file.readAsString().then((text) => {
        if (mine !== latest)
          return;
        applyPlan(parsePlan(text));
        call.setParamValue('plan', text);
      }).catch((e) => error.value = String(e?.message ?? e));
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
      if (el.classList.contains('u2-badge')) {
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
    noPlan.hidden = planned;
    (noPlan.lastElementChild as HTMLElement).innerHTML = '<b>No plan linked.</b> To compare each group\'s effect ' +
      'with the one assumed, save the plan in Plan (<i>Save plan</i>) and pick it in the Plan field, or drop its ' +
      `.pmxplan file there. Below, each group against the control, and the controls and the end of the study ` +
      `against the ${v('line')} atlas.`;

    // each treated group: the test against the control first, then the effect assumed, the one observed and its
    // interval, on one axis, and the verdict against the plan
    const g = r.versusPlan;
    const colors = colorsOf(r);
    const head = planRow(['Group', 'Result', 'Assumed', 'Observed (95% CI)', planned ? percentAxis(-50) : '',
      planned ? 'Against the plan' : ''], 'pmx-plan-head');
    const LO = -50;
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
      const metric = x('metric') === 'TC' ? 'T/C' : 'TGI';
      const [obs, lo, hi, assumed] = [nx('observed'), nx('lower'), nx('upper'), nx('assumed')];
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
        obs === null ? '—' : observedCell(`${metric} ${Math.round(obs)}%`,
          lo === null || hi === null ? '' : `(${Math.round(lo)}–${Math.round(hi)}%)`),
        obs === null || lo === null || hi === null || assumed === null || eLo === null ? span('') :
          planBar(obs, lo, hi, assumed, [eLo, eHi], LO, colors[x('group')], metric),
        planned ? verdictBadge(x('verdict')) : span('')],
      'pmx-plan-row'));
    }
    const legend = span('Each bar: the shaded range is what a study run as planned shows (95 of 100 labs), the dark ' +
      'tick the effect the plan assumed, the dot the observed effect with its 95% interval. The effect is as planned ' +
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
      `${source === 'the plan' ? 'planned' : `${source} predicts`} grows ${daily(sLo)} to ` +
      `${sHi === null ? 'more' : daily(sHi)} a day.`));

    // the day the median control reached the endpoint, against a study run as planned
    // (the range only: its median falls on a later measurement day than Plan's end of study, the median lab's day,
    // whenever that day is the last measurement before a gap, and would read as a contradiction)
    const [eLo, eHi, obs] = [n('end_q025'), n('end_q975'), n('observed_end')];
    const range = eLo === null ? '' : eHi === null ? `from day ${Math.round(eLo)}` :
      `between day ${Math.round(eLo)} and day ${Math.round(eHi)}`;
    // an open top: more than 2.5% of such studies never get there, said with their share
    const never = n('end_not_reached');
    const expected = `a study run as ${source === 'the plan' ? 'planned' : `${source} predicts`} reaches it ` +
      `${range}` + (eLo !== null && eHi === null && never !== null && never > 0 ?
      `, and ${pct(never)} of such studies never do` : '');
    endLine.replaceChildren(verdictBadge(v('end_verdict') || null), said(obs === null ? 'The median control never ' +
      `reached ${v('endpoint')} mm³ in the study; ${expected}.` :
      `The median control reached ${v('endpoint')} mm³ on day ${obs}; ${expected}.`));

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
    planSummary.textContent = [
      `${slows} of ${g.rowCount} treated group${g.rowCount > 1 ? 's' : ''} slow${g.rowCount > 1 ? '' : 's'} growth`,
      v('control_verdict') ? `Controls grew ${lower(v('control_verdict'))}` : null,
      planned && compared ? `${held} of ${compared} group${compared > 1 ? 's' : ''} as planned` +
        (unclear ? `, ${unclear} inconclusive` : '') : null,
      obs === null ? 'the endpoint was not reached' : `the endpoint was reached on day ${obs}` +
        (range ? ` (expected ${range})` : ''),
    ].filter((x) => x).join(' · ') + (planned ? ` (plan: ${v('plan')})` : '');
  });
  const predictedAxis = volumeAxis();
  const predictedGrowth = new BandChart(computed(() => result.value?.predicted), {x: 'day', lower: 'lower',
    upper: 'upper', line: 'median', split: 'series', xTitle: 'Day after randomization', yTitle: 'Tumor volume',
    log: predictedAxis.log, markers: true});

  return {own: [predictedAxis.input, predictedGrowth], content: [planSummary, planNotes, noPlan, effectsBlock,
    block('Control growth against the prediction', 'The median control, against the range of a study\'s median ' +
      'control run as planned (95 of 100 labs)', controlLine, predictedAxis.row, sized(predictedGrowth, 340)),
    block('End of the study', 'The day the median control reached the endpoint, against a study run as planned',
      endLine),
    block('Contribute to the atlas', 'This study\'s controls could strengthen the model\'s next predictions for ' +
      'the line', atlasLine)]};
}
