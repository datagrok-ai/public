/* What the Plan and Analyze views share, so the two read alike: the figures' formats (whole percents, three
   significant digits, whole days), the titled blocks, a group's swatch,
   the engine's notes, the Linear / Log volume axis of a growth chart, and the panes with the loader over the
   results. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {buttonWithIcon, computed, div, divH, divV, dot, icon, loader, RadioInput, Section, signal, span, StatCard,
  TabStrip, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';

/** A fraction as a whole percent: 0.734 → "73%". */
export const pct = (x: unknown) => typeof x === 'number' && Number.isFinite(x) ? `${Math.round(x * 100)}%` : '—';

/** A number to `digits` significant digits (an integer as it is). */
export const num = (v: unknown, digits = 3) =>
  typeof v === 'number' && !Number.isInteger(v) ? Number(v.toPrecision(digits)) : v;

/** A p-value: two significant digits, below 0.001 as "< 0.001". */
export const pValue = (p: unknown) => typeof p !== 'number' || !Number.isFinite(p) ? '—' :
  p < 0.001 ? '< 0.001' : String(num(p, 2));

/** A change in whole percent with a true minus sign: -4.2 → "−4%". */
export const signed = (x: number) => pct(x / 100).replace('-', '−');

/** A study day: whole, as the calendar counts them. */
export const dayText = (d: unknown) => typeof d === 'number' && Number.isFinite(d) ? `Day ${Math.round(d)}` : '—';

/** Narrowest the two panes get; below it the results pane scrolls instead of squeezing its charts. */
export const DECLARATION_MIN = 300;
export const RESULTS_MIN = 560;

/** Where x (%, from `lo` to 100) falls on an effect axis, as a CSS length. */
export const axisAt = (x: number, lo: number) =>
  `${(100 * (Math.max(lo, Math.min(100, x)) - lo) / (100 - lo)).toFixed(1)}%`;

/** A vertical guide on an effect axis at x: the twin's line at 0, light ones elsewhere. */
export function guide(x: number, lo: number): HTMLElement {
  const g = div([]);
  Object.assign(g.style, {position: 'absolute', left: axisAt(x, lo), top: '0', bottom: '0', width: '1px',
    background: x === 0 ? 'var(--dg-text-color-muted, #9aa3ad)' : 'var(--dg-border-color, #e1e6ec)'});
  return g;
}

/** A figure of the results with what it covers in its tooltip, owned by `owner` (the view's form). */
export function statCard(owner: {runInScope: <T>(f: () => T) => T}, label: string,
  value: ReadonlySignal<StatValue | undefined>, tip: () => string, format?: (x: StatValue) => string): HTMLElement {
  return owner.runInScope(() => {
    const c = new StatCard({label, value, format});
    Tooltip.bind(c.root, tip);
    return c.root;
  });
}

/** The figures in a grid that keeps them one size and wraps them as the pane narrows (u2's row does not wrap). */
export function cardGrid(cards: HTMLElement[], cls: string, minWidth: number): HTMLElement {
  const grid = div(cards, cls);
  Object.assign(grid.style, {display: 'grid', gridTemplateColumns: `repeat(auto-fit, minmax(${minWidth}px, 1fr))`,
    gap: 'var(--dg-space-m)'});
  return grid;
}

/** A Copy button for a text ready to paste, on its own row at the right. */
export function copyRow(text: () => string, what: string): HTMLElement {
  const copy = () => void navigator.clipboard.writeText(text())
    .then(() => grok.shell.info(`${what} copied`), () => grok.shell.warning(`${what}: the clipboard is not available`));
  const row = divH([buttonWithIcon('Copy', 'copy', copy)]);
  Object.assign(row.style, {justifyContent: 'flex-end', marginBottom: 'var(--dg-space-s)'});
  return row;
}

/** Whether a u2 funcForm is running its function (it says so in aria-busy), as a signal the form owns. */
export function formRunning(form: {root: HTMLElement, own: (f: () => void) => unknown}): ReadonlySignal<boolean> {
  const running = signal(false);
  const watch = new MutationObserver(() => running.value = form.root.getAttribute('aria-busy') === 'true');
  watch.observe(form.root, {attributes: true, attributeFilter: ['aria-busy']});
  form.own(() => watch.disconnect());
  return running;
}

/** A finished call's outputs by name. */
export const outputsOf = <T>(c: DG.FuncCall): T =>
  Object.fromEntries(c.func.outputs.map((p) => [p.name, c.outputs[p.name]])) as T;

/** A titled block of the results column. */
export function block(title: string, hint: HTMLElement | string, ...content: HTMLElement[]): HTMLElement {
  const s = new Section({title, hint, collapsible: false});
  content.forEach((c) => s.body.append(c));
  return s.root;
}

/** A viewer control given a fixed height (a viewer fills its host; the results column scrolls). */
export function sized<T extends {root: HTMLElement}>(v: T, height: number): HTMLElement {
  v.root.style.height = `${height}px`;
  return v.root;
}

/** A u2 table whose cells wrap instead of cutting their text with an ellipsis: on a narrow pane a column takes
 *  its widest word and the text goes on to the next line, nothing hidden. */
export function wrapCells(table: HTMLElement): void {
  for (const cell of table.querySelectorAll<HTMLElement>('th, td')) {
    Object.assign(cell.style, {whiteSpace: 'normal', overflow: 'visible', textOverflow: 'clip', height: 'auto',
      paddingTop: '3px', paddingBottom: '3px'});
  }
}

/** A group's name after a swatch of its colour, the colour every chart gives it. */
export function swatch(name: string, color: string): HTMLElement {
  const box = span('');
  Object.assign(box.style, {display: 'inline-block', width: '10px', height: '10px', borderRadius: '2px',
    background: color, marginRight: '6px', verticalAlign: '-1px', flex: '0 0 auto'});
  // the box stays beside the name's first line when a narrow column wraps it
  const el = span('');
  Object.assign(el.style, {display: 'inline-flex', alignItems: 'baseline'});
  el.append(box, span(name));
  return el;
}

/** A note on one line, the dot on the text's first line (u2 has no banner). */
export function noteLine(text: string, kind: 'warning' | 'success' | 'error' = 'warning'): HTMLElement {
  const mark = dot(kind);
  Object.assign(mark.style, {flex: '0 0 auto', marginTop: '0.45em'});
  const line = divH([mark, span(text)], 'pmx-note');
  Object.assign(line.style, {alignItems: 'flex-start', gap: 'var(--dg-space-s)', lineHeight: '1.45'});
  return line;
}

/** The engine's notes (" | "-separated) as note lines. */
export const notes = (text: unknown) => String(text ?? '').split(' | ').filter((n) => n).map((n) => noteLine(n));

/** A box for the notes under the block they are about. */
export function notesBox(): HTMLElement {
  const box = divV([], 'pmx-notes');
  Object.assign(box.style, {gap: 'var(--dg-space-s)', margin: 'var(--dg-space-l) 0 0'});
  return box;
}

/** An informative note in a tinted box (the link between the app's stages), `html` its text. */
export function infoNote(html: string): HTMLElement {
  const text = span('');
  text.innerHTML = html;
  const note = divH([icon('info-circle'), text], 'pmx-info-note');
  Object.assign(note.style, {gap: 'var(--dg-space-m)', alignItems: 'flex-start', marginTop: 'var(--dg-space-l)',
    padding: 'var(--dg-space-m) var(--dg-space-l)',
    borderRadius: '4px', lineHeight: '1.45', background: 'var(--dg-accent-light, #eaf3fb)',
    color: 'var(--dg-text-color, #1f2a37)'});
  return note;
}

/** A growth chart's volume axis, Linear or Log (Log by default: exponential growth reads as a line): the row that
 *  offers it above the chart, and whether Log is chosen (the chart's `log`). */
export function volumeAxis(): {input: RadioInput, row: HTMLElement, log: ReadonlySignal<boolean>} {
  const input = new RadioInput({items: ['Linear', 'Log'], buttons: true, value: 'Log',
    itemTooltips: {Linear: 'Linear volume axis', Log: 'Logarithmic volume axis'}});
  const row = divH([span('Volume axis'), input.root]);
  Object.assign(row.style, {justifyContent: 'flex-end', alignItems: 'center', gap: 'var(--dg-space-m)'});
  return {input, row, log: computed(() => input.value.value === 'Log')};
}

export type Stage = 'plan' | 'monitor' | 'analyze';
const STAGES: {stage: Stage, label: string, tip: string}[] = [
  {stage: 'plan', label: 'Plan', tip: 'Design the study before it runs: animals, end of study, calendar and power'},
  {stage: 'monitor', label: 'Monitor', tip: 'Coming soon: follow the study live against its plan, from implantation ' +
    'to the end'},
  {stage: 'analyze', label: 'Analyze', tip: 'Analyze the finished study as planned, each treated animal against its ' +
    'untreated twin'},
];

/** The app's stages in order, the current one marked; `open` brings another one up (Monitor is not built yet). A
 *  stage moves the whole page, so the strip goes in the view's ribbon, above the declaration and the results. */
export function stageStrip(current: Stage, open: (stage: Stage) => void): HTMLElement {
  const strip = divH([], 'pmx-stages');
  strip.setAttribute('role', 'navigation');
  strip.setAttribute('aria-label', 'Study stages');
  Object.assign(strip.style, {alignItems: 'center', gap: 'var(--dg-space-s)', fontSize: '13px', flexWrap: 'wrap'});
  STAGES.forEach(({stage, label, tip}, i) => {
    if (i > 0)
      strip.append(span('›'));
    const on = stage === current;
    const ready = stage !== 'monitor';
    const step = document.createElement(on || !ready ? 'span' : 'button');
    step.textContent = `${i + 1} ${label}`;
    Object.assign(step.style, {border: '0', background: 'none', padding: '2px 4px', font: 'inherit',
      borderBottom: on ? '2px solid var(--dg-accent, #2083d5)' : '2px solid transparent',
      color: on ? 'var(--dg-text-color, #1f2a37)' : ready ? 'var(--dg-accent, #2083d5)' :
        'var(--dg-text-color-muted, #9aa3ad)', fontWeight: on ? '600' : '400', cursor: on || !ready ? 'default' :
        'pointer'});
    if (on)
      step.setAttribute('aria-current', 'step');
    if (step instanceof HTMLButtonElement)
      step.addEventListener('click', () => open(stage));
    step.title = tip;
    strip.append(step);
  });
  return strip;
}

/** A tab of a stage's results: what it is for in one plain line (`intro`), its blocks, or for one not built yet
 *  (`soon`) what it will show. */
export interface ResultTab {id: string, label: string, intro: string, content?: HTMLElement[], soon?: string}

/** The results' tabs, one per question, under the figures that answer the main one; a view opens on its first. */
export function resultTabs(view: string, tabs: ResultTab[]): TabStrip {
  const strip = new TabStrip({tabs: tabs.map((t) => {
    const intro = span(t.intro, 'pmx-tab-intro');
    Object.assign(intro.style, {color: 'var(--dg-text-color-muted, #68727f)', lineHeight: '1.45'});
    const body = t.soon !== undefined ? [infoNote(`<b>Coming soon.</b> ${t.soon}`)] : t.content ?? [];
    const panel = divV([intro, ...body], `pmx-${view}-tab-${t.id}`);
    Object.assign(panel.style, {gap: 'var(--dg-space-l)', paddingTop: 'var(--dg-space-m)'});
    return {id: t.id, label: t.label, content: panel, tooltip: t.soon !== undefined ? `${t.label}: coming soon` :
      t.label};
  })});
  return strip;
}

/** An informative line (no warning): what the view shows no longer matches its form. */
export function staleNote(text: string): HTMLElement {
  const line = divH([icon('info-circle'), span(text)], 'pmx-stale');
  Object.assign(line.style, {gap: 'var(--dg-space-s)', alignItems: 'center', marginBottom: 'var(--dg-space-m)',
    color: 'var(--dg-text-color-muted, #68727f)'});
  return line;
}

/** A declaration form once built: its Advanced category folded (what most studies leave as it is), and each
 *  `auto` field's empty value worded as what the run does then. */
export function settleForm(form: {settled: Promise<unknown>, root: HTMLElement,
  getInput: (name: string) => {root: HTMLElement} | null | undefined}, auto: {[field: string]: string}): void {
  void form.settled.then(() => {
    // ponytail: funcForm opens every category; Advanced is folded by its header until u2 takes a
    // folded-categories option
    const advanced = [...form.root.querySelectorAll<HTMLElement>('.u2-form-category')]
      .find((h) => h.textContent?.trim().startsWith('Advanced'));
    if (advanced?.getAttribute('aria-expanded') === 'true')
      advanced.click();
    // ponytail: u2's NumberInput takes no placeholder; written on its input until u2 takes one
    for (const [name, text] of Object.entries(auto)) {
      const box = form.getInput(name)?.root.querySelector('input');
      if (box)
        box.placeholder = text;
    }
  });
}

/** The loader over the whole visible results pane while a run computes (`text`: what it does). */
export function busyOverlay(text: string): HTMLElement {
  const spin = loader(text);
  Object.assign(spin.style, {flexDirection: 'column', gap: 'var(--dg-space-l)', fontSize: '14px'});
  // u2's spinner at a size that reads as the page loading (its 16px suits an inline wait)
  Object.assign((spin.querySelector('.u2-loader-spinner') as HTMLElement).style, {width: '40px', height: '40px',
    borderWidth: '4px'});
  const overlay = div([spin], 'pmx-busy');
  Object.assign(overlay.style, {position: 'absolute', inset: '0', display: 'flex', alignItems: 'center',
    justifyContent: 'center', background: 'color-mix(in srgb, var(--dg-background, #fff) 75%, transparent)',
    zIndex: '2'});
  return overlay;
}

/** The results pane, scrolling, in a frame that does not scroll: the loader always covers the visible part. */
export function resultsFrame(content: HTMLElement[], overlay: HTMLElement): HTMLElement {
  const results = div(content, 'pmx-results');
  Object.assign(results.style, {height: '100%', overflow: 'auto', boxSizing: 'border-box',
    padding: 'var(--dg-space-m) var(--dg-space-l) var(--dg-space-xl, 24px)'});
  const frame = div([results, overlay], 'pmx-results-frame');
  Object.assign(frame.style, {position: 'relative', height: '100%', overflow: 'hidden',
    borderLeft: '1px solid var(--dg-border-color, #e1e6ec)'});
  return frame;
}

/** The declaration in a pane of its own: it scrolls by itself and clips what it holds, so nothing of it runs under
 *  the results. */
export function declarationPane(form: HTMLElement): HTMLElement {
  const pane = div([form], 'pmx-declaration');
  Object.assign(pane.style, {height: '100%', overflowY: 'auto', overflowX: 'hidden', boxSizing: 'border-box',
    padding: 'var(--dg-space-m) var(--dg-space-l) var(--dg-space-l)', minWidth: `${DECLARATION_MIN}px`});
  return pane;
}
