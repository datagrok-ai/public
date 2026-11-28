/* What the Plan and Analyze views share, so the two read alike: the figures' formats (whole percents, three
   significant digits, whole days), the titled blocks, a group's swatch,
   the engine's notes, the Linear / Log volume axis of a growth chart, and the panes with the loader over the
   results. */
import {maxOf, minOf} from './format';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {buttonWithIcon, computed, div, divH, divV, dot, effect, icon, loader, Section, signal, span,
  link, Splitter, StatCard, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import type {Dialog} from '@datagrok-libraries/u2';
import type {FuncCallForm, ViewHost} from '@datagrok-libraries/u2/dg';

import {_package} from './package';
import {cruxWake} from './crux';
import {addColorItem} from './colors';
import {PmxFigure} from './pmx-objects';
import {isFiniteNumber} from './engine';
import {minus} from './format';

/** A fraction as a whole percent, a true minus when negative: 0.734 → "73%", -0.64 → "−64%". */
export const pct = (x: unknown) => isFiniteNumber(x) ? minus(`${Math.round(x * 100) || 0}%`) : '—';

/** A TGI (percent) as every view and document shows it: a whole percent ("90%"), the tiles, the groups table and the
 *  CTD alike (round 4, analyst: "TGI 90.2%" beside "90%"). */
export const tgiText = (x: unknown) => pct(isFiniteNumber(x) ? x / 100 : x);

/** A number to `digits` significant digits (an integer as it is). */
export const num = (v: unknown, digits = 3) =>
  typeof v === 'number' && !Number.isInteger(v) ? Number(v.toPrecision(digits)) : v;

/** A p-value: two significant digits, below 0.001 as "< 0.001". */
export const pValue = (p: unknown) => !isFiniteNumber(p) ? '—' : p < 0.001 ? '< 0.001' : String(num(p, 2));

/** An adjusted p's hover, the unadjusted p: "Unadjusted p = 0.03", "Unadjusted p < 0.001" (never "= < 0.001");
 *  null: none known. Said even where the correction left it unchanged (QA 0.1.3 r4). */
export const unadjustedP = (raw: unknown): string | null => !isFiniteNumber(raw) ? null :
  `Unadjusted p ${raw < 0.001 ? '' : '= '}${pValue(raw)}`;

/** A change in whole percent with a true minus sign: -4.2 → "−4%". */
export const signed = (x: number) => pct(x / 100);

/** An interval in whole percent: "12–42%", and with a negative bound "−33 to 42%" (a hyphen and a dash read as
 *  "33 to 42"). */
export function percentRange(lo: number, hi: number): string {
  const [a, b] = [Math.round(lo), Math.round(hi)];
  const m = (x: number) => minus(String(x));
  return a < 0 || b < 0 ? `${m(a)} to ${m(b)}%` : `${a}–${b}%`;
}

/** A study day: whole, as the calendar counts them. */
export const dayText = (d: unknown) => isFiniteNumber(d) ? `Day ${Math.round(d)}` : '—';

/** An effect's metric as the app writes it: "T/C", else "TGI". */
export const metricWord = (metric: unknown) => metric === 'TC' ? 'T/C' : 'TGI';

/** Today's date (UTC): "2026-11-26". */
export const isoToday = () => new Date().toISOString().slice(0, 10);

/** Under a document's title: the study's ID (when it has one), where the document comes from and today's date. */
export const docSubtitle = (id: string | undefined, from: string) =>
  [id, from, `prepared ${isoToday()}`].filter((x) => x).join(' · ');

/** One of the user's own settings of the app (grok.userSettings, kept per user). */
export const userSetting = (key: string) => grok.userSettings.getValue('PmxStudio', key);
/** Sets it, sent at once: a reload before the 10 s sync would read the old value. */
export function setUserSetting(key: string, value: string): void {
  grok.userSettings.add('PmxStudio', key, value);
  void grok.userSettings.flush();
}

/** Narrowest the declaration pane gets. */
export const DECLARATION_MIN = 260;

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

/** A figure of the results with what it covers in its tooltip, owned by `owner` (the view's form); `subject`: what
 *  the figure is of (the study), a change from the previous result never read across two of them. */
export function statCard(owner: {build: <T>(f: () => T) => T}, label: string | ReadonlySignal<string>,
  value: ReadonlySignal<StatValue | undefined>, tip: () => string, format?: (x: StatValue) => string,
  compared = true, subject?: () => unknown): HTMLElement {
  // a figure no result follows (Compare's): no change, nothing to pin
  if (!compared) {
    return owner.build(() => {
      const c = new StatCard({label, value, format});
      c.root.classList.add('pmx-stat');
      Tooltip.bind(c.root, tip);
      return c.root;
    });
  }
  return owner.build(() => {
    // the change from the previous result, or from the figures pinned (cardGrid's right click): u2's delta, as a
    // number (a percent when the card shows one), in one colour (more animals is not better, more power is)
    const pinned = signal(false);
    const delta = signal<number | undefined>(undefined);
    let last: number | undefined;
    let base: number | undefined;
    let wasPinned = false;
    let change: [number, number] | null = null;
    const percent = !!format && /%$/.test(format(1));
    // a change as the card says it; one that rounds to nothing is no change (no "↓ −0%")
    const said = (d: number) => percent ? format!(Math.abs(d)) : String(+Math.abs(d).toFixed(2));
    let of: unknown;
    effect(() => {
      const v = value.value;
      const now = subject?.();
      if (now !== of) {
        of = now;
        last = undefined;
      }
      const pin = pinned.value;
      const x = isFiniteNumber(v) ? v : undefined;
      if (pin !== wasPinned) {
        wasPinned = pin;
        base = pin ? x ?? last : undefined;
      }
      if (v === undefined)
        return;
      const ref = pin ? base : last;
      delta.value = x !== undefined && ref !== undefined && said(x - ref) !== said(0) ? x - ref : undefined;
      if (delta.value !== undefined)
        change = [ref!, x!];
      last = x;
    });
    const c = new StatCard({label, value, format, delta, deltaFormat: (d) => `${d > 0 ? '+' : '−'}${said(d)}`});
    c.root.classList.add('pmx-stat');
    // the change itself, on the arrow: from what, to what, against which result
    const chip = c.root.querySelector<HTMLElement>('.u2-stat-delta');
    const shown = (x: number) => format ? format(x) : String(+x.toFixed(2));
    if (chip) {
      Tooltip.bind(chip, () => change ? `Now ${shown(change[1])}; ${pinned.value ? 'the figures pinned' :
        'the previous result of this study'}: ${shown(change[0])}` : '');
    }
    c.root.addEventListener('pmx-pin', (e) => pinned.value = (e as CustomEvent<boolean>).detail);
    // (the card's own sentence ends with its full stop, given or not: the tooltip runs on after it)
    Tooltip.bind(c.root, () => `${tip().replace(/([^.!?…])\s*$/, '$1.')}\n\nThe arrow: the change from ` +
      `${pinned.value ? 'the figures pinned' :
        'the previous result'}. Right-click to ${pinned.value ? 'unpin them' : 'pin these figures and compare'}.`);
    return c.root;
  });
}

/** The figures in a grid that keeps them one size and wraps them as the pane narrows (u2's row does not wrap). */
export function cardGrid(cards: HTMLElement[], cls: string, minWidth: number, pinnable = true): HTMLElement {
  const grid = div(cards, cls);
  Object.assign(grid.style, {display: 'grid', gridTemplateColumns: `repeat(auto-fit, minmax(${minWidth}px, 1fr))`,
    gap: 'var(--dg-space-m)'});
  // pin and compare: the figures now become what every next result is compared with, until unpinned
  const pin = (on: boolean) => {
    grid.classList.toggle('pmx-pinned', on);
    for (const c of cards)
      c.dispatchEvent(new CustomEvent('pmx-pin', {detail: on}));
  };
  if (!pinnable)
    return grid;
  // its right click: Pin or Unpin, the figures' actions (handlers.ts)
  ui.bind(new PmxFigure('The figures', {pin: {on: () => grid.classList.contains('pmx-pinned'), set: pin}}), grid,
    {contextMenu: true});
  return grid;
}

/** `text` to the clipboard, said when `what` names it; a browser that refuses the clipboard is said, never an
 *  unhandled rejection. */
export const copyText = (text: string, what?: string): void => void navigator.clipboard.writeText(text)
  .then(() => what && grok.shell.info(`${what} copied`),
    () => grok.shell.warning(`${what ?? 'Copy'}: the browser did not allow the clipboard`));

/** A Copy button for a text ready to paste, on its own row at the right. */
export function copyRow(text: () => string, what: string): HTMLElement {
  const copy = () => copyText(text(), what);
  const row = divH([buttonWithIcon('Copy', 'copy', copy)]);
  Object.assign(row.style, {justifyContent: 'flex-end', marginBottom: 'var(--dg-space-s)'});
  return row;
}

/** `fn` once the form is built, unless its view closed first: an effect registered on a disposed form throws a u2
 *  usage error (a test, or a user, can close a view before it settles). */
export function whenSettled(form: FuncCallForm, fn: () => void): void {
  void form.settled.then(() => {
    if (!form.scope.isDisposed)
      fn();
  });
}

/** The view has no Run: each change of the form recomputes it once the changes stop for a moment (`DEBOUNCE` ms)
 *  and nothing is missing or invalid. `changed` goes in the form's onInputChanged, and with `TABLE_EDITS` in the study
 *  table's data events (values typed one by one: a longer pause); `cancel` drops a pending recompute (a view that
 *  runs by itself, after filling the form, calls it). Not u2's showRun, which also saves the call,
 *  uploading its tables (the user's study) to the server. */
export function autoRun(form: () => FuncCallForm, run: () => unknown):
  {changed: (wait?: number) => void, cancel: () => void} {
  let timer = 0;
  const cancel = () => window.clearTimeout(timer);
  const changed = (wait = DEBOUNCE) => {
    cruxWake(); // (the engine at rest reopens during the pause, not after it)
    cancel();
    timer = window.setTimeout(() => {
      const f = form();
      if (f.missingRequired.peek().length === 0 && f.invalidFields.peek().length === 0)
        void run();
    }, wait);
  };
  return {changed, cancel};
}

const DEBOUNCE = 150;
/** The pause after an edit of the study table before the view recomputes (ms): one value typed, then the next. */
export const TABLE_EDITS = 1000;

/** Calls `edited` on each edit of the table the last result read (a value typed, rows pasted, added or removed): the
 *  table is the same input, so the form does not see it change. Returns what stops listening. */
export function onTableEdits(form: FuncCallForm, result: ReadonlySignal<unknown>, table: () => unknown,
  edited: () => void): () => void {
  let subs: {unsubscribe(): void}[] = [];
  const stop = () => {
    subs.forEach((s) => s.unsubscribe());
    subs = [];
  };
  form.effect(() => {
    void result.value;
    stop();
    const df = table();
    if (df instanceof DG.DataFrame)
      subs = [df.onDataChanged, df.onRowsAdded, df.onRowsRemoved].map((e) => e.subscribe(edited));
  });
  return stop;
}

/** What a call is asked, to tell whether a run still matches its form when it ends; `compareSchedules`,
 *  `detectableEffect` and `withTwin` say only how much the run computes, not what it is about. `only`: just those
 *  inputs. */
export const askedOf = (call: DG.FuncCall, only?: string[]) => JSON.stringify(call.func.inputs
  .filter((p) => only ? only.includes(p.name) :
    !['compareSchedules', 'detectableEffect', 'withTwin'].includes(p.name)).map((p) => {
    const v = call.inputs[p.name];
    // (a column param holds the platform's own object, which JSON cannot write)
    return v instanceof DG.DataFrame ? v.toCsv() : v instanceof DG.Column ? v.name : v;
  }));

/** A new result with each of its tables that holds what the previous one's held replaced by that very table: a panel
 *  that reads only unchanged tables is given the same objects, and does not redraw. */
export function keepSame<R extends object>(prev: R | undefined, next: R): R {
  if (!prev)
    return next;
  const out = {...next};
  for (const k of Object.keys(next) as (keyof R)[]) {
    const [a, b] = [prev[k], next[k]];
    if (a instanceof DG.DataFrame && b instanceof DG.DataFrame && a.rowCount === b.rowCount &&
      a.columns.names().join() === b.columns.names().join() && a.toCsv() === b.toCsv())
      out[k] = a as R[keyof R];
  }
  return out;
}

/** `result` as a panel that reads only its tables `keys` (and what `extra` gives) sees it: a new value only when one of
 *  those changed (keepSame keeps an unchanged table the same object), so the panel redraws only then. */
export function only<R extends object>(result: ReadonlySignal<R | undefined>, keys: (keyof R)[],
  extra: (r: R) => unknown[] = () => []): ReadonlySignal<R | undefined> {
  let last: R | undefined;
  let lastKey: unknown[] = [];
  return computed(() => {
    const r = result.value;
    const key = r ? [...keys.map((k) => r[k]), ...extra(r)] : [];
    if (r && last && key.length === lastKey.length && key.every((x, i) => x === lastKey[i]))
      return last;
    [last, lastKey] = [r, key];
    return r;
  });
}

/** A call's inputs by name as its function's own code receives them (a column the call names, its table's column),
 *  for a view to run the engine directly: the platform's FuncCall.call() then spends 100–200 ms of the main thread on
 *  the call's table outputs, after every run. The function stays registered for Grokky, the API and the tests. */
export function callInputs<T>(call: DG.FuncCall, extra: object = {}): T {
  const tables = call.func.inputs.map((p) => call.getParamValue(p.name)).filter((v) => v instanceof DG.DataFrame);
  return {...Object.fromEntries(call.func.inputs.map((p) => {
    const v = call.getParamValue(p.name);
    return [p.name, typeof v === 'string' && p.propertyType === DG.TYPE.COLUMN ?
      tables.map((t: DG.DataFrame) => t.col(v)).find((c) => c) ?? null : v];
  })), ...extra} as T;
}

/** A titled block of the results column; what it shows (`hint`) on the (i) after its title, not a line under it. */
export function block(title: string, hint: HTMLElement | string, ...content: HTMLElement[]): HTMLElement {
  const s = new Section({title, collapsible: false});
  s.header.classList.add('pmx-block-header');
  const said = () => typeof hint === 'string' ? hint : hint.textContent ?? '';
  if (said()) {
    const info = span('i', 'pmx-info');
    info.setAttribute('aria-label', said());
    ui.tooltip.bind(info, said);
    s.header.append(info);
  }
  s.root.classList.add('pmx-block');
  s.body.classList.add('pmx-block-body');
  content.forEach((c) => c.classList.contains('pmx-axis') ? s.header.append(c) : s.body.append(c));
  return s.root;
}

/** A viewer control given a fixed height (a viewer fills its host; the results column scrolls) and its host's whole
 *  width: its chart is absolutely placed, so the control has no width of its own, and in a flex host (a context
 *  panel's accordion pane) it shrank to 0 px (QA 0.1.3 A3/T10). The chart follows the width as it comes (opened, the
 *  panel shown, resized): the control's ResizeObserver. */
export function sized<T extends {root: HTMLElement}>(v: T, height: number): HTMLElement {
  Object.assign(v.root.style, {height: `${height}px`, width: '100%', minWidth: '0', boxSizing: 'border-box'});
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

/** A group's name after a swatch of its colour, the colour every chart gives it; a right click on the swatch offers
 *  Change color… as a chart's legend does. `group`: the group, when `name` says more. */
export function swatch(name: string, color: string, group: string = name): HTMLElement {
  const box = span('');
  Object.assign(box.style, {display: 'inline-block', width: '10px', height: '10px', borderRadius: '2px',
    background: color, marginRight: '6px', verticalAlign: '-1px', flex: '0 0 auto'});
  box.addEventListener('contextmenu', (e) => {
    const menu = DG.Menu.popup();
    if (!addColorItem(menu, group, color))
      return;
    e.preventDefault();
    e.stopPropagation();
    menu.show({causedBy: e});
  });
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
  // (what follows a tab is said on hover: engine.ts NOTE_TIP)
  const [shown, tip] = text.split('\t');
  const said = span(shown);
  if (tip)
    said.title = tip;
  const line = divH([mark, said], 'pmx-note');
  Object.assign(line.style, {alignItems: 'flex-start', gap: 'var(--dg-space-s)', lineHeight: '1.45'});
  return line;
}

/** What a failure tells the user: PMX Studio's own errors (a plain Error, or a sentence) are written to be read; any
 *  other (the engine, a bug) is `fallback`. The technical text goes to the console. */
export function plainError(e: any, fallback: string): string {
  // a message written for the user (a string, a plain Error) is said, not logged as a console error; anything else is
  // a failure, logged (QA 0.1.3 round 3, technician: "The plan needs 60 animals…" in the console as an error)
  const said = typeof e === 'string' || (e instanceof Error && e.constructor === Error && !!e.message);
  if (!said)
    console.error('PMX Studio:', e);
  return typeof e === 'string' ? inWords(e) : said ? inWords((e as Error).message) : fallback;
}

/** An engine's message in the user's words: without the engine function it came from ("stats.tumor_study_design: "),
 *  its bounds in symbols ("(> 0)") or British spelling. */
export const inWords = (message: string): string => message
  .replace(/\b(?:stats|tgi|pk|crux|nca)\.[\w.]+:\s*/g, '')
  .replace(/\s*\((?:[<>]=?|=)\s*-?[\d.]+\)/g, '')
  .replace(/randomis/g, 'randomiz').replace(/tumour/g, 'tumor').replace(/colour/g, 'color');

/** The engine's notes (" | "-separated) as note lines. */
export const notes = (text: unknown) => String(text ?? '').split(' | ').filter((n) => n).map((n) => noteLine(n));

/** A box for the notes under the block they are about. */
export function notesBox(): HTMLElement {
  const box = divV([], 'pmx-notes');
  Object.assign(box.style, {gap: 'var(--dg-space-s)', margin: 'var(--dg-space-l) 0 0'});
  return box;
}

/** The example: a published MC38 study (Yu et al., Nat Commun 2024, CC BY 4.0), its citation and changes in
 *  files/example/README.md. */
export const EXAMPLE = {file: 'example/mc38-ld4172-anti-pd1.csv', saveAs: 'pmx-example-mc38-yu2024-ccby.csv'};
export const EXAMPLE_SOURCE = 'A published MC38 study: Yu X. et al., Nature Communications 15 (2024), ' +
  'doi:10.1038/s41467-024-55006-2, CC BY 4.0';

/** Downloads the example study (a static file of the package). */
export const downloadExample = () => void _package.files.readAsText(EXAMPLE.file)
  .then((text) => DG.Utils.download(EXAMPLE.saveAs, text, 'text/csv'))
  .catch((e) => grok.shell.error(plainError(e, 'The example could not be downloaded. Try again in a moment.')));

/** The empty Monitor's and Analyze's one quiet line: a study to try (a static file of the package, no engine) and
 *  where one's own starts. */
export function exampleLine(): HTMLElement {
  const get = link('Download an example study', downloadExample);
  get.title = EXAMPLE_SOURCE;
  const line = span('No study at hand? ', 'pmx-muted');
  line.append(get, '. Your own starts in Plan.');
  return line;
}

/** A quiet line with the grey mark: something the user should know but need not act on. */
export function markLine(text: string): HTMLElement {
  const line = divH([icon('info-circle'), span(text)], 'ui-label-details pmx-mark-line');
  Object.assign(line.style, {gap: 'var(--dg-space-s)', alignItems: 'flex-start', lineHeight: '1.45'});
  return line;
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

/** A quiet line under a field, in the platform's detail style (ui.labels.details): what the app read there or can do,
 *  while `content` gives something (text, or text with links in it); null hides it. Call it in a scope (it follows
 *  `content` there): every feature the app keeps out of sight is offered this one way. */
export function contextHint(content: () => string | (string | HTMLElement)[] | null): HTMLElement {
  const e = div([], 'ui-label-details pmx-context-hint');
  effect(() => {
    const c = content();
    e.replaceChildren(...c === null ? [] : typeof c === 'string' ? [c] : c);
    e.style.display = c === null || c.length === 0 ? 'none' : '';
  });
  return e;
}

/** A value typed for a parameter a card sets, checked as its field did: a whole number from `min` to `max` (the
 *  value), or null (back to the computed one), else what is wrong with it. */
export function fixedValue(raw: unknown, min: number, max: number): {value: number | null} | {error: string} {
  if (raw === null || raw === undefined || raw === '')
    return {value: null};
  const x = typeof raw === 'number' ? raw : Number(raw);
  if (!Number.isFinite(x))
    return {error: 'Enter a number'};
  if (!Number.isInteger(x))
    return {error: 'Enter a whole number'};
  if (x < min || x > max)
    return {error: `Enter a number from ${min} to ${max}`};
  return {value: x};
}

/** A word in a contextHint that does something on a click (opens a choice, applies a value). */
export function hintLink(text: string, onClick: (e: MouseEvent | KeyboardEvent) => void, tooltip?: string):
  HTMLElement {
  // u2's link: its look, and Enter as a click
  const a = link(text, onClick);
  a.classList.add('pmx-hint-link');
  if (tooltip)
    a.title = tooltip;
  return a;
}

/** A growth chart's axis switch, plain words in the panel's header ("Linear · Log"), the one in use darker and
 *  underlined; `modes`: each word and its tooltip. Placed in a panel's content, the panel moves it to its header. */
export function axisSwitch(modes: [string, string][], now: {value: string, peek: () => string},
  onPick?: (m: string) => void): HTMLElement {
  const row = span('', 'pmx-axis pmx-noprint');
  const words = modes.map(([m, tip]) => {
    const w = span(m, 'pmx-axis-mode');
    w.title = tip;
    w.tabIndex = 0;
    w.setAttribute('role', 'button');
    const pick = () => {
      now.value = m;
      onPick?.(m);
      draw();
    };
    w.addEventListener('click', (e) => { e.stopPropagation(); pick(); });
    w.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        e.stopPropagation();
        pick();
      }
    });
    return w;
  });
  const draw = () => words.forEach((w, i) => {
    const on = modes[i][0] === now.peek();
    w.classList.toggle('is-on', on);
    w.setAttribute('aria-pressed', String(on));
  });
  words.forEach((w, i) => row.append(...(i ? [span('·', 'pmx-axis-sep'), w] : [w])));
  draw();
  return row;
}

/** A growth chart's volume axis, Linear or Log: Linear by default (as a lab's figures are drawn), else the user's
 *  last choice, which every chart then opens on. `row`: the switch for the panel's header; `log`: the chart's. */
export function volumeAxis(): {input: {dispose: () => void}, row: HTMLElement, log: ReadonlySignal<boolean>} {
  const mode = signal(userSetting('volumeAxis') === 'Log' ? 'Log' : 'Linear');
  const row = axisSwitch([['Linear', 'Volume on a linear axis'], ['Log', 'Volume on a logarithmic axis']], mode,
    (m) => setUserSetting('volumeAxis', m));
  return {input: {dispose: () => {}}, row, log: computed(() => mode.value === 'Log')};
}

export type Stage = 'home' | 'plan' | 'monitor' | 'analyze';
/** Marks `view` as `stage`'s: every stage's view is named after the app, so the app finds its own by this mark. */
export function markStage<V extends DG.ViewBase>(view: V, stage: Stage): V {
  view.root.dataset.pmxStage = stage;
  view.helpUrl = helpPage();
  return view;
}

/** The app's help, the package's README as published with it: the help panel of every view of the app. */
export const helpPage = () => `${_package.webRoot.replace(/\/?$/, '/')}README.md`;
/** The app's help in the platform's help panel, opened if hidden (the README of the package installed, not
 *  GitHub's). */
export function showHelp(): void {
  grok.shell.windows.help.visible = true;
  grok.shell.windows.help.showHelp(helpPage());
}
export const stageOf = (view: DG.ViewBase | null | undefined) => view?.root?.dataset.pmxStage as Stage | undefined;
export const STAGES: {stage: Stage, label: string, tip: string}[] = [
  {stage: 'home', label: 'Overview', tip: 'Every stage, and your studies'},
  {stage: 'plan', label: 'Plan', tip: 'Design the study before it runs: animals, end of study, calendar and power'},
  {stage: 'monitor', label: 'Monitor', tip: 'Follow the running study: who needs attention, and how it is going'},
  {stage: 'analyze', label: 'Analyze', tip: 'Analyze the finished study as planned, each treated animal against ' +
    'what it would measure untreated'},
];

/** An informative line (no warning): what the view shows no longer matches its form. */
export function staleNote(text: string): HTMLElement {
  const line = divH([icon('info-circle'), span(text)], 'pmx-stale');
  Object.assign(line.style, {gap: 'var(--dg-space-s)', alignItems: 'center', marginBottom: 'var(--dg-space-m)',
    color: 'var(--dg-text-color-muted, #68727f)'});
  return line;
}

/** A declaration form once built: its Advanced category folded (what most studies leave as it is), and each
 *  `auto` field's empty value worded as what the run does then. */
export function settleForm(form: FuncCallForm, auto: {[field: string]: string}): void {
  void form.settled.then(() => {
    // (a view closed before its form settled: nothing to fold)
    if (form.scope.isDisposed)
      return;
    foldColumns(form);
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

/** The study's column pickers (Animal, Group, Day, Volume) found by themselves: one quiet line, "change" opens them;
 *  they open by themselves only while one is empty. Nothing happens in a form without them (Plan). */
function foldColumns(form: FuncCallForm): void {
  const inputs = ['animal', 'group', 'day', 'volume'].map((n) => form.getInput(n))
    .filter((i): i is NonNullable<typeof i> => !!i);
  if (!inputs.length)
    return;
  let asked = false;
  const line = span('');
  const show = () => {
    const missing = inputs.some((i) => i.value === null || i.value === undefined);
    for (const i of inputs)
      i.root.style.display = asked || missing ? '' : 'none';
    line.hidden = missing;
  };
  line.append('Columns: found automatically · ', hintLink('change', () => {
    asked = !asked;
    show();
  }));
  inputs[0].root.parentElement?.insertBefore(line, inputs[0].root);
  form.effect(show);
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

/** A thin moving bar along the top of the results while they are recomputed, what was shown left in view. */
export function busyBar(): HTMLElement {
  return div([], 'pmx-busy-bar');
}

/** The results pane, scrolling down (never sideways: it fits its width, and a table or calendar too wide for it
 *  scrolls by itself, studio.css), in a frame that does not scroll: the loader always covers the visible part, the
 *  bar sits on its top edge. */
export function resultsFrame(content: HTMLElement[], ...overlays: HTMLElement[]): HTMLElement {
  const results = div(content, 'pmx-results');
  Object.assign(results.style, {height: '100%', overflowX: 'hidden', overflowY: 'auto', boxSizing: 'border-box',
    padding: 'var(--dg-space-m) var(--dg-space-l) var(--dg-space-xl, 24px)'});
  const frame = div([results, ...overlays], 'pmx-results-frame');
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

/** A stage's two panes: its declaration, and its results under the run's error and, while a field is missing or
 *  invalid, a note that they are behind (dimmed); the first run under the loader, every later one under the thin
 *  bar, what was shown left in view. `before`: shown between the note and the results (an empty state). */
export function stagePanes(form: FuncCallForm, content: HTMLElement, o: {busy: ReadonlySignal<boolean>,
  result: ReadonlySignal<unknown>, error: ReadonlySignal<string | null>, stale: ReadonlySignal<boolean>,
  busyText: string, staleText: string, before?: HTMLElement[], empty?: ReadonlySignal<boolean>}): Splitter {
  const overlay = busyOverlay(o.busyText);
  const bar = busyBar();
  form.effect(() => {
    overlay.hidden = !o.busy.value || o.result.value !== undefined;
    bar.hidden = !o.busy.value || o.result.value === undefined;
  });
  const failure = div([], 'pmx-stage-error');
  form.effect(() => {
    // an empty state (o.empty) is said plainly; a failure in the error's colour
    const empty = o.empty?.value ?? false;
    failure.style.color = empty ? '' : 'var(--dg-error-color, #c62828)';
    failure.hidden = o.error.value === null;
    failure.replaceChildren(...(o.error.value === null ? [] : [empty ? span(o.error.value) :
      noteLine(o.error.value, 'error')]));
    // the results of the inputs before are not this failure's: hidden, not left readable under it
    content.style.display = o.error.value !== null && !empty ? 'none' : '';
  });
  const outdated = staleNote(`A field is missing or invalid: ${o.staleText}`);
  form.effect(() => {
    outdated.hidden = !o.stale.value || form.missingRequired.value.length + form.invalidFields.value.length === 0;
    content.style.opacity = o.stale.value ? '0.5' : '';
  });
  return new Splitter([declarationPane(form.root),
    resultsFrame([failure, outdated, ...o.before ?? [], content], overlay, bar)],
  {direction: 'horizontal', sizes: [0.2, 0.8], minSize: DECLARATION_MIN});
}

/** u2 auto-picks a column for every column param, the optional ones too (the nearest name): an optional column keeps
 *  the pick only when its name says what it is, so a study without body weights does not read its volumes as
 *  weights. */
export function keepNamedColumns(form: FuncCallForm, optional: {[input: string]: (column: string) => boolean}): void {
  form.effect(() => {
    void form.getInput('study')?.value;
    whenSettled(form, () => {
      for (const [name, named] of Object.entries(optional)) {
        const input = form.getInput(name);
        const column = input?.peek();
        if (input && typeof column === 'string' && !named(column))
          input.value = null;
      }
    });
  });
}

/** The name of a body weight column, for the optional column above. */
export const isBodyWeight = (c: string) => /weight|\bbw\b/i.test(c) && !/tum/i.test(c);

/** The dialogs open now, by purpose ('preparing': its launcher is still reading what it will show). */
const modals = new Map<string, DG.Dialog | 'preparing'>();

/** The number of purposes held: for the tests. */
export const modalsHeld = () => modals.size;

/** Whether a dialog is on screen: one closed has left the page (or is hidden). */
const onScreen = (d: DG.Dialog) => d.root.isConnected && getComputedStyle(d.root).display !== 'none';

/** The dialog open for `purpose` (undefined: none); a closed one is forgotten here. */
export function modalOf(purpose: string): DG.Dialog | undefined {
  const d = modals.get(purpose);
  if (!d || d === 'preparing')
    return undefined;
  if (onScreen(d))
    return d;
  modals.delete(purpose);
  return undefined;
}

/** Frees `purpose` at once, its dialog just closed (a launcher opening another in its place). */
export const releaseModal = (purpose: string) => void modals.delete(purpose);

/** Puts an open dialog in front of the others. */
function raise(d: DG.Dialog): void {
  const parent = d.root.parentElement;
  if (parent && parent.lastElementChild !== d.root)
    parent.append(d.root);
}

/** A new dialog for `purpose`, or null when one is open for it already (brought to the front instead): every dialog
 *  the app opens is made here, so a second click never opens a second dialog over the first. The purpose is free
 *  again once its dialog closes (leaves the page). */
export function modal(purpose: string, title: Parameters<typeof ui.dialog>[0] = purpose):
  DG.Dialog | null {
  const open = modalOf(purpose);
  if (open) {
    raise(open);
    return null;
  }
  // (while 'preparing', only its own launcher gets here: a second click stopped in preparing)
  // (a dialog made and never shown, which reports no close, is forgotten here; every launcher shows its dialog in the
  // turn it makes it)
  for (const [k, v] of [...modals]) {
    if (v !== 'preparing' && !onScreen(v))
      modals.delete(k);
  }
  // (with the app's help: a dialog without one turned the help panel to "Help unavailable" while open, QA 0.1.3
  // designer D20)
  const d = ui.dialog(typeof title === 'object' ? {...title, helpUrl: title.helpUrl ?? helpPage()} :
    {title: title ?? '', helpUrl: helpPage()});
  modals.set(purpose, d);
  // free as soon as it closes (its root may stay on the page a moment while it fades)
  d.onClose.subscribe(() => {
    if (modals.get(purpose) === d)
      modals.delete(purpose);
  });
  // (every dialog of the app shown with its buttons in the page as they show: shownOrder)
  const show = d.show.bind(d);
  d.show = (options?: Parameters<typeof show>[0]) => {
    const shown = show(options);
    shownOrder(d.root);
    keepFocus(d.root);
    return shown;
  };
  return d;
}

/** What Tab stops at in a dialog, in page order. */
const focusables = (root: HTMLElement): HTMLElement[] => [...root.querySelectorAll<HTMLElement>(
  'button, input, select, textarea, a[href], [tabindex]:not([tabindex="-1"])')]
  .filter((e) => e !== root && !(e as HTMLButtonElement).disabled && e.offsetParent !== null);

/** A dialog shown takes the focus, its first field, else its first button (Cancel: nothing happens on Enter by
 *  mistake), unless its launcher put it somewhere in it already; Tab and Shift+Tab go round inside it (QA 0.1.3 r5:
 *  the focus sat on the dialog itself and the first Tab went to the page behind). */
export function keepFocus(root: HTMLElement): void {
  setTimeout(() => {
    const at = document.activeElement;
    if (at && at !== root && root.contains(at))
      return;
    const all = focusables(root);
    (all.find((e) => ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.tagName)) ??
      root.querySelector<HTMLElement>('.d4-dialog-footer button:not([disabled])') ?? all[0])?.focus();
  });
  root.addEventListener('keydown', (e) => {
    if (e.key !== 'Tab')
      return;
    const all = focusables(root);
    if (!all.length)
      return;
    const i = all.indexOf(document.activeElement as HTMLElement);
    const to = i < 0 ? all[0] : e.shiftKey && i === 0 ? all[all.length - 1] :
      !e.shiftKey && i === all.length - 1 ? all[0] : null;
    if (to) {
      e.preventDefault();
      to.focus();
    }
  });
}

/** The u2 dialogs on screen, by purpose. */
const shownDialogs = new Map<string, Dialog>();

/** Closes every dialog of the app on screen as Cancel does (its close handlers run, nothing is applied): the app moved
 *  to another page under it (QA 0.1.3: Back left Compare with other studies open over the next stage). */
export function closeModals(): void {
  for (const [purpose, d] of [...modals]) {
    if (d === 'preparing')
      continue;
    modals.delete(purpose);
    if (onScreen(d))
      d.close();
  }
  for (const [purpose, d] of [...shownDialogs]) {
    shownDialogs.delete(purpose);
    if (d.isOpen.peek())
      d.close('cancel');
  }
}

/** A u2 dialog for `purpose`, made by `make` and shown pinned to `host`'s view, unless one is open for it already (as
 *  modal() does for the platform's): a second click never stacks a second dialog over the first. */
export function showOnce(host: ViewHost, purpose: string, make: () => Dialog): void {
  if (shownDialogs.get(purpose)?.isOpen.peek())
    return;
  shownDialogs.set(purpose, host.show(host.build(make)));
}

/** For a launcher that reads something before its dialog shows: a second click meanwhile, or while its dialog is
 *  open, does nothing (the open dialog comes to the front). Undefined when it did nothing. */
export async function preparing<T>(purpose: string, prepare: () => Promise<T>): Promise<T | undefined> {
  const open = modalOf(purpose);
  if (open || modals.get(purpose) === 'preparing') {
    if (open)
      raise(open);
    return undefined;
  }
  modals.set(purpose, 'preparing');
  try {
    return await prepare();
  } finally {
    if (modals.get(purpose) === 'preparing')
      modals.delete(purpose);
  }
}

/** A dialog's OK button named for what it does (Export, Exclude, Save): returned, as getButton finds it only while it
 *  still says OK. Call it after onOK. */
export function okNamed(d: DG.Dialog, verb: string): HTMLButtonElement {
  const ok = d.getButton('OK');
  ok.textContent = verb;
  return ok;
}

/** A dialog that only shows something: its one way out named Close (Datagrok names it Cancel). */
export function closeNamed(d: DG.Dialog): DG.Dialog {
  const cancel = d.getButton('CANCEL');
  if (cancel)
    cancel.textContent = 'Close';
  return d;
}

/** Shows `d` in the middle of the window, still draggable by its title and with its ✕: the platform's `center` puts
 *  it in a fixed frame that has neither, and without it the dialog opens by the menu that asked for it. */
export function showCentred(d: DG.Dialog, options: {width?: number} = {}): DG.Dialog {
  // a tooltip open as the dialog opens (a card's, under the pointer as its menu was used) stayed over the dialog
  ui.tooltip.hide();
  d.show(options);
  const r = d.root;
  r.style.left = `${Math.max(0, (window.innerWidth - r.offsetWidth) / 2)}px`;
  r.style.top = `${Math.max(0, (window.innerHeight - r.offsetHeight) / 3)}px`;
  return d;
}

/** A dialog's buttons put in the page in the order they show, so Tab goes as the eye does: the platform lays its
 *  footer out reversed, and Tab went Discard, Save, Cancel over Cancel · Save · Discard (QA 0.1.3 browser check). */
export function shownOrder(root: HTMLElement): void {
  const first = root.querySelector<HTMLElement>('.d4-dialog-footer button, .d4-command-bar button');
  const bar = first?.parentElement;
  if (!bar)
    return;
  const buttons = [...bar.children].filter((e): e is HTMLElement => e.tagName === 'BUTTON');
  const left = (e: HTMLElement) => e.getBoundingClientRect().left;
  const shown = [...buttons].sort((a, b) => left(a) - left(b));
  if (shown.every((b, i) => b === buttons[i]))
    return;
  const css = getComputedStyle(bar);
  if (css.flexDirection.endsWith('-reverse')) {
    // (the same side, read the other way)
    const flip: {[k: string]: string} = {'flex-start': 'flex-end', 'start': 'end', 'normal': 'flex-end',
      'flex-end': 'flex-start', 'end': 'start'};
    bar.style.justifyContent = flip[css.justifyContent] ?? css.justifyContent;
    bar.style.flexDirection = css.flexDirection.replace('-reverse', '');
  }
  for (const b of buttons)
    b.style.order = '0';
  const at = document.createComment('');
  buttons[0].before(at);
  at.replaceWith(...shown);
}

/** A dialog's one short line, quiet, its detail (if any) on hover. */
export function dialogLine(text: string, tip?: string): HTMLElement {
  const e = ui.divText(text);
  Object.assign(e.style, {color: 'var(--grey-5, #68727f)', marginBottom: '8px', maxWidth: '520px'});
  if (tip)
    ui.tooltip.bind(e, tip);
  return e;
}

/** A line of an object's context panel: a label and its value, the detail in a tooltip. */
export interface PanelFact {label: string, value: string | HTMLElement, tip?: string}

/** An object's answer, first in its panel: `text` (its result in the words of the table it comes from), `sub` (a
 *  quiet line: the p-value, the range), `tone` the mark beside it (ok: an effect; warn: against; neutral: none). */
export interface PanelVerdict {text: string, sub?: string, tone: 'ok' | 'warn' | 'neutral'}

/** A folded section of an object's context panel, built when opened; `summary` a few quiet words after its title. */
export interface PanelSection {title: string, content: () => HTMLElement, open?: boolean, summary?: string}

/** What the context panel shows for an object (a group, an animal, a day), alike for each: its name after its
 *  colour, a quiet line under it, its answer (none where there is none: never made up), a few facts, then folded
 *  sections for the detail. */
export function objectPanel(o: {title: string, color?: string, subtitle?: string, subtitleTip?: string,
  verdict?: PanelVerdict | null, /** right under the answer (an animal's effect bar) */ chart?: HTMLElement,
  facts: PanelFact[], sections?: PanelSection[], /** its sections' accordion, once built */
  accordion?: (acc: DG.Accordion) => void}): HTMLElement {
  const head = document.createElement('h3');
  head.className = 'pmx-object-title';
  head.style.margin = '0';
  head.append(o.color ? swatch(o.title, o.color) : span(o.title));
  const facts = factList(o.facts);
  facts.className = 'pmx-object-facts';
  const el = div([head], 'pmx-object-panel');
  el.style.padding = '8px 12px 10px';
  if (o.subtitle) {
    const sub = div([o.subtitle], 'pmx-object-subtitle');
    sub.style.color = 'var(--dg-text-color-muted, #7d8792)';
    if (o.subtitleTip)
      sub.title = o.subtitleTip;
    el.append(sub);
  }
  if (o.verdict) {
    const v = div([span(o.verdict.text)], `pmx-object-verdict pmx-${o.verdict.tone}`);
    if (o.verdict.sub)
      v.append(div([o.verdict.sub], 'pmx-object-verdict-sub'));
    el.append(v);
  }
  if (o.chart)
    el.append(o.chart);
  if (o.facts.length)
    el.append(facts);
  if (o.sections?.length) {
    const acc = ui.accordion();
    for (const s of o.sections) {
      const pane = acc.addPane(s.title, s.content, !!s.open);
      // the summary beside the title, quiet (no header found: in the title)
      const header = s.summary ? pane.root.querySelector('.d4-accordion-pane-header') : null;
      if (header) {
        const sum = span(` ${s.summary}`, 'pmx-pane-summary');
        Object.assign(sum.style, {color: 'var(--dg-text-color-muted, #7d8792)', fontWeight: '400',
          marginLeft: '6px'});
        header.append(sum);
      } else if (s.summary)
        pane.name = `${s.title} · ${s.summary}`;
    }
    acc.root.classList.add('pmx-object-sections');
    el.append(acc.root);
    o.accordion?.(acc);
  }
  return el;
}

/** Labels and values, two columns (an object's facts, or a panel section's). */
export function factList(list: PanelFact[]): HTMLElement {
  const facts = document.createElement('dl');
  Object.assign(facts.style, {display: 'grid', gridTemplateColumns: 'max-content 1fr', gap: '4px 12px',
    margin: '8px 0'});
  for (const f of list) {
    const dt = document.createElement('dt');
    dt.textContent = f.label;
    dt.style.color = 'var(--dg-text-color-muted, #7d8792)';
    const dd = document.createElement('dd');
    dd.style.margin = '0';
    dd.append(f.value);
    if (f.tip)
      ui.tooltip.bind(dd, f.tip); // the panel is built outside any u2 owner: u2's Tooltip would throw
    facts.append(dt, dd);
  }
  return facts;
}

/** Plain lines, one under the other (a panel section's text). */
export function textLines(text: (string | HTMLElement)[]): HTMLElement {
  const el = div(text.map((t) => div([t])));
  Object.assign(el.style, {display: 'flex', flexDirection: 'column', gap: '4px'});
  return el;
}

/** Shows `x` in the context panel, opening it: the user asked for depth (a click on a group or an animal). An animal
 *  or a group is given as its object (PmxAnimal, PmxGroup): its handler renders its details and its actions. A tick
 *  later, as NCA does, so the platform's own pick on the same click does not replace it. Returns `x`. */
export function showDepth<T extends object>(x: T): T {
  setTimeout(() => {
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = x;
  }, 0);
  return x;
}

/** A tiny line of `ys` over `xs` (a context panel's curve): no axes, the last point dotted. */
export function sparkline(xs: number[], ys: number[], color = 'var(--dg-accent, #2083d5)', w = 120, h = 28):
HTMLElement {
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('width', String(w));
  svg.setAttribute('height', String(h));
  const box = span('', 'pmx-sparkline');
  box.style.display = 'inline-flex';
  box.append(svg);
  if (xs.length === 0)
    return box;
  const [x0, x1, y0, y1] = [minOf(xs), maxOf(xs), minOf(ys), maxOf(ys)];
  const px = (x: number) => 2 + (x1 > x0 ? (x - x0) / (x1 - x0) : 0.5) * (w - 4);
  const py = (y: number) => h - 2 - (y1 > y0 ? (y - y0) / (y1 - y0) : 0.5) * (h - 4);
  const line = document.createElementNS(NS, 'polyline');
  line.setAttribute('points', xs.map((x, i) => `${px(x)},${py(ys[i])}`).join(' '));
  line.setAttribute('fill', 'none');
  line.setAttribute('stroke', color);
  line.setAttribute('stroke-width', '1.5');
  const end = document.createElementNS(NS, 'circle');
  end.setAttribute('cx', String(px(xs[xs.length - 1])));
  end.setAttribute('cy', String(py(ys[ys.length - 1])));
  end.setAttribute('r', '2.5');
  end.setAttribute('fill', color);
  svg.append(line, end);
  return box;
}
