/* PMX Studio's home: the stages as cards joined the way Flow joins its steps (each edge names the file that passes
   between two stages), the growth of a study behind them, the user's studies, and the way in for a study under way:
   a file dropped anywhere on the page, or picked in the Open a study dialog, is read and opened in its stage. */
import '../css/hub.css';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable, computed, signal} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {markStage, stageOf} from './ui';
import {steps} from './study-header';
import type {Stage} from './ui';
import {openStage, showStage} from './stages';
import {planView} from './plan-view';
import {analyzeView} from './analyze-view';
import {monitorView, openStudyFile} from './monitor-view';
import {newStudyTable, parseStudy, tableOf} from './study-file';
import {PLAN_TAG, myPlans, parsePlan} from './plan-file';
import {LibraryEntry, libraryReader, matches} from './library';
import {dayToday} from './calendar';
import {Brought, OPENED_EXTS, OpenStage, bring, readFile, suggestedStage} from './open-study';

/** The home view. */
export async function hubView(): Promise<DG.ViewBase> {
  return markStage(ViewHost.create({name: 'PMX Studio', content: (self) => hubContent(self)}).view, 'home');
}

// ---- small builders: the page is static markup with a few listeners, so plain elements do -------------------------

function h<K extends keyof HTMLElementTagNameMap>(tag: K, cls?: string, children: (Node | string)[] = [],
  attrs: {[name: string]: string} = {}): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  if (cls)
    e.className = cls;
  for (const [k, v] of Object.entries(attrs))
    e.setAttribute(k, v);
  e.append(...children);
  return e;
}

/** An inline icon in the current colour (the paths are constants of this file). */
function icon(paths: string, size = 22, width = 1.6): SVGSVGElement {
  const t = document.createElement('template');
  t.innerHTML = `<svg width="${size}" height="${size}" viewBox="0 0 24 24" fill="none" stroke="currentColor" ` +
    `stroke-width="${width}" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${paths}</svg>`;
  return t.content.firstElementChild as SVGSVGElement;
}

const ICON = {
  plan: '<rect x="3.5" y="5" width="17" height="15" rx="2"/>' +
    '<path d="M3.5 10h17M8 3v4M16 3v4M8 14h2M12 14h2M16 14h.5M8 17h2"/>',
  monitor: '<path d="M3 12h4l2.5-6 4 12 2.5-6H21"/>',
  analyze: '<path d="M4 4v16h16"/><path d="M7 16c3-1 5-4 7-7s3-4 5-5"/><path d="M7 17c3 0 6-1 8-3s3-3 4-4" ' +
    'stroke-dasharray="2 2.5"/>',
  open: '<path d="M12 15V4M7.5 8.5L12 4l4.5 4.5"/><path d="M4 15v3a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2v-3"/>',
  drop: '<path d="M12 4v11M7.5 10.5L12 15l4.5-4.5"/><path d="M4 15v3a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2v-3"/>',
  file: '<path d="M6 3h8l4 4v14H6z"/><path d="M14 3v4h4M9 12h6M9 15h6M9 18h4"/>',
  ok: '<path d="M5 12l5 5 9-10"/>',
  none: '<path d="M6 12h12"/>',
  warn: '<path d="M12 3l9.5 17h-19z"/><path d="M12 10v4M12 17h0"/>',
  info: '<circle cx="12" cy="12" r="9"/><path d="M12 11v5M12 8h0"/>',
  add: '<path d="M12 5v14M5 12h14"/>',
};

// ---- the growth behind the page ------------------------------------------------------------------------------------

/** A study as the atlas draws one, for the page's background: a control on a reduced Gompertz curve and two treated
 *  arms (kill rates from day 2), each with the band of 90% of its animals widening over 28 days and its median.
 *  Illustrative parameters, fixed: the same picture every time. */
function growthArt(): SVGSVGElement {
  const xAt = (t: number) => (t / 28 * 1440).toFixed(1);
  const yAt = (v: number) => (400 - v * 70).toFixed(1);
  // Gompertz on a linear axis, as tumors grow: slow from the implant, fastest midway, levelling off at each group's
  // plateau (the treated lower and later); the 90% band log-normal around the median
  const gompertz = (vMax: number, b: number) => (t: number) =>
    vMax * Math.exp(Math.log(0.15 / vMax) * Math.exp(-b * t));
  const steps = Array.from({length: 57}, (_, i) => i / 2);
  const arms = [{cls: 'arm2', m: gompertz(1.9, 0.11)}, {cls: 'arm1', m: gompertz(3.1, 0.13)},
    {cls: 'ctrl', m: gompertz(4.4, 0.16)}];
  const groups = arms.map(({cls, m}) => {
    const s = (t: number) => 0.1 + 0.004 * t;
    const band = (z: number) => 'M' + steps.map((t) => `${xAt(t)},${yAt(m(t) * Math.exp(z * s(t)))}`).join(' L') +
      ' L' + steps.slice().reverse().map((t) => `${xAt(t)},${yAt(m(t) * Math.exp(-z * s(t)))}`).join(' L') + ' Z';
    const median = 'M' + steps.map((t) => `${xAt(t)},${yAt(m(t))}`).join(' L');
    return `<g class="${cls}"><path class="b90" d="${band(1.645)}"/>` +
      `<path class="med" pathLength="1" d="${median}"/></g>`;
  });
  const t = document.createElement('template');
  t.innerHTML = `<svg class="pmx-bg" viewBox="0 0 1440 440" preserveAspectRatio="xMaxYMin slice" ` +
    `aria-hidden="true">${groups.join('')}</svg>`;
  return t.content.firstElementChild as SVGSVGElement;
}

// ---- the stages and the edges between them -------------------------------------------------------------------------

/** An edge between two neighbouring cards, from the first's output port to the next one's input, named by the file
 *  that passes between them. */
function edge(chip: string): HTMLElement {
  const wire = 'M0,70 C38,70 34,30 72,30';
  const t = document.createElement('template');
  t.innerHTML = `<svg class="pmx-wire-h" width="72" height="100" viewBox="0 0 72 100"><path class="pmx-glow" ` +
    `d="${wire}"/><path class="pmx-wire" d="${wire}"/><circle class="pmx-port" cx="0" cy="70" r="4.5"/>` +
    `<circle class="pmx-port" cx="72" cy="30" r="4.5"/></svg>` +
    '<svg class="pmx-wire-v" width="10" height="64" viewBox="0 0 10 64"><path class="pmx-wire" d="M5,0 L5,64"/></svg>';
  return h('div', 'pmx-link', [t.content, h('span', 'pmx-chip', [chip])], {'aria-hidden': 'true'});
}

/** The edge under the cards, from Plan to Analyze: the plan Analyze reads the study against. */
function arc(chip: string): HTMLElement {
  const t = document.createElement('template');
  t.innerHTML = '<svg width="40" height="48" style="left:-2px"><path class="pmx-glow" d="M2,0 C2,28 14,38 40,38"/>' +
    '<path class="pmx-wire" d="M2,0 C2,28 14,38 40,38"/><circle class="pmx-port" cx="2" cy="0" r="4.5"/></svg>' +
    '<svg height="48" style="left:38px;width:calc(100% - 76px)"><line class="pmx-glow" x1="0" y1="38" x2="100%" ' +
    'y2="38"/><line class="pmx-wire" x1="0" y1="38" x2="100%" y2="38"/></svg>' +
    '<svg width="40" height="48" style="right:-2px"><path class="pmx-glow" d="M0,38 C26,38 38,28 38,0"/>' +
    '<path class="pmx-wire" d="M0,38 C26,38 38,28 38,0"/><circle class="pmx-port" cx="38" cy="0" r="4.5"/></svg>';
  const a = h('div', 'pmx-arc', [t.content, h('span', 'pmx-chip', [chip])], {'aria-hidden': 'true'});
  // from the first card's centre to the third's: three cards share the width left by the two 72px edges
  Object.assign(a.style, {marginLeft: 'calc((100% - 144px) / 6)', marginRight: 'calc((100% - 144px) / 6)'});
  return a;
}

interface Card {stage: Stage, n: string, name: string, what: string, icon: string}
const CARDS: Card[] = [
  {stage: 'plan', n: '01', name: 'Plan', what: 'Animals, duration, calendar and power, before the study runs.',
    icon: ICON.plan},
  {stage: 'monitor', n: '02', name: 'Monitor', what: 'A running study against its plan: new measurements, animals ' +
    'that need attention, planned reviews.', icon: ICON.monitor},
  {stage: 'analyze', n: '03', name: 'Analyze', what: 'A finished study as planned: each arm against the control, ' +
    'each animal against its untreated twin.', icon: ICON.analyze},
];

function card(c: Card, onOpen: () => void, onHover: (on: boolean) => void): HTMLButtonElement {
  const b = h('button', 'pmx-card pmx-rise', [
    h('span', 'pmx-card-top', [h('span', 'pmx-ico', [icon(c.icon)]), h('span', 'pmx-num', [c.n])]),
    h('span', 'pmx-name', [c.name + ' ', h('span', 'pmx-go', ['→'], {'aria-hidden': 'true'})]),
    h('span', 'pmx-what', [c.what]),
  ], {'type': 'button', 'data-stage': c.stage});
  b.addEventListener('click', onOpen);
  for (const [ev, on] of [['mouseenter', true], ['focus', true], ['mouseleave', false], ['blur', false]] as const)
    b.addEventListener(ev, () => onHover(on));
  return b;
}

// ---- the page ------------------------------------------------------------------------------------------------------

function hubContent(self: ViewHost): ViewContent {
  const open = (stage: Stage) => stage === 'plan' ? void openStage('plan') : openStudyDialog();
  // the card under the pointer lights the edges that touch it and dims the others
  const edges: {el: HTMLElement, ends: Stage[]}[] = [];
  const hover = (stage: Stage | null) => {
    for (const e of edges) {
      const on = stage !== null && e.ends.includes(stage);
      e.el.classList.toggle('is-on', on);
      e.el.classList.toggle('is-dim', stage !== null && !on);
    }
  };
  const cards = CARDS.map((c, i) => {
    const el = card(c, () => open(c.stage), (on) => hover(on ? c.stage : null));
    el.style.animationDelay = `${0.16 + i * 0.08}s`;
    return el;
  });
  const between = [edge('.pmxplan'), edge('.pmxstudy')];
  between.forEach((el, i) => edges.push({el, ends: [CARDS[i].stage, CARDS[i + 1].stage]}));
  const under = arc('.pmxplan');
  edges.push({el: under, ends: ['plan', 'analyze']});

  const flow = h('div', 'pmx-flow', [cards[0], between[0], cards[1], between[1], cards[2]]);
  flow.style.gridTemplateColumns = 'minmax(0,1fr) 72px minmax(0,1fr) 72px minmax(0,1fr)';

  const strip = h('button', 'pmx-strip pmx-rise', [
    h('span', 'pmx-ico', [icon(ICON.open, 20)]),
    h('span', '', [h('span', 'pmx-strip-title', ['Open a study']),
      h('span', 'pmx-strip-text', ['A plan, a followed study or a table of measurements. Or drop it anywhere.'])]),
  ], {type: 'button'});
  strip.style.animationDelay = '.44s';
  strip.addEventListener('click', () => openStudyDialog());

  const studies = yourStudies();
  studies.root.style.animationDelay = '.52s';
  // back on the Overview after saving a plan or a study elsewhere: the list is read again
  self.subscribe(grok.events.onCurrentViewChanged, () => {
    if (stageOf(grok.shell.v) === 'home')
      studies.reload();
  });

  const title = h('h1', 'pmx-hub-title pmx-rise', ['PMX Studio']);
  const lede = h('p', 'pmx-hub-lede pmx-rise', ['Xenograft and syngeneic efficacy studies, from design to analysis.']);
  lede.style.animationDelay = '.08s';
  const page = h('div', 'pmx-hub-page', [title, lede, flow, under, strip, studies.root]);

  // a file dragged anywhere on the page: the whole page takes it, and the dialog opens with it read
  const dropLayer = h('div', 'pmx-drop', [h('div', 'pmx-drop-frame'),
    h('div', 'pmx-drop-text', [h('span', 'pmx-ico', [icon(ICON.drop, 26, 1.7)]), 'Drop to open'])],
  {'aria-hidden': 'true'});
  const root = h('div', 'pmx-hub', [growthArt(), page, dropLayer]);
  takeDrops(root, (on) => dropLayer.classList.toggle('is-open', on), (files) => openStudyDialog(files));

  // no path of its own: opened as the app, the platform gives it the app's; opened later, showStage does
  return {root, ribbon: [[h('div', 'pmx-head', [steps('home', (s) => void openStage(s))])]]};
}

/** How far a study has gone: Plan — Monitor — Analyze, the stages behind it filled. */
function stageMarker(e: LibraryEntry): HTMLElement {
  const marker = h('span', 'pmx-stage-marker', [1, 2, 3].map((k) => h('span', k <= e.stage ? 'on' : '')));
  marker.title = e.stage === 1 ? 'Planned' : e.stage === 2 ? 'Followed in Monitor' : 'Finished: analyzed';
  return h('span', 'pmx-stage-cell', [marker, ['Plan', 'Monitor', 'Analyze'][e.stage - 1]]);
}

/** The study day today is, from its day 0 (null: day 0 not fixed): the calendar's own count. */
const todayOf = (e: LibraryEntry) => dayToday(e.start);
const shortDate = (iso: string) => new Date(`${iso}T00:00`).toLocaleDateString('en-GB',
  {weekday: 'short', day: 'numeric', month: 'short'});

/** How far a study is: a followed one, its last day measured of the planned end (a bar); a plan, when it starts. */
function progressOf(e: LibraryEntry): HTMLElement {
  if (e.kind === 'plan')
    return h('span', 'pmx-muted', [e.start ? `starts ${shortDate(e.start)}` : 'not started']);
  if (e.stage === 3)
    return h('span', 'pmx-muted', [`finished, day ${e.day}`]);
  const share = e.day !== null && e.end ? Math.min(1, e.day / e.end) : 0;
  const bar = h('span', 'pmx-progress', [h('i', '')]);
  (bar.firstChild as HTMLElement).style.width = `${Math.round(share * 100)}%`;
  return h('span', 'pmx-progress-cell', [bar, `day ${e.day ?? 0}${e.end ? ` of ~${e.end}` : ''}`]);
}

/** Today across the followed studies: how many, and which are due a measurement today (day 0 known). */
function todayLine(entries: LibraryEntry[]): string {
  const followed = entries.filter((e) => e.stage === 2);
  if (!followed.length)
    return '';
  const due = followed.filter((e) => todayOf(e) !== null && e.next === todayOf(e)).length;
  const measured = followed.filter((e) => todayOf(e) !== null && e.day === todayOf(e)).length;
  const date = new Date().toLocaleDateString('en-GB', {weekday: 'short', day: 'numeric', month: 'short'});
  return [`Today, ${date}: ${followed.length} ${followed.length === 1 ? 'study' : 'studies'} followed`,
    due ? `${due} to measure today` : '', measured ? `${measured} measured today` : ''].filter((s) => s).join(' · ');
}

/** The user's plans and studies, newest first, found by any word of their ID, name, line, compounds, owner or tags;
 *  a plan opens in Plan, a study in Monitor; while there are none, a line that says so. Built with the page (the view
 *  owns the table), filled once the folders are read; `reload` reads them again. */
function yourStudies(): {root: HTMLElement, reload: () => void} {
  const all = signal<LibraryEntry[]>([]);
  const query = signal('');
  const stage = signal('All');
  const inStage = (e: LibraryEntry, s: string) => s === 'All' || ['Planned', 'Followed', 'Finished'][e.stage - 1] === s;
  const items = computed(() => all.value.filter((e) => matches(e, query.value) && inStage(e, stage.value)));
  const host = h('div', 'pmx-studies pmx-rise');
  host.hidden = true;
  const none = h('div', 'pmx-muted', ['No study matches']);
  const showNone = () => none.hidden = items.peek().length > 0;
  const today = h('div', 'pmx-library-today');
  // the stage filters, each with its count
  const chips = h('div', 'pmx-filters');
  const drawChips = () => chips.replaceChildren(...['All', 'Planned', 'Followed', 'Finished'].map((s) => {
    const b = h('button', `pmx-filter${stage.peek() === s ? ' is-on' : ''}`,
      [`${s} · ${all.peek().filter((e) => inStage(e, s)).length}`]);
    b.setAttribute('aria-pressed', String(stage.peek() === s));
    b.onclick = () => {
      stage.value = s;
      drawChips();
      showNone();
    };
    return b;
  }));
  // the user's own as soon as they are read (nothing is computed for the list: the page opens without the engine)
  const show = (entries: LibraryEntry[]) => {
    all.value = entries;
    host.hidden = false;
    const empty = entries.length === 0;
    for (const el of [search.root, chips, table.root]) // (styles that set display would beat `hidden`)
      el.style.display = empty ? 'none' : '';
    none.textContent = empty ? 'No studies yet' : 'No study matches';
    today.textContent = todayLine(entries);
    today.hidden = !today.textContent;
    drawChips();
    showNone();
  };
  // (asked as the Overview is built and again as it becomes the current view: one reading)
  const reload = libraryReader(show);
  reload();
  // a row opens its file as it is now, not as it was when the list was read (the study may have been saved since),
  // and once: a second click while it opens would open it twice
  let opening = false;
  const openEntry = async (e: LibraryEntry): Promise<unknown> => {
    const text = e.file ? await grok.dapi.files.readAsText(e.file) : e.text;
    return e.kind === 'plan' ? showStage('plan', await planView(parsePlan(text))) : openStudyFile(text);
  };
  const search = ui.input.string('', {placeholder: 'Search studies',
    tooltipText: 'Any word of the ID, name, line, compound, owner or tags',
    onValueChanged: (v) => {
      query.value = v ?? '';
      showNone();
    }});
  const muted = (s: string) => h('span', 'pmx-muted', [s]);
  // fixed columns, so the table keeps its shape whatever the filter shows: a longer cell ends in an ellipsis and
  // shows its whole text on hover (u2-table-clip: the table titles it when cut)
  const clip = (c: Node | string, cls = '') => h('div', `u2-table-clip pmx-clip ${cls}`, [c]);
  const table = new BasicTable<LibraryEntry>({items, keyOf: (e) => e.file?.fullPath ?? e.id + e.name, fixedLayout: true,
    columns: [
      {header: 'Study', width: '22%', render: (e) => h('div', 'pmx-study-cell', [clip(e.name), clip(muted(e.id))])},
      {header: 'Compound and schedule', width: '20%', render: (e) => clip(e.treatment)},
      {header: 'Line', width: '7%', render: (e) => clip(e.line)},
      {header: 'Stage', width: '11%', render: (e) => clip(stageMarker(e))},
      {header: 'Progress', width: '15%', render: (e) => clip(progressOf(e))},
      {header: 'Owner', width: '8%', render: (e) => clip(muted(e.owner))},
      {header: 'Updated', width: '9%', render: (e) => clip(muted(e.updated ? new Date(e.updated)
        .toLocaleDateString('en-GB', {day: 'numeric', month: 'short', year: 'numeric'}) : ''))},
      {header: 'Tags', width: '8%', render: (e) => {
        const tags = clip(h('span', '', e.tags.map((t) => h('span', 'pmx-tag', [t]))));
        tags.title = e.tags.join(', ');
        return tags;
      }},
    ],
    onRowClick: (e) => {
      if (opening)
        return;
      opening = true;
      void openEntry(e).catch((err) => grok.shell.error(String(err?.message ?? err))).finally(() => opening = false);
    }});
  table.root.style.cursor = 'pointer';
  host.append(h('div', 'pmx-library-head', [h('h3', '', ['Your studies']), search.root, chips]), today,
    h('div', 'pmx-library', [table.root, none]));
  return {root: host, reload};
}

/** Takes the files dragged over `target` (and only files): `over` says when some are over it, `got` receives them.
 *  The platform's own drop handling (which would open them as tables) does not see them. */
function takeDrops(target: HTMLElement, over: (on: boolean) => void, got: (files: File[]) => void): void {
  let depth = 0;
  const hasFiles = (e: DragEvent) => !!e.dataTransfer && [...e.dataTransfer.types].includes('Files');
  target.addEventListener('dragenter', (e) => {
    if (!hasFiles(e))
      return;
    e.preventDefault();
    if (depth++ === 0)
      over(true);
  });
  target.addEventListener('dragover', (e) => {
    if (!hasFiles(e))
      return;
    e.preventDefault();
    e.stopPropagation();
    e.dataTransfer!.dropEffect = 'copy';
  });
  target.addEventListener('dragleave', (e) => {
    if (hasFiles(e) && --depth <= 0) {
      depth = 0;
      over(false);
    }
  });
  target.addEventListener('drop', (e) => {
    if (!hasFiles(e))
      return;
    e.preventDefault();
    e.stopPropagation();
    depth = 0;
    over(false);
    got([...e.dataTransfer!.files]);
  });
}

// ---- Open a study --------------------------------------------------------------------------------------------------

type Phase = 'drop' | 'loading' | 'loaded' | 'error';
interface Option {stage: OpenStage, title: string, desc: string, note: string, disabled: boolean}

/** The Open a study dialog: a file dropped on it, browsed for, or one of the user's plans; once read, what it holds
 *  and the stage it opens in (suggested from its dates). `first`: files dropped on the page, read at once. */
export function openStudyDialog(first?: File[]): void {
  let phase: Phase = 'drop';
  let brought: Brought = {};
  let chosen: OpenStage | null = null;
  let error = '';
  let loadingName = '';
  let latest = 0;

  const body = h('div', 'pmx-open');
  const dialog = ui.dialog({title: 'Open a study'});
  dialog.add(body);
  dialog.onOK(() => void go().catch((e) => grok.shell.error(String(e?.message ?? e))), {closeOnEnter: false});
  // the OK button says where the study opens, and waits until there is one
  const ok: HTMLButtonElement = dialog.getButton('OK') ?? document.createElement('button');

  const options = (): {question: string, list: Option[], hint: string} => {
    const {table, plan, study} = brought;
    if (study) {
      // a followed study is Monitor's file: nothing to choose
      return {question: 'A followed study', hint: '', list: [
        {stage: 'monitor', title: 'Open it', desc: 'Monitor', note: '', disabled: false}]};
    }
    if (table) {
      const dated = !!plan && table.lastDay !== null && plan.end !== null;
      return {question: 'Has the study finished?',
        hint: dated ? `Day ${table.lastDay} of about ${plan!.end}, from the dates alone.` : '',
        list: [
          {stage: 'monitor', title: 'Still running', desc: 'Monitor', disabled: !plan,
            note: !plan ? 'Add the plan it is run with to follow it' : ''},
          {stage: 'analyze', title: 'Finished', note: '', disabled: false,
            desc: plan ? 'Analyze, against its plan' : 'Analyze, against the historical controls'},
        ]};
    }
    return {question: 'Has the study started?', hint: '', list: [
      {stage: 'plan', title: 'Not started', desc: 'Plan', note: '', disabled: false},
      {stage: 'monitor', title: 'Starting now', desc: 'Monitor', note: '', disabled: false}]};
  };

  const selected = (): Option | null => {
    if (phase !== 'loaded')
      return null;
    const list = options().list;
    const pick = list.find((o) => o.stage === chosen && !o.disabled) ??
      list.find((o) => o.stage === suggestedStage(brought) && !o.disabled) ??
      list.find((o) => !o.disabled);
    return pick ?? null;
  };

  const dropZone = (small: boolean, text: string): HTMLElement => {
    const input = h('input', 'pmx-file-input', [], {type: 'file', multiple: '',
      accept: OPENED_EXTS.map((x) => '.' + x).join(',')});
    input.addEventListener('change', () => {
      void take([...input.files ?? []]);
      input.value = '';
    });
    const zone = h('label', 'pmx-dz' + (small ? ' pmx-dz-small' : ''), small ?
      [h('span', 'pmx-ico', [icon(ICON.add, 14, 2)]), text, input] :
      [h('span', 'pmx-ico', [icon(ICON.open)]),
        h('span', 'pmx-dz-main', ['Drag a file here, or ', h('b', '', ['browse'])]),
        h('span', 'pmx-exts', [OPENED_EXTS.filter((x) => x !== 'txt').map((x) => '.' + x).join(' · ')]), input]);
    return zone;
  };

  const render = () => {
    const parts: HTMLElement[] = [];
    if (phase === 'drop' || phase === 'error') {
      if (phase === 'error')
        parts.push(h('div', 'pmx-error', [icon(ICON.warn, 16, 1.8), h('span', '', [error])]));
      const saved = ui.link('Saved plans', () => void savedPlans(), 'Your plans saved in Plan');
      parts.push(dropZone(false, ''), h('div', 'pmx-alts', [saved]));
    } else if (phase === 'loading') {
      parts.push(h('div', 'pmx-loading', [
        h('div', 'pmx-file-row', [icon(ICON.file, 20), h('div', '', [h('b', '', [loadingName])])]),
        h('div', 'pmx-track', [h('div')])]));
    } else {
      const {table, plan, study} = brought;
      const rows = [...(table ? [[table.file, 'Measurements']] : []), ...(plan ? [[plan.file, 'Plan']] : []),
        ...(study ? [[study.file, 'Followed study']] : [])]
        .map(([name, kind]) => h('div', 'pmx-file-row', [icon(ICON.file, 20),
          h('div', '', [h('b', '', [name]), h('span', '', [kind])])]));
      parts.push(h('div', 'pmx-opts', rows));
      const holds = holdsOf(brought);
      if (holds.length)
        parts.push(h('ul', 'pmx-holds', holds));
      const {question, list, hint} = options();
      const pick = selected();
      parts.push(h('div', 'pmx-ask', [h('div', 'pmx-q', [question]), h('div', 'pmx-opts', list.map((o) => {
        const radio = h('input', '', [], {type: 'radio', name: 'pmx-stage'});
        radio.checked = pick?.stage === o.stage;
        radio.disabled = o.disabled;
        radio.addEventListener('change', () => {
          chosen = o.stage;
          render();
        });
        return h('label', 'pmx-opt' + (o.disabled ? ' is-off' : pick?.stage === o.stage ? ' is-on' : ''), [radio,
          h('span', '', [h('b', '', [o.title]), o.desc, ...(o.note ? [h('em', '', [o.note])] : [])])]);
      }))]));
      if (hint)
        parts.push(h('div', 'pmx-hint', [icon(ICON.info, 16, 1.8), h('span', '', [hint])]));
      if (error)
        parts.push(h('div', 'pmx-error', [icon(ICON.warn, 16, 1.8), h('span', '', [error])]));
      const missing = table && !plan ? 'Add its .pmxplan' : plan && !table && !study ? 'Add its measurements' : '';
      if (missing)
        parts.push(dropZone(true, missing));
    }
    body.replaceChildren(h('div', 'pmx-in', parts));
    const pick = selected();
    ok.textContent = pick ? `Open in ${pick.desc.split(',')[0]}` : 'Open';
    ok.disabled = !pick;
    ok.classList.toggle('disabled', !pick);
  };

  /** Reads the files and adds what they hold to what the dialog already has; a later read overtakes an earlier. */
  const take = async (files: File[]) => {
    if (files.length === 0)
      return;
    const mine = ++latest;
    phase = 'loading';
    loadingName = files.map((f) => f.name).join(', ');
    render();
    const started = Date.now();
    const reads = await Promise.all(files.map((f) => f.text().then((t) => readFile(f.name, t),
      () => readFile(f.name, ''))));
    // the bar shows for a moment even for a small file: the change of state reads as the file being taken in
    await new Promise((r) => setTimeout(r, Math.max(0, 450 - (Date.now() - started))));
    if (mine !== latest)
      return;
    const {brought: next, bad} = bring(brought, reads);
    error = bad.length ? `Not a study file: ${bad.join(', ')}` : '';
    brought = next;
    chosen = null;
    phase = next.table || next.plan || next.study ? 'loaded' : 'error';
    render();
  };

  const savedPlans = async () => {
    const files = await myPlans().catch(() => [] as DG.FileInfo[]);
    if (files.length === 0) {
      grok.shell.info('No saved plans yet');
      return;
    }
    DG.Menu.popup().items(files.map((f) => f.fileName), (name: string) => {
      const f = files.find((x) => x.fileName === name)!;
      void grok.dapi.files.readAsText(f).then((text) => take([new File([text], f.fileName)]));
    }).show();
  };

  const go = async () => {
    const pick = selected();
    const {table, plan, study} = brought;
    if (!pick)
      return;
    if (pick.stage === 'plan' && plan)
      showStage('plan', await planView(parsePlan(plan.text)));
    else if (pick.stage === 'analyze' && table) {
      if (plan)
        table.df.setTag(PLAN_TAG, plan.text);
      showStage('analyze', await analyzeView({study: table.df}));
    } else if (pick.stage === 'monitor') {
      // a followed study as saved; a table under way with its plan; a plan starting now, with an empty table
      const df = study ? tableOf(parseStudy(study.text)) : table?.df ?? (plan ? newStudyTable(parsePlan(plan.text)) :
        null);
      if (!df)
        return;
      if (table && plan)
        df.setTag(PLAN_TAG, plan.text);
      showStage('monitor', await monitorView(df));
    }
  };

  takeDrops(dialog.root, (on) => body.querySelectorAll('.pmx-dz').forEach((z) => z.classList.toggle('is-over', on)),
    (files) => void take(files));
  render();
  dialog.show();
  if (first?.length)
    void take(first);
}

/** What the files brought hold, one line each. */
function holdsOf(b: Brought): HTMLElement[] {
  const line = (kind: 'ok' | 'none' | 'warn', text: string, sub = '') => h('li', '', [
    h('span', kind, [icon(kind === 'ok' ? ICON.ok : kind === 'warn' ? ICON.warn : ICON.none, 16, 2)]),
    h('span', '', [text, ...(sub ? [h('small', '', [sub])] : [])])]);
  const out: HTMLElement[] = [];
  if (b.table) {
    const t = b.table;
    const facts = [t.animals !== null ? `${t.animals} animals` : '', t.groups !== null ? `${t.groups} groups` : '',
      t.lastDay !== null ? `day 0–${t.lastDay}` : ''].filter(Boolean);
    const what = facts.join(' · ') || 'Measurements';
    out.push(t.missing.length ? line('warn', what, `Not found: ${t.missing.join(', ')}. ` +
      'You pick it next.') : line('ok', what, t.found.join(' · ')));
    out.push(t.weight ? line('ok', 'Body weight') : line('none', 'No body weight'));
    out.push(b.plan ? line('ok', `Plan: ${b.plan.name}`) : line('none', 'No plan'));
  } else if (b.plan && !b.study) {
    const facts = [b.plan.line, b.plan.end !== null ? `about ${b.plan.end} days` : ''].filter(Boolean);
    out.push(line('ok', b.plan.name, facts.join(' · ')));
  }
  return out;
}
