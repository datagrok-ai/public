/* PMX Studio's home: the stages as cards joined the way Flow joins its steps (each edge names the file that passes
   between two stages), the growth of a study behind them, the user's studies, and the way in for a study under way:
   a file dropped anywhere on the page, or picked in the Open a study dialog, is read and opened in its stage. */
import {maxOf} from './format';
import '../css/hub.css';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable, computed, effect, signal} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {STAGES, isoToday, markLine, markStage, modal, modalOf, noteLine, plainError, preparing, releaseModal,
  setUserSetting, stageOf, userSetting, downloadExample, EXAMPLE_SOURCE, showHelp, objectPanel, okNamed} from './ui';
import {PmxSavedStudy} from './handlers';
import {studySummary as savedSummary} from './file-preview';
import {esc} from './sheet';
import type {Stage} from './ui';
import {OPENED_AS, bringUp, openStage, showStage} from './stages';
import {_package} from './package';
import {cruxWarm} from './crux';
import {planView} from './plan-view';
import {OPTIONAL_COLUMNS, analyzeView, keptEffect} from './analyze-view';
import {combinationSteps, notCombined} from './synergy';
import {followPlan, monitorView, openStudyFile} from './monitor-view';
import {EXPOSURE_TAG, FINISHED_TAG, LINE_TAG, newStudyTable, parseStudy, roleColumns, saveStudy,
  studyOf as savedStudyOf, studyNameOf, tableOf, type SavedStudy} from './study-file';
import {concentrationMapping, exposureJson, ncaByArm, ncaInstalled, runNca, type NcaByArm} from './exposure';
import {analysisDeclarationOf, analysisFrames} from './analyze';
import {withoutExcluded} from './study-table';
import {PLAN_TAG, parsePlan, planTag} from './plan-file';
import {domainStorageOn, entryText, firstOpen, firstOpenNotice, importFolders, keptRows, keptRowsLine, list,
  markFromFile, markTableFromFile, myEntries, open, ownerOf, remove, restore, save, savedPlans, saveError, share}
  from './store';
import {emptyMeta, metaOf, setMetaOn} from './study-meta';
import {linePrior, tumorModels} from './priors';

/** The Overview table's column widths in percent of the room after the tick's column (PICK_WIDTH), adding up to 100
 *  (fixed layout: more would push the last column out of view). With domain storage, a last column for Share and
 *  Delete (Restore among the deleted). */
/** The Overview's tick column: the checkbox and its padding, whatever the page's width. */
export const PICK_WIDTH = '32px';
export const libraryWidths = (rows: boolean) => ({study: 19, drugs: 10, treatment: rows ? 14 : 15, model: 8,
  // (the action column holds Restore among the deleted: at 5 it read "Rest…")
  stage: rows ? 10 : 11, progress: rows ? 11 : 13, owner: 8, updated: 9, tags: rows ? 4 : 7,
  ...rows ? {action: 7} : {}});

/** Shares a plan or study of the user's, to read, with the groups picked (a user is a group of one). */
function shareDialog(e: LibraryEntry): void {
  const who = ui.input.userGroups('Share with');
  const d = modal('share', `Share "${e.name}"`);
  if (!d)
    return;
  d.add(who).add(ui.divText('They can open it and save a copy of their own; only you ' +
    'can change it.', 'pmx-muted')).onOK(() => {
    const groups = who.value ?? [];
    if (!groups.length)
      return;
    void share(e.row!.kind, e.row!.id, groups).then(
      // (a toast reads HTML: the name is the user's, shared with others, and never markup)
      () => grok.shell.info(`"${esc(e.name)}" is shared to read with ` +
        esc(groups.map((g) => g.friendlyName).join(', '))),
      (err) => grok.shell.error(esc(plainError(err, `"${e.name}" could not be shared just now. Try again.`))));
  }).show();
}

/** The user setting that says the user's files were brought in (domain storage). */
const FILES_IN = 'filesBroughtIn';
import {LibraryEntry, entryOf, libraryReader, matches} from './library';
import {compareTexts, compareView} from './compare-view';
import {checklistDialog} from './arrive';
import {CTD_TITLE, M4S_URL, commonDrug, ctdBlocks, ctdOfSaved, m4sOrder, summaryItems, type CtdStudy} from './ctd';
import {DOCX_TYPE, docx} from './docx';
import {ToolLink, download, toolbox} from './report';

/** A study Compare can read: one followed, with measurements. */
const comparable = (e: LibraryEntry) => e.kind === 'study' && e.day !== null;
import {dayToday} from './calendar';
import {Brought, MeasuredTable, OPENED_EXTS, OpenStage, bring, fileText, readFile, routedNote, suggestedStage}
  from './open-study';
import {Question, studySummary, waiting} from './study-reader';
import {plural} from './wide-table';
import {mergedWarning} from './xlsx-raw';
import {APP_TABLE_TAG} from './native';

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
  play: '<path d="M8 5.5v13l10.5-6.5z"/>',
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
function edge(): HTMLElement {
  const wire = 'M0,70 C38,70 34,30 72,30';
  const t = document.createElement('template');
  t.innerHTML = `<svg class="pmx-wire-h" width="72" height="100" viewBox="0 0 72 100"><path class="pmx-glow" ` +
    `d="${wire}"/><path class="pmx-wire" d="${wire}"/><circle class="pmx-port" cx="0" cy="70" r="4.5"/>` +
    `<circle class="pmx-port" cx="72" cy="30" r="4.5"/></svg>` +
    '<svg class="pmx-wire-v" width="10" height="64" viewBox="0 0 10 64"><path class="pmx-wire" d="M5,0 L5,64"/></svg>';
  return h('div', 'pmx-link', [t.content], {'aria-hidden': 'true'});
}

/** The edge under the cards, from Plan to Analyze: the plan Analyze reads the study against. */
function arc(): HTMLElement {
  const t = document.createElement('template');
  t.innerHTML = '<svg width="40" height="48" style="left:-2px"><path class="pmx-glow" d="M2,0 C2,28 14,38 40,38"/>' +
    '<path class="pmx-wire" d="M2,0 C2,28 14,38 40,38"/><circle class="pmx-port" cx="2" cy="0" r="4.5"/></svg>' +
    '<svg height="48" style="left:38px;width:calc(100% - 76px)"><line class="pmx-glow" x1="0" y1="38" x2="100%" ' +
    'y2="38"/><line class="pmx-wire" x1="0" y1="38" x2="100%" y2="38"/></svg>' +
    '<svg width="40" height="48" style="right:-2px"><path class="pmx-glow" d="M0,38 C26,38 38,28 38,0"/>' +
    '<path class="pmx-wire" d="M0,38 C26,38 38,28 38,0"/><circle class="pmx-port" cx="38" cy="0" r="4.5"/></svg>';
  const a = h('div', 'pmx-arc', [t.content], {'aria-hidden': 'true'});
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
    'each animal against what it would measure untreated.', icon: ICON.analyze},
];

/** A secondary way in, at a card's foot or under the empty list: ▷, what it does, and a line on what it opens. */
function demoButton(title: string, what: string, run: () => void): HTMLButtonElement {
  const b = h('button', 'pmx-demo', [h('span', 'pmx-demo-ico', [icon(ICON.play, 12, 1.8)]),
    h('span', '', [h('span', 'pmx-demo-title', [title]), h('span', 'pmx-demo-what', [what])])], {type: 'button'});
  b.addEventListener('click', run);
  return b;
}

/** A stage's card: the upper area one button that opens the stage, its demo a second button at the foot (no button
 *  inside a button: each is reached and read on its own). */
function card(c: Card, onOpen: () => void, onDemo: () => void, onHover: (on: boolean) => void): HTMLElement {
  const main = h('button', 'pmx-card-main', [
    h('span', 'pmx-card-top', [h('span', 'pmx-ico', [icon(c.icon)]), h('span', 'pmx-num', [c.n])]),
    h('span', 'pmx-name', [c.name + ' ', h('span', 'pmx-go', ['→'], {'aria-hidden': 'true'})]),
    h('span', 'pmx-what', [c.what]),
  ], {'type': 'button'});
  main.addEventListener('click', onOpen);
  const demo = demoButton('Try with demo data', DEMO[c.stage as OpenStage].what, onDemo);
  demo.setAttribute('aria-label', `Try ${c.name} with demo data: ${DEMO[c.stage as OpenStage].what}`);
  const el = h('div', 'pmx-card pmx-rise', [main, h('div', 'pmx-card-foot', [demo])], {'data-stage': c.stage});
  for (const [ev, on] of [['mouseenter', true], ['focusin', true], ['mouseleave', false], ['focusout', false]] as const)
    el.addEventListener(ev, () => onHover(on));
  // a stage about to be opened: the engine made ready, if the idle Overview has not done it yet
  el.addEventListener('mouseenter', cruxWarm, {once: true});
  return el;
}

// ---- the page ------------------------------------------------------------------------------------------------------

function hubContent(self: ViewHost): ViewContent {
  const open = (stage: Stage) => stage === 'plan' ? void openStage('plan') :
    void openStudyDialog(undefined, {stage: stage as OpenStage});
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
    const el = card(c, () => open(c.stage), () => void openDemo(c.stage as OpenStage),
      (on) => hover(on ? c.stage : null));
    el.style.animationDelay = `${0.16 + i * 0.08}s`;
    return el;
  });
  const between = [edge(), edge()];
  between.forEach((el, i) => edges.push({el, ends: [CARDS[i].stage, CARDS[i + 1].stage]}));
  const under = arc();
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

  // the logo is decoration beside the name it would repeat: empty alt; width and height set, so the title keeps its
  // size before the image arrives
  const logo = h('img', 'pmx-hub-logo', [], {src: `${_package.webRoot.replace(/\/?$/, '/')}img/pmx-studio-logo.svg`,
    alt: '', width: '32', height: '32', draggable: 'false',
    style: 'margin-right: 12px; vertical-align: -5px'}); // 32 px inside the title's 41 px line
  const title = h('h1', 'pmx-hub-title pmx-rise', [logo, 'PMX Studio']);
  const lede = h('p', 'pmx-hub-lede pmx-rise', ['Xenograft and syngeneic efficacy studies, from design to analysis.']);
  lede.style.animationDelay = '.08s';
  const page = h('div', 'pmx-hub-page', [title, lede, flow, under, strip, studies.root]);

  // a file dragged anywhere on the page: the whole page takes it, and the dialog opens with it read
  const layer = dropLayer();
  const root = h('div', 'pmx-hub', [growthArt(), page, layer]);
  takeDrops(root, (on) => layer.classList.toggle('is-open', on), (files) => openStudyDialog(files));
  // Ctrl+O (Cmd+O) here too, as in every stage: the Open dialog (QA 0.1.3 round 3 chaos 20)
  const onKey = (e: KeyboardEvent) => {
    if ((e.ctrlKey || e.metaKey) && !e.altKey && !e.shiftKey && e.key.toLowerCase() === 'o' &&
      !!grok.shell.v?.root?.contains(root)) {
      e.preventDefault();
      e.stopPropagation();
      openStudyDialog();
    }
  };
  document.addEventListener('keydown', onKey, true);
  self.own(() => document.removeEventListener('keydown', onKey, true));
  // the Overview drawn, the browser idle: the engine made ready for the stage the user opens next (once a session)
  const idle = (window as any).requestIdleCallback ?? ((f: () => void) => setTimeout(f, 1500));
  idle(() => cruxWarm(), {timeout: 5000});

  // no path of its own: opened as the app, the platform gives it the app's; opened later, showStage does
  return {root, ribbon: [],
    toolbox: hubToolbox(() => void compareSaved(studies.root))};
}

/** The studies ticked in the Overview's list, by entryKey: kept while the app is open (the Overview opened again
 *  shows them ticked), and ticked in Analyze's Compare with other studies. */
let pickedInOverview: ReadonlySet<string> = new Set();
const entryKey = (e: LibraryEntry) => e.file?.fullPath ?? e.row?.id ?? e.id + e.name;

/** Compare with the user's studies that have measurements (those ticked in the Overview when two or more are), when
 *  there are two or more; else the list (`studies`) shown, with what Compare needs. */
async function compareSaved(studies?: HTMLElement): Promise<void> {
  try {
    const all = (await myEntries()).filter(comparable);
    const ticked = all.filter((e) => pickedInOverview.has(entryKey(e)));
    const eligible = ticked.length >= 2 ? ticked : all;
    if (eligible.length >= 2) {
      await showStage('home', await compareView(eligible));
      return;
    }
  } catch (e) {
    grok.shell.error(plainError(e, 'Your studies could not be read. Try again in a moment.'));
    return;
  }
  studies?.scrollIntoView({behavior: 'smooth', block: 'start'});
  grok.shell.info('Compare needs two studies with measurements or more: save them from Monitor or Analyze, or ' +
    'try the demo studies.');
}

/** Compare opened on `es` (two or more), in the app's view: Back comes back here. */
function openCompare(es: LibraryEntry[]): Promise<void> {
  return compareView(es).then((c) => showStage('home', c)).then(() => {},
    (e) => grok.shell.error(plainError(e, 'Compare could not be opened.')));
}

/** Analyze's Compare with other studies: the user's studies with measurements, the one in Analyze (its saved row
 *  `rowId`) and those ticked in the Overview ticked; Compare opens on two or more. */
export async function compareWith(rowId: string | undefined): Promise<void> {
  let all: LibraryEntry[];
  try {
    all = (await myEntries()).filter(comparable);
  } catch (e) {
    grok.shell.error(plainError(e, 'Your studies could not be read. Try again in a moment.'));
    return;
  }
  const here = all.find((e) => !!rowId && e.row?.id === rowId);
  if (all.length < 2) {
    grok.shell.info('Compare needs two saved studies with measurements or more: save this one, and another, first.');
    return;
  }
  const d = modal('compare-with', 'Compare with other studies');
  if (!d)
    return;
  const boxes = all.map((e): [LibraryEntry, DG.InputBase<boolean>] => [e, ui.input.bool(
    `${e.name}${e.id ? ` · ${e.id}` : ''}${e === here ? ' (this study)' : ''}`,
    {value: e === here || pickedInOverview.has(entryKey(e)), onValueChanged: () => allow()})]);
  const picks = () => boxes.filter(([, b]) => b.value).map(([e]) => e);
  // (the ticks here are this comparison's: the Overview's own stay as the user left them, QA 0.1.3 analyst #10)
  // why Compare waits, said under the list (a disabled button gave no reason: QA 0.1.3 round 3, chaos 23)
  const why = h('div', 'pmx-muted', []);
  d.add(ui.divV(boxes.map(([, b]) => b.root), {style: {maxHeight: '60vh', overflow: 'auto'}}))
    .add(why)
    .onOK(() => void openCompare(picks()));
  // (the OK button exists once onOK is set)
  const ok = okNamed(d, 'Compare');
  const allow = () => {
    const n = picks().length;
    ok.disabled = n < 2;
    why.textContent = n >= 2 ? '' : (here ? '' : 'This study is not in the list: save it (File › Save) to compare ' +
      'it. ') + (n === 1 ? 'Tick one more study to compare.' : 'Tick two studies or more to compare.');
    ok.title = why.textContent;
  };
  allow();
  d.show();
}

/** Brings the open view opened for `key` (OPENED_AS: a row id, a file path, the demo Compare) to the front; false when
 *  none is open. (In a click: a switch made at once was undone by the platform, and the row looked dead.) */
export function showOpened(key: string): boolean {
  const v = [...grok.shell.views].find((x) => x.root.dataset[OPENED_AS] === key);
  if (v)
    bringUp(v);
  return !!v;
}

/** A row of the Overview, or a study picked in Ctrl+O's Saved studies, opened: one function (Ctrl+O's read it as a
 *  file, "Not saved yet · copy of …" with no ?study=, QA 0.1.3 browser check). Open already (a double click, or
 *  opened before): that view comes to the front. `stage`: a study's, picked (else Analyze when finished, Monitor). */
export async function openEntryOnce(e: LibraryEntry, stage?: 'monitor' | 'analyze'): Promise<void> {
  if (showOpened(entryKey(e)))
    return;
  await openedAs(entryKey(e), () => openEntry(e, stage)).catch((err) => grok.shell.error(esc(plainError(err,
    `"${e.name}" could not be opened. Try again; if it fails again, reload the page.`))));
}

/** An entry as it is now, not as it was when the list was read (the study may have been saved since): a study row
 *  with its row, so that its Save updates it; a plan row remembered as it is read. */
async function openEntry(e: LibraryEntry, stage?: 'monitor' | 'analyze'): Promise<unknown> {
  if (e.row?.kind === 'study') {
    const {doc, saved} = await open('study', e.row.id);
    return openStudyFile(JSON.stringify(doc), saved, stage);
  }
  const text = e.row ? await entryText(e) : e.file ? await grok.dapi.files.readAsText(e.file).catch((err) => {
    console.warn('PMX Studio:', err);
    throw new Error(`"${e.name}" could not be read from your files. Check your connection and try again.`);
  }) : e.text;
  if (e.kind !== 'plan')
    return openStudyFile(text, null, stage);
  // (a file, not a row: its first save is the user's copy when its owner is someone else)
  const p = parsePlan(text);
  if (!e.row)
    markFromFile(p, ownerOf('plan', p));
  return showStage('plan', await planView(p));
}

/** Runs `open`, and marks the app's views it opened as opened for `key`. */
async function openedAs(key: string, open: () => Promise<unknown>): Promise<void> {
  const had = new Set([...grok.shell.views].map((v) => v.root));
  await open();
  for (const v of grok.shell.views) {
    if (!had.has(v.root) && stageOf(v))
      v.root.dataset[OPENED_AS] = key;
  }
}

/** The app's help, in the platform's help panel (as F1 shows it for every view of the app). */
const HELP_LINK: ToolLink = ['PMX Studio help', showHelp,
  'What each stage does, the files it opens and saves'];

/** The Overview's toolbox: the page's ways in, as links (PMX Studio, the demo data, help). */
function hubToolbox(compareSaved: () => void): HTMLElement {
  return toolbox([
    ['PMX Studio', [
      ['Open a study…', () => void openStudyDialog(), 'A plan, a followed study or a table of measurements'],
      ['New plan', () => void planView().then((v) => showStage('plan', v)), 'Plan a new study'],
      ['Compare saved studies', compareSaved, 'Your studies with measurements, side by side']]],
    ['Demo', [
      ['Plan: ' + DEMO.plan.what, () => void openDemo('plan'), 'Opens the demo plan in Plan'],
      ['Monitor: ' + DEMO.monitor.what, () => void openDemo('monitor'), 'Opens the demo study in Monitor'],
      ['Analyze: ' + DEMO.analyze.what, () => void openDemo('analyze'), 'Opens the finished demo study in Analyze'],
      ['Compare two demo studies', () => void compareDemo(), 'CT26 and MC38 trametinib, saved for you and opened ' +
        'in Compare']]],
    ['Help', [HELP_LINK]]]);
}

/** A stage's toolbox links of its own: in Study, Open… (for this stage) and the study's details; in Help, its demo,
 *  the example study (Monitor, Analyze) and the app's help. */
export function stageLinks(at: OpenAt, details: () => void): {study: ToolLink[], help: ToolLink[]} {
  return {study: [['Open…', () => void openStudyDialog(undefined, at), `Open a study in ${labelOf(at.stage)}. ` +
    `${TAKES[at.stage]}.`], ['Study details…', details, 'Its ID, dates, people, groups and notes']],
  help: [[`Demo: ${DEMO[at.stage].what}`, () => void openDemo(at.stage), `Opens the demo in ${labelOf(at.stage)}`],
    ...at.stage === 'plan' ? [] : [['Download an example study', downloadExample, EXAMPLE_SOURCE] as ToolLink],
    HELP_LINK]};
}

/** Compare's help links: its demo and the app's help. */
export const compareHelp = (): ToolLink[] => [['Compare two demo studies', () => void compareDemo(),
  'CT26 and MC38 trametinib, saved for you and opened in Compare'], HELP_LINK];

/** What shows over a page while files are dragged over it. */
function dropLayer(): HTMLElement {
  return h('div', 'pmx-drop', [h('div', 'pmx-drop-frame'),
    h('div', 'pmx-drop-text', [h('span', 'pmx-ico', [icon(ICON.drop, 26, 1.7)]), 'Drop to open'])],
  {'aria-hidden': 'true'});
}

/** A saved study in the context panel: what it is and how far it got, its answer only when the study kept its
 *  analysis (read from what it saved: no analysis is run for the panel), its model, groups and animals, progress and
 *  date (not its owner's login). `text` reads its saved record. */
export function savedStudyPanel(e: LibraryEntry, text: () => Promise<string>): HTMLElement {
  const updated = e.updated ? new Date(e.updated).toLocaleDateString('en-GB', {day: 'numeric', month: 'short',
    year: 'numeric'}) : '';
  const build = (s: ReturnType<typeof parseStudy> | null) => {
    const ctd = s ? ctdOfSaved(s) : null;
    const tgi = s?.exposure?.effect?.tgi ?? {};
    const treated = ctd ? ctd.groups.filter((g) => !g.control) : [];
    // (a group that grew more than the control is never the best: then none slows growth)
    const best = Object.entries(tgi).filter(([, x]) => Number.isFinite(x) && x > 0).sort((a, b) => b[1] - a[1])[0];
    const slow = treated.filter((g) => g.activity === 'Slows growth').length;
    // the groups as Analyze counts them (its kept analysis), else the named groups of its details, each once
    const groups = ctd?.groups.length || new Set((s?.meta?.groups ?? []).map((g) => g.group.trim()).filter((g) => g))
      .size;
    const animals = ctd ? ctd.groups.reduce((t, g) => t + (g.n || 0), 0) : 0;
    const sub = treated.length ? `${slow} of ${treated.length} treated groups slow growth` : undefined;
    return objectPanel({title: e.name, subtitle: `${e.id ? `${e.id} · ` : ''}${e.stage === 3 ? 'Finished' :
      e.stage === 1 ? 'Planned' : `Running, day ${e.day ?? 0}`}`,
    verdict: best ? {text: `Best group: ${best[0]}, TGI ${Math.round(best[1])}%`, sub,
      tone: slow > 0 ? 'ok' : 'neutral'} :
      Object.keys(tgi).length ? {text: 'No group slows growth', sub, tone: 'warn'} : null,
    facts: [{label: 'Model', value: e.model || e.line || ctd?.line || '—'},
      ...groups ? [{label: 'Groups', value: `${groups}${animals ? `, ${animals} animals` : ''}`}] : [],
      // (never "day 21 of 14": a finished study, or one run past its planned end, says so)
      {label: 'Progress', value: e.stage === 1 ? 'Not started' : e.stage === 3 ? `Day ${e.day ?? 0}, finished` :
        !e.end ? `Day ${e.day ?? 0}` : (e.day ?? 0) > e.end ? `Day ${e.day}, past the planned ${e.end}` :
          `Day ${e.day ?? 0} of ${e.end}`},
      ...updated ? [{label: 'Updated', value: updated}] : []]});
  };
  const el = ui.div([build(null)]);
  void text().then((t) => el.replaceChildren(build(parseStudy(t))), (err) => console.warn('PMX Studio:', err));
  return el;
}

/** How far a study has gone: Plan — Monitor — Analyze, the stages behind it filled. */
function stageMarker(e: LibraryEntry): HTMLElement {
  const marker = h('span', 'pmx-stage-marker', [1, 2, 3].map((k) => h('span', k <= e.stage ? 'on' : '')));
  marker.title = e.stage === 1 ? 'Planned' : e.stage === 2 ? 'Followed in Monitor' : 'Finished: analyzed';
  return h('span', 'pmx-stage-cell', [marker, ['Plan', 'Monitor', 'Analyze'][e.stage - 1]]);
}

/** The study day today is, from its day 0 (null: day 0 not fixed): the calendar's own count. */
const todayOf = (e: LibraryEntry) => dayToday(e.start);
const shortDate = (d: Date) => d.toLocaleDateString('en-GB', {weekday: 'short', day: 'numeric', month: 'short'});

/** How far a study is: a followed one, its last day measured of the planned end (a bar); a plan, when it starts. */
function progressOf(e: LibraryEntry): HTMLElement {
  if (e.kind === 'plan')
    return h('span', 'pmx-muted', [e.start ? `starts ${shortDate(new Date(`${e.start}T00:00`))}` : 'not started']);
  if (e.stage === 3)
    return h('span', 'pmx-muted', [`finished, day ${e.day}`]);
  const share = e.day !== null && e.end ? Math.min(1, e.day / e.end) : 0;
  const bar = h('span', 'pmx-progress', [h('i', '')]);
  (bar.firstChild as HTMLElement).style.width = `${Math.round(share * 100)}%`;
  // past the planned end: said so, not "day 21 of ~14"
  const of = !e.end ? '' : (e.day ?? 0) > e.end ? `, past the planned ${e.end}` : ` of ~${e.end}`;
  return h('span', 'pmx-progress-cell', [bar, `day ${e.day ?? 0}${of}`]);
}

/** Today across the followed studies: how many, and which are due a measurement today (day 0 known). */
function todayLine(entries: LibraryEntry[]): string {
  const followed = entries.filter((e) => e.stage === 2);
  if (!followed.length)
    return '';
  const due = followed.filter((e) => todayOf(e) !== null && e.next === todayOf(e)).length;
  const measured = followed.filter((e) => todayOf(e) !== null && e.day === todayOf(e)).length;
  return [`Today, ${shortDate(new Date())}: ${followed.length} ${followed.length === 1 ? 'study' : 'studies'} followed`,
    due ? `${due} to measure today` : '', measured ? `${measured} measured today` : ''].filter((s) => s).join(' · ');
}

/** A drug's regulatory summary from the picked studies: each one's rows computed from its data as Analyze computes
 *  them (ctdOfStudy, in M4S order); the dialog says which are in, with Download (.docx) beside Copy. */
async function ctdSummary(entries: LibraryEntry[], drug: string): Promise<void> {
  const {texts, skipped} = await compareTexts(entries);
  const progress = DG.TaskBarProgressIndicator.create(`Analyzing ${texts.size === 1 ? '1 study' :
    `${texts.size} studies`}...`);
  const kept = await Promise.all([...texts].map(async ([e, t]) => ({e, s: await ctdOfStudy(parseStudy(t))})))
    .finally(() => progress.close());
  const included = m4sOrder(kept.flatMap((x) => x.s ? [x.s] : []));
  const left = [...kept.filter((x) => !x.s).map((x) => x.e.name || x.e.id), ...skipped];
  const sub = `${drug} · ${included.length} stud${included.length === 1 ? 'y' : 'ies'} · prepared ${isoToday()}`;
  checklistDialog({title: `CTD summary for ${drug}`, guide: 'ICH M4S(R2)', url: M4S_URL, asker: 'The summary',
    download: !included.length ? undefined : async () => download(`${drug} - CTD summary.docx`,
      await docx(ctdBlocks(included, {sub, drug}), {title: CTD_TITLE, landscape: true}), DOCX_TYPE)},
  summaryItems(included, left));
}

/** The user's plans and studies, newest first, found by any word of their ID, name, line, compounds, owner or tags;
 *  a plan opens in Plan, a study in Monitor (a finished one in Analyze); while there are none, a line that says so.
 *  Built with the page (the view owns the table), filled once the folders are read; `reload` reads them again. */
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
  // no studies yet: two demo studies to compare, saved as the user's own
  const tryCompare = demoButton('Or compare two demo studies', 'CT26 and MC38 trametinib, saved for you and opened ' +
    'in Compare', () => void compareDemo());
  tryCompare.hidden = true;
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
    none.textContent = empty ? trash ? 'Nothing deleted' : 'No studies yet' : 'No study matches';
    tryCompare.hidden = !empty || trash;
    // (nothing is compared from the trash: no hint to pick)
    compareHint.hidden = trash || chosen.peek().length >= 2;
    today.textContent = todayLine(entries);
    today.hidden = !today.textContent;
    drawChips();
    showNone();
  };
  // (asked as the Overview is built and again as it becomes the current view: one reading)
  // (with domain storage: the rows the user can see, theirs and those shared with them)
  // (and, asked for, those in the trash: a click brings one back)
  let trash = false;
  // the first time with rows on, the user's files are brought in before the list is read (only read, once: the
  // list shows the studies they had); if that fails, the Bring your files in link stays to try again
  let bringLink: HTMLElement | null = null;
  let firstRead = domainStorageOn() && !userSetting(FILES_IN) ? firstOpen().then((r) => {
    setUserSetting(FILES_IN, 'yes');
    const {info, warning} = firstOpenNotice(r.report);
    // (they name the user's files: a toast reads HTML)
    if (info)
      grok.shell.info(esc(info));
    if (warning)
      grok.shell.warning(esc(warning));
    return r.entries;
  }, () => {
    if (bringLink)
      bringLink.hidden = false;
    return null;
  }) : null;
  const rows = async () => {
    const first = firstRead ? await firstRead : null;
    firstRead = null;
    return first ?? list('', 1000, undefined, trash);
  };
  const reload = domainStorageOn() ? () => void rows().then(show, () =>
    grok.shell.warning('Your studies could not be listed just now: open the Overview again in a moment.')) :
    libraryReader(show);
  reload();
  // a row opens once (openEntryOnce): a second click while it opens would open it twice
  let opening = false;
  const search = ui.input.string('', {placeholder: 'Search studies',
    tooltipText: 'Any word of the ID, name, line, compound, owner, run by or tags',
    onValueChanged: (v) => {
      query.value = v ?? '';
      showNone();
    }});
  const muted = (s: string) => h('span', 'pmx-muted', [s]);
  // fixed columns, so the table keeps its shape whatever the filter shows: a longer cell ends in an ellipsis and
  // shows its whole text on hover (u2-table-clip: the table titles it when cut)
  const clip = (c: Node | string, cls = '') => h('div', `u2-table-clip pmx-clip ${cls}`, [c]);
  // studies picked to compare: those with measurements
  const keyOf = entryKey;
  // (the picks made before, kept: back from Compare they are still ticked)
  const picked = signal<Set<string>>(new Set(pickedInOverview));
  const openRow = (e: LibraryEntry) => {
    // (a deleted one does not open: Restore brings it back first)
    if (opening || trash)
      return;
    opening = true;
    void openEntryOnce(e).finally(() => opening = false);
  };
  const pick = (e: LibraryEntry) => {
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.disabled = !comparable(e);
    box.checked = picked.peek().has(keyOf(e));
    box.title = comparable(e) ? 'Pick to compare' : 'Only studies with measurements can be compared';
    box.setAttribute('aria-label', `Compare ${e.name}`);
    box.onclick = (ev) => ev.stopPropagation();
    box.onchange = () => {
      const next = new Set(picked.peek());
      if (box.checked)
        next.add(keyOf(e));
      else
        next.delete(keyOf(e));
      picked.value = pickedInOverview = next;
      // the study last ticked in the context panel, shown there when the user has it open: a tick does not open it
      // (QA 0.1.3: opening it narrowed the page, the cards above re-wrapped and the list jumped under the pointer);
      // a click on its row still opens it; unticked, it leaves
      if (box.checked) {
        grok.shell.o = new PmxSavedStudy(e.name, {key: keyOf(e), stage: ['Plan', 'Monitor', 'Analyze'][e.stage - 1],
          details: () => savedStudyPanel(e, () => e.text ? Promise.resolve(e.text) : e.file ?
            grok.dapi.files.readAsText(e.file) : entryText(e)),
          open: () => openRow(e), compare: () => chosen.peek().length >= 2 ? () => compare.click() : null});
      } else if (grok.shell.o instanceof PmxSavedStudy && grok.shell.o.o.key === keyOf(e))
        grok.shell.o = null;
    };
    return box;
  };
  const chosen = computed(() => all.value.filter((e) => comparable(e) && picked.value.has(keyOf(e))));
  const compare = ui.button('Compare', () => void openCompare(chosen.peek()));
  const compareHint = h('span', 'pmx-muted', ['Pick two studies with measurements or more to compare them']);
  // studies of one drug picked: its regulatory summary, the rows each study's last analysis kept
  const drug = computed(() => chosen.value.length >= 2 ? commonDrug(chosen.value.map((e) => e.drugs)) : null);
  const summary = ui.button('CTD summary', () => {
    const [entries, name] = [chosen.peek(), drug.peek()];
    if (name)
      void ctdSummary(entries, name).catch((e) => grok.shell.error(plainError(e, 'The summary could not be made.')));
  });
  summary.title = 'The studies as rows of the nonclinical dossier (CTD), with Download (.docx)';
  effect(() => {
    const n = chosen.value.length;
    compare.textContent = `Compare (${n})`;
    compare.disabled = n < 2;
    compare.hidden = n < 2;
    compareHint.hidden = n >= 2 || trash;
    summary.hidden = !drug.value;
    summary.textContent = `CTD summary for ${drug.value ?? ''}`;
  });
  const W = libraryWidths(domainStorageOn());
  const pc = (n: number) => `${n}%`;
  const table = new BasicTable<LibraryEntry>({items, keyOf, fixedLayout: true,
    columns: [
      // (the tick's column a fixed width, outside the shares: as a share it fell to 9 px at 900 px with the context
      // panel open and cut the checkbox)
      {header: '', width: PICK_WIDTH, render: (e) => trash ? '' : pick(e)},
      {header: 'Study', width: pc(W.study), render: (e) => h('div', 'pmx-study-cell', [clip(e.name),
        clip(muted(e.copyOf ? `${e.id} · copy of ${e.copyOf.owner || 'a shared study'}` : e.id))])},
      {header: 'Drug(s)', width: pc(W.drugs), render: (e) => clip(e.drugs)},
      {header: 'Compound and schedule', width: pc(W.treatment), render: (e) => clip(e.treatment)},
      {header: 'Model', width: pc(W.model), render: (e) => clip(e.line)},
      {header: 'Stage', width: pc(W.stage), render: (e) => clip(stageMarker(e))},
      {header: 'Progress', width: pc(W.progress), render: (e) => clip(progressOf(e))},
      {header: 'Owner · Run by', width: pc(W.owner), render: (e) =>
        clip(muted(e.runBy ? `${e.owner} · ${e.runBy}` : e.owner))},
      {header: 'Updated', width: pc(W.updated), render: (e) => clip(muted(e.updated ? new Date(e.updated)
        .toLocaleDateString('en-GB', {day: 'numeric', month: 'short', year: 'numeric'}) : ''))},
      {header: 'Tags', width: pc(W.tags), render: (e) => {
        const tags = clip(h('span', '', e.tags.map((t) => h('span', 'pmx-tag', [t]))));
        tags.title = e.tags.join(', ');
        return tags;
      }},
      ...domainStorageOn() ? [{header: '', width: pc(W.action!), render: (e: LibraryEntry) =>
        !e.row || e.readOnly ? '' :
          trash ? ui.link('Restore', (ev?: MouseEvent) => {
            ev?.stopPropagation(); // ui.link calls back with no event
            void restore(e.row!.kind, e.row!.id).then(() => {
              grok.shell.info(`"${esc(e.name)}" is back in your studies`);
              reload();
            }, (err) => grok.shell.error(saveError(err, e.row!.kind)));
          }, 'Bring it back to your studies') :
            h('span', 'pmx-row-actions', [ui.iconFA('share-alt', (ev: MouseEvent) => {
              ev.stopPropagation();
              shareDialog(e);
            }, 'Share to read'), ui.iconFA('trash-alt', (ev: MouseEvent) => {
              ev.stopPropagation();
              modal('delete', `Delete "${e.name}"?`)
                ?.add(ui.divText('It goes to the trash: Show deleted brings it back.'))
                .onOK(() => void remove(e.row!.kind, e.row!.id).then(reload,
                  (err) => grok.shell.error(saveError(err, e.row!.kind)))).show();
            }, 'Delete (to the trash)')])}] : [],
    ],
    onRowClick: (e) => openRow(e)});
  table.root.style.cursor = 'pointer';
  const title = h('h3', '', ['Your studies']);
  const deleted = ui.link('Show deleted', () => {
    trash = !trash;
    // (nothing deleted is compared: the picks are dropped, and the boxes are not shown there)
    picked.value = pickedInOverview = new Set();
    title.textContent = trash ? 'Deleted' : 'Your studies';
    deleted.textContent = trash ? 'Back to your studies' : 'Show deleted';
    reload();
  }, 'Plans and studies you deleted');
  deleted.hidden = !domainStorageOn();
  // with domain storage, when the first opening could not bring the files in: try again (the files stay)
  const bring = ui.link('Bring your files in', () => {
    bring.hidden = true;
    void importFolders().then((r) => {
      setUserSetting(FILES_IN, 'yes');
      const parts = [`${r.imported} brought in`, ...r.already ? [`${r.already} already in`] : [],
        ...r.renumbered.length ? [`${r.renumbered.length} under a new ID, one with theirs being there already ` +
          `(${r.renumbered.join(', ')})`] : []];
      grok.shell.info(esc(`Your files: ${parts.join(', ')}. The files are left as they were.`));
      if (r.damaged.length)
        grok.shell.warning(`Could not be read, so not brought in: ${esc(r.damaged.join(', '))}`);
      reload();
    }, (e) => {
      bring.hidden = false;
      grok.shell.error(plainError(e, 'Your files could not be brought in just now. Try again in a moment.'));
    });
  }, 'Your plans and studies saved as files before: each becomes a row, once (run again, nothing doubles)');
  bring.hidden = true;
  bringLink = bring;
  const head = h('div', 'pmx-library-head', [title, search.root, chips, compare, summary, compareHint, bring,
    deleted]);
  // storage off: rows saved while it was on are not listed, so the Overview says they are there (not lost)
  const kept = h('div', 'pmx-muted');
  kept.hidden = true;
  if (!domainStorageOn()) {
    void keptRows().then((n) => {
      kept.textContent = keptRowsLine(n);
      kept.hidden = !kept.textContent;
    }, () => {}); // (no domain tables to read: nothing to say)
  }
  host.append(head, today, kept, h('div', 'pmx-library', [table.root, none, tryCompare]));
  return {root: host, reload};
}

/** Takes the files dragged over `target` (and only files): `over` says when some are over it, `got` receives them.
 *  The platform's own drop handling (which would open them as tables) does not see them. */
function takeDrops(target: HTMLElement, over: (on: boolean) => void, got: (files: File[]) => void): void {
  let depth = 0;
  const hasFiles = (e: DragEvent) => !!e.dataTransfer && [...e.dataTransfer.types].includes('Files');
  // every drag event over the page stays the page's: the platform shows its own drop overlay ("Incoming!") on a
  // dragenter it sees, and hides it on a drop it sees, so it was left on screen by a drop taken here
  target.addEventListener('dragenter', (e) => {
    if (!hasFiles(e))
      return;
    e.preventDefault();
    e.stopPropagation();
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
    if (hasFiles(e))
      e.stopPropagation();
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

/** Where the dialog was opened from: a stage's Open…, or a file dropped on its page. `history`: Plan's own use of a
 *  table, as the growth history of a model outside the atlas. */
export interface OpenAt {stage: OpenStage, history?: (t: MeasuredTable) => void}

const labelOf = (s: OpenStage) => STAGES.find((x) => x.stage === s)!.label;
/** What a stage's Open a study takes, said under its drop zone. */
const TAKES: {[s in OpenStage]: string} = {
  plan: 'A saved plan, or the growth history of a model outside the atlas',
  monitor: 'A followed study, or the measurements of a running study with its plan',
  analyze: 'The measurements of a finished study, with its plan. Mouse PK is added in Analyze, once it is open',
};

/** A stage's page taking files dropped on it, as the Overview does: the Open a study dialog opens with them read.
 *  The layer goes first: the page's last child is the panes, which take the height left. */
export function takesFiles(page: HTMLElement, at: OpenAt): HTMLElement {
  const layer = dropLayer();
  page.style.position = 'relative';
  page.prepend(layer);
  takeDrops(page, (on) => layer.classList.toggle('is-open', on), (files) => void openStudyDialog(files, at));
  return page;
}

const OPEN_DIALOG = 'open-study';
/** What the Open a study dialog on screen does with files dropped on the page meanwhile: reads them in it. */
let takeInOpen: (files: File[]) => Promise<void> = async () => {};
/** Where the Open a study dialog on screen was opened from (its stage, and Plan's history use). */
let openAt: OpenAt | undefined;
const sameAt = (a?: OpenAt, b?: OpenAt) => a?.stage === b?.stage && !!a?.history === !!b?.history;

/** The Open a study dialog: a file dropped on it, browsed for, or one of the user's plans; once read, what it holds
 *  and the stage it opens in (suggested from its dates). `first`: files dropped on the page, read at once. One at a
 *  time: while one is open, it is the one answered (with `first` read in it). `at`: the stage it was opened from,
 *  which a file opens in when it can; a file of another stage opens in its own, and the dialog says so. */
export function openStudyDialog(first?: File[], at?: OpenAt): DG.Dialog {
  let phase: Phase = 'drop';
  let brought: Brought = {};
  let chosen: OpenStage | null = at?.stage ?? null;
  let error = '';
  let loadingName = '';
  let latest = 0;
  // a study picked in Saved studies, alone in the dialog: it opens as its Overview row does (openEntryOnce)
  let savedPick: LibraryEntry | null = null;

  // one open from elsewhere (the Overview's, another stage's) gives way: a file dropped here is read for this stage
  const other = modalOf(OPEN_DIALOG);
  if (other && !sameAt(openAt, at)) {
    other.close();
    releaseModal(OPEN_DIALOG);
  }
  const dialog = modal(OPEN_DIALOG, {title: at ? `Open a study in ${labelOf(at.stage)}` : 'Open a study'});
  if (!dialog) {
    if (first?.length)
      void takeInOpen(first);
    return modalOf(OPEN_DIALOG)!;
  }
  openAt = at;
  const body = h('div', 'pmx-open');
  dialog.add(body);
  dialog.onOK(() => void go().catch((e) => grok.shell.error(plainError(e, 'The study could not be opened. Try ' +
    'again; if it fails again, reload the page.'))), {closeOnEnter: false});
  // the OK button says where the study opens, and waits until there is one
  const ok: HTMLButtonElement = dialog.getButton('OK') ?? document.createElement('button');
  ok.classList.add('pmx-primary', 'ui-btn-raised');
  // and up to two secondary buttons, the other stages it may open in (named and shown by render)
  const altStages: (OpenStage | null)[] = [null, null];
  const alts = altStages.map((_, i) => {
    const text = `pmx-alt-${i}`;
    dialog.addButton(text, () => {
      const stage = altStages[i];
      if (!stage)
        return;
      chosen = stage;
      dialog.close();
      void go().catch((e) => grok.shell.error(plainError(e, 'The study could not be opened. Try again; if it ' +
        'fails again, reload the page.')));
    });
    const b: HTMLButtonElement = dialog.getButton(text) ?? document.createElement('button');
    b.style.display = 'none';
    return b;
  });

  const options = (): {list: Option[], hint: string} => {
    const {table, plan, study} = brought;
    if (study) {
      // a followed study is Monitor's file: nothing to choose
      return {hint: '', list: [
        {stage: 'monitor', title: 'Open it', desc: 'Monitor', note: '', disabled: false}]};
    }
    if (table) {
      const dated = !!plan && table.lastDay !== null && plan.end !== null;
      // in Plan, a table may be a model's growth history (a model outside the atlas)
      const history = !!at?.history && !plan;
      const asHistory = 'A table of measurements opens in Analyze. To plan with it as the growth history of a model ' +
        'outside the atlas, pick that below.';
      return {
        hint: history ? asHistory : !dated ? '' : table.lastDay! > plan!.end! ?
          // (from the last day measured: a table may have days and no dates)
          `Measured to day ${table.lastDay}, past the planned end (day ${plan!.end}).` :
          `Measured to day ${table.lastDay} of about ${plan!.end}.`,
        list: [
          ...history ? [{stage: 'plan' as const, title: 'A model\'s growth history', note: '', disabled: false,
            desc: 'Plan, with a model outside the atlas'}] : [],
          {stage: 'monitor', title: 'Still running', desc: 'Monitor', disabled: !plan,
            note: !plan ? 'Add the plan it is run with to follow it' : ''},
          {stage: 'analyze', title: 'Finished', note: '', disabled: false,
            desc: plan ? 'Analyze, against its plan' : 'Analyze, against the historical controls'},
        ]};
    }
    return {hint: '', list: [
      {stage: 'plan', title: 'Not started', desc: 'Plan', note: '', disabled: false},
      {stage: 'monitor', title: 'Starting now', desc: 'Monitor', note: '', disabled: false}]};
  };

  const selected = (): Option | null => {
    // a table not read, or waiting on an answer, opens nowhere
    const read = brought.table?.read;
    if (phase !== 'loaded' || read && (read.error || waiting(read)))
      return null;
    const list = options().list;
    const pick = list.find((o) => o.stage === chosen && !o.disabled) ??
      list.find((o) => o.stage === suggestedStage(brought) && !o.disabled) ??
      list.find((o) => !o.disabled);
    return pick ?? null;
  };

  // the tumor model of a table without its plan: none assumed (not the one picked last, not MC38) until the user
  // picks it here (QA 0.1.3 round 3: analyst 3, chaos 5, technician N2)
  let line: string | null = null;
  let lines: string[] | null = null;
  const modelPicker = (): HTMLElement => {
    if (!lines) {
      lines = [];
      void tumorModels().then((m) => {
        lines = Object.keys(m).sort();
        render();
      }, () => {});
    }
    const input = ui.input.choice<string>('Tumor model', {value: line ?? undefined, items: lines, nullable: true,
      onValueChanged: (v: string | null) => {
        line = v || null;
        render();
      }});
    ui.tooltip.bind(input.root, 'The table names no model: its untreated growth is predicted with the one picked ' +
      'here. Pick the study’s own.');
    return h('div', 'pmx-ask', [...line ? [] : [h('div', 'pmx-q', ['Pick the tumor model'])], input.root]);
  };

  const askTable = (q: Question): HTMLElement => {
    const t = brought.table!;
    const now = t.answers[q.id] ?? q.preset;
    return h('div', 'pmx-ask', [h('div', 'pmx-q', [q.text]), h('div', 'pmx-opts', q.options.map(([value, label]) => {
      const radio = h('input', '', [], {type: 'radio', name: `pmx-read-${q.id}`});
      radio.checked = now === value;
      radio.addEventListener('change', () => {
        const again = readFile(t.file, t.text, {...t.answers, [q.id]: value}, t.merged);
        if (again.kind === 'table')
          brought = {...brought, table: again};
        render();
      });
      return h('label', 'pmx-opt' + (now === value ? ' is-on' : ''), [radio, h('span', '', [label])]);
    }))]);
  };

  const dropZone = (small: boolean, text: string): HTMLElement => {
    const input = h('input', 'pmx-file-input', [], {type: 'file', multiple: '',
      accept: OPENED_EXTS.map((x) => '.' + x).join(',')});
    input.addEventListener('change', () => {
      void take([...input.files ?? []]);
      input.value = '';
    });
    return h('label', 'pmx-dz' + (small ? ' pmx-dz-small' : ''), small ?
      [h('span', 'pmx-ico', [icon(ICON.add, 14, 2)]), text, input] :
      [h('span', 'pmx-ico', [icon(ICON.open)]),
        h('span', 'pmx-dz-main', ['Drag a file here, or ', h('b', '', ['browse'])]),
        ...at ? [h('span', 'pmx-muted', [TAKES[at.stage]])] : [],
        h('span', 'pmx-exts', [OPENED_EXTS.filter((x) => x !== 'txt').map((x) => '.' + x).join(' · ')]), input]);
  };
  const errorLine = () => h('div', 'pmx-error', [icon(ICON.warn, 16, 1.8), h('span', '', [error])]);
  // (in the Overview and in Monitor, the studies saved too: the same picker as Ctrl+O there, QA 0.1.3 r5)
  const savedPlansLink = (studies = false) => h('div', 'pmx-alts', [ui.link('Saved plans', () => void pickPlan(),
    'Your plans saved in Plan'), ...studies ? [ui.link('Saved studies', () => void pickPlan('study'),
    'Your studies saved in Monitor or Analyze')] : []]);

  const render = () => {
    const parts: HTMLElement[] = [];
    if (phase === 'drop' || phase === 'error') {
      if (phase === 'error')
        parts.push(errorLine());
      parts.push(dropZone(false, ''), savedPlansLink(!at || at.stage === 'monitor'));
    } else if (phase === 'loading') {
      parts.push(h('div', 'pmx-loading', [
        h('div', 'pmx-file-row', [icon(ICON.file, 20), h('div', '', [h('b', '', [loadingName])])]),
        h('div', 'pmx-track', [h('div')])]));
    } else {
      const {table, plan, study} = brought;
      parts.push(readSummary(brought));
      const notes = readNotes(brought);
      if (notes.length)
        parts.push(h('div', 'pmx-read-notes', notes));
      // what the table's reader could not tell: one choice each, the table read again with every answer
      for (const q of table?.read.questions ?? [])
        parts.push(askTable(q));
      // a table without its plan has no tumor model: the one it will be analyzed with, said and picked here (#8)
      // (only for a table read: one that was not has no growth to predict, QA 0.1.3 round 4 chaos 24)
      if (table && !plan && table.read.rows.length > 0 && !table.read.error)
        parts.push(modelPicker());
      const {hint} = options();
      const pick = selected();
      if (hint)
        parts.push(h('div', 'pmx-hint', [icon(ICON.info, 16, 1.8), h('span', '', [hint])]));
      // opened from a stage, a file of another opens in its own: said, never refused
      // (the history offered in Plan says it in its hint)
      const away = at && pick && !(at.history && table && !plan) ? routedNote(brought, at.stage, pick.stage) : null;
      if (away)
        parts.push(h('div', 'pmx-hint', [icon(ICON.info, 16, 1.8), h('span', '', [away])]));
      if (error)
        parts.push(errorLine());
      const missing = table && !plan ? 'Add its .pmxplan' : plan && !table && !study ? 'Add its measurements' : '';
      if (missing)
        parts.push(dropZone(true, missing));
      // measurements without their plan: a plan saved in PMX Studio links as well as a .pmxplan file
      if (table && !plan)
        parts.push(savedPlansLink());
    }
    body.replaceChildren(h('div', 'pmx-in', parts));
    const pick = selected();
    ok.textContent = pick ? `Open in ${pick.desc.split(',')[0]}` : 'Open';
    ok.title = pick ? `${pick.title}: ${pick.desc}` : '';
    ok.disabled = !pick;
    ok.classList.toggle('disabled', !pick);
    // the other stages it may open in, as secondary buttons beside the primary (a table that could not be read, or
    // waiting on an answer, none); one that cannot yet says why on hover
    const others = !pick ? [] : options().list.filter((o) => o.stage !== pick.stage);
    alts.forEach((b, i) => {
      const o = others[i];
      altStages[i] = o?.stage ?? null;
      b.style.display = o ? '' : 'none';
      if (!o)
        return;
      b.textContent = `Open in ${o.desc.split(',')[0]}`;
      b.disabled = o.disabled;
      b.classList.toggle('disabled', o.disabled);
      b.title = o.note || `${o.title}: ${o.desc}`;
    });
  };

  /** Reads the files and adds what they hold to what the dialog already has; a later read overtakes an earlier. */
  const take = async (files: File[]) => {
    if (files.length === 0)
      return;
    savedPick = null;
    const mine = ++latest;
    phase = 'loading';
    loadingName = files.map((f) => f.name).join(', ');
    render();
    const started = Date.now();
    const reads = await Promise.all(files.map((f) => fileText(f).then((t) => readFile(f.name, t.text, {}, t.merged),
      (e) => ({kind: 'bad' as const, file: f.name, why: plainError(e, 'The file could not be read.')}))));
    // the bar shows for a moment even for a small file: the change of state reads as the file being taken in
    await new Promise((r) => setTimeout(r, Math.max(0, 450 - (Date.now() - started))));
    if (mine !== latest)
      return;
    const {brought: next} = bring(brought, reads);
    // each file set aside, with what is wrong with it
    error = reads.map((r) => r.kind === 'bad' ? `"${r.file}" could not be opened. ${r.why}` : '').filter((s) => s)
      .join(' ');
    brought = next;
    // a table in Plan: a study, opened in Analyze, unless the user picks it as a growth history
    chosen = next.table && !next.plan && at?.history ? null : at?.stage ?? null;
    phase = next.table || next.plan || next.study ? 'loaded' : 'error';
    render();
  };

  const pickPlan = async (kind: 'plan' | 'study' = 'plan') => {
    const plans = await savedPlans(undefined, kind).catch((e) => {
      grok.shell.error(plainError(e, `Your saved ${kind === 'plan' ? 'plans' : 'studies'} could not be listed. Try ` +
        'again in a moment.'));
      return null;
    });
    if (plans?.length === 0) {
      grok.shell.info(kind === 'plan' ? 'No saved plans yet: save one in Plan, or drop a .pmxplan file here.' :
        'No saved studies yet: follow a plan in Monitor and save it.');
    }
    if (!plans?.length)
      return;
    DG.Menu.popup().items(plans.map((p) => p.name), (name: string) => {
      const p = plans.find((x) => x.name === name)!;
      void p.read().then(async (f) => {
        await take([new File([await f.readAsString()], f.fileName)]);
        savedPick = kind === 'study' && brought.study && !brought.table && !brought.plan ? p.entry : null;
      }, (e) =>
        grok.shell.error(esc(plainError(e, `"${name}" could not be read. Try again in a moment.`))));
    }).show();
  };

  const go = async () => {
    const pick = selected();
    const {table, plan} = brought;
    if (!pick)
      return;
    if (pick.stage === 'plan' && table && !plan && at?.history)
      at.history(table);
    else if (savedPick && pick.stage !== 'plan')
      await openEntryOnce(savedPick, pick.stage);
    else {
      // the model said in the dialog is the study's (no longer one Analyze assumes)
      if (table && !plan && line)
        table.df.setTag(LINE_TAG, line);
      await openBrought(brought, pick.stage);
    }
  };

  takeDrops(dialog.root, (on) => body.querySelectorAll('.pmx-dz').forEach((z) => z.classList.toggle('is-over', on)),
    (files) => void take(files));
  render();
  takeInOpen = take;
  dialog.show();
  if (first?.length)
    void take(first);
  return dialog;
}

/** Opens what the files brought hold in `stage`; `exposure`, a concentration table Analyze reads with it (NCA). All
 *  of it read from files: a first save makes the user's copy of what someone else owns (as read now). */
export async function openBrought(b: Brought, stage: OpenStage, exposure?: DG.DataFrame): Promise<void> {
  const {table, plan, study} = b;
  if (stage === 'plan' && plan) {
    markFromFile(plan.saved, ownerOf('plan', plan.saved));
    await showStage('plan', await planView(plan.saved));
  } else if (stage === 'analyze' && table) {
    if (plan)
      table.df.setTag(PLAN_TAG, planTag(plan.text));
    markTableFromFile(table.df, metaOf(table.df)?.owner ?? '');
    await showStage('analyze', await analyzeView({study: table.df, exposure}));
  } else if (stage === 'monitor' && plan && !table && !study) {
    // a plan starting now: its study already followed opens (Follow's way), else a new empty table
    await followPlan(plan.saved, () => {
      const df = newStudyTable(plan.saved);
      markTableFromFile(df, ownerOf('plan', plan.saved));
      return df;
    });
  } else if (stage === 'monitor') {
    // a followed study as saved; a table under way with its plan
    const df = study ? tableOf(study.saved) : table?.df ?? (plan ? newStudyTable(plan.saved) : null);
    if (!df)
      return;
    if (table && plan)
      df.setTag(PLAN_TAG, planTag(plan.text));
    markTableFromFile(df, study ? ownerOf('study', study.saved) : table ? metaOf(df)?.owner ?? '' :
      ownerOf('plan', plan!.saved));
    await showStage('monitor', await monitorView(df));
  }
}

// ---- the demo data ------------------------------------------------------------------------------------------------

/** The demo studies (files/demo, with the package: synthetic, see its README), one per stage and two to compare. */
export const DEMO = {
  plan: {plan: 'ct26-trametinib-antipd1.pmxplan', what: 'CT26, trametinib + anti-PD-1'},
  monitor: {study: 'ct26-running-day9.pmxstudy', what: 'The same study under way, animals flagged'},
  analyze: {table: 'ct26-trametinib-antipd1-finished.csv', plan: 'ct26-trametinib-antipd1.pmxplan',
    // (read on day 16, the day its median control reached the endpoint: "Day 21" beside "evaluated on day 16" read as
    // a contradiction, round 3 designer 6)
    pk: 'ct26-trametinib-pk-mouse.csv', what: 'Finished study (read on day 16), with mouse PK'},
  compare: [{table: 'ct26-trametinib-antipd1-finished.csv', plan: 'ct26-trametinib-antipd1.pmxplan',
    pk: 'ct26-trametinib-pk-mouse.csv', line: 'CT26'},
  {table: 'mc38-trametinib-finished.csv', pk: 'mc38-trametinib-pk-mouse.csv', line: 'MC38'}],
} as const;

const demoText = (name: string) => _package.files.readAsText(`demo/${name}`);

/** The demo's day 0, a Monday as in its plan (measured Mon / Wed / Fri: days 0, 2, 4, 7, 9…), so that today is the
 *  first day from `from` on: Monitor's (from 7) is at day 7 to 13, day 9 on a Wednesday; Analyze's (from 21) at
 *  day 21 to 27, ended. The weekdays hold whatever today is (QA 0.1.3: dated so that today was day 9, a Friday
 *  measured the study on Wed / Fri / Sun). */
export function demoStart(from: number, now = new Date()): string {
  const weekday = (now.getDay() + 6) % 7;
  const today = from + ((weekday - from) % 7 + 7) % 7;
  return new Date(Date.UTC(now.getFullYear(), now.getMonth(), now.getDate() - today)).toISOString().slice(0, 10);
}

/** A demo file's text dated by demoStart; a study's measurements after today left out (Monday is day 7: day 9's
 *  values are not taken yet), so the header, the calendar and "Measure today" agree on the day the demo is at. */
export function demoDated(text: string, from: number, now = new Date()): string {
  const start = demoStart(from, now);
  const dated = text.replace(/("start"\s*:\s*)"\d{4}-\d{2}-\d{2}"/g, `$1"${start}"`);
  const today = Math.round((Date.UTC(now.getFullYear(), now.getMonth(), now.getDate()) -
    Date.parse(`${start}T00:00:00Z`)) / 86400000);
  if (!dated.trimStart().startsWith('{'))
    return dated;
  const s = JSON.parse(dated);
  if (!Array.isArray(s.rows) || s.extra || s.rows.every((r: {day?: number}) => !(Number(r.day) > today)))
    return dated;
  return JSON.stringify({...s, rows: s.rows.filter((r: {day?: number}) => !(Number(r.day) > today))});
}

/** The first day each stage's demo may be at today (demoStart). */
const DEMO_LAST_DAY: {[stage: string]: number} = {monitor: 7, analyze: 21};

/** A demo file as Open a study reads it (a damaged one would be the package's fault: said as an error), its day 0
 *  dated for today when `lastDay` is given (the files carry one fixed date). */
async function demoRead(name: string, lastDay?: number): Promise<Brought> {
  const text = await demoText(name);
  const {brought, bad} = bring({}, [readFile(name, lastDay === undefined ? text : demoDated(text, lastDay))]);
  if (bad.length)
    throw new Error(`The demo file ${name} could not be read.`);
  return brought;
}

/** What a stage's demo opens with. */
export async function demoBrought(stage: OpenStage): Promise<Brought> {
  const d = DEMO[stage];
  const parts = await Promise.all(Object.entries(d).filter(([k]) => k === 'plan' || k === 'study' || k === 'table')
    .map(([, name]) => demoRead(name, DEMO_LAST_DAY[stage])));
  return Object.assign({}, ...parts);
}

/** A concentration table of the demo (null without the NCA package: nothing reads it). */
async function demoPk(name: string): Promise<DG.DataFrame | null> {
  if (!ncaInstalled())
    return null;
  const df = DG.DataFrame.fromCsv(await demoText(name));
  df.name = name.replace(/\.csv$/, '');
  // the app's own table: closing or reloading the untouched demo asks nothing (native.ts appTable)
  df.setTag(APP_TABLE_TAG, 'true');
  return df;
}

/** Opens a stage's demo there (Analyze with its mouse PK, when the NCA package is there). Once: a second click while
 *  it opens does nothing. */
export function openDemo(stage: OpenStage): Promise<void> {
  return preparing(`demo-${stage}`, async () => {
    const d = DEMO[stage];
    const exposure = 'pk' in d ? await demoPk(d.pk) : null;
    await openBrought(await demoBrought(stage), stage, exposure ?? undefined);
  }).then(() => {}, (e) => grok.shell.error(plainError(e, 'The demo could not be opened. Try again in a moment.')));
}

/** The demo plan's ID: its earlier study is the demo's own (demoFinishedTable), never one of the user's library. */
export const DEMO_PLAN_ID = 'PMX-DEMO-CT26';

let demoFinished: Promise<DG.DataFrame> | null = null;
/** The demo's finished study (Analyze's demo table, read as the demo reads it), marked finished, of its model: what
 *  the demo plan scales its arms' effects from. Read once; a failed read is not kept. */
export function demoFinishedTable(): Promise<DG.DataFrame> {
  demoFinished ??= demoRead(DEMO.analyze.table).then(({table}) => {
    const df = table!.df;
    df.setTag(FINISHED_TAG, 'true');
    df.setTag(LINE_TAG, 'CT26');
    df.name = 'CT26 demo study, finished';
    return df;
  }).catch((e) => {
    demoFinished = null;
    throw e;
  });
  return demoFinished;
}

/** The tag the demo studies are saved with: how a later try finds them among the user's studies. */
export const DEMO_TAG = 'PMX Studio demo';
/** LibraryEntry.stage of a finished study. */
const FINISHED_STAGE = 3;

/** A finished demo study as Analyze saves it: its plan, the model and the analysis's effect (each arm's TGI on its
 *  evaluation day, which Compare reads), NCA's exposure from its mouse PK when the NCA package is there. Run headless,
 *  as the Analyze view runs it, without the twin (Compare does not read it). */
async function demoStudy(d: typeof DEMO.compare[number], df: DG.DataFrame, plan?: Brought['plan']):
  Promise<SavedStudy> {
  const pk = await demoPk(d.pk);
  const mapping = pk ? await concentrationMapping(pk) : null;
  let nca: NcaByArm | null = null;
  const effect = await headlessEffect(df, {line: d.line, plan: plan?.text ?? null,
    endpoint: plan?.saved.predicted.endpoint ?? 1000}, async (arms) =>
    nca = pk && mapping ? ncaByArm(arms, await runNca(pk, mapping), pk.name, metaOf(df)) : null);
  df.setTag(EXPOSURE_TAG, JSON.stringify({...exposureJson(new Map(), nca) ?? {}, effect}));
  return savedStudyOf(df, null);
}

/** What Analyze keeps of a study's analysis (keptEffect: the comparisons, body weight as weighed, each combination
 *  against its drugs), run headless as the Analyze view runs it on the study, without the twin: the excluded rows
 *  left out, the plan's model and endpoint. `nca`: the exposure of the arms, when there is one. */
async function headlessEffect(df: DG.DataFrame, o: {line: string, plan: string | null, endpoint: number,
  humane?: number | null, weightLimit?: number | null}, nca: (arms: string[]) => Promise<NcaByArm | null> =
  async () => null): Promise<NonNullable<NonNullable<SavedStudy['exposure']>['effect']>> {
  // the body weight and removals Analyze takes by name, so the regulatory table says what Analyze shows (round 4,
  // analyst A5: "BW nadir; deaths" — beside Analyze's 2 of 10)
  const named = (input: string) => df.columns.toList().find((c) => OPTIONAL_COLUMNS[input](c.name)) ?? null;
  const roles = roleColumns(df);
  const x = withoutExcluded({study: df, animal: roles.animal, group: roles.group, day: roles.day,
    volume: roles.volume, plan: o.plan, line: o.line, control: '', endpoint: o.endpoint, correction: 'bonferroni',
    withTwin: false, bodyWeight: named('bodyWeight'), tumorWeight: named('tumorWeight'),
    removalReason: named('removalReason'), ...o.humane ? {humane: o.humane} : {},
    ...o.weightLimit ? {weightLimit: o.weightLimit} : {}});
  const frames = await analysisFrames(analysisDeclarationOf(x));
  const arms = frames.comparisons.col('arm')!.toList() as string[];
  const exposure = await nca(arms);
  // each combination against its drugs, as Analyze's ladder (the same function)
  const lp = o.line ? await linePrior(o.line).catch(() => null) : null;
  const logBeta = lp?.study_prior?.mean_log_beta;
  const {found} = await combinationSteps({study: x.study, animal: x.animal.name, group: x.group.name,
    day: x.day.name, volume: x.volume.name, bodyWeight: x.bodyWeight?.name ?? null}, arms,
  notCombined(o.plan), frames.analysis.get('control', 0), o.line,
  frames.analysis.get('evaluation_day', 0) ?? Infinity, logBeta === undefined ? undefined : Math.exp(logBeta),
  () => '');
  return keptEffect(frames, metaOf(df), exposure, found);
}

/** A saved study's rows of the regulatory table, computed from its data as Analyze computes them (an analysis kept
 *  before the table had its growth estimate, body weight or combinations said only "p = 0.088", QA 0.1.3 r5); each
 *  group's exposure as last analyzed (its concentrations are not in the study); its kept rows if it cannot be run. */
export async function ctdOfStudy(s: SavedStudy): Promise<CtdStudy | null> {
  const kept = ctdOfSaved(s);
  try {
    const df = tableOf(s);
    const plan = df.getTag(PLAN_TAG) || null;
    const {ctd} = await headlessEffect(df, {line: s.plan?.predicted.line || s.line || df.getTag(LINE_TAG) || '',
      plan, endpoint: s.plan?.predicted.endpoint ?? 1000, humane: s.humane, weightLimit: s.weightLimit});
    if (!ctd)
      return kept;
    const pk = new Map((kept?.groups ?? []).filter((g) => g.pk).map((g) => [g.name, g.pk!]));
    return {...ctd, id: s.meta?.id || ctd.id,
      groups: ctd.groups.map((g) => pk.has(g.name) ? {...g, pk: pk.get(g.name)} : g)};
  } catch (e) {
    console.warn('PMX Studio:', e);
    return kept;
  }
}

/** The finished demo studies among the user's (`mine`), each saved with `saveDoc` unless there already (by its tag
 *  and name; the deleted are not listed): the entries Compare opens. */
export async function demoStudies(mine: LibraryEntry[], saveDoc: (doc: SavedStudy) => Promise<unknown>):
  Promise<LibraryEntry[]> {
  const entries: LibraryEntry[] = [];
  for (const d of DEMO.compare) {
    // read and tagged first (quick): its name, given by its plan, finds it among the saved
    // (dated as Analyze's demo is: a finished study's day 0 three weeks back, never a date after its data)
    const {table, plan} = {...await demoRead(d.table, DEMO_LAST_DAY.analyze),
      ...'plan' in d ? await demoRead(d.plan, DEMO_LAST_DAY.analyze) : {}};
    const df = table!.df;
    if (plan)
      df.setTag(PLAN_TAG, planTag(plan.text));
    df.setTag(FINISHED_TAG, 'true');
    df.setTag(LINE_TAG, d.line);
    // (an ID for each: the one of its plan, else the demo's)
    const meta = metaOf(df) ?? emptyMeta();
    setMetaOn(df, {...meta, id: meta.id || plan?.saved.study?.id || `PMX-DEMO-${d.line}`, tags: [DEMO_TAG]});
    const name = studyNameOf(df);
    const there = mine.find((e) => e.kind === 'study' && e.stage === FINISHED_STAGE && e.tags.includes(DEMO_TAG) &&
      e.name === name);
    if (there) {
      entries.push(there);
      continue;
    }
    const doc = await demoStudy(d, df, plan);
    try {
      await saveDoc(doc);
    } catch (e) {
      // the first saved, the second not: said, so the user knows one is among their studies
      const before = entries.length ? `"${entries.map((x) => x.name).join('", "')}" is saved; ` : '';
      throw new Error(`${before}"${doc.name}" could not be saved: ${plainError(e, 'try again in a moment')}`);
    }
    entries.push(entryOf(JSON.stringify(doc), Date.now())!);
  }
  return entries;
}

/** The two finished demo studies saved among the user's studies (unless saved already), then opened in Compare
 *  (the Compare it opened before, while open). */
export function compareDemo(): Promise<void> {
  if (showOpened('demo-compare'))
    return Promise.resolve();
  return preparing('demo-compare', async () => {
    const progress = DG.TaskBarProgressIndicator.create('Saving the demo studies...');
    try {
      const entries = await demoStudies(await myEntries(),
        (doc) => domainStorageOn() ? save('study', doc) : saveStudy(doc));
      const view = await showStage('home', await compareView(entries));
      view.root.dataset[OPENED_AS] = 'demo-compare';
    } finally {
      progress.close();
    }
  }).then(() => {}, (e) => grok.shell.error(plainError(e, 'The demo studies could not be compared. Try again in ' +
    'a moment.')));
}

/** What the files brought are, as the stages' header names a study: its name, then groups, animals and days, then the
 *  files it was read from. */
function readSummary(b: Brought): HTMLElement {
  const {table: t, plan, study} = b;
  let name = '';
  let facts: string[] = [];
  if (t) {
    name = t.df.name || t.file;
    facts = t.read.error ? ['Not read'] : [t.groups !== null ? plural(t.groups, 'group') : '',
      t.animals !== null ? plural(t.animals, 'animal') : '',
      t.lastDay !== null ? `days ${t.firstDay}–${t.lastDay}` : '',
      t.weight ? 'body weight' : 'no body weight'];
  } else if (study) {
    const s = study.saved;
    const animals = new Set(s.rows.map((r) => r.animal)).size;
    const days = s.rows.map((r) => r.day).filter((d): d is number => typeof d === 'number');
    name = s.name || study.file;
    facts = ['Followed study', plural(animals, 'animal'), days.length ? `measured to day ${maxOf(days)}` : ''];
  } else if (plan) {
    name = plan.name;
    facts = ['Plan', plan.line, plan.end !== null ? `about ${plan.end} days` : ''];
  }
  const files = [t?.file, plan ? (t ? `plan "${plan.name}"` : plan.file) : '', study?.file].filter(Boolean)
    .join(' · ');
  return h('div', 'pmx-read', [h('span', 'pmx-ico', [icon(ICON.file, 20)]), h('div', 'pmx-read-text', [
    h('div', 'pmx-read-name', [name]), h('div', 'pmx-read-facts', [facts.filter(Boolean).join(' · ')]),
    h('div', 'pmx-read-files', [files])])]);
}

/** What the reading found, in the stages' notice style: what was read and how (the check that nothing was lost or
 *  guessed) quiet, what was not read or read as one a warning. */
function readNotes(b: Brought): HTMLElement[] {
  const t = b.table;
  if (t) {
    if (t.read.error)
      return [noteLine(t.read.error, 'error')];
    if (waiting(t.read))
      return [markLine(`${t.read.questions.length > 1 ? 'Questions' : 'A question'} before the table is read`)];
    const lines = studySummary(t.read, !!b.plan?.blinded);
    const merged = mergedWarning(t.merged, !!b.plan?.blinded);
    return [...marks(lines.filter((l) => !l.startsWith('Not read'))),
      ...lines.filter((l) => l.startsWith('Not read')).map((l) => noteLine(l)), ...merged ? [noteLine(merged)] : []];
  }
  if (b.study)
    return marks(savedSummary(b.study.saved).slice(1));
  return [];
}

/** Lines of what was read under one info mark (one per line read as noise: a mark per group). */
function marks(lines: string[]): HTMLElement[] {
  if (!lines.length)
    return [];
  const block = markLine('');
  block.lastElementChild!.replaceWith(h('div', '', lines.map((l) => h('div', '', [l]))));
  return [block];
}
