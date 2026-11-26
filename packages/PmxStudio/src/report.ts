/* What a stage exports, from the view's own menu (Export), as other Datagrok apps do: the report, a page to print or
   save as PDF from the browser's print dialog (the stage's panels as shown, each chart as an image, then the study's
   details, notes and log), and the stage's files. Datagrok's viewers save themselves as PNG from their own menu; the
   echarts charts get the same two items on a right click. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as echarts from 'echarts/core';

import {startMismatch} from './calendar';
import {addColorItem} from './colors';
import {entered, esc} from './sheet';
import {provenance} from './plan-file';
import type {LogEntry} from './study-file';
import {StudyMeta, animalsLine, dateOf, dateText, ofGroups} from './study-meta';
import {isoToday, plainError} from './ui';
import {PmxFigure, popupFor} from './handlers';

/** A menu item: its name, what it does, its group (Export unless named) and, for one not always available, why
 *  it is not now (null: it is). */
export type Exports = [string, () => void | Promise<void>, string?, (() => string | null)?][];

export const failed = (e: any) =>
  grok.shell.error(plainError(e, 'The export could not be made. Try again; if it fails again, reload the page.'));

/** Runs a menu item as the menu and the toolbox do: one not available now says why (the ribbon's menu shows it
 *  enabled, and with isEnabled it swallowed the click: QA pass 4, Reset colors did nothing, said nothing). */
function runItem(x: Exports[number]): void {
  const why = x[3]?.();
  if (why)
    return void grok.shell.info(esc(`${x[0]}: ${why}.`));
  // clicks while it runs, and for a moment after (a download is done at once), are one click: one file, not three
  if (running.has(x[1]))
    return;
  running.add(x[1]);
  const started = Date.now();
  void Promise.resolve().then(x[1]).catch(failed).finally(() =>
    setTimeout(() => running.delete(x[1]), Math.max(0, ONE_CLICK_MS - (Date.now() - started))));
}

/** The items running now (their actions). */
const running = new Set<Exports[number][1]>();
/** How long after it starts an item ignores more clicks. */
export const ONE_CLICK_MS = 1000;

/** Runs `x` as a click on it in the menu or the toolbox does (for tests). */
export const clickItem = (x: Exports[number]) => runItem(x);

// a group 'Data > Open as table' is a submenu of Data
const where = (x: Exports[number]) => (x[2] ?? 'Export').split(' > ');

/** The view's menu: the stage's exports in an Export group, its other items in theirs (View: the layout). */
export function exportMenu(view: DG.ViewBase, items: Exports): void {
  const menu = DG.Menu.create();
  // one submenu with a single item is that item, named by the submenu
  for (const name of [...new Set(['Export', ...items.map((x) => where(x)[0])])]) {
    const group = menu.group(name);
    const subs = new Map<string, DG.Menu>();
    const mine = items.filter((x) => where(x)[0] === name);
    for (const x of mine) {
      const sub = where(x)[1];
      if (sub && mine.filter((y) => where(y)[1] === sub).length === 1) {
        group.item(sub, () => runItem(x));
        continue;
      }
      if (sub && !subs.has(sub))
        subs.set(sub, group.group(sub));
      (sub ? subs.get(sub)! : group).item(x[0], () => runItem(x));
    }
    subs.forEach((s) => s.endGroup());
    group.endGroup();
  }
  view.ribbonMenu = menu;
}

/** A toolbox link: its text, what it does, its tooltip. */
export type ToolLink = [string, () => void, string?];

/** A toolbox: sections of links, each open (the Overview's and every page's); a section with no link is left out. */
export function toolbox(sections: [string, ToolLink[]][]): HTMLElement {
  const acc = ui.accordion('pmx-toolbox');
  for (const [title, links] of sections.filter(([, l]) => l.length)) {
    acc.addPane(title, () => ui.divV(links.map(([text, run, tip]) => ui.link(text, () => run(), tip)),
      {style: {gap: '6px'}}), true);
  }
  return acc.root;
}

/** The toolbox sections a menu item goes to, by its group: Study, Export, Documents, View (Data's tables too). */
export const TOOLBOX_SECTIONS = ['Study', 'Export', 'Documents', 'View'];

/** A page's toolbox from the items of its view menu (the same items, run the same way): `study` first in Study
 *  (Open…, the study's details), then the menu's by section; `help` last. Read when the page is built: every page
 *  adds its items before it returns. */
export function stageToolbox(items: Exports, study: ToolLink[], help: ToolLink[]): HTMLElement {
  const section = (x: Exports[number]) => where(x)[0] === 'Data' ? 'View' : where(x)[0];
  const label = (x: Exports[number]) => where(x)[1] ? `${where(x)[1]}: ${x[0]}` : x[0];
  const of = (s: string): ToolLink[] => items.filter((x) => section(x) === s).map((x) => [label(x), () => runItem(x)]);
  return toolbox([...TOOLBOX_SECTIONS.map((s): [string, ToolLink[]] => [s, [...s === 'Study' ? study : [], ...of(s)]]),
    ['Help', help]]);
}

/** Saves a text, a file's bytes or a data URL as a file, and says so, naming it: the browser may refuse a download
 *  without a word the page can see (a page that asks for several in a row is blocked until the user allows it), so
 *  the toast also says what to do when the file does not arrive; one that could not be made is said as an error. */
export function download(name: string, content: string | Uint8Array, type = 'text/plain',
  said: string[] = [name]): void {
  try {
    const url = typeof content === 'string' && content.startsWith('data:') ? content :
      URL.createObjectURL(new Blob([typeof content === 'string' ? content : content.slice()], {type}));
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    a.click();
    if (url !== content) // (let go once the browser has taken it)
      setTimeout(() => URL.revokeObjectURL(url), 60_000);
    // (a toast reads HTML: the file is named after the study; several files at once: one toast, the last names all)
    if (said.length) {
      grok.shell.info(`Downloading ${said.map((n) => `"${esc(n)}"`).join(' and ')}. Not in your downloads? The ` +
        'browser blocked it: allow downloads from this site (the icon in the address bar) and export again.');
    }
  } catch (e) {
    grok.shell.error(esc(`"${name}" could not be downloaded: ${e instanceof Error ? e.message : e}`));
  }
}

export const pngOf = (chart: echarts.ECharts) =>
  chart.getDataURL({type: 'png', pixelRatio: 2, backgroundColor: '#fff'});

/** The legend item under a mouse event on `chart`: its name and its series' colour (null: not on a legend item). */
function legendItemAt(chart: echarts.ECharts, e: MouseEvent): {name: string, color: string} | null {
  const r = chart.getDom().getBoundingClientRect();
  // ponytail: echarts marks a legend item's group only with its private __legendDataIndex; if a version drops it,
  // the menu simply has no colour item
  let el: any = chart.getZr().handler.findHover(e.clientX - r.left, e.clientY - r.top)?.target;
  while (el && el.__legendDataIndex === undefined)
    el = el.parent;
  const name = el?.children?.().find((c: any) => c.type === 'text')?.style?.text;
  // a series' colour: its own (the KM, spider and waterfall charts) or its items' (the band charts)
  const series = ((chart.getOption().series ?? []) as {name?: string, color?: unknown, itemStyle?: {color?: unknown}}[])
    .find((s) => s.name === name);
  const color = series?.itemStyle?.color ?? series?.color;
  return typeof name === 'string' && typeof color === 'string' ? {name, color} : null;
}

/** The height drawn of `chart`'s first legend, in pixels (0: none, or not measurable). */
function legendHeight(chart: echarts.ECharts): number {
  // ponytail: echarts' getModel / getViewOfComponentModel are internal (not in its typings); if a version drops
  // them, the height reads 0 and the plot keeps its screen margin
  const c = chart as any;
  const model = c.getModel?.()?.getComponent?.('legend');
  const h = model ? c.getViewOfComponentModel?.(model)?.group?.getBoundingRect?.()?.height : 0;
  return Number.isFinite(h) ? h : 0;
}

/** An echarts chart's right-click menu: the figure's actions (handlers.ts: Save as PNG, Copy as image), named after
 *  its panel; on a group's legend item, Change color… first, the platform's colour picker, the colour then the
 *  study's (a canvas cannot be bound: the menu is built from the same actions). */
export function chartMenu(host: HTMLElement): void {
  host.addEventListener('contextmenu', (e) => {
    const chart = echarts.getInstanceByDom(host);
    // a right click the chart has taken for itself (an animal's menu on its bar) is not the figure's
    if (!chart || e.defaultPrevented)
      return;
    e.preventDefault();
    const name = host.closest<HTMLElement>('.pmx-panel')?.dataset.panel ?? 'figure';
    const item = legendItemAt(chart, e);
    popupFor(new PmxFigure(name, {chart: host}), {causedBy: e}, (menu) => {
      if (item && addColorItem(menu, item.name, item.color))
        menu.separator();
    });
  });
}

/** A chart as an image: echarts' own, a Datagrok viewer drawn by the platform's html2canvas; null: not drawn. */
async function imageOf(el: HTMLElement): Promise<string | null> {
  const chart = echarts.getInstanceByDom(el);
  if (chart) {
    // a scrolled legend ("1/2" and arrows) would leave groups off the paper: every group shown for the image; a
    // legend in more than one row (that one, or one wrapped at the sheet's width) pushes the plot down by its extra
    // rows, else they print over the plot's top labels and gridline (QA pass 4); the screen's put back after
    const option = chart.getOption();
    const legends = (option.legend ?? []) as {type?: string, itemHeight?: number}[];
    const grids = (option.grid ?? []) as {top?: unknown}[];
    const scrolled = legends.some((l) => l.type === 'scroll');
    if (scrolled)
      chart.setOption({legend: legends.map((l) => ({...l, type: 'plain'}))}, {replaceMerge: ['legend']});
    // (one row: the item's height and the legend's padding, 5 above and below, echarts' defaults)
    const extra = legends.length ? Math.max(0, legendHeight(chart) - (legends[0].itemHeight ?? 14) - 10) : 0;
    if (extra > 0)
      chart.setOption({grid: grids.map((g) => typeof g.top === 'number' ? {top: g.top + extra} : {})});
    try {
      return pngOf(chart);
    } finally {
      if (extra > 0)
        chart.setOption({grid: grids.map((g) => typeof g.top === 'number' ? {top: g.top} : {})});
      if (scrolled)
        chart.setOption({legend: legends}, {replaceMerge: ['legend']});
    }
  }
  try {
    await DG.Utils.loadJsCss(['/js/common/html2canvas.min.js']);
    const h2c = (window as any).html2canvas;
    return (await h2c(el, {backgroundColor: '#fff', scale: 2, logging: false})).toDataURL('image/png');
  } catch {
    return null;
  }
}

/** A table as HTML: its header row, then its rows, every cell as text. */
export const htmlTable = (head: string[], rows: unknown[][]) =>
  `<table><tr>${head.map((h) => `<th>${esc(h)}</th>`).join('')}` +
  `</tr>${rows.map((r) => `<tr>${r.map((c) => `<td>${esc(c)}</td>`).join('')}</tr>`).join('')}</table>`;

/** The report: `content` (the stage's results) as shown, its charts as images, after a title; then the study's
 *  details, its notes by day and its log. */
export async function reportHtml(o: {title: string, stage: string, meta: StudyMeta | null, content: HTMLElement,
  log?: LogEntry[], annex?: [string, HTMLElement][],
  /** the study's groups: details kept for another name (a group renamed, a plan's "Control") are not printed */
  groups?: string[],
  /** the atlas of the plan the report rests on: a hash, null when the plan records none; undefined: no plan */
  planAtlas?: string | null,
  /** the table's last measured day: day 0's date is not printed when the table runs past today's study day */
  lastDay?: number | null}): Promise<string> {
  const charts = [...o.content.querySelectorAll<HTMLElement>('[_echarts_instance_], .d4-viewer')]
    .filter((el) => el.offsetParent !== null && !el.parentElement?.closest('[_echarts_instance_], .d4-viewer'));
  const images = await Promise.all(charts.map(imageOf));
  charts.forEach((el, i) => el.dataset.pmxFigure = String(i));
  // a choice inside the text ("Randomize on a Monday"): printed as the words it shows (read before the copy, which
  // does not keep what was picked); the form's own inputs (.ui-input-root) still go below
  const choices = [...o.content.querySelectorAll<HTMLSelectElement | HTMLInputElement>(
    'select, input[type="text"], input[type="number"], input:not([type])')];
  choices.forEach((el) => el.dataset.pmxValue = el instanceof HTMLSelectElement ?
    el.selectedOptions[0]?.text ?? '' : el.value);
  const page = o.content.cloneNode(true) as HTMLElement;
  charts.forEach((el) => delete el.dataset.pmxFigure);
  choices.forEach((el) => delete el.dataset.pmxValue);
  for (const el of page.querySelectorAll<HTMLElement>('[data-pmx-value]'))
    el.replaceWith(Object.assign(document.createElement('b'), {textContent: el.dataset.pmxValue}));
  // a panel's hint says what it shows; its sentences on how to work it (double-click, right-click) are not printed
  for (const el of page.querySelectorAll<HTMLElement>('.u2-section-hint')) {
    el.textContent = (el.textContent ?? '').split(/(?<=\.)\s+/)
      .filter((s) => !/\b(double-click|right-click|click|drag|hover)\b/i.test(s)).join(' ');
  }
  for (const el of page.querySelectorAll<HTMLElement>('[data-pmx-figure]')) {
    const src = images[Number(el.dataset.pmxFigure)];
    const img = document.createElement(src ? 'img' : 'p');
    if (src)
      Object.assign(img, {src, alt: 'figure'});
    else
      img.textContent = '(this figure could not be drawn)';
    Object.assign(img.style, {width: '100%', height: 'auto'});
    // the screen's fixed heights around a chart (ui.ts sized) would let the image, scaled to the sheet's width, run
    // over what follows it (the growth chart over the calendar's heading)
    for (let p = el.parentElement; p && p !== page.parentElement; p = p.parentElement)
      p.style.height = '';
    el.replaceWith(img);
  }
  // what is for working, not for reading: the change badges on the figures and the links that act included
  page.querySelectorAll('button, input, select, textarea, [hidden], .pmx-noprint, .ui-input-root, .u2-stat-delta, ' +
    '.pmx-hint-link').forEach((e) => e.remove());

  const m = o.meta;
  const prov = await provenance().catch(() => undefined);
  const planSaid = o.planAtlas === undefined || o.planAtlas === prov?.atlas ? '' :
    `; its plan: atlas ${o.planAtlas === null ? 'not recorded' : o.planAtlas.slice(0, 12)}`;
  // a day 0 dated after the measurements it should precede is not this table's: said, and no day dated from it
  const mismatch = startMismatch(m?.start ?? null, o.lastDay);
  const day0 = mismatch ? null : dateOf(m?.start ?? null, 0);
  // the host strain in the header, a report's first question (said when the study's details lack it)
  const given = m?.animals?.strain?.trim() ?? '';
  const strain = !m ? '' : given ? `${m.animals?.species || 'mouse'}, ${given}` :
    'Strain not recorded in the study\'s details';
  const sub = [m?.id, o.stage, strain, mismatch ?? (day0 ? `Day 0: ${dateText(day0, true)}` : ''),
    prov ? `PMX Studio ${prov.app}, atlas ${prov.atlas.slice(0, 12)}${planSaid}` : '',
    `Exported ${isoToday()} by ${grok.shell.user?.login ?? ''}`].filter((x) => x);
  const groups = m ? ofGroups(m, o.groups) : [];
  const details = !m ? '' : '<h2>Study details</h2>' + [m.description ? `<p>${esc(m.description)}</p>` : '',
    htmlTable(['Owner', 'Project', 'Tags', 'Blinded entry'], [[m.owner, m.project, m.tags.join(', '),
      m.blinded ? 'yes' : 'no']]),
    // the animals as the details give them (species, strain, sex, age, weight, source)
    animalsLine(m.animals).filled ? htmlTable(['Animals'], [[animalsLine(m.animals).filled]]) : '',
    groups.length ? htmlTable(['Group', 'Compound', 'Dose', 'Route', 'Schedule'], groups.map((g) =>
      [g.group, g.compound, g.dose === null ? '' : `${entered(g.dose)} ${g.unit}`, g.route, g.schedule])) :
      ''].join('');
  const notes = !m?.notes.length ? '' : '<h2>Notebook</h2>' + htmlTable(['Day', 'Kind', 'Animal', 'Note', 'Who'],
    [...m.notes].sort((a, b) => a.day - b.day).map((n) => {
      const d = day0 ? dateOf(m.start, n.day) : null;
      return [`${n.day}${d ? ` (${dateText(d)})` : ''}`, n.kind, n.animal, n.text, n.who];
    }));
  // what the stage keeps off screen (text and tables, no charts): under its title
  const annex = (o.annex ?? []).map(([h, el]) => `<h2>${esc(h)}</h2>${el.outerHTML}`).join('');
  const log = !o.log?.length ? '' : '<h2>Log</h2>' + htmlTable(['When', 'Who', 'What'],
    o.log.map((l) => [l.when.slice(0, 16).replace('T', ' '), l.who, l.what]));
  // the app's own styles, so the panels read as on screen
  const styles = [...document.querySelectorAll('style, link[rel="stylesheet"]')].map((s) => s.outerHTML).join('\n');
  // (a cell breaks a word only when it alone overflows, break-word: "anywhere" let the table size a label column to
  // a few letters, "Randomi-zation" in the protocol)
  return `<!doctype html><html><head><meta charset="utf-8"><title>${esc(o.title)}</title>${styles}
<style>
html,body{height:auto;overflow:visible;background:#fff}body{margin:24px;font:13px/1.45 Roboto,system-ui,sans-serif}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:24px 0 8px}.pmx-report-sub{color:#555;margin-bottom:16px}
table{border-collapse:collapse;margin:8px 0}th,td{border:1px solid #bbb;padding:3px 6px;text-align:left}
th{background:#f0f0f0}.pmx-panel{break-inside:avoid}img{break-inside:avoid}
body *{max-width:100%!important;box-sizing:border-box}
html,body{width:auto!important;min-width:0!important}
body *{min-width:0!important}
td,th{white-space:normal!important;overflow-wrap:break-word;hyphens:manual}
table{table-layout:auto!important}
@media print{body{margin:0}}
</style></head><body>
<h1>${esc(o.title)}</h1><div class="pmx-report-sub">${sub.map(esc).join(' · ')}</div>
${page.outerHTML}${annex}${details}${notes}${log}
</body></html>`;
}
