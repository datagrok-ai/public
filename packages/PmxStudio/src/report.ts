/* What a stage exports, from the view's own menu (Export), as other Datagrok apps do: the report, a page to print or
   save as PDF from the browser's print dialog (the stage's panels as shown, each chart as an image, then the study's
   details, notes and log), and the stage's files. Datagrok's viewers save themselves as PNG from their own menu; the
   echarts charts get the same two items on a right click. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as echarts from 'echarts/core';

import type {LogEntry} from './study-file';
import {StudyMeta, dateOf, dateText} from './study-meta';

/** A menu item: its name, what it does and its group (Export unless named). */
export type Exports = [string, () => void | Promise<void>, string?][];

const failed = (e: any) => grok.shell.error(String(e?.message ?? e));

/** The view's menu: the stage's exports in an Export group, its other items in theirs (View: the layout). */
export function exportMenu(view: DG.ViewBase, items: Exports): void {
  const menu = DG.Menu.create();
  for (const name of [...new Set(['Export', ...items.map((x) => x[2] ?? 'Export')])]) {
    const group = menu.group(name);
    for (const [label, run] of items.filter((x) => (x[2] ?? 'Export') === name))
      group.item(label, () => void Promise.resolve().then(run).catch(failed));
    group.endGroup();
  }
  view.ribbonMenu = menu;
}

/** Saves a text or a data URL as a file. */
export function download(name: string, content: string, type = 'text/plain'): void {
  const a = document.createElement('a');
  a.href = content.startsWith('data:') ? content : URL.createObjectURL(new Blob([content], {type}));
  a.download = name;
  a.click();
}

const pngOf = (chart: echarts.ECharts) => chart.getDataURL({type: 'png', pixelRatio: 2, backgroundColor: '#fff'});

/** An echarts chart's right-click menu: Save as PNG, Copy as image, named after its panel. */
export function chartMenu(host: HTMLElement): void {
  host.addEventListener('contextmenu', (e) => {
    const chart = echarts.getInstanceByDom(host);
    if (!chart)
      return;
    e.preventDefault();
    const name = host.closest<HTMLElement>('.pmx-panel')?.dataset.panel ?? 'figure';
    DG.Menu.popup()
      .item('Save as PNG', () => download(`${name}.png`, pngOf(chart)))
      .item('Copy as image', () => void fetch(pngOf(chart)).then((r) => r.blob())
        .then((b) => navigator.clipboard.write([new ClipboardItem({'image/png': b})])).catch(failed))
      .show({causedBy: e});
  });
}

/** A chart as an image: echarts' own, a Datagrok viewer drawn by the platform's html2canvas; null: not drawn. */
async function imageOf(el: HTMLElement): Promise<string | null> {
  const chart = echarts.getInstanceByDom(el);
  if (chart)
    return pngOf(chart);
  try {
    await DG.Utils.loadJsCss(['/js/common/html2canvas.min.js']);
    const h2c = (window as any).html2canvas;
    return (await h2c(el, {backgroundColor: '#fff', scale: 2, logging: false})).toDataURL('image/png');
  } catch {
    return null;
  }
}

const esc = (s: unknown) => String(s ?? '').replace(/[&<>"]/g, (c) =>
  ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c]!);
const table = (head: string[], rows: unknown[][]) => `<table><tr>${head.map((h) => `<th>${esc(h)}</th>`).join('')}` +
  `</tr>${rows.map((r) => `<tr>${r.map((c) => `<td>${esc(c)}</td>`).join('')}</tr>`).join('')}</table>`;

/** The report: `content` (the stage's results) as shown, its charts as images, after a title; then the study's
 *  details, its notes by day and its log. */
export async function reportHtml(o: {title: string, stage: string, meta: StudyMeta | null, content: HTMLElement,
  log?: LogEntry[], annex?: [string, HTMLElement][]}): Promise<string> {
  const charts = [...o.content.querySelectorAll<HTMLElement>('[_echarts_instance_], .d4-viewer')]
    .filter((el) => el.offsetParent !== null && !el.parentElement?.closest('[_echarts_instance_], .d4-viewer'));
  const images = await Promise.all(charts.map(imageOf));
  charts.forEach((el, i) => el.dataset.pmxFigure = String(i));
  const page = o.content.cloneNode(true) as HTMLElement;
  charts.forEach((el) => delete el.dataset.pmxFigure);
  for (const el of page.querySelectorAll<HTMLElement>('[data-pmx-figure]')) {
    const src = images[Number(el.dataset.pmxFigure)];
    const img = document.createElement(src ? 'img' : 'p');
    if (src)
      Object.assign(img, {src, alt: 'figure'});
    else
      img.textContent = '(this figure could not be drawn)';
    img.style.width = '100%';
    el.replaceWith(img);
  }
  // what is for working, not for reading
  page.querySelectorAll('button, input, select, textarea, [hidden], .pmx-noprint, .ui-input-root')
    .forEach((e) => e.remove());

  const m = o.meta;
  const day0 = dateOf(m?.start ?? null, 0);
  const sub = [m?.id, o.stage, day0 ? `Day 0: ${dateText(day0, true)}` : '',
    `Exported ${new Date().toISOString().slice(0, 10)} by ${grok.shell.user?.login ?? ''}`].filter((x) => x);
  const details = !m ? '' : '<h2>Study details</h2>' + [m.description ? `<p>${esc(m.description)}</p>` : '',
    table(['Owner', 'Project', 'Tags', 'Blinded entry'], [[m.owner, m.project, m.tags.join(', '),
      m.blinded ? 'yes' : 'no']]),
    m.groups.length ? table(['Group', 'Compound', 'Dose', 'Route', 'Schedule'], m.groups.map((g) =>
      [g.group, g.compound, g.dose === null ? '' : `${g.dose} ${g.unit}`, g.route, g.schedule])) : ''].join('');
  const notes = !m?.notes.length ? '' : '<h2>Notebook</h2>' + table(['Day', 'Kind', 'Animal', 'Note', 'Who'],
    [...m.notes].sort((a, b) => a.day - b.day).map((n) => {
      const d = dateOf(m.start, n.day);
      return [`${n.day}${d ? ` (${dateText(d)})` : ''}`, n.kind, n.animal, n.text, n.who];
    }));
  // what the stage keeps off screen (text and tables, no charts): under its title
  const annex = (o.annex ?? []).map(([h, el]) => `<h2>${esc(h)}</h2>${el.outerHTML}`).join('');
  const log = !o.log?.length ? '' : '<h2>Log</h2>' + table(['When', 'Who', 'What'],
    o.log.map((l) => [l.when.slice(0, 16).replace('T', ' '), l.who, l.what]));
  // the app's own styles, so the panels read as on screen
  const styles = [...document.querySelectorAll('style, link[rel="stylesheet"]')].map((s) => s.outerHTML).join('\n');
  return `<!doctype html><html><head><meta charset="utf-8"><title>${esc(o.title)}</title>${styles}
<style>
html,body{height:auto;overflow:visible;background:#fff}body{margin:24px;font:13px/1.45 Roboto,system-ui,sans-serif}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:24px 0 8px}.pmx-report-sub{color:#555;margin-bottom:16px}
table{border-collapse:collapse;margin:8px 0}th,td{border:1px solid #bbb;padding:3px 6px;text-align:left}
th{background:#f0f0f0}.pmx-panel{break-inside:avoid}img{break-inside:avoid}
@media print{body{margin:0}}
</style></head><body>
<h1>${esc(o.title)}</h1><div class="pmx-report-sub">${sub.map(esc).join(' · ')}</div>
${page.outerHTML}${annex}${details}${notes}${log}
</body></html>`;
}
