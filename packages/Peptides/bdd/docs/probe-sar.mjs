import {createRequire} from 'node:module';
const require = createRequire('C:/Users/rizhi/Desktop/GROK-CORE/reddata/public/libraries/bdd/package.json');
const {chromium} = require('@playwright/test');

const storageState = 'C:/Users/rizhi/Desktop/GROK-CORE/reddata/public/packages/Peptides/bdd/e2e/.auth.json';
const browser = await chromium.launch({args: ['--disable-accelerated-2d-canvas']});
const page = await (await browser.newContext({storageState, viewport: {width: 1600, height: 1000}})).newPage();
page.on('pageerror', (e) => console.log('pageerror:', String(e).slice(0, 200)));
await page.goto('http://localhost:8888/', {waitUntil: 'domcontentloaded', timeout: 180000});
await page.locator('[name="Browse"]').waitFor({timeout: 180000});
const t0 = Date.now();
const lap = (m) => console.log(`${((Date.now() - t0) / 1000).toFixed(1)}s ${m}`);

await page.evaluate(async () => {
  grok.shell.windows.simpleMode = true;
  const df = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
  grok.shell.addTableView(df);
  await new Promise((r) => { const s = df.onSemanticTypeDetected.subscribe(() => { s.unsubscribe(); r(); }); setTimeout(r, 5000); });
});
lap('table open');
console.log(await page.evaluate(() => {
  const df = grok.shell.t;
  return {rows: df.rowCount, cols: df.columns.names(), semType: df.col('AlignedSequence').semType, units: df.col('AlignedSequence').getTag('units'), renderer: df.col('AlignedSequence').getTag('cell.renderer')};
}));
await page.evaluate(() => { grok.shell.windows.showContextPanel = true; grok.shell.o = grok.shell.t.col('AlignedSequence'); });
await page.locator('[name="pane-Peptides"]').waitFor({timeout: 60000});
lap('Peptides pane');
console.log('panes:', await page.evaluate(() => Array.from(document.querySelectorAll('.grok-prop-panel [name^="pane-"]')).map((e) => e.getAttribute('name')).join(' | ')));
const header = page.locator('[name="pane-Peptides"] .d4-accordion-pane-header');
if (!(await header.evaluate((e) => e.classList.contains('expanded'))))
  await header.click();
await page.locator('[name="button-Launch-SAR"]').waitFor({timeout: 60000});
lap('pane expanded');
console.log('pane names:', await page.evaluate(() => Array.from(document.querySelectorAll('[name="pane-Peptides"] [name]')).map((e) => `${e.tagName.toLowerCase()}[${e.getAttribute('name')}]`).join(' ')));
console.log('pane weblogo:', await page.evaluate(() => {
  const host = document.querySelector('[name="pane-Peptides"] .bio-wl-host');
  const root = host?.closest('[name^="viewer-"]') ?? host?.parentElement;
  const w = root ? DG.Widget.find(root) : null;
  return {host: !!host, rootName: root?.getAttribute('name'), widget: w?.constructor?.name, hasStatus: typeof w?.getWidgetStatus === 'function', canvases: document.querySelectorAll('[name="pane-Peptides"] canvas').length};
}));
await page.locator('[name="button-Launch-SAR"]').click();
await page.waitForFunction(() => Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']), null, {timeout: 120000});
lap('model');
await page.waitForFunction(() => {
  const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
  const types = Array.from(tv.viewers).map((v) => v.type);
  return ['Sequence Variability Map', 'Most Potent Residues', 'MCL', 'Logo Summary Table'].every((t) => types.includes(t));
}, null, {timeout: 180000});
lap('4 viewers');
await page.waitForTimeout(3000);
console.log(await page.evaluate(() => {
  const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
  const out = {view: tv.name, viewers: [], roots: [], ribbon: [], grids: []};
  for (const v of Array.from(tv.viewers)) {
    out.viewers.push({type: v.type, ctor: v.constructor?.name, rootName: v.root?.getAttribute('name'), status: typeof v.getWidgetStatus, pending: typeof v.isRenderPending, onRendered: typeof v.onRendered, found: DG.Widget.find(v.root)?.constructor?.name});
  }
  for (const e of document.querySelectorAll('[name^="viewer-"]')) {
    const r = e.getBoundingClientRect();
    out.roots.push(`${e.getAttribute('name')} ${Math.round(r.width)}x${Math.round(r.height)} in=${e.closest('[name^="viewer-"]:not([name="' + e.getAttribute('name') + '"])')?.getAttribute('name') ?? '-'}`);
  }
  for (const e of document.querySelectorAll('.d4-ribbon i, .d4-ribbon [name], .d4-ribbon-item'))
    out.ribbon.push(`${e.tagName.toLowerCase()} name=${e.getAttribute('name')} aria=${e.getAttribute('aria-label')} cls=${e.className}`);
  for (const e of document.querySelectorAll('.d4-grid')) {
    const w = DG.Widget.find(e);
    const st = w?.getWidgetStatus?.();
    out.grids.push({name: e.getAttribute('name'), inside: e.parentElement?.closest('[name^="viewer-"]')?.getAttribute('name'), widget: w?.constructor?.name, type: w?.type, areas: st ? Object.keys(st.hitAreas ?? {}).slice(0, 6) : null, values: st ? Object.keys(st.values ?? {}).slice(0, 12) : null});
  }
  return JSON.stringify(out, null, 1);
}));
console.log('svm inputs:', await page.evaluate(() => Array.from(document.querySelectorAll('[name="viewer-Sequence-Variability-Map"] [name], [name="viewer-Most-Potent-Residues"] [name], [name="viewer-Logo-Summary-Table"] [name]')).map((e) => `${e.tagName.toLowerCase()}[${e.getAttribute('name')}]`).join(' ')));
console.log('dock tabs:', await page.evaluate(() => Array.from(document.querySelectorAll('.dock-spawn-tab-handle-text, .tab-handle-text, .d4-tab-header, .panel-titlebar-text')).map((e) => e.textContent.trim()).filter(Boolean).join(' | ')));
console.log('model settings:', await page.evaluate(() => {
  const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
  const m = tv.dataFrame.temp['peptidesModel'];
  return JSON.stringify({settings: m.settings, cols: tv.dataFrame.columns.names(), mcl: m._mclCols, seqSpace: m._sequenceSpaceCols});
}));
// the SVM inner grid: which cells exist, a click on a cell through the grid's own status
console.log('svm grid status sample:', await page.evaluate(() => {
  const root = document.querySelector('[name="viewer-Sequence-Variability-Map"] .d4-grid');
  const w = DG.Widget.find(root);
  const st = w.getWidgetStatus();
  const areas = Object.entries(st.hitAreas).filter(([k]) => k.startsWith('cell ')).slice(0, 5);
  const vals = Object.entries(st.values).filter(([k]) => !k.startsWith('color of') && !k.startsWith('text of') && !k.startsWith('cell type') && !k.startsWith('column width'));
  return JSON.stringify({areas, vals, textSample: Object.entries(st.values).filter(([k]) => k.startsWith('text of')).slice(0, 8)});
}));
await page.screenshot({path: 'C:/Users/rizhi/AppData/Local/Temp/claude/c--Users-rizhi-Desktop-GROK-CORE-reddata/5064346e-d8d5-4cff-a425-aab38edc88c6/scratchpad/sar.png'});
await browser.close();
