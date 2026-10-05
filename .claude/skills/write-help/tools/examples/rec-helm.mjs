import {film} from '../rec-lib.mjs';

const setup = `
const df = await grok.data.files.openTable('System:AppData/Helm/samples/helm-showcase.csv');
await grok.data.detectSemanticTypes(df);
const tv = grok.shell.addTableView(df);
await new Promise(r => setTimeout(r, 2500));
tv.grid.col('HELM').width = 330;
`;
const pre = async (page) => {
  await page.evaluate(() => grok.shell.tv.grid.scrollToCell('HELM', 6));
  await page.waitForTimeout(1200);
  await page.addStyleTag({content: '.d4-balloon, .d4-balloon-container { display: none !important; }'});
};
const cellValue = `grok.shell.t.col('HELM').get(7)`;

await film('helm-editor', setup, async (page, m) => {
  const tid = (s) => page.locator(`[data-testid="${s}"]`).first();
  await page.waitForTimeout(600);
  const r = await page.evaluate(() => { const g = grok.shell.tv.grid; const b = g.cell('HELM', 7).bounds; const rc = g.canvas.getBoundingClientRect(); return [rc.left + b.midX, rc.top + b.midY]; });
  console.log('cell at', r);
  await m.dblclick(r[0], r[1]);
  await page.waitForSelector('[data-testid="editor-svg"]', {timeout: 30000});
  await page.waitForTimeout(2500);

  await m.clickEl(page.locator('[data-testid^="canvas-atom-"]', {hasText: /^E$/}).first(), {after: 900});
  await m.clickEl(tid('palette-tab-PEPTIDE'), {after: 700});
  await m.clickEl(tid('palette-search'), {after: 200});
  await m.type('W', 150);
  await page.waitForTimeout(800);
  await m.clickEl(tid('palette-tile-W'), {after: 1600});

  for (const [tab, wait] of [['tab-helm', 2000], ['tab-properties', 2200], ['tab-extra-molecular-structure', 5500], ['tab-extra-composition-analysis', 2500]]) {
    const t = tid(tab);
    if (!(await t.count())) { console.log('no tab', tab); continue; }
    await m.clickEl(t, {after: wait});
  }
  const ok = page.locator('[name="button-OK"]:visible').last();
  await m.clickEl(ok, {after: 1500});
  console.log('cell', await page.evaluate(cellValue));
  await m.moveTo(r[0] + 150, r[1] + 150);
  await page.waitForTimeout(2200);
}, {pre, start: [600, 400], thumbAt: 0.62, noTooltips: true, colors: 128, out: (process.env.REC_OUT ?? 'out') + '/helm-editor'});
