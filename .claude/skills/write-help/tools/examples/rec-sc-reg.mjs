import {film, walkMenu} from '../rec-lib.mjs';

const setup = `
const full = await grok.data.getDemoTable('demog.csv');
const df = full.clone(DG.BitSet.create(full.rowCount, (i) => i % 8 === 0));
df.name = 'demog';
const tv = grok.shell.addTableView(df);
await new Promise(r => setTimeout(r, 500));
tv.dockManager.close(tv.grid.root);
tv.addViewer('Scatter plot', {xColumnName: 'HEIGHT', yColumnName: 'WEIGHT'});
`;
const prop = (p) => `grok.shell.tv.viewers.find(v => v.type === 'Scatter plot').props.${p}`;
const pre = async (page) => {
  await page.addStyleTag({content: '.d4-balloon, .d4-balloon-container { display: none !important; }'});
};

await film('scatter-plot-regression-line', setup, async (page, m) => {
  await page.waitForTimeout(700);
  await m.click(200, 170, {button: 'right', after: 700});
  await walkMenu(m, ['Tools', 'Show Regression Line']);
  for (let i = 0; i < 3; i++) { await page.keyboard.press('Escape'); await page.waitForTimeout(120); }
  await page.waitForTimeout(1800);
  console.log('regression', await page.evaluate(prop('showRegressionLine')));

  await m.clickEl(page.locator('[name="div-column-combobox-color"]'), {after: 900});
  await m.type('SEX', 130);
  await page.waitForTimeout(600);
  const pick = await page.evaluate(() => { const r = document.activeElement.getBoundingClientRect(); return [r.x, r.y]; });
  await m.click(pick[0] + 40, pick[1] + 64, {after: 800});
  await page.waitForTimeout(1800);
  console.log('color', await page.evaluate(prop('colorColumnName')));

  const item = (t) => page.locator('.d4-legend-item:visible', {hasText: new RegExp('^' + t)}).first();
  await m.clickEl(item('M'), {after: 1500});
  await m.clickEl(item('M'), {after: 1500});

  const vb = await page.locator('.d4-viewer').first().boundingBox();
  await m.moveTo(vb.x + vb.width - 200, vb.y + 20);
  await page.waitForTimeout(500);
  await m.clickEl(page.locator('.panel-titlebar [name="icon-font-icon-settings"]:visible').first(), {after: 1200});
  const bars = await page.evaluate(() => [...document.querySelectorAll('.splitbar-vertical')].filter(e => e.offsetParent).map(e => { const r = e.getBoundingClientRect(); return [r.x + r.width / 2, r.y + r.height / 2]; }));
  const sb = bars.sort((p, q) => q[0] - p[0])[0];
  await m.drag(sb[0], 330, sb[0] - 80, 330, {steps: 14, stepMs: 40});
  const search = page.locator('.grok-prop-panel input.d4-search-input:visible, input.d4-search-input:visible').last();
  await m.clickEl(search, {after: 200});
  await m.type('equation', 90);
  await page.waitForTimeout(600);
  await m.clickEl(page.locator('[name="prop-view-show-regression-line-equation"]:visible'), {after: 1600});
  console.log('equation', await page.evaluate(prop('showRegressionLineEquation')));
  await m.clickEl(search, {after: 200});
  await page.keyboard.press('Control+A');
  await m.type('correlation', 90);
  await page.waitForTimeout(600);
  await m.clickEl(page.locator('[name="prop-view-show-pearson-correlation"]:visible'), {after: 1600});
  console.log('pearson', await page.evaluate(prop('showPearsonCorrelation')));
  const close = page.locator('.grok-prop-panel').locator('xpath=ancestor::*[contains(@class,"panel-base")]').locator('[name="icon-times"], .grok-font-icon-close, i.fa-times').first();
  await m.clickEl(close, {after: 800});
  await m.moveTo(vb.x + vb.width / 2, vb.y + vb.height - 60);
  await page.waitForTimeout(2200);
}, {pre, start: [500, 300], thumbAt: 0.98, noTooltips: true, colors: 64, out: (process.env.REC_OUT ?? 'out') + '/scatter-plot-regression-line'});
