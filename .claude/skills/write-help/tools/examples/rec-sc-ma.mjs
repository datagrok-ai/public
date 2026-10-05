import {film} from '../rec-lib.mjs';

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

await film('scatter-plot-moving-average', setup, async (page, m) => {
  const row = (re) => page.locator('.property-grid-item:visible', {hasText: re}).first();
  await page.waitForTimeout(500);
  const vb = await page.locator('.d4-viewer').first().boundingBox();
  await m.moveTo(vb.x + vb.width - 200, vb.y + 20);
  await page.waitForTimeout(500);
  await m.clickEl(page.locator('.panel-titlebar [name="icon-font-icon-settings"]:visible').first(), {after: 1200});
  const bars = await page.evaluate(() => [...document.querySelectorAll('.splitbar-vertical')].filter(e => e.offsetParent).map(e => { const r = e.getBoundingClientRect(); return [r.x + r.width / 2, r.y + r.height / 2]; }));
  const sbar = bars.sort((p, q) => q[0] - p[0])[0];
  await m.drag(sbar[0], 330, sbar[0] - 80, 330, {steps: 14, stepMs: 40});
  await m.clickEl(page.locator('input.d4-search-input:visible').last(), {after: 200});
  await m.type('moving', 90);
  await page.waitForTimeout(600);
  await m.clickEl(page.locator('[name="prop-view-show-moving-average-line"]:visible'), {after: 1500});

  const sb = await row(/^Moving Average Window$/).locator('input[type="range"]').boundingBox();
  const at = (v) => sb.x + 7 + (sb.width - 14) * (v - 2) / 98;
  const sy = sb.y + sb.height / 2;
  await m.drag(at(10), sy, at(100), sy, {steps: 34, stepMs: 60});
  await page.waitForTimeout(600);
  await m.drag(at(100), sy, at(30), sy, {steps: 20, stepMs: 60});
  await page.waitForTimeout(900);
  console.log('window', await page.evaluate(prop('movingAverageWindow')));

  await m.clickEl(page.locator('[name="div-column-combobox-x"]'), {after: 900});
  await m.type('STARTED', 110);
  await page.waitForTimeout(600);
  const xp = await page.evaluate(() => { const r = document.activeElement.getBoundingClientRect(); return [r.x, r.y]; });
  console.log('x picker', xp);
  await m.click(xp[0] + 40, xp[1] + 48, {after: 1500});
  console.log('x', await page.evaluate(prop('xColumnName')));

  const sb2 = await row(/^Moving Average Window$/).locator('input[type="range"]').boundingBox();
  const at2 = (v) => sb2.x + 7 + (sb2.width - 14) * (v - 2) / 98;
  await m.drag(at2(30), sy, at2(5), sy, {steps: 18, stepMs: 60});
  await page.waitForTimeout(900);
  console.log('window2', await page.evaluate(prop('movingAverageWindow')));
  const valueCell = async (re) => { const b = await row(re).boundingBox(); return [b.x + b.width - 60, b.y + b.height / 2]; };
  let [ux, uy] = await valueCell(/^Moving Average Window Unit/);
  await m.click(ux, uy, {after: 400});
  for (const unit of ['Days', 'Weeks', 'Months']) {
    const steps = unit === 'Days' ? 2 : 1;
    for (let i = 0; i < steps; i++) { await page.keyboard.press('ArrowDown'); await page.waitForTimeout(150); }
    await page.keyboard.press('Enter');
    await page.waitForTimeout(1800);
    console.log(unit, '->', await page.evaluate(prop('movingAverageWindowUnit')));
    if (unit !== 'Months') { [ux, uy] = await valueCell(/^Moving Average Window Unit/); await m.click(ux, uy, {after: 400}); }
  }
  await m.clickEl(page.locator('[name="prop-view-show-moving-average-deviation"]:visible'), {after: 1500});
  console.log('deviation', await page.evaluate(prop('showMovingAverageDeviation')));
  const close = page.locator('.grok-prop-panel').locator('xpath=ancestor::*[contains(@class,"panel-base")]').locator('[name="icon-times"], .grok-font-icon-close, i.fa-times').first();
  await m.clickEl(close, {after: 1000});

  await m.clickEl(page.locator('[name="div-column-combobox-color"]'), {after: 900});
  await m.type('SEX', 130);
  await page.waitForTimeout(600);
  const pick = await page.evaluate(() => { const r = document.activeElement.getBoundingClientRect(); return [r.x, r.y]; });
  await m.click(pick[0] + 40, pick[1] + 64, {after: 800});
  await page.waitForTimeout(1500);
  console.log('color', await page.evaluate(prop('colorColumnName')));

  const item = (t) => page.locator('.d4-legend-item:visible', {hasText: new RegExp('^' + t)}).first();
  await m.clickEl(item('M'), {after: 1500});
  await m.clickEl(item('M'), {after: 1500});
  await m.moveTo(vb.x + vb.width / 2, vb.y + vb.height - 60);
  await page.waitForTimeout(2200);
}, {pre, start: [500, 300], thumbAt: 0.98, noTooltips: true, colors: 64, out: (process.env.REC_OUT ?? 'out') + '/scatter-plot-moving-average'});
