import {film} from '../rec-lib.mjs';

const setup = `
const df = await grok.data.getDemoTable('demog.csv');
const tv = grok.shell.addTableView(df);
await new Promise(r => setTimeout(r, 500));
const v = DG.Viewer.lineChart(df, {xColumnName: 'STARTED', yColumnNames: ['AGE', 'HEIGHT', 'WEIGHT']});
tv.dockManager.dock(v, DG.DOCK_TYPE.RIGHT, null, 'Line chart', 0.72);
`;
const prop = (p) => `grok.shell.tv.viewers.find(v => v.type === 'Line chart').props.${p}`;
const pre = async (page) => {
  await page.addStyleTag({content: '.d4-balloon, .d4-balloon-container { display: none !important; }'});
};

await film('line-chart-small', setup, async (page, m) => {
  const vb = await page.locator('[name="viewer-Line-chart"]').first().boundingBox();
  const axisY = vb.y + vb.height - 37;
  const left = vb.x + 46, right = vb.x + vb.width - 6;
  await m.moveTo((left + right) / 2, axisY - 2);
  await page.waitForTimeout(900);
  await m.drag(right - 2, axisY, right - 300, axisY, {steps: 30, stepMs: 45});
  await page.waitForTimeout(700);
  await m.drag(right - 300, axisY, right - 2, axisY, {steps: 30, stepMs: 45});
  await page.waitForTimeout(700);

  await m.dblclick(vb.x + vb.width / 2, vb.y + vb.height / 2);
  await page.waitForTimeout(1500);
  console.log('y cols', await page.evaluate(prop('yColumnNames')), await page.evaluate(prop('multiAxis')));

  await m.clickEl(page.locator('[name="add-split"]'), {after: 500});
  await m.type('SEX', 130);
  await page.waitForTimeout(600);
  const sp = await page.evaluate(() => { const r = document.activeElement.getBoundingClientRect(); return [r.x, r.y]; });
  await m.click(sp[0] + 40, sp[1] + 64, {after: 1500});
  console.log('split', await page.evaluate(prop('splitColumnNames')));

  const item = (t) => page.locator('.d4-legend-item:visible', {hasText: new RegExp('^' + t)}).first();
  await m.clickEl(item('F'), {after: 1500});
  await m.clickEl(item('M'), {after: 1500});
  await m.clickEl(item('M'), {after: 1500});

  const xb = await page.locator('[name^="div-column-combobox"]', {hasText: 'STARTED'}).last().boundingBox();
  await m.click(xb.x + xb.width - 30, xb.y + xb.height / 2, {after: 900});
  await m.type('AGE', 130);
  await page.waitForTimeout(600);
  const xp = await page.evaluate(() => { const r = document.activeElement.getBoundingClientRect(); return [r.x, r.y]; });
  await m.click(xp[0] + 40, xp[1] + 48, {after: 1200});
  console.log('x', await page.evaluate(prop('xColumnName')));
  await m.moveTo(vb.x + vb.width / 2, vb.y + vb.height / 2 + 40);
  await page.waitForTimeout(2200);
}, {pre, start: [420, 300], thumbAt: 0.98, noTooltips: false, colors: 64, fps: 10, size: [400, 300],
  crop: async (page) => { const b = await page.locator('[name="viewer-Line-chart"]').first().boundingBox(); const h = b.width * 3 / 4; console.log('viewer', JSON.stringify(b)); return [b.x, b.y + b.height - h, b.width, h]; },
  out: (process.env.REC_OUT ?? 'out') + '/line-chart-small'});
