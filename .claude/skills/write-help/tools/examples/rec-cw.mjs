import {film, W, H} from '../rec-lib.mjs';

// Never clicks SAVE AND APPLY: the sync block exists only in the page and is gone when the browser closes.
const setup = ``;
const pre = async (page) => {
  await page.addStyleTag({content: '.d4-balloon, .d4-balloon-container { display: none !important; }'});
};

await film('log-export-cw', setup, async (page, m) => {
  const bringIntoView = async (loc) => {
    for (let i = 0; i < 12; i++) {
      const b = await loc.boundingBox();
      if (b && b.y + b.height < H - 60 && b.y > 60) return;
      await page.mouse.wheel(0, b && b.y > 60 ? 120 : -120);
      await page.waitForTimeout(250);
    }
  };
  await page.waitForTimeout(600);
  const gear = page.locator('[name="Settings"]').first();
  await m.clickEl(gear, {after: 2500});
  const logger = page.locator('.grok-view :text-is("Logger")').first();
  await logger.waitFor({timeout: 20000});
  await m.moveTo(W - 200, H * 0.6);
  await bringIntoView(logger);
  await m.clickEl(logger, {after: 1000});
  const hdr = page.locator('.d4-accordion-pane-header:visible', {hasText: /^Log sync/}).first();
  await hdr.waitFor({timeout: 30000});
  await page.waitForTimeout(1200);
  await m.clickEl(hdr.locator('.d4-accordion-pane-header-title, span, div').first().or(hdr), {after: 1500});
  const add = page.locator(':text-is("ADD NEW SYNC BLOCK"), :text-is("Add new sync block")').first();
  await add.waitFor({timeout: 30000});
  await page.waitForTimeout(1500);
  await m.clickEl(add, {after: 1800});

  const row = (label) => page.locator('.grok-view .ui-input-root:visible', {has: page.locator('.ui-input-label', {hasText: new RegExp('^' + label + '$')})}).first();
  const cloud = row('Cloud').locator('select');
  await bringIntoView(cloud);
  await m.clickEl(cloud, {after: 400});
  for (let i = 0; i < 2; i++) { await page.keyboard.press('ArrowUp'); await page.waitForTimeout(250); }
  await page.keyboard.press('Enter');
  await page.waitForTimeout(1500);
  console.log('cloud', await cloud.inputValue());

  const levels = row('Levels');
  await bringIntoView(levels);
  for (const lvl of ['error', 'warning']) {
    const cb = levels.locator('label, .ui-input-bool, div', {hasText: new RegExp('^' + lvl + '$')}).first();
    await m.clickEl(cb, {after: 500});
  }
  const lg = row('Log Group').locator('input');
  await bringIntoView(lg);
  await m.clickEl(lg, {after: 200});
  await m.type('datagrok-logs', 80);
  const st = row('Stream').locator('input');
  await bringIntoView(st);
  await m.clickEl(st, {after: 200});
  await m.type('datlas', 80);
  await page.waitForTimeout(800);
  await m.hoverEl(page.locator(':text-is("SAVE AND APPLY"), :text-is("Save and apply")').first());
  await page.waitForTimeout(2200);
}, {pre, start: [500, 300], thumbAt: 0.98, noTooltips: true, colors: 128,
  out: (process.env.REC_OUT ?? 'out') + '/log-export-cw'});
