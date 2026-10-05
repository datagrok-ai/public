import {film, W, H} from '../rec-lib.mjs';

// Saves one project named NAME on the stand; run 'node dd-projects.mjs delete' after every take.
// Point TREE and QUERY at a parameterized query the stand has (Browse path to its connection, and the query label).
const TREE = ['Databases', 'Postgres', 'Northwind'];
const QUERY = 'customers in @country';
const NAME = 'Customers by country';
const setup = `grok.shell.windows.showBrowse = true;`;
const node = (page, t) => page.locator('.d4-tree-view-group-label:visible, .d4-tree-view-item-label:visible', {hasText: new RegExp('^' + t + '$')}).first();
const query = (page) => page.locator('.d4-tree-view-item-label:visible, .d4-tree-view-group-label:visible', {hasText: QUERY}).first();

// Not recorded: open the tree down to the query.
const pre = async (page) => {
  await page.addStyleTag({content: '.d4-balloon, .d4-balloon-container { display: none !important; }'});
  for (const t of TREE) {
    await node(page, t).locator('xpath=..').locator('.d4-tree-view-tri').first().click();
    await page.waitForTimeout(4000);
  }
  await query(page).scrollIntoViewIfNeeded();
  await page.mouse.move(700, 300);
  await page.waitForTimeout(1000);
};

await film('dynamic-dashboards', setup, async (page, m) => {
  const wheelTo = async (loc) => {
    for (let i = 0; i < 20; i++) {
      const b = await loc.boundingBox();
      if (b && b.y > 70 && b.y + b.height < H - 60) return;
      await m.moveTo(200, 350);
      await page.mouse.wheel(0, !b || b.y > 70 ? 100 : -100);
      await page.waitForTimeout(200);
    }
  };
  await page.waitForTimeout(1000);
  await m.clickEl(query(page), {button: 'right', after: 900});
  await m.clickEl(page.locator('.d4-menu-item-label:visible').getByText('Run', {exact: true}).first(), {after: 2500});
  const pin = page.locator('.d4-dialog input').first();
  await pin.waitFor({timeout: 20000});
  await m.clickEl(pin, {after: 200});
  await page.keyboard.press('Control+A');
  await m.type('USA', 120);
  await m.clickEl(page.locator('.d4-dialog [name="button-OK"]').first(), {after: 1000});
  await page.locator('.ui-btn:has-text("REFRESH")').first().waitFor({timeout: 30000});
  await page.waitForTimeout(2500);
  const country = page.locator('.d4-accordion-pane:has(.d4-accordion-pane-header:has-text("Source")) input').first();
  const refresh = page.locator('.ui-btn:has-text("REFRESH")').first();

  await m.clickEl(page.locator('[name="icon-pie-chart"]:visible').first(), {after: 2500});

  await m.clickEl(page.locator('[name="button-Save"]').first(), {after: 2500});
  const title = page.locator('.d4-dialog [contenteditable="true"], .d4-dialog input:visible').first();
  await m.clickEl(title, {after: 200});
  await page.keyboard.press('Control+A');
  await m.type(NAME, 70);
  await page.waitForTimeout(600);
  await m.hoverEl(page.locator('.d4-dialog :text-is("Data sync")').first());
  await page.waitForTimeout(1500);
  await m.clickEl(page.locator('.d4-dialog [name="button-OK"]').first(), {after: 4000});
  console.log('saved');

  const share = page.locator('.d4-dialog [name="button-CANCEL"]').first();
  if (await share.count()) { await page.waitForTimeout(1200); await m.clickEl(share, {after: 1200}); }
  const tab = page.locator('.tab-handle:visible:has-text("USA")').first();
  await m.hoverEl(tab);
  await page.waitForTimeout(400);
  await m.clickEl(tab.locator('[class*=close]').first(), {after: 1500});

  await wheelTo(node(page, 'Dashboards'));
  await m.clickEl(node(page, 'Dashboards'), {after: 4000});
  const card = page.locator('.grok-gallery-grid-item:has-text("' + NAME + '"), .d4-gallery-card:has-text("' + NAME + '"), :text-is("' + NAME + '")').first();
  const [cx, cy] = await m.center(card);
  await m.moveTo(cx, cy);
  await page.waitForTimeout(300);
  await page.mouse.dblclick(cx, cy, {delay: 60});
  await page.locator('.ui-btn:has-text("REFRESH")').first().waitFor({timeout: 30000});
  await page.waitForTimeout(3000);
  await m.clickEl(country, {after: 200});
  await page.keyboard.press('Control+A');
  await m.type('France', 120);
  await m.clickEl(refresh, {after: 3500});
  await m.moveTo(700, 420);
  await page.waitForTimeout(2200);
}, {pre, start: [700, 300], thumbAt: 0.98, noTooltips: true, colors: 128,
  out: (process.env.REC_OUT ?? 'out') + '/dynamic-dashboards'});
