/* Steps the "Work with data" how-to clips (features/howto/05-data) need and the library's vocabulary lacks.
   Each one is a gesture a person makes, written so that it films well: the pointer never flies across the
   screen for a test's sake. */
import type {Page} from '@playwright/test';
import {dataset, Given, When} from '@datagrok-libraries/bdd';
import {el, expect, gestures, menuNames, openTopMenu, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

// Two Northwind demo files that share a key, for joining: every order line names its product by ProductID.
// (categories.csv, customers.csv, orders.csv and suppliers.csv have unquoted commas and open with stray columns.)
dataset('northwind-products', {path: 'System:DemoFiles/northwind/products.csv',
  description: '77 Northwind products: ProductID, ProductName, CategoryID, UnitPrice, stock; the view is "products"'});
dataset('northwind-order-details', {path: 'System:DemoFiles/northwind/order-details.csv',
  description: '2155 Northwind order lines: OrderID, ProductID, UnitPrice, Quantity, Discount (view "order-details")'});

/** The `+` of a pivot table's tag row opens the column picker; the pointer then steps just beside the `+`,
 * clear of the picker (a row it rests on is previewed onto the chip), where the library's step sends it to
 * the page's corner — on film a flight across the video that pops the sidebar's toolbox out on the way. */
export const addToPivotRow = When('user picks {string} with the plus of the {string} row of pivot table viewer',
  async (page: Page, column: string, row: string) => {
    const target = el('pivot table viewer');
    const box = await viewers.hitArea(page, target, `add ${row}`, true);
    const c = viewers.centerOf(box);
    await page.mouse.click(c.x, c.y);
    const popup = page.locator('.d4-column-grid').filter({visible: true}).last();
    await popup.waitFor({state: 'visible', timeout: 5000}).catch(() => undefined);
    const picker = await popup.boundingBox().catch(() => null);
    const clear = (p: {x: number; y: number}): boolean => p.x > 0 && p.y > 0 && (!picker ||
      p.x < picker.x - 4 || p.x > picker.x + picker.width + 4 ||
      p.y < picker.y - 4 || p.y > picker.y + picker.height + 4);
    const beside = [{x: box.x - 16, y: c.y}, {x: c.x, y: box.y - 16}, {x: box.x + box.width + 16, y: c.y}];
    const away = beside.find(clear) ?? {x: 2, y: 2};
    await page.mouse.move(away.x, away.y);
    await gestures.pickInColumnGrid(page, column, `the "${row}" row`);
    await viewers.settle(page, target, 3000);
  }, {tier: 'ui', description: 'the + of a tag row of the pivot table, the column typed into its picker and taken, ' +
    'the pointer resting beside the +'});

/** A tour of a top menu group: the group opens as the library's walk opens it (left open), and the pointer
 * moves down onto the item's label, which the menu lights as it would for a person reading it. Nothing is
 * picked; the next item of the same group is reached by a move within the open dropdown, so it stays open. */
export const pointAtTopMenuItem = When('user points at {string} in the top menu', async (page: Page, path: string) => {
  const {segments, names} = menuNames(path);
  const leaf = page.locator(`[name="${names[names.length - 1]}"]`).filter({visible: true}).first();
  if (!await leaf.isVisible())
    await openTopMenu(page, segments.slice(0, -1).join(' > '), false);
  await leaf.waitFor({state: 'visible', timeout: pollMs(5000)});
  const label = leaf.locator(':scope > .d4-menu-item-label').first();
  const box = (await label.boundingBox()) ?? (await leaf.boundingBox())!;
  await page.mouse.move(box.x + Math.min(40, box.width / 2), box.y + box.height / 2, {steps: 6});
}, {tier: 'ui', description: 'opens the group (left open) and rests the pointer on the item\'s label; picks nothing'});

/** Find and Replace runs ReplaceValues through the platform's default result handling, and its scalar output
 * (the count of replaced values) opens the Variables pane on the right of the shell — not what the clip is about,
 * so a how-to closes it off camera (the defect is in the report of the 05-data clips). */
export const closeVariablesPane = Given('the variables pane is closed', async (page: Page) => {
  await page.evaluate(() => grok.shell.windows.showVariables = false);
  await expect.poll(() => page.evaluate(() => grok.shell.windows.showVariables),
    {message: 'the Variables pane closed'}).toBe(false);
}, {tier: 'api', description: 'View > Variables switched off through the JS API'});
