/* The gestures only the scatter plot needs, and only because the library has no phrase for them
   yet: the lasso polygon the Lasso Tool selects with, the drag of a range-slider handle over an
   axis strip, and a pick in an on-viewer column selector (the library's `selects` takes the
   selector's hidden time-unit `<select>` first). Everything else in the scatter plot features is
   the library's `viewers` tier and the platform's data steps, read from the hit areas and readings
   the plot reports (`marker of row <n>`, `empty space`, `rows shown`, `x axis min` and the rest —
   see `core/client/d4/lib/src/viewers/scatterplot/CLAUDE.md`). */
import {expect, Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import {el, viewers} from '@datagrok-libraries/bdd/runtime';

const PLOT = () => el('scatter plot viewer');

const settle = async (page: Page): Promise<void> => {
  const loc = await viewers.viewerLocator(page, PLOT());
  await loc.evaluate((e) => (window as any).__bdd.settle(e, 300));
};

/** A Shift-drag along a closed polygon inside the area — what the platform reads as a lasso while
 * the Lasso Tool is on. The pointer walks each leg in several moves: the browser delivers pointer
 * moves frame-aligned, and the selector builds the polygon from the moves it sees. */
export const dragLasso = When('user drags a lasso over the {string} area of scatter plot viewer',
  async (page: Page, area: string) => {
    const b = await viewers.hitArea(page, PLOT(), area, true);
    const points = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [0.25, 0.25]]
      .map(([fx, fy]) => ({x: b.x + b.width * fx, y: b.y + b.height * fy}));
    await page.keyboard.down('Shift');
    await page.mouse.move(points[0].x, points[0].y);
    await page.mouse.down();
    for (const p of points.slice(1))
      await page.mouse.move(p.x, p.y, {steps: 8});
    await page.mouse.up();
    await page.keyboard.up('Shift');
  }, {tier: 'ui', description: 'the Lasso Tool must be on; the polygon covers the middle half of the area'});

/** The axis range slider is hover-revealed over its axis strip, so the pointer goes there first and
 * the step waits for the slider to show: a DOM read issued right after the move can otherwise run
 * before the move's handler, with the handles still at their unrendered positions. */
export const dragRangeHandle = When('user drags the {word} handle of the {word} range slider of scatter plot viewer by {int} pixels',
  async (page: Page, handle: string, axis: string, px: number) => {
    const target = PLOT();
    if (handle !== 'min' && handle !== 'max')
      throw new Error(`a range slider has a min and a max handle, not a "${handle}" one`);
    if (axis !== 'x' && axis !== 'y')
      throw new Error(`the scatter plot has an x and a y range slider, not a "${axis}" one`);
    const strip = viewers.centerOf(await viewers.hitArea(page, target, `${axis} axis`, true));
    await page.mouse.move(strip.x - 3, strip.y);
    await page.mouse.move(strip.x, strip.y);
    const loc = await viewers.viewerLocator(page, target);
    const slider = loc.locator(`svg[name="${axis}-slider"]`);
    await slider.waitFor({state: 'visible'});
    const box = await slider.locator(`[name="${handle}-handle"]`).boundingBox();
    if (!box)
      throw new Error(`${target.phrase}: the ${axis} range slider shows no ${handle} handle`);
    const x = box.x + box.width / 2;
    const y = box.y + box.height / 2;
    await page.mouse.move(x, y);
    await page.mouse.down();
    if (axis === 'x')
      await page.mouse.move(x + px, y, {steps: 4});
    else
      await page.mouse.move(x, y + px, {steps: 4});
    await page.mouse.up();
    await settle(page);
  }, {tier: 'ui', description: 'a positive count moves the handle right (x) or down (y)'});

/** The on-viewer column selector (`div-column-combobox-<property>`): a mouse-down on its caption
 * opens the column grid, typing opens the grid's search, Enter takes the column typed. The Color
 * and Size selectors are hover-revealed — hover the viewer first. */
export const pickColumn = When('user picks {string} in the {word} column selector of scatter plot viewer',
  async (page: Page, column: string, which: string) => {
    const target = PLOT();
    const loc = await viewers.viewerLocator(page, target);
    const selector = loc.locator(`[name="div-column-combobox-${which.toLowerCase()}"]`);
    await selector.waitFor({state: 'visible', timeout: 5000});
    const box = await selector.boundingBox();
    if (!box)
      throw new Error(`${target.phrase}: the ${which} column selector has no box`);
    await page.mouse.move(box.x + Math.min(10, box.width / 2), box.y + box.height / 2);
    await page.mouse.down();
    await page.mouse.up();
    await page.locator('.d4-column-grid').last().waitFor({state: 'visible', timeout: 5000});
    await page.keyboard.type(column);
    await page.keyboard.press('Enter');
    await expect(selector.locator('.d4-column-selector-column')).toHaveText(column, {timeout: 5000});
    await settle(page);
  }, {tier: 'ui', description: 'X, Y, Color or Size — the selector is named by the property it binds'});
