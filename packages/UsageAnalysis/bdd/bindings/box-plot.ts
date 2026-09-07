/* The steps only the box plot needs: its value-axis zoom through the vertical range slider, the
   viewport it reports, and a double-click on plot space with no marker under it. Everything else
   in the box plot features is the library's `viewers` tier (`grok-bdd list-steps`). */
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

const Y_SLIDER = 'svg[type="range-slider"][name="y-slider"]';

async function viewportHeight(page: Page, target: ElementRef): Promise<number> {
  const loc = await viewers.viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.viewerOf(el).viewport.height);
}

/** Drags the top handle of the value-axis slider down by 40% of the slider, and remembers the
 * full range for the range assertions below. */
export const zoomValueAxis = When('user zooms into the value axis of {widget}', async (page: Page, target: ElementRef) => {
  const loc = await viewers.viewerLocator(page, target);
  const full = await viewportHeight(page, target);
  await loc.evaluate((el, h) => { (window as any).__bdd.viewerOf(el).root.__fullRange = h; }, full);
  const slider = loc.locator(Y_SLIDER);
  const track = await slider.boundingBox();
  const top = await slider.locator('[name="max-handle"]').boundingBox();
  if (!track || !top)
    throw new Error(`${target.phrase}: no value-axis range slider`);
  const x = top.x + top.width / 2;
  const y = top.y + top.height / 2;
  await page.mouse.move(x, y);
  await page.mouse.down();
  await page.mouse.move(x, y + track.height * 0.4, {steps: 4});
  await page.mouse.up();
}, {tier: 'ui'});

export const narrowedRange = Then('{widget} should show a narrowed value range', async (page: Page, target: ElementRef) => {
  const loc = await viewers.viewerLocator(page, target);
  const full: number = await loc.evaluate((el) => (window as any).__bdd.viewerOf(el).root.__fullRange);
  await expect.poll(() => viewportHeight(page, target), {message: 'the value range did not narrow'}).toBeLessThan(full * 0.95);
});

export const fullRange = Then('{widget} should show the full value range again', async (page: Page, target: ElementRef) => {
  const loc = await viewers.viewerLocator(page, target);
  const full: number = await loc.evaluate((el) => (window as any).__bdd.viewerOf(el).root.__fullRange);
  await expect.poll(async () => Math.abs(await viewportHeight(page, target) - full) / full, {message: 'the value range is not back to full'}).toBeLessThan(0.01);
});

/** The top of the view area: values that high are rare, so nothing is under the pointer. */
export const doubleClickEmptySpace = When('user double-clicks on empty plot space of {widget}', async (page: Page, target: ElementRef) => {
  const view = await viewers.hitArea(page, target, 'view');
  await page.mouse.dblclick(view.x + view.width / 2, view.y + view.height * 0.05);
}, {tier: 'ui'});
