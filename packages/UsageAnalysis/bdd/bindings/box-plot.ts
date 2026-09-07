/* The steps only the box plot needs: its value-axis zoom through the vertical range slider, and
   clicks on plot space with no marker under them. Everything else in the box plot features is the
   library's `viewers` tier and the platform's data steps (`grok-bdd list-steps`). */
import {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import {ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

const Y_SLIDER = 'svg[type="range-slider"][name="y-slider"]';

/** Drags the top handle of the value-axis slider down by 40% of the slider; the range before is
 * the snapshot "a narrower value range than before" reads. The slider lays itself out when the
 * pointer enters the axis strip, so the pointer goes there first. */
export const zoomValueAxis = When('user zooms into the value axis of {widget}', async (page: Page, target: ElementRef) => {
  const axis = viewers.centerOf(await viewers.hitArea(page, target, 'y axis', true));
  await page.mouse.move(axis.x - 3, axis.y);
  await page.mouse.move(axis.x, axis.y);
  const loc = await viewers.viewerLocator(page, target);
  const slider = loc.locator(Y_SLIDER);
  const track = await slider.boundingBox();
  // the handle on top is the max handle, or the min handle on an inverted axis
  const handles = [await slider.locator('[name="max-handle"]').boundingBox(), await slider.locator('[name="min-handle"]').boundingBox()];
  const top = handles.filter((h) => h !== null).sort((a, b) => a!.y - b!.y)[0];
  if (!track || !top)
    throw new Error(`${target.phrase}: no value-axis range slider`);
  const x = top.x + top.width / 2;
  const y = top.y + top.height / 2;
  await page.mouse.move(x, y);
  await page.mouse.down();
  await page.mouse.move(x, y + track.height * 0.4, {steps: 4});
  await page.mouse.up();
}, {tier: 'ui'});

/** The top of the view area: values that high are rare, so nothing is under the pointer. */
const emptySpace = async (page: Page, target: ElementRef): Promise<{x: number; y: number}> => {
  const view = await viewers.hitArea(page, target, 'view', true);
  return {x: view.x + view.width / 2, y: view.y + view.height * 0.05};
};

export const clickEmptySpace = When('user clicks on empty plot space of {widget}', async (page: Page, target: ElementRef) => {
  const p = await emptySpace(page, target);
  await page.mouse.click(p.x, p.y);
}, {tier: 'ui', description: 'plot space with no marker under it — a click there clears the selection'});

export const doubleClickEmptySpace = When('user double-clicks on empty plot space of {widget}', async (page: Page, target: ElementRef) => {
  const p = await emptySpace(page, target);
  await page.mouse.dblclick(p.x, p.y);
}, {tier: 'ui', description: 'resets the view'});
