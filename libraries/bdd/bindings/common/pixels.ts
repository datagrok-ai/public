/* Pixel claims that do not go through a viewer. The `viewers` tier reads a viewer's own canvas and
   the areas it reports, which is the better answer whenever the thing on screen IS a viewer of an
   open table view. These two steps are for what it cannot reach: the small multiples of a facet
   plot, which are not viewers, and a chart inside a Vue-rendered view, whose canvas answers
   `toDataURL` with a blank image and whose element a locator does not consider visible.

   Both work at the page level — every canvas inside the phrase for a colour count, and a screenshot
   of the phrase's rectangle for a change — so neither needs the element to be registered anywhere. */
import {Page} from '@playwright/test';
import {expect} from '../../src/runtime/patience.js';
import {Then, When} from '../../src/registry.js';
import type {ElementRef} from '../../src/runtime/args.js';
import {locate} from '../../src/runtime/locate.js';

/** Distinct colours over every canvas inside the phrase: alpha under half ignored, near-white and
 * near-black dropped (paper and axis text), greys dropped (gridlines), and the rest rounded to four
 * bits per channel so antialiasing does not read as a colour of its own. The same reduction the
 * hand-written spec used, as a step. */
export const canvasColors = Then('the canvases of {element} should be painted in at least {int} colors',
  async (page: Page, target: ElementRef, min: number) => {
    const loc = await locate(page, target);
    const count = (): Promise<number> => loc.evaluateAll((els) => {
      const buckets = new Set<string>();
      for (const el of els) {
        const canvases = el instanceof HTMLCanvasElement ? [el] : Array.from(el.querySelectorAll('canvas'));
        for (const c of canvases as HTMLCanvasElement[]) {
          const ctx = c.getContext('2d');
          if (!ctx || c.width === 0 || c.height === 0)
            continue;
          let img: ImageData;
          try {
            img = ctx.getImageData(0, 0, c.width, c.height);
          }
          catch {
            continue;
          }
          const d = img.data;
          for (let i = 0; i < d.length; i += 4) {
            const r = d[i], g = d[i + 1], b = d[i + 2];
            if (d[i + 3] < 128 || (r > 240 && g > 240 && b > 240) || (r < 30 && g < 30 && b < 30))
              continue;
            if (Math.abs(r - g) < 15 && Math.abs(g - b) < 15)
              continue;
            buckets.add(`${r >> 4}|${g >> 4}|${b >> 4}`);
          }
        }
      }
      return buckets.size;
    });
    await expect.poll(count, {message: `distinct colours over the canvases of ${target.phrase}`})
      .toBeGreaterThanOrEqual(min);
  }, {tier: 'ui', description: 'over every canvas the phrase contains, so a facet of small multiples answers as one'});

const pictures = new WeakMap<Page, Map<string, Buffer>>();

/** A screenshot of the element's rectangle, taken through the page rather than the element: a chart
 * inside a Vue view is painted where the locator says it is, even when the locator does not think
 * the element visible enough to screenshot itself. */
async function picture(page: Page, target: ElementRef): Promise<Buffer> {
  const loc = (await locate(page, target)).first();
  await loc.scrollIntoViewIfNeeded().catch(() => undefined);
  const box = await loc.boundingBox();
  if (!box || box.width < 1 || box.height < 1)
    throw new Error(`${target.phrase} has no rectangle on the page to picture`);
  return page.screenshot({clip: box});
}

export const takePicture = When('user takes a picture of {element}', async (page: Page, target: ElementRef) => {
  const taken = pictures.get(page) ?? new Map<string, Buffer>();
  pictures.set(page, taken);
  taken.set(target.phrase, await picture(page, target));
}, {tier: 'ui', description: 'the pixels of its rectangle now, to compare with after the next change'});

export const lookedDifferent = Then('{element} should look different', async (page: Page, target: ElementRef) => {
  const before = pictures.get(page)?.get(target.phrase);
  if (!before)
    throw new Error(`no picture of ${target.phrase} to compare with — "user takes a picture of it" first`);
  await expect.poll(async () => (await picture(page, target)).equals(before),
    {message: `the pixels of ${target.phrase} against the picture taken before`}).toBe(false);
});

export const lookedSame = Then('{element} should look the same', async (page: Page, target: ElementRef) => {
  const before = pictures.get(page)?.get(target.phrase);
  if (!before)
    throw new Error(`no picture of ${target.phrase} to compare with — "user takes a picture of it" first`);
  expect((await picture(page, target)).equals(before), `the pixels of ${target.phrase}`).toBe(true);
});
