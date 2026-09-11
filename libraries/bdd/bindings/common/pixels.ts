/* A pixel claim that does not go through one viewer. The `viewers` tier reads a viewer's own canvas
   and the areas it reports, which is the answer whenever the thing on screen IS a viewer — of a
   table view or, through `DG.Widget.find`, of any view. This step is for a claim no single viewer
   can make: the colours across the small multiples of a facet, which are several viewers. It works
   at the page level, over every canvas inside the phrase, so the element needs no registration. */
import {Page} from '@playwright/test';
import {expect} from '../../src/runtime/patience.js';
import {Then} from '../../src/registry.js';
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
