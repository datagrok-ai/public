/* The steps only the PC plot needs. Everything else in the pc-plot features is the library's
   `viewers` tier and the platform's data steps, read from the hit areas and readings the plot
   reports (`axis "AGE"`, `band "AGE" - "HEIGHT"`, `line of row <n>`, `range max handle "AGE"`,
   `axis order`, `lines drawn`, `range max of "AGE"` and the rest — see
   `core/client/d4/lib/src/viewers/pc_plot/CLAUDE.md`).
   The axis-order claim stays here until a second viewer has named axes; its slider drag and its
   row-membership claim are in the library now (`bindings/tiers/viewers/widgets.ts`). */
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {el, ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

const PLOT = (): ElementRef => el('pc plot viewer');

const settle = async (page: Page, target: ElementRef = PLOT()): Promise<void> => {
  const loc = await viewers.viewerLocator(page, target);
  await loc.evaluate((e) => (window as any).__bdd.settle(e, 300));
};

/** The axes in drawing order: what the plot reports (`axis order`) AND where it drew them — the
 * `axis <col>` strips left to right. A name list alone would echo the property. */
export const axesShouldBe = Then('the axes of {widget} should be {string}', async (page: Page, target: ElementRef, expected: string) => {
  const want = expected.split(',').map((s) => s.trim());
  let shown = '';
  await expect.poll(async () => {
    const order = String(await viewers.readValue(page, target, 'axis order'));
    if (order !== want.join(', ')) {
      shown = `the plot reports "${order}"`;
      return false;
    }
    const areas = await viewers.hitAreas(page, target);
    const xs = want.map((name) => areas[`axis "${name}"`]);
    const missing = want.filter((name, i) => xs[i] === undefined);
    if (missing.length > 0) {
      shown = `it draws no axis for ${missing.join(', ')}; it has: ${Object.keys(areas).filter((k) => k.startsWith('axis "')).join(', ') || 'none'}`;
      return false;
    }
    shown = want.map((name, i) => `${name} at ${Math.round(xs[i].x)}`).join(', ');
    return xs.every((b, i) => i === 0 || b.x > xs[i - 1].x);
  }, {timeout: 5000, message: `the axes of ${target.phrase} are not "${expected}": ${shown}`}).toBe(true);
}, {description: 'the reported drawing order, and the axis strips actually left to right in it'});

/** A reorder drag of a column label. The plot starts the platform's drag-and-drop from a debounced
 * mouse-down, so the step waits for the drag to have begun (`body.d4-drag`, what
 * `AppEvents.beginDragDrop` sets) instead of holding the button for a fixed time. */
export const dragAxisLabel = When('user drags the {string} axis label of pc plot viewer onto the {string} axis label',
  async (page: Page, from: string, to: string) => {
    const target = PLOT();
    const a = viewers.centerOf(await viewers.hitArea(page, target, `axis label "${from}"`, true));
    const b = viewers.centerOf(await viewers.hitArea(page, target, `axis label "${to}"`));
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move(a.x + (b.x - a.x) / 4, a.y);
    await page.locator('body.d4-drag').waitFor({state: 'attached', timeout: 5000});
    await page.mouse.move(b.x, b.y, {steps: 4});
    await page.mouse.up();
    await settle(page);
  }, {tier: 'ui', description: 'the dragged column takes the slot the other label sits in'});


