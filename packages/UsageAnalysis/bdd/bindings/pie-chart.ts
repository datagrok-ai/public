/* The steps only the pie chart needs, all read from what the chart reports about the frame it drew
   (`getWidgetStatus`): the wedges' shares and start angles for a sort claim, the distance between a
   wedge and the centre of the disc for Shift, and the swatch of the colour dialog a legend item
   opens. Everything else in the pie chart features is the library's `viewers` tier and the
   platform's data steps (`grok-bdd list-steps`).

   Two of these are not pie chart business and should be promoted to the library:
*/
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {el, ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

interface Slice {
  name: string;
  share: number;
  startAngle: number;
  outerRadius: number;
}

/** Every reading the viewer reports right now, keyed as it names them. */
async function readings(page: Page, target: ElementRef): Promise<Record<string, unknown>> {
  await viewers.installViewerRuntime(page);
  const loc = await viewers.viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.viewerOf(el).getWidgetStatus()?.values ?? {});
}

/** The wedges of the last frame, in the order they were drawn (start angles accumulate, so they
 * increase with the drawing order and never wrap). */
async function slicesOf(page: Page, target: ElementRef): Promise<Slice[]> {
  const values = await readings(page, target);
  const of = (prefix: string, name: string): number => Number(values[`${prefix} of "${name}"`]);
  const names = Object.keys(values).filter((k) => k.startsWith('share of "')).map((k) => k.slice('share of "'.length, -1));
  if (names.length === 0)
    throw new Error(`${target.phrase} reports no slices; it reports: ${Object.keys(values).join(', ') || 'nothing'}`);
  const slices = names.map((name) => ({name, share: of('share', name), startAngle: of('start angle', name), outerRadius: of('outer radius', name)}));
  slices.sort((a, b) => a.startAngle - b.startAngle);
  return slices;
}

async function expectOrder(page: Page, target: ElementRef, holds: (s: Slice[]) => boolean, what: string): Promise<void> {
  let shown = '';
  await expect.poll(async () => {
    const slices = await slicesOf(page, target);
    shown = slices.map((s) => `${s.name || '(empty)'} ${s.share.toFixed(1)}% at ${s.startAngle.toFixed(1)}°`).join(', ');
    return holds(slices);
  }, {timeout: 5000, message: `the slices of ${target.phrase} are not ordered by ${what}; in drawing order they are: ${shown}`}).toBe(true);
}

export const slicesByShareDescending = Then('the slices of {widget} should be ordered by share descending', (page: Page, target: ElementRef) =>
  expectOrder(page, target, (s) => s.every((x, i) => i === 0 || x.share <= s[i - 1].share + 1e-9), 'share, largest first'),
{description: 'the biggest wedge drawn first, every later one no bigger than the one before it'});

export const slicesByShareAscending = Then('the slices of {widget} should be ordered by share ascending', (page: Page, target: ElementRef) =>
  expectOrder(page, target, (s) => s.every((x, i) => i === 0 || x.share >= s[i - 1].share - 1e-9), 'share, smallest first'));

export const slicesByCategory = Then('the slices of {widget} should be ordered by category', (page: Page, target: ElementRef) =>
  expectOrder(page, target, (s) => s.every((x, i) => i === 0 || x.name.localeCompare(s[i - 1].name) >= 0), 'category name'),
{description: 'the wedges drawn in alphabetical order of their categories — which on demog-1000 is a different order from by-value'});

/** Shift moves every wedge out along its own mid-angle: the point the chart reports inside a wedge
 * sits at `shift + 0.7 x outer radius` from the centre of the disc, against `0.7 x outer radius`
 * on an unexploded pie. Pie mode only — a donut's reported point sits between the two radii. */
export const slicesOffCentre = Then('the slices of {widget} should sit {int} pixels off the centre', async (page: Page, target: ElementRef, shift: number) => {
  let shown = '';
  await expect.poll(async () => {
    const areas = await viewers.hitAreas(page, target);
    const values = await readings(page, target);
    const disc = areas['pie'];
    if (!disc)
      throw new Error(`${target.phrase} has no "pie" area — it painted no disc`);
    const centre = viewers.centerOf(disc);
    const offsets = Object.keys(areas).filter((k) => k.startsWith('slice ')).map((k) => {
      const name = k.slice('slice '.length).replace(/^"|"$/g, '');
      const p = viewers.centerOf(areas[k]);
      return {name, off: Math.hypot(p.x - centre.x, p.y - centre.y) - 0.7 * Number(values[`outer radius of "${name}"`])};
    });
    if (offsets.length === 0)
      throw new Error(`${target.phrase} reports no slices; it has: ${Object.keys(areas).join(', ')}`);
    shown = offsets.map((o) => `${o.name} ${o.off.toFixed(1)}`).join(', ');
    return offsets.every((o) => Math.abs(o.off - shift) <= 2);
  }, {timeout: 5000, message: `the slices of ${target.phrase} do not sit ${shift} pixels off the centre; they sit at: ${shown}`}).toBe(true);
}, {description: 'the exploded-pie claim: every wedge that far out of the disc\'s centre (within a pixel of rounding)'});


export const pickCategoryColumn = When('user picks {string} in the category selector of pie chart viewer',
  async (page: Page, column: string) => {
    const target = el('pie chart viewer');
    const loc = await viewers.viewerLocator(page, target);
    const selector = loc.locator('[name="div-column-combobox-category"]');
    await selector.waitFor({state: 'visible', timeout: 5000});
    const box = await selector.boundingBox();
    if (!box)
      throw new Error(`${target.phrase}: the category selector has no box`);
    await page.mouse.move(box.x + Math.min(10, box.width / 2), box.y + box.height / 2);
    await page.mouse.down();
    await page.mouse.up();
    await page.locator('.d4-column-grid').last().waitFor({state: 'visible', timeout: 5000});
    await page.keyboard.type(column);
    await page.keyboard.press('Enter');
    await expect(selector.locator('.d4-column-selector-column')).toHaveText(column, {timeout: 5000});
    await loc.evaluate((e) => (window as any).__bdd.settle(e, 300));
  }, {tier: 'ui', description: 'the category column re-picked on the chart itself, as a user re-picks it'});

/** The palette of the colour dialog a legend item opens: every swatch carries its own hex as a
 * name (`color_picker.dart`), so the pick is by colour and not by position. */
