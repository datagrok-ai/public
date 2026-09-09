/* The steps only the bar chart needs, all read from the hit areas the chart reports (`bar <category>`
   for every bar of the last frame, stack segments as `A | B`): a click on plot space no bar covers,
   the order of the bars along the category axis, their lengths, and how high the tallest reaches.
   Everything else in the bar chart features is the library's `viewers` tier and the platform's
   data steps (`grok-bdd list-steps`). */
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {el, ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

interface Bar {
  name: string;
  box: viewers.Box;
}

const spread = (xs: number[]): number => Math.max(...xs) - Math.min(...xs);

/** The category bars of the last frame, stack segments aside. */
async function categoryBars(page: Page, target: ElementRef): Promise<{bars: Bar[]; view: viewers.Box}> {
  const areas = await viewers.hitAreas(page, target);
  const bars = Object.keys(areas).filter((k) => k.startsWith('bar ') && !k.includes(' | ')).map((k) => ({name: k.slice(4), box: areas[k]}));
  if (bars.length === 0)
    throw new Error(`${target.phrase} reports no bars; it has: ${Object.keys(areas).join(', ') || 'no areas'}`);
  return {bars, view: areas['view']};
}

/** Bars side by side (vertical bars: their spans along x do not overlap) rather than one under
 * another (horizontal bars all start at the value axis). */
function sideBySide(bars: Bar[]): boolean {
  const sorted = [...bars].sort((a, b) => a.box.x - b.box.x);
  return sorted.every((b, i) => i === 0 || b.box.x >= sorted[i - 1].box.x + sorted[i - 1].box.width - 1);
}

/** Vertical bars left to right, their tops in that order: descending values have tops going
 * down (a negative bar's top is the baseline, below every positive bar's). */
async function expectOrder(page: Page, target: ElementRef, order: 'descend' | 'ascend'): Promise<void> {
  let shown = '';
  await expect.poll(async () => {
    const {bars} = await categoryBars(page, target);
    if (!sideBySide(bars))
      throw new Error(`${target.phrase} lays its bars out one under another; "from left to right" needs vertical bars`);
    bars.sort((a, b) => a.box.x - b.box.x);
    shown = bars.map((b) => `${b.name || '(empty)'} at ${Math.round(b.box.y)}`).join(', ');
    return bars.every((b, i) => i === 0 || (order === 'descend' ? b.box.y >= bars[i - 1].box.y - 0.5 : b.box.y <= bars[i - 1].box.y + 0.5));
  }, {timeout: 5000, message: `the bars of ${target.phrase} do not ${order} from left to right; their tops: ${shown}`}).toBe(true);
}

export const barsDescend = Then('the bars of {widget} should descend from left to right', (page: Page, target: ElementRef) =>
  expectOrder(page, target, 'descend'), {description: 'vertical bars, each no taller than the one to its left (a negative bar counts as the shortest)'});

export const barsAscend = Then('the bars of {widget} should ascend from left to right', (page: Page, target: ElementRef) =>
  expectOrder(page, target, 'ascend'));

export const barsStacked = Then('the bars of {widget} should lie one under another', async (page: Page, target: ElementRef) => {
  let shown = '';
  await expect.poll(async () => {
    const {bars} = await categoryBars(page, target);
    shown = bars.map((b) => `${b.name || '(empty)'} at x ${Math.round(b.box.x)}`).join(', ');
    return !sideBySide(bars);
  }, {timeout: 5000, message: `the bars of ${target.phrase} lie side by side: ${shown}`}).toBe(true);
}, {description: 'horizontal bars, one row per category'});

/** A negative bar hangs from the baseline: its top is where a positive bar's bottom is, and it has
 * a height of its own. */
export const barHangsBelow = Then('the {string} area of {widget} should hang below the {string} area', async (page: Page, below: string, target: ElementRef, above: string) => {
  let shown = '';
  await expect.poll(async () => {
    const areas = await viewers.hitAreas(page, target);
    const a = areas[below];
    const b = areas[above];
    if (!a || !b)
      throw new Error(`${target.phrase} has no "${a ? above : below}" area; it has: ${Object.keys(areas).join(', ')}`);
    shown = `"${below}" spans y ${Math.round(a.y)}..${Math.round(a.y + a.height)}, "${above}" y ${Math.round(b.y)}..${Math.round(b.y + b.height)}`;
    return a.y >= b.y + b.height - 1 && a.height > 1;
  }, {timeout: 5000, message: `the "${below}" area of ${target.phrase} does not hang below the "${above}" area: ${shown}`}).toBe(true);
}, {description: 'a negative bar under a positive one: it starts at the baseline the other ends on, and reaches down'});

export const zoomCategories = When('user zooms into the categories from the {string} area to the {string} area of {widget}',
  async (page: Page, from: string, to: string, target: ElementRef) => {
    const a = viewers.centerOf(await viewers.hitArea(page, target, from, true));
    const b = viewers.centerOf(await viewers.hitArea(page, target, to));
    await page.keyboard.down('Alt');
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move(b.x, b.y, {steps: 3});
    await page.mouse.up();
    await page.keyboard.up('Alt');
  }, {tier: 'ui', description: 'an Alt-drag between two bars — the chart zooms its category axis to the bars the drag spans'});

/** The bars' lengths along the value axis: equal within a pixel, or not. */
async function expectLengths(page: Page, target: ElementRef, equal: boolean): Promise<void> {
  let shown = '';
  await expect.poll(async () => {
    const {bars} = await categoryBars(page, target);
    const lengths = bars.map((b) => sideBySide(bars) ? b.box.height : b.box.width);
    shown = bars.map((b, i) => `${b.name || '(empty)'} ${Math.round(lengths[i])}`).join(', ');
    return spread(lengths) <= 1 === equal;
  }, {timeout: 5000, message: `the bars of ${target.phrase} ${equal ? 'differ in length' : 'are of equal length'}: ${shown}`}).toBe(true);
}

export const barsEqual = Then('the bars of {widget} should be of equal length', (page: Page, target: ElementRef) => expectLengths(page, target, true),
  {description: 'every category bar as long as the others — Relative Values with a Stack column'});

export const barsDiffer = Then('the bars of {widget} should differ in length', (page: Page, target: ElementRef) => expectLengths(page, target, false));

export const tallestReaches = Then('the tallest bar of {widget} should start within the top {int}% of the view', async (page: Page, target: ElementRef, percent: number) => {
  let shown = '';
  await expect.poll(async () => {
    const {bars, view} = await categoryBars(page, target);
    if (!sideBySide(bars))
      throw new Error(`${target.phrase} lays its bars out one under another; the top of the view is for vertical bars`);
    const top = Math.min(...bars.map((b) => b.box.y));
    const fraction = (top - view.y) / view.height * 100;
    shown = `${Math.round(fraction)}%`;
    return fraction <= percent;
  }, {timeout: 5000, message: `the tallest bar of ${target.phrase} starts ${shown} down the view, not within the top ${percent}%`}).toBe(true);
}, {description: 'no blank band above the tallest bar wider than that share of the view'});

/** A point of the view no bar covers: of a coarse grid over the view, the one farthest from every
 * bar; the baseline is taken so "should have repainted" can follow. */
async function emptySpace(page: Page): Promise<{x: number; y: number}> {
  const target = el('bar chart viewer');
  const view = await viewers.hitArea(page, target, 'view', true);
  const areas = await viewers.hitAreas(page, target);
  const bars = Object.keys(areas).filter((k) => k.startsWith('bar ')).map((k) => areas[k]);
  const distance = (x: number, y: number, b: viewers.Box): number =>
    Math.hypot(Math.max(b.x - x, 0, x - b.x - b.width), Math.max(b.y - y, 0, y - b.y - b.height));
  let best: {x: number; y: number; d: number} | undefined;
  for (let i = 1; i < 10; i++) {
    for (let j = 1; j < 10; j++) {
      const x = view.x + view.width * i / 10;
      const y = view.y + view.height * j / 10;
      const d = bars.length === 0 ? Infinity : Math.min(...bars.map((b) => distance(x, y, b)));
      if (!best || d > best.d)
        best = {x, y, d};
    }
  }
  if (!best || best.d < 3)
    throw new Error('bar chart viewer has no empty plot space: the bars cover the view');
  return best;
}

export const clickEmptySpace = When('user clicks on empty plot space of bar chart viewer', async (page: Page) => {
  const p = await emptySpace(page);
  await page.mouse.click(p.x, p.y);
}, {tier: 'ui', description: 'plot space with no bar under it — a click there releases the bar filter'});

export const doubleClickEmptySpace = When('user double-clicks on empty plot space of bar chart viewer', async (page: Page) => {
  const p = await emptySpace(page);
  await page.mouse.dblclick(p.x, p.y);
}, {tier: 'ui', description: 'resets the view'});
