/* The `viewers` tier, part two: steps that were written for one viewer and turned out to say
   nothing about it. A phrase belongs here as soon as a second viewer wants it — the tile viewer's
   description placement is every viewer's Description Position, the pie chart's category selector
   is the scatter plot's X selector under another name, the trellis plot's cells are the matrix
   plot's. Each one reads a hit area or a reading the viewer reports, so a viewer joins simply by
   reporting the same key; where a name is a convention rather than a rectangle, the convention is
   written next to the step.
   What stays in a package binding is a gesture that only that viewer's DOM has (the tile viewer's
   sketch designer, the pc plot's axis-label drag) or an arithmetic only it can check (the pivot's
   aggregation against a groupBy). */
import {Page} from '@playwright/test';
import {expect, pollMs} from '../../../src/runtime/patience.js';
import {Then, When} from '../../../src/registry.js';
import type {ElementRef} from '../../../src/runtime/args.js';
import {exactText} from '../../../src/runtime/locate.js';
import * as g from '../../../src/runtime/gestures.js';
import * as v from '../../../src/runtime/viewers.js';

async function settle(page: Page, target: ElementRef, capMs = 300): Promise<void> {
  const loc = await v.viewerLocator(page, target);
  await loc.evaluate((e, ms) => (window as any).__bdd.settle(e, ms), capMs);
}

// --- readings that are lists ----------------------------------------------------------------------

/* A reading whose value is a comma-separated list: the tile viewer's `fields` and `lane names`, the
   pc plot's `axes`, the trellis plot's `x columns`. */
async function expectReadingContains(page: Page, target: ElementRef, name: string, item: string, negate: boolean): Promise<void> {
  let last = '';
  const poll = expect.poll(async () => {
    last = String(await v.readValue(page, target, name));
    return last.split(/\s*,\s*/).includes(item);
  }, {timeout: pollMs(5000), message: `"${name}" of ${target.phrase} is "${last}"`});
  await (negate ? poll.not : poll).toBe(true);
}

export const readingContains = Then('the {string} reading of {widget} should contain {string}',
  (page: Page, name: string, target: ElementRef, item: string) => expectReadingContains(page, target, name, item, false),
{description: 'a reading the viewer reports as a comma-separated list holds that member'});

export const readingNotContains = Then('the {string} reading of {widget} should not contain {string}',
  (page: Page, name: string, target: ElementRef, item: string) => expectReadingContains(page, target, name, item, true),
{description: 'the negative; a feature pairs it with a positive one that proves the member was there to lose'});

// --- the viewer's own context menu ------------------------------------------------------------------

/** The region the widget declares as its context-menu target (`getWidgetStatus().shortcuts`
 * `ContextMenu` → a hit-area name), or nothing when it declares none and the centre will do. A
 * viewer whose middle is covered by something with a menu of its own — the tile viewer's cards,
 * the trellis plot's inner viewers, the pivot's grid — declares the strip where its own menu is. */
async function menuArea(page: Page, target: ElementRef): Promise<string | undefined> {
  await v.installViewerRuntime(page);
  const loc = await v.viewerLocator(page, target);
  const name: string | null = await loc.evaluate((e) =>
    (window as any).__bdd.viewerOf(e).getWidgetStatus()?.shortcuts?.['ContextMenu'] ?? null);
  return name ?? undefined;
}

export const openViewerMenu = When('user opens the viewer menu of {widget}', async (page: Page, target: ElementRef) => {
  await v.openContextMenuOf(page, target, await menuArea(page, target));
}, {tier: 'ui', description: 'right-clicks the region the widget declares as its ContextMenu shortcut, its centre otherwise'});

export const pickFromViewerMenu = When('user picks {string} from the viewer menu of {widget}', async (page: Page, path: string, target: ElementRef) => {
  await v.snapshot(page, target);
  await v.openContextMenuOf(page, target, await menuArea(page, target));
  await v.pickMenuPath(page, path);
}, {tier: 'ui', description: 'the same, then the path picked; the baseline is taken first so a repaint check can follow'});

// --- where the description sits -----------------------------------------------------------------

async function expectDescriptionPlace(page: Page, target: ElementRef, above: boolean): Promise<void> {
  const loc = await v.viewerLocator(page, target);
  const view = await v.hitArea(page, target, 'view');
  const box = await loc.locator('.d4-viewer-description').filter({visible: true}).first().boundingBox();
  if (box === null)
    throw new Error(`${target.phrase} shows no description`);
  const ok = above ? box.y + box.height <= view.y + 1 : box.y + 1 >= view.y + view.height;
  expect(ok, `the description of ${target.phrase} is not ${above ? 'above' : 'below'} its content: the description spans ` +
    `${Math.round(box.y)}..${Math.round(box.y + box.height)} and the content ${Math.round(view.y)}..${Math.round(view.y + view.height)}`).toBe(true);
}

export const descriptionAbove = Then('the description of {widget} should be above its content',
  (page: Page, target: ElementRef) => expectDescriptionPlace(page, target, true),
{description: 'Description Position Top: the description box ends where the "view" area begins'});

export const descriptionBelow = Then('the description of {widget} should be below its content',
  (page: Page, target: ElementRef) => expectDescriptionPlace(page, target, false),
{description: 'Description Position Bottom: the description box starts where the "view" area ends'});

// --- plot space nothing covers --------------------------------------------------------------------

/** A point inside the `view` area that no other hit area is near: of a coarse grid over the view,
 * the one farthest from every drawn thing. A viewer that knows better says so by reporting an
 * `empty space` area, and then that is used instead. The baseline is taken so a repaint check can
 * follow the click. */
async function emptySpace(page: Page, target: ElementRef): Promise<{x: number; y: number}> {
  const areas = await v.hitAreas(page, target);
  await v.snapshot(page, target);
  if (areas['empty space'] !== undefined)
    return v.centerOf(areas['empty space']);
  const view = areas['view'];
  if (view === undefined)
    throw new Error(`${target.phrase} reports no "view" area to look for empty space in; it has: ${Object.keys(areas).join(', ')}`);
  const drawn = Object.keys(areas).filter((k) => k !== 'view' && k !== 'plot').map((k) => areas[k]);
  const distance = (x: number, y: number, b: v.Box): number =>
    Math.hypot(Math.max(b.x - x, 0, x - b.x - b.width), Math.max(b.y - y, 0, y - b.y - b.height));
  let best: {x: number; y: number; d: number} | undefined;
  for (let i = 1; i < 10; i++)
    for (let j = 1; j < 10; j++) {
      const x = view.x + view.width * i / 10;
      const y = view.y + view.height * j / 10;
      const d = drawn.length === 0 ? Infinity : Math.min(...drawn.map((b) => distance(x, y, b)));
      if (!best || d > best.d)
        best = {x, y, d};
    }
  if (!best || best.d < 3)
    throw new Error(`${target.phrase} has no empty plot space: what it draws covers the view`);
  return best;
}

export const clickEmptySpace = When('user clicks on empty plot space of {widget}', async (page: Page, target: ElementRef) => {
  const p = await emptySpace(page, target);
  await page.mouse.click(p.x, p.y);
}, {tier: 'ui', description: 'plot space with nothing under it — a click there clears the selection or releases the filter'});

export const doubleClickEmptySpace = When('user double-clicks on empty plot space of {widget}', async (page: Page, target: ElementRef) => {
  const p = await emptySpace(page, target);
  await page.mouse.dblclick(p.x, p.y);
}, {tier: 'ui', description: 'resets the view'});

// --- range sliders --------------------------------------------------------------------------------

/** One handle of one range slider, dragged by a pixel count. The sliders are revealed on
 * `mouseenter`, so the pointer goes over the viewer first and the handle is then located either by
 * the hit area the viewer reports for it (`range max handle "AGE"`) or, for a viewer that reports
 * none, by the slider's own name in the DOM (`svg[name="x-slider"]`). A positive count moves the
 * handle right on a horizontal slider and down on a vertical one. */
export const dragRangeHandle = When('user drags the {word} handle of the {string} range slider of {widget} by {int} pixels',
  async (page: Page, handle: string, slider: string, target: ElementRef, px: number) => {
    if (handle !== 'min' && handle !== 'max')
      throw new Error(`a range slider has a min and a max handle, not a "${handle}" one`);
    // a viewer added or closed next to this one is still being laid out: the pointer would go to
    // where the plot was, and the sliders it reveals on enter would never appear
    await settle(page, target);
    await page.mouse.move(0, 0);
    const centre = v.centerOf(await v.hitArea(page, target, 'view'));
    await page.mouse.move(centre.x, centre.y);
    const name = `range ${handle} handle "${slider}"`;
    let box: {x: number; y: number; width: number; height: number} | null = null;
    let horizontal = slider === 'x';
    // the slider lays itself out a frame after the pointer reveals it, so whether the viewer
    // reports the handle is polled for and not read once — a single read that came too early sent
    // this down the DOM fallback, which then failed on a viewer that does report it
    let areas = await v.hitAreas(page, target);
    if (areas[name] === undefined) {
      await expect.poll(async () => (areas = await v.hitAreas(page, target))[name] !== undefined,
        {timeout: pollMs(2000)}).toBe(true).catch(() => undefined);
    }
    if (areas[name] !== undefined)
      box = await v.hitArea(page, target, name, true);
    else {
      // one axis per plot names its strip "x axis", one axis per column names it `axis "AGE"`
      const stripName = Object.keys(areas).find((k) => k === `${slider} axis` || k === `axis "${slider}"`);
      if (stripName === undefined)
        throw new Error(`${target.phrase} has no "${name}" area and no axis strip for "${slider}"; it has: ${Object.keys(areas).join(', ') || 'none'}`);
      const strip = await v.hitArea(page, target, stripName, true);
      horizontal = strip.width >= strip.height;
      const s = v.centerOf(strip);
      await page.mouse.move(s.x - 3, s.y);
      await page.mouse.move(s.x, s.y);
      const loc = await v.viewerLocator(page, target);
      const svg = loc.locator(`svg[name="${slider}-slider"]`);
      await svg.waitFor({state: 'visible', timeout: 5000}).catch(() => {
        throw new Error(`${target.phrase} has neither a "${name}" area nor a "${slider}-slider"`);
      });
      box = await svg.locator(`[name="${handle}-handle"]`).boundingBox();
    }
    if (box === null)
      throw new Error(`${target.phrase}: the "${slider}" range slider shows no ${handle} handle`);
    const x = box.x + box.width / 2;
    const y = box.y + box.height / 2;
    await page.mouse.move(x, y);
    await page.mouse.down();
    await page.mouse.move(horizontal ? x + px / 2 : x, horizontal ? y : y + px / 2);
    await page.mouse.move(horizontal ? x + px : x, horizontal ? y : y + px);
    await page.mouse.up();
    await settle(page, target);
  }, {tier: 'ui', description: 'the slider is named "x"/"y" on a plot with axis sliders and by the column on a plot with one per axis'});

/** A slider handle dragged to the start or the end of its OWN track, both taken from the regions
 * the viewer reports (`<axis> scroll slider`, `<axis> scroll <min|max> handle`). The by-pixels step
 * above needs a count that depends on how wide the viewer happened to be laid out; this one does
 * not, which is what a "narrow the viewport as far as it goes" claim needs. */
export const dragHandleToEnd = When('user drags the {word} handle of the {string} scroll slider of {widget} to its {word}',
  async (page: Page, handle: string, axis: string, target: ElementRef, where: string) => {
    if (handle !== 'min' && handle !== 'max')
      throw new Error(`a slider has a min and a max handle, not a "${handle}" one`);
    if (where !== 'start' && where !== 'end')
      throw new Error(`a track has a start and an end, not a "${where}"`);
    const view = v.centerOf(await v.hitArea(page, target, 'view'));
    await page.mouse.move(view.x, view.y);
    const track = await v.hitArea(page, target, `${axis} scroll slider`, true);
    const from = v.centerOf(await v.hitArea(page, target, `${axis} scroll ${handle} handle`, true));
    const horizontal = track.width >= track.height;
    const end = where === 'end'
      ? (horizontal ? track.x + track.width : track.y + track.height)
      : (horizontal ? track.x : track.y);
    const to = horizontal ? {x: end, y: from.y} : {x: from.x, y: end};
    await page.mouse.move(from.x, from.y);
    await page.mouse.down();
    await page.mouse.move((from.x + to.x) / 2, (from.y + to.y) / 2);
    await page.mouse.move(to.x, to.y);
    await page.mouse.up();
    await settle(page, target, 3000);
  }, {tier: 'ui', description: 'no pixels-per-unit anywhere — the track and the handle are both reported regions'});

/** The value axis zoomed to its upper part: the top handle of the vertical slider dragged down by
 * two fifths of the track. */
export const zoomValueAxis = When('user zooms into the value axis of {widget}', async (page: Page, target: ElementRef) => {
  const axis = v.centerOf(await v.hitArea(page, target, 'y axis', true));
  await page.mouse.move(axis.x - 3, axis.y);
  await page.mouse.move(axis.x, axis.y);
  const loc = await v.viewerLocator(page, target);
  const slider = loc.locator('svg[type="range-slider"][name="y-slider"]');
  await slider.waitFor({state: 'visible'});
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
}, {tier: 'ui', description: 'the "y axis" area must be reported; pair it with "should show a narrower value range than before"'});

// --- the selectors drawn on the viewer -------------------------------------------------------------

/** An on-viewer column selector (`div-column-combobox-<property>`): a mouse-down on its caption
 * opens the column grid, typing opens the grid's search and Enter takes the column typed. The name
 * is the property the selector binds, lowercased — `x`, `y`, `color`, `size`, `category`, `value`.
 * Some selectors are hover-revealed, so the pointer goes over the viewer first. */
export const pickInColumnSelector = When('user picks {string} in the {string} column selector of {widget}',
  async (page: Page, column: string, which: string, target: ElementRef) => {
    const loc = await v.viewerLocator(page, target);
    const centre = v.centerOf(await v.hitArea(page, target, 'view'));
    await page.mouse.move(centre.x, centre.y);
    const selector = loc.locator(`[name="div-column-combobox-${which.toLowerCase()}"]`);
    await selector.waitFor({state: 'visible', timeout: 5000});
    const box = await selector.boundingBox();
    if (!box)
      throw new Error(`${target.phrase}: the "${which}" column selector has no box`);
    await page.mouse.move(box.x + Math.min(10, box.width / 2), box.y + box.height / 2);
    await page.mouse.down();
    await page.mouse.up();
    // the pointer stays where it is: these selectors are revealed by the hover, and leaving the
    // viewer takes the popup with them. A row it rests on is previewed onto the selector, which a
    // pick that lands overwrites anyway
    await g.pickInColumnGrid(page, column, `the "${which}" column selector of ${target.phrase}`, selector);
    await expect(selector.locator('.d4-column-selector-column')).toHaveText(column, {timeout: pollMs(5000)});
    await settle(page, target);
  }, {tier: 'ui', description: 'the column re-picked on the chart itself, the way a user re-picks it'});

/** The same selector, typed into and committed, with no claim about the result: a selector that
 * offers only some columns leaves the column it had, and the feature reads it afterwards. */
export const typeInColumnSelector = When('user types {string} into the {string} column selector of {widget}',
  async (page: Page, column: string, which: string, target: ElementRef) => {
    const loc = await v.viewerLocator(page, target);
    const centre = v.centerOf(await v.hitArea(page, target, 'view'));
    await page.mouse.move(centre.x, centre.y);
    const selector = loc.locator(`[name="div-column-combobox-${which.toLowerCase()}"]`);
    await selector.waitFor({state: 'visible', timeout: 5000});
    const box = await selector.boundingBox();
    if (!box)
      throw new Error(`${target.phrase}: the "${which}" column selector has no box`);
    await page.mouse.move(box.x + Math.min(10, box.width / 2), box.y + box.height / 2);
    await page.mouse.down();
    await page.mouse.up();
    await g.typeInColumnGrid(page, column, `the "${which}" column selector of ${target.phrase}`, selector);
    await settle(page, target);
  }, {tier: 'ui', description: 'for the negative: the selector takes the name typed only when it offers that column'});

/** The type selector of a viewer that hosts another viewer (a `ComboPopup` named "viewer
 * selector"): a click opens the list, and the pick goes through the host's `setViewerType`, which
 * is what announces the type change — writing the property does not. */
export const pickInnerViewer = When('user picks {string} in the viewer selector of {widget}',
  async (page: Page, type: string, target: ElementRef) => {
    const loc = await v.viewerLocator(page, target);
    const selector = loc.locator('[name="viewer selector"]').first();
    await selector.waitFor({state: 'visible', timeout: 5000});
    await v.snapshot(page, target);
    await selector.click();
    const popup = page.locator('.d4-combo-popup-expanded').last();
    await popup.waitFor({state: 'visible', timeout: 5000});
    const item = popup.locator('.d4-list-item').filter({hasText: exactText(type)}).first();
    if (await item.count() === 0)
      throw new Error(`no "${type}" in the viewer selector; it offers: ${(await popup.locator('.d4-list-item').allTextContents()).join(', ')}`);
    await item.click();
    await settle(page, target, 3000);
  }, {tier: 'ui', description: 'the control panel must be shown; the pick is the gesture that fires the type-changed event'});

/** A setting of the viewer inside every cell: `setOptions({innerViewerLook: {...}})`, the only
 * channel the cells read back — a mutation of the look object itself is not. The name is the
 * inner look's own field (`allowZoom`, `xColumnName`, `colorColumnName`). */
export const setInnerProperty = When('user sets {string} inner property of {widget} to {string}',
  async (page: Page, name: string, target: ElementRef, value: string) => {
    const loc = await v.viewerLocator(page, target);
    await v.snapshot(page, target);
    const coerced: unknown = /^(true|false)$/i.test(value) ? /^true$/i.test(value)
      : value !== '' && !Number.isNaN(Number(value)) ? Number(value) : value;
    await loc.evaluate((e, [n, val]) => {
      (window as any).__bdd.viewerOf(e).setOptions({innerViewerLook: {[n as string]: val}});
    }, [name, coerced] as [string, unknown]);
    await settle(page, target, 3000);
  }, {tier: 'api', description: 'by the field name the inner look serializes, not by a property caption'});

// --- geometry of a viewer laid out in cells ----------------------------------------------------------

export const cellsWideTall = Then('the cells of {widget} should be {int} wide and {int} tall',
  async (page: Page, target: ElementRef, wide: number, tall: number) => {
    await v.expectReading(page, target, 'columns', 'equal', wide);
    await v.expectReading(page, target, 'rows', 'equal', tall);
  }, {description: 'the viewport in cells, from the "columns" and "rows" readings'});

export const areaHangsBelow = Then('the {string} area of {widget} should hang below the {string} area',
  async (page: Page, below: string, target: ElementRef, above: string) => {
    let shown = '';
    await expect.poll(async () => {
      const areas = await v.hitAreas(page, target);
      const a = areas[below];
      const b = areas[above];
      if (!a || !b) {
        shown = `it has: ${Object.keys(areas).join(', ') || 'none'}`;
        return false;
      }
      shown = `"${below}" spans y ${Math.round(a.y)}..${Math.round(a.y + a.height)}, "${above}" y ${Math.round(b.y)}..${Math.round(b.y + b.height)}`;
      return a.y >= b.y + b.height - 1 && a.height > 1;
    }, {timeout: pollMs(5000), message: `the "${below}" area of ${target.phrase} does not hang below the "${above}" area: ${shown}`}).toBe(true);
  }, {description: 'one area starts where the other ends and reaches further down — a negative bar under a positive one'});

// --- one row of what the viewer drew -----------------------------------------------------------------

/** Rows as the `line of row <n>` / `tile of row <n>` / `marker of row <n>` hit areas count them,
 * from 1: the selection is read off the frame the viewer is bound to. */
async function rowSelected(page: Page, target: ElementRef, row: number): Promise<boolean> {
  await v.installViewerRuntime(page);
  const loc = await v.viewerLocator(page, target);
  return loc.evaluate((e, n) => (window as any).__bdd.viewerOf(e).dataFrame.selection.get(n - 1) as boolean, row);
}

export const lineSelected = Then('the line of row {int} of {widget} should be selected', async (page: Page, row: number, target: ElementRef) => {
  await expect.poll(() => rowSelected(page, target, row), {message: `row ${row} of the table ${target.phrase} draws`}).toBe(true);
}, {description: 'the row behind a line the viewer drew, by the number its hit areas use'});

export const lineNotSelected = Then('the line of row {int} of {widget} should not be selected', async (page: Page, row: number, target: ElementRef) => {
  await expect.poll(() => rowSelected(page, target, row), {message: `row ${row} of the table ${target.phrase} draws`}).toBe(false);
});

// --- viewers that lay a card out per row ---------------------------------------------------------

/* The convention: a card-laying viewer reports one `<COLUMN> of row <r>` reading per field it drew,
   whose value is the grid's display string for that cell, and a `fields` reading listing the
   columns the card is composed of. The tile viewer and the forms viewer both do. */
async function shownColumn(page: Page, target: ElementRef, column: string): Promise<{row: number; text: string}[]> {
  await v.installViewerRuntime(page);
  const loc = await v.viewerLocator(page, target);
  const values: Record<string, unknown> = await loc.evaluate((e) =>
    (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {});
  const pattern = new RegExp(`^${column.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')} of row (\\d+)$`);
  const out: {row: number; text: string}[] = [];
  for (const key of Object.keys(values)) {
    const m = pattern.exec(key);
    if (m)
      out.push({row: Number(m[1]), text: String(values[key])});
  }
  return out;
}

async function expectEveryCard(page: Page, target: ElementRef, column: string, holds: (text: string) => boolean, what: string): Promise<void> {
  let bad = '';
  let seen = 0;
  await expect.poll(async () => {
    const cells = await shownColumn(page, target, column);
    seen = cells.length;
    if (seen === 0)
      return false;
    const wrong = cells.filter((c) => !holds(c.text));
    bad = wrong.map((c) => `row ${c.row} shows "${c.text}"`).join(', ');
    return wrong.length === 0;
  }, {timeout: pollMs(5000), message: seen === 0
    ? `${target.phrase} has laid out no card with a "${column}" field`
    : `not every card of ${target.phrase} shows ${what} in "${column}": ${bad}`}).toBe(true);
}

export const everyTileShows = Then('every tile of {widget} should show {string} in {string}',
  (page: Page, target: ElementRef, value: string, column: string) =>
    expectEveryCard(page, target, column, (t) => t === value, `"${value}"`),
{description: 'every card the viewer has laid out shows that text for the column — the "all tiles are M" claim'});

export const everyTileBetween = Then('every tile of {widget} should show a value between {float} and {float} in {string}',
  (page: Page, target: ElementRef, lo: number, hi: number, column: string) =>
    expectEveryCard(page, target, column, (t) => Number(t) >= lo && Number(t) <= hi, `a value between ${lo} and ${hi}`),
{description: 'the "every age is over 50" claim, read off the cards and not off the table'});

let rememberedFields: string[] | null = null;

async function fieldsOf(page: Page, target: ElementRef): Promise<string[]> {
  const text = String(await v.readValue(page, target, 'fields'));
  return text === '' ? [] : text.split(/\s*,\s*/);
}

export const rememberFields = When('user remembers the fields of {widget}', async (page: Page, target: ElementRef) => {
  rememberedFields = await fieldsOf(page, target);
  if (rememberedFields.length === 0)
    throw new Error(`${target.phrase} shows no fields to remember`);
}, {tier: 'api', description: 'the card\'s composition, for the contrast that follows'});

/* Neither half of the contrast names the column that left or the one that took its place: an
   auto-generated card orders its fields by relevance, and the excluded one is whatever the score
   left over. */
async function expectComposition(page: Page, target: ElementRef, refilled: boolean): Promise<void> {
  if (rememberedFields === null)
    throw new Error('nothing was remembered: put "user remembers the fields of <widget>" before the change');
  const before = rememberedFields;
  let now: string[] = [];
  await expect.poll(async () => {
    now = await fieldsOf(page, target);
    const gone = before.filter((f) => !now.includes(f));
    const gained = now.filter((f) => !before.includes(f));
    return now.length === before.length && (refilled
      ? gone.length === 1 && gained.length === 1
      : gone.length === 0 && gained.length === 0);
  }, {timeout: pollMs(8000), message: refilled
    ? `${target.phrase} did not refill the freed slot: it showed ${before.join(', ')} and now shows ${now.join(', ')}`
    : `${target.phrase} did not keep its composition: it showed ${before.join(', ')} and now shows ${now.join(', ')}`}).toBe(true);
}

export const fieldsRefilled = Then('the fields of {widget} should have refilled the freed slot',
  (page: Page, target: ElementRef) => expectComposition(page, target, true),
{description: 'the auto-generated card: the column that left took a field with it and a column that had none took its place'});

export const fieldsAsRemembered = Then('the fields of {widget} should be as remembered',
  (page: Page, target: ElementRef) => expectComposition(page, target, false),
{description: 'the designed card: the same fields as before, so the field of a column that left stays, empty'});

/** A Shift-drag along a closed polygon inside the area — what the platform reads as a lasso while the
 * Lasso Tool is on. The pointer walks each leg in several moves: the browser delivers pointer moves
 * frame-aligned, and the selector builds the polygon from the moves it sees. */
export const dragLasso = When('user drags a lasso over the {string} area of {widget}',
  async (page: Page, area: string, target: ElementRef) => {
    const b = await v.hitArea(page, target, area, true);
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

// --- readings that are prose, and readings that must not be there ----------------------------

/* `should contain` splits a reading on commas and asks for membership — the right question for a
   list (`fields`, `node names`, `stats`), the wrong one for a paragraph. */
async function expectTextIn(page: Page, name: string, target: ElementRef, text: string, negate: boolean): Promise<void> {
  let last = '';
  const holds = async (): Promise<boolean> => {
    last = String(await v.readValue(page, target, name));
    return last.includes(text);
  };
  try {
    const poll = expect.poll(holds, {timeout: pollMs(5000)});
    await (negate ? poll.not : poll).toBe(true);
  }
  catch {
    // built after the poll, so it names the text that was actually there
    throw new Error(`"${name}" of ${target.phrase} is ${JSON.stringify(last)}`);
  }
}

export const readingIncludes = Then('the {string} reading of {widget} should include the text {string}',
  (page: Page, name: string, target: ElementRef, text: string) => expectTextIn(page, name, target, text, false),
{description: 'a reading the viewer reports as prose — a markup viewer\'s rendered text — holds that phrase'});

export const readingExcludes = Then('the {string} reading of {widget} should not include the text {string}',
  (page: Page, name: string, target: ElementRef, text: string) => expectTextIn(page, name, target, text, true),
{description: 'the negative; pair it with a positive that proves the phrase was there to lose'});

/** The mirror of `should not have a "x" area`: the viewer refuses the reading in this state. A
 * viewer in a message state keeps no counts, and its previous chart object may still be alive —
 * "no `words` reading" is then a real claim about what the status will not answer. */
export const noSuchReading = Then('{widget} should not report a {string} reading',
  async (page: Page, target: ElementRef, name: string) => {
    await v.installViewerRuntime(page);
    const loc = await v.viewerLocator(page, target);
    let has: string[] = [];
    await expect.poll(async () => {
      const r: {now?: unknown; has: string[]} =
        await loc.evaluate((el, n) => (window as any).__bdd.valueChange(el, n), name);
      has = r.has;
      return r.now !== undefined && r.now !== null;
    }, {timeout: pollMs(5000), message: `${target.phrase} reports: ${has.join(', ') || 'no readings'}`}).toBe(false);
  }, {description: 'the reading is absent, not merely empty'});

// --- two of the viewer's own areas, compared as they stand -------------------------------------

type Dim = 'taller' | 'shorter' | 'wider' | 'narrower';

const AREA_DIM: Record<Dim, (b: v.Box) => number> = {
  taller: (b) => b.height, shorter: (b) => b.height,
  wider: (b) => b.width, narrower: (b) => b.width,
};

export const areaBiggerThanArea = Then('the {string} area of {widget} should be {word} than the {string} area',
  async (page: Page, a: string, target: ElementRef, comparison: string, b: string) => {
    if (!(comparison in AREA_DIM))
      throw new Error(`an area is taller, shorter, wider or narrower than another, not "${comparison}"`);
    const dim = comparison as Dim;
    const unit = dim === 'taller' || dim === 'shorter' ? 'tall' : 'wide';
    let shown = '';
    // both areas are read together and polled for: a viewer laying itself out reports neither, and
    // a single read of a layout in progress is a race, not a claim about the layout
    const holds = async (): Promise<boolean> => {
      const areas = await v.hitAreas(page, target);
      if (areas[a] === undefined || areas[b] === undefined) {
        shown = `${target.phrase} has no "${areas[a] ? b : a}" area; it has: ${Object.keys(areas).join(', ') || 'none'}`;
        return false;
      }
      const [x, y] = [AREA_DIM[dim](areas[a]), AREA_DIM[dim](areas[b])];
      shown = `the "${a}" area is ${Math.round(x)} ${unit} and "${b}" is ${Math.round(y)}`;
      return dim === 'taller' || dim === 'wider' ? x > y : x < y;
    };
    try {
      await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
    }
    catch {
      throw new Error(shown);
    }
  }, {description: 'the library can say "taller than before" and "the same height"; this says one is bigger than the other now'});

// --- the place a widget occupied ----------------------------------------------------------------

/* Any claim about what happens after a widget is closed needs the box it HAD — by then the element
   is gone, and no viewer owns the rectangle. */
const places = new WeakMap<Page, v.Box>();

export const rememberPlace = When('user remembers the place of {widget}', async (page: Page, target: ElementRef) => {
  const loc = await v.viewerLocator(page, target);
  const box = await loc.boundingBox();
  if (box === null)
    throw new Error(`${target.phrase} has no box on screen`);
  places.set(page, box);
}, {tier: 'ui', description: 'the rectangle the widget occupies now, for a claim made after it is gone'});

export const moveAcrossPlace = When('user moves the pointer across the remembered place', async (page: Page) => {
  const box = places.get(page);
  if (box === undefined)
    throw new Error('no place was remembered — "user remembers the place of <widget>" comes first');
  const y = box.y + box.height / 2;
  for (const fraction of [0.3, 0.5, 0.7])
    await page.mouse.move(box.x + box.width * fraction, y, {steps: 3});
  // a pointer move is delivered frame-aligned and a viewer's pointer handlers run inside that
  // frame, so two animation frames is the browser's own signal, not a timer
  await page.evaluate(() => new Promise<void>((done) =>
    requestAnimationFrame(() => requestAnimationFrame(() => done()))));
}, {tier: 'ui', description: 'three real moves over the region and two animation frames — what a disposed viewer would throw on'});

// --- a drag between two widgets ---------------------------------------------------------------------

/** The platform's drag starts on the first move after the button goes down, so the pointer travels
 * in steps; a single jump drops nothing. */
export const dragAreaOntoWidget = When('user drags the {string} area of {widget} onto the {string} area of {widget}',
  async (page: Page, from: string, source: ElementRef, to: string, target: ElementRef) => {
    const a = v.centerOf(await v.hitArea(page, source, from, true));
    const b = v.centerOf(await v.hitArea(page, target, to, true));
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    for (let i = 1; i <= 10; i++)
      await page.mouse.move(a.x + (b.x - a.x) * i / 10, a.y + (b.y - a.y) * i / 10);
    await page.mouse.up();
    await settle(page, target);
  }, {tier: 'ui', description: 'a grid column header onto a pivot row, a card onto a lane — drag and drop across two widgets'});
