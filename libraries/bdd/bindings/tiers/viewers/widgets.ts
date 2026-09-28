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
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from '../../../src/runtime/patience.js';
import {Then, When} from '../../../src/registry.js';
import type {ElementRef} from '../../../src/runtime/args.js';
import {el} from '../../../src/runtime/args.js';
import {exactText, locate} from '../../../src/runtime/locate.js';
import * as g from '../../../src/runtime/gestures.js';
import * as guide from '../../../src/runtime/guide.js';
import * as v from '../../../src/runtime/viewers.js';

const settle = v.settle;

// --- readings that are lists ----------------------------------------------------------------------

/* A reading whose value is a comma-separated list: the tile viewer's `fields` and `lane names`, the
   pc plot's `axes`, the trellis plot's `x columns`. */
async function expectReadingContains(page: Page, target: ElementRef, name: string, item: string, negate: boolean): Promise<void> {
  if (item.includes(','))
    throw new Error(`"${item}" holds a comma, so it is never one member of the "${name}" list and a negative claim on it cannot fail; name one member`);
  let last = '';
  const poll = expect.poll(async () => {
    const r = await v.readingOf(page, target, name);
    last = String(r);
    return r instanceof v.MissingReading ? negate : last.split(/\s*,\s*/).includes(item);
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

export const hoverEmptySpace = When('user hovers over empty plot space of {widget}', async (page: Page, target: ElementRef) => {
  const p = await emptySpace(page, target);
  await page.mouse.move(p.x - 3, p.y - 3);
  await page.mouse.move(p.x, p.y);
}, {tier: 'ui', description: 'the pointer on the plot with nothing under it — what a viewer\'s mouse-over group falls back to when no category is hovered'});

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

/** An on-viewer column selector (`div-column-combobox-<property>`), opened: a mouse-down on its
 * caption opens the column grid. The name is the property the selector binds, lowercased — `x`,
 * `y`, `color`, `size`, `category`, `value`. The selectors are hover-revealed, so the pointer goes
 * over the viewer first and stays there: leaving the viewer takes the popup with them. */
async function openColumnSelector(page: Page, target: ElementRef, which: string): Promise<Locator> {
  const loc = await v.viewerLocator(page, target);
  const centre = v.centerOf(await v.hitArea(page, target, 'view'));
  await page.mouse.move(centre.x, centre.y);
  // the selector is named by its property without spaces: "Category 1" is category1
  const selector = loc.locator(`[name="div-column-combobox-${which.toLowerCase().replace(/\s+/g, '')}"]`);
  await selector.waitFor({state: 'visible', timeout: 5000});
  // the guide lights the selector, not the whole viewer the phrase named
  await guide.located(page, selector);
  await g.openColumnSelector(page, selector, false);
  return selector;
}

export const pickInColumnSelector = When('user picks {string} in the {string} column selector of {widget}',
  async (page: Page, column: string, which: string, target: ElementRef) => {
    const selector = await openColumnSelector(page, target, which);
    await g.pickInColumnGrid(page, column, `the "${which}" column selector of ${target.phrase}`, selector);
    await expect(selector.locator('.d4-column-selector-column')).toHaveText(column, {timeout: pollMs(5000)});
    await settle(page, target);
  }, {tier: 'ui', description: 'the column re-picked on the chart itself, the way a user re-picks it'});

/** The same selector, typed into and committed, with no claim about the result: a selector that
 * offers only some columns leaves the column it had, and the feature reads it afterwards. */
export const typeInColumnSelector = When('user types {string} into the {string} column selector of {widget}',
  async (page: Page, column: string, which: string, target: ElementRef) => {
    const selector = await openColumnSelector(page, target, which);
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
      const viewer = (window as any).__bdd.viewerOf(e);
      const list = Array.isArray(viewer.getOptions(true)?.look?.innerViewerLook?.[n as string]) && typeof val === 'string';
      viewer.setOptions({innerViewerLook: {[n as string]: list ? (val === '' ? [] : (val as string).split(/\s*,\s*/)) : val}});
    }, [name, coerced] as [string, unknown]);
    await settle(page, target, 3000);
  }, {tier: 'api', description: 'by the field name the inner look serializes, not by a property caption; a list field takes "A, B, C"'});

export const innerPropertyShouldBe = Then('{string} inner property of {widget} should be {string}',
  async (page: Page, name: string, target: ElementRef, value: string) => {
    const loc = await v.viewerLocator(page, target);
    await expect.poll(() => loc.evaluate((e, n) => {
      const look = (window as any).__bdd.viewerOf(e).getOptions(true)?.look?.innerViewerLook;
      const x = look?.[n as string];
      return x === undefined || x === null ? '' : String(x);
    }, name), {timeout: pollMs(5000), message: `"${name}" of the inner viewer of ${target.phrase}`}).toBe(value);
  }, {description: 'the inner look as the viewer serializes it (`getOptions().look.innerViewerLook`), "" when the field is not written'});

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
    await g.withKeys(page, ['Shift'], async () => {
      await page.mouse.move(points[0].x, points[0].y);
      await page.mouse.down();
      for (const p of points.slice(1))
        await page.mouse.move(p.x, p.y, {steps: 8});
      await page.mouse.up();
    });
  }, {tier: 'ui', description: 'the Lasso Tool must be on; the polygon covers the middle half of the area'});

// --- readings that are prose, and readings that must not be there ----------------------------

/* `should contain` splits a reading on commas and asks for membership — the right question for a
   list (`fields`, `node names`, `stats`), the wrong one for a paragraph. */
async function expectTextIn(page: Page, name: string, target: ElementRef, text: string, negate: boolean): Promise<void> {
  let last = '';
  const holds = async (): Promise<boolean> => {
    const r = await v.readingOf(page, target, name);
    last = String(r);
    return r instanceof v.MissingReading ? negate : last.includes(text);
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

// --- the column list of a "Select columns..." dialog -------------------------------------------------

/* The dialog a column-list property opens (the "..." of Y, Group By, Aggregate): a search box, All /
   None, and a grid whose rows are the table's columns — `__name` holds the column name and `x` the
   check box. The grid reports `cell <r> of <col>` and `text of cell <r> of <col>` like any grid, so a
   row is found by its name and its box is clicked where the grid says it is. */
interface ColumnListRow {row: number; name: string; checked: boolean; x: number; y: number}

async function columnListRows(page: Page, target: ElementRef): Promise<ColumnListRow[]> {
  await v.installViewerRuntime(page);
  const host = (await locate(page, target)).filter({visible: true}).first();
  const grid = host.locator('[name="viewer-Grid"]').filter({visible: true}).first();
  await grid.waitFor({timeout: pollMs(5000)});
  return grid.evaluate((el) => {
    const w = (window as any).__bdd.viewerOf(el);
    const s = w.getWidgetStatus();
    const canvas = (s.parts?.canvas ?? el.querySelector('canvas')) as Element;
    const c = canvas.getBoundingClientRect();
    const rows: ColumnListRow[] = [];
    for (const key of Object.keys(s.hitAreas)) {
      const m = /^cell (\d+) of __name$/.exec(key);
      const box = m ? s.hitAreas[`cell ${m[1]} of x`] : undefined;
      if (!m || !box)
        continue;
      rows.push({row: Number(m[1]), name: String(s.values[`text of cell ${m[1]} of __name`] ?? ''),
        checked: s.values[`text of cell ${m[1]} of x`] === 'true', x: c.x + box.x + box.width / 2, y: c.y + box.y + box.height / 2});
    }
    return rows.sort((a, b) => a.row - b.row);
  });
}

const rowNamed = (rows: ColumnListRow[], name: string): ColumnListRow | undefined => rows.find((r) => r.name === name);
const namesOf = (rows: ColumnListRow[]): string => rows.map((r) => r.name).join(', ');

export const toggleInColumnList = When('user toggles the {string} column in the column list of {element}',
  async (page: Page, column: string, target: ElementRef) => {
    let rows: ColumnListRow[] = [];
    await expect.poll(async () => rowNamed(rows = await columnListRows(page, target), column) !== undefined,
      {timeout: pollMs(5000), message: `no "${column}" row in the column list of ${target.phrase}; it shows: ${namesOf(rows)}`}).toBe(true);
    let row = rowNamed(rows, column)!;
    // a list in a dialog that just opened is still finding its place: the box is clicked where
    // two reads in a row agree it is
    await expect.poll(async () => {
      const again = rowNamed(await columnListRows(page, target), column);
      const still = again !== undefined && again.x === row.x && again.y === row.y;
      row = again ?? row;
      return still;
    }, {timeout: pollMs(3000), message: `the "${column}" row of the column list of ${target.phrase} keeps moving`}).toBe(true);
    await page.mouse.move(row.x, row.y);
    await page.mouse.click(row.x, row.y);
    const flipped = async (ms: number): Promise<boolean> => expect.poll(async () => rowNamed(await columnListRows(page, target), column)?.checked,
      {timeout: pollMs(ms)}).toBe(!row.checked).then(() => true, () => false);
    // the first click on a list that just opened is now and then spent on giving it the focus
    if (!(await flipped(1500)))
      await page.mouse.click(row.x, row.y);
    await expect.poll(async () => rowNamed(await columnListRows(page, target), column)?.checked,
      {timeout: pollMs(5000), message: `the box of "${column}" in the column list of ${target.phrase} did not flip`}).toBe(!row.checked);
  }, {tier: 'ui', description: 'a click on the check box of the row that names the column; done once the box has flipped'});

async function expectChecked(page: Page, column: string, target: ElementRef, checked: boolean): Promise<void> {
  let rows: ColumnListRow[] = [];
  await expect.poll(async () => rowNamed(rows = await columnListRows(page, target), column)?.checked,
    {timeout: pollMs(5000), message: `"${column}" in the column list of ${target.phrase} (it shows: ${namesOf(rows)})`}).toBe(checked);
}

export const checkedInColumnList = Then('the {string} column should be checked in the column list of {element}',
  (page: Page, column: string, target: ElementRef) => expectChecked(page, column, target, true));

export const uncheckedInColumnList = Then('the {string} column should not be checked in the column list of {element}',
  (page: Page, column: string, target: ElementRef) => expectChecked(page, column, target, false));

export const columnListStartsWith = Then('the column list of {element} should start with {string}',
  async (page: Page, target: ElementRef, column: string) => {
    let rows: ColumnListRow[] = [];
    await expect.poll(async () => namesOf((rows = await columnListRows(page, target)).slice(0, 1)),
      {timeout: pollMs(5000), message: `the first row of the column list of ${target.phrase} (it shows: ${namesOf(rows)})`}).toBe(column);
  }, {description: 'the row the search puts first'});

// --- one element inside another -------------------------------------------------------------------

export const liesWithin = Then('{element} should lie within {element}', async (page: Page, inner: ElementRef, outer: ElementRef) => {
  const a = await (await locate(page, inner)).filter({visible: true}).first().boundingBox();
  const b = await (await locate(page, outer)).filter({visible: true}).first().boundingBox();
  if (!a || !b)
    throw new Error(`${!a ? inner.phrase : outer.phrase} has no box on the page`);
  const inside = a.x >= b.x - 0.5 && a.y >= b.y - 0.5 && a.x + a.width <= b.x + b.width + 0.5 && a.y + a.height <= b.y + b.height + 0.5;
  expect(inside, `${inner.phrase} spans ${Math.round(a.x)},${Math.round(a.y)}..${Math.round(a.x + a.width)},${Math.round(a.y + a.height)}, ` +
    `outside ${outer.phrase} at ${Math.round(b.x)},${Math.round(b.y)}..${Math.round(b.x + b.width)},${Math.round(b.y + b.height)}`).toBe(true);
}, {description: 'the first box wholly inside the second (half a pixel of rounding allowed)'});

// --- the text an element shows ----------------------------------------------------------------------

/* `should contain text` reads `textContent`, which keeps the text of children that are hidden — the
   status bar keeps its "Filtered: N" panel in the DOM and only hides it. These read `innerText`: what
   is rendered. */
async function expectShownText(page: Page, target: ElementRef, text: string, negate: boolean): Promise<void> {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  let last = '';
  const poll = expect.poll(async () => (last = (await loc.innerText()).replace(/\s+/g, ' ')).includes(text), {timeout: pollMs(5000)});
  try {
    await (negate ? poll.not : poll).toBe(true);
  }
  catch {
    throw new Error(`${target.phrase} ${negate ? 'still shows' : 'does not show'} "${text}"; it shows "${last}"`);
  }
}

export const showsText = Then('{element} should show the text {string}', (page: Page, target: ElementRef, text: string) =>
  expectShownText(page, target, text, false), {description: 'the rendered text (innerText) contains it; hidden children do not count'});

export const notShowsText = Then('{element} should not show the text {string}', (page: Page, target: ElementRef, text: string) =>
  expectShownText(page, target, text, true), {description: 'the rendered text (innerText) does not contain it'});

// --- docking a viewer by its title bar ---------------------------------------------------------------

/* Dragging a viewer by its title bar shows dock-spawn's wheels: a compass of `.dock-wheel-left|right|
   top|down|fill` items over the panel under the pointer, and a second `.dock-wheel-base` with one item
   at each edge of the view. Both carry the same side classes; only the compass holds the fill item.
   dock-spawn marks the item under the pointer on mouseover, and the drop waits for that mark. */
const DOCK_SIDES: Record<string, string> = {left: 'left', right: 'right', top: 'top', bottom: 'down'};
type Rect = {x: number; y: number; width: number; height: number};
const centreOf = (r: Rect) => ({x: r.x + r.width / 2, y: r.y + r.height / 2});

async function panelOf(page: Page, target: ElementRef): Promise<Locator> {
  const loc = await v.viewerLocator(page, target);
  return loc.locator('xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]');
}

async function dockViewer(page: Page, target: ElementRef, side: string, over: ElementRef | null): Promise<void> {
  const wheelSide = DOCK_SIDES[side];
  if (!wheelSide)
    throw new Error(`a viewer docks to the left, right, top or bottom — not "${side}"`);
  const title = (await panelOf(page, target)).locator('.panel-titlebar').first();
  const t = await title.boundingBox();
  if (!t)
    throw new Error(`${target.phrase} has no title bar to drag`);
  const view = await (await locate(page, over ?? el('open tableview'))).filter({visible: true}).first().boundingBox();
  if (!view)
    throw new Error(`${over ? over.phrase : 'the table view'} has no box to dock against`);
  const aim = centreOf(view);
  const what = over ? `the "${side}" item of the dock compass over ${over.phrase}` : `the "${side}" edge item of the view`;
  await page.mouse.move(t.x + Math.min(30, t.width / 4), t.y + t.height / 2);
  await page.mouse.down();
  await page.mouse.move(t.x + Math.min(30, t.width / 4) + 20, t.y + t.height / 2 + 20, {steps: 4});
  await page.mouse.move(aim.x, aim.y, {steps: 8});
  const wheel = over ? '.dock-wheel-base:has(.dock-wheel-fill)' : '.dock-wheel-base:not(:has(.dock-wheel-fill))';
  const item = page.locator(`${wheel} .dock-wheel-item.dock-wheel-${wheelSide}`).filter({visible: true}).first();
  try {
    await item.waitFor({timeout: pollMs(5000)}).catch(() => {
      throw new Error(`dragging ${target.phrase} showed no ${what}`);
    });
    const box = await item.boundingBox();
    if (!box)
      throw new Error(`${what} has no box`);
    const p = centreOf(box);
    await page.mouse.move(p.x, p.y, {steps: 6});
    await expect(item, `${what} under the pointer`).toHaveClass(new RegExp(`dock-wheel-${wheelSide}-icon-hover`), {timeout: pollMs(3000)});
  }
  catch (e) {
    await page.mouse.up();
    throw e;
  }
  await page.mouse.up();
  await v.settleAll(page);
}

export const dockToViewEdge = When('user docks {widget} to the {word} edge of the view', (page: Page, target: ElementRef, side: string) =>
  dockViewer(page, target, side, null), {tier: 'ui', description: 'drags the viewer by its title bar onto the dock wheel item at that edge of the table view'});

export const dockBesideViewer = When('user docks {widget} to the {word} side of {widget}', (page: Page, target: ElementRef, side: string, other: ElementRef) =>
  dockViewer(page, target, side, other), {tier: 'ui', description: 'drags the viewer by its title bar onto the dock wheel that appears over the other viewer, on the item of that side'});

async function panelBox(page: Page, target: ElementRef): Promise<Rect> {
  const b = await (await panelOf(page, target)).boundingBox();
  if (!b)
    throw new Error(`${target.phrase} is not on screen`);
  return b;
}

const EDGE = 6;

/** Whether the viewer's panel runs along an edge of the area the docked panels cover; the reason
 * names both boxes. Read in a poll: dock-spawn relays the panels out after the drop. */
async function dockedArea(page: Page): Promise<{x: number; y: number; right: number; bottom: number}> {
  const view = (await locate(page, el('open tableview'))).filter({visible: true}).first();
  const root = await view.boundingBox();
  const views = await view.locator('.panel-base').evaluateAll((els) =>
    els.map((e) => e.getBoundingClientRect()).filter((r) => r.width > 0).map((r) => ({x: r.x, y: r.y, width: r.width, height: r.height})));
  if (!root)
    throw new Error('the table view has no box');
  // the grid is the view's document, not a dock panel: the view's own box gives the left, right and bottom
  return {x: Math.min(root.x, ...views.map((r) => r.x)), y: Math.min(...views.map((r) => r.y)),
    right: Math.max(root.x + root.width, ...views.map((r) => r.x + r.width)), bottom: Math.max(root.y + root.height, ...views.map((r) => r.y + r.height))};
}

async function alongEdge(page: Page, target: ElementRef, side: string): Promise<{ok: boolean; why: string}> {
  const a = await panelBox(page, target);
  const all = await dockedArea(page);
  const ok = side === 'right' ? Math.abs(a.x + a.width - all.right) <= EDGE && a.y <= all.y + EDGE && a.y + a.height >= all.bottom - EDGE :
    side === 'left' ? Math.abs(a.x - all.x) <= EDGE && a.y <= all.y + EDGE && a.y + a.height >= all.bottom - EDGE :
      side === 'bottom' ? Math.abs(a.y + a.height - all.bottom) <= EDGE && a.x <= all.x + EDGE && a.x + a.width >= all.right - EDGE :
        side === 'top' ? Math.abs(a.y - all.y) <= EDGE && a.x <= all.x + EDGE && a.x + a.width >= all.right - EDGE : undefined;
  if (ok === undefined)
    throw new Error(`a view has a left, right, top or bottom edge — not "${side}"`);
  return {ok, why: `${target.phrase} spans ${Math.round(a.x)},${Math.round(a.y)}..${Math.round(a.x + a.width)},${Math.round(a.y + a.height)}; ` +
    `the docked viewers span ${Math.round(all.x)},${Math.round(all.y)}..${Math.round(all.right)},${Math.round(all.bottom)}`};
}

/** Whether the viewer's panel touches both edges of a corner of the docked area, whatever its length. */
async function inCorner(page: Page, target: ElementRef, vertical: string, horizontal: string): Promise<{ok: boolean; why: string}> {
  if (!['top', 'bottom'].includes(vertical) || !['left', 'right'].includes(horizontal))
    throw new Error(`a view has a top or bottom, left or right corner — not "${vertical} ${horizontal}"`);
  const a = await panelBox(page, target);
  const all = await dockedArea(page);
  const ok = (vertical === 'top' ? Math.abs(a.y - all.y) : Math.abs(a.y + a.height - all.bottom)) <= EDGE &&
    (horizontal === 'left' ? Math.abs(a.x - all.x) : Math.abs(a.x + a.width - all.right)) <= EDGE;
  return {ok, why: `${target.phrase} spans ${Math.round(a.x)},${Math.round(a.y)}..${Math.round(a.x + a.width)},${Math.round(a.y + a.height)}; ` +
    `the docked viewers span ${Math.round(all.x)},${Math.round(all.y)}..${Math.round(all.right)},${Math.round(all.bottom)}`};
}

export const dockedInCorner = Then('{widget} should be docked in the {word} {word} corner of the view', (page: Page, target: ElementRef, vertical: string, horizontal: string) =>
  expectDocked(page, () => inCorner(page, target, vertical, horizontal), false),
{description: 'the viewer\'s panel touches both edges of that corner of the docked area, whatever its length (6 px slack)'});

export const notDockedInCorner = Then('{widget} should not be docked in the {word} {word} corner of the view', (page: Page, target: ElementRef, vertical: string, horizontal: string) =>
  expectDocked(page, () => inCorner(page, target, vertical, horizontal), true),
{description: 'the state before a drag'});

async function besides(page: Page, target: ElementRef, where: string, other: ElementRef): Promise<{ok: boolean; why: string}> {
  const a = await panelBox(page, target);
  const b = await panelBox(page, other);
  const ok = where === 'below' ? Math.abs(a.y - (b.y + b.height)) <= EDGE && Math.abs(a.x - b.x) <= EDGE && Math.abs(a.width - b.width) <= EDGE :
    where === 'above' ? Math.abs(a.y + a.height - b.y) <= EDGE && Math.abs(a.x - b.x) <= EDGE && Math.abs(a.width - b.width) <= EDGE :
      where === 'left-of' ? Math.abs(a.x + a.width - b.x) <= EDGE && Math.abs(a.y - b.y) <= EDGE && Math.abs(a.height - b.height) <= EDGE :
        where === 'right-of' ? Math.abs(a.x - (b.x + b.width)) <= EDGE && Math.abs(a.y - b.y) <= EDGE && Math.abs(a.height - b.height) <= EDGE : undefined;
  if (ok === undefined)
    throw new Error(`a viewer is docked below, above, left-of or right-of another — not "${where}"`);
  return {ok, why: `${target.phrase} spans ${Math.round(a.x)},${Math.round(a.y)} ${Math.round(a.width)}x${Math.round(a.height)}; ` +
    `${other.phrase} ${Math.round(b.x)},${Math.round(b.y)} ${Math.round(b.width)}x${Math.round(b.height)}`};
}

async function expectDocked(page: Page, read: () => Promise<{ok: boolean; why: string}>, negate: boolean): Promise<void> {
  await v.settleAll(page);
  let why = '';
  try {
    await expect.poll(async () => { const r = await read(); why = r.why; return r.ok; }, {timeout: pollMs(5000)}).toBe(!negate);
  }
  catch {
    throw new Error(`${negate ? 'still docked there' : 'not docked there'}: ${why}`);
  }
}

export const dockedAtViewEdge = Then('{widget} should be docked along the {word} edge of the view', (page: Page, target: ElementRef, side: string) =>
  expectDocked(page, () => alongEdge(page, target, side), false),
{description: 'the viewer\'s panel runs the whole length of that edge of the docked area and touches it (6 px slack)'});

export const notDockedAtViewEdge = Then('{widget} should not be docked along the {word} edge of the view', (page: Page, target: ElementRef, side: string) =>
  expectDocked(page, () => alongEdge(page, target, side), true),
{description: 'the state before a drag, so that the docking after it is the drag\'s doing'});

export const dockedBeside = Then('{widget} should be docked {word} {widget}', (page: Page, target: ElementRef, where: string, other: ElementRef) =>
  expectDocked(page, () => besides(page, target, where, other), false),
{description: 'below / above: the two panels touch and share their width; left-of / right-of: they touch and share their height (6 px slack)'});

export const notDockedBeside = Then('{widget} should not be docked {word} {widget}', (page: Page, target: ElementRef, where: string, other: ElementRef) =>
  expectDocked(page, () => besides(page, target, where, other), true),
{description: 'the state before a drag'});

export const openViewerHelp = When('user opens the help of {widget}', async (page: Page, target: ElementRef) => {
  const icon = (await panelOf(page, target)).locator('.panel-titlebar [name="icon-font-icon-help"]').first();
  await (await panelOf(page, target)).locator('.panel-titlebar').first().hover();
  await icon.click();
}, {tier: 'ui', description: 'the "?" icon of the viewer\'s title bar (shown on hover)'});
// --- where an area sits -------------------------------------------------------------------------------

/* Hit areas are read in the viewer's own canvas coordinates here, so a rectangle remembered before
   a gesture, a property change or a resize-and-restore compares with the one read after it whatever
   the page did around the viewer. An annotation region's title is the case: it must stay where the
   viewer put it when the title is clicked, when the plot style changes, when the viewer regrows. */
const areaPlaces = new WeakMap<Page, Map<string, v.Box>>();
const placeKey = (target: ElementRef, area: string): string => `${target.phrase}|${area.toLowerCase().replace(/[^a-z0-9]/g, '')}`;
const fmtBox = (b: v.Box | undefined): string => b ? `${Math.round(b.x)},${Math.round(b.y)} ${Math.round(b.width)}x${Math.round(b.height)}` : 'none';
const sameBox = (a: v.Box, b: v.Box): boolean => [a.x - b.x, a.y - b.y, a.width - b.width, a.height - b.height].every((d) => Math.abs(d) <= 1);

const areaRect = async (page: Page, target: ElementRef, area: string): Promise<v.Box> => (await v.areaRects(page, target, [area]))[0];

export const rememberAreaPlace = When('user remembers the place of the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  if (!areaPlaces.has(page))
    areaPlaces.set(page, new Map());
  areaPlaces.get(page)!.set(placeKey(target, area), await areaRect(page, target, area));
}, {tier: 'api', description: 'the rectangle the viewer reports for the area now, kept for "placed as remembered" after any number of changes'});

export const areaPlacedAsRemembered = Then('the {string} area of {widget} should be placed as remembered', async (page: Page, area: string, target: ElementRef) => {
  const remembered = areaPlaces.get(page)?.get(placeKey(target, area));
  if (remembered === undefined)
    throw new Error(`no place was remembered for the "${area}" area of ${target.phrase} — "user remembers the place of the <area> area of <widget>" comes first`);
  const now = await areaRect(page, target, area);
  expect(sameBox(remembered, now), `the "${area}" area of ${target.phrase} was remembered at ${fmtBox(remembered)} and is now at ${fmtBox(now)}`).toBe(true);
}, {description: 'the same rectangle within a pixel, read once the viewer is quiet — a region title that a click or a resize-and-restore must not move'});

type Relation = 'inside' | 'outside' | 'above' | 'below' | 'left' | 'right';
const RELATION: Record<Relation, (a: v.Box, b: v.Box) => boolean> = {
  inside: (a, b) => a.x >= b.x - 1 && a.y >= b.y - 1 && a.x + a.width <= b.x + b.width + 1 && a.y + a.height <= b.y + b.height + 1,
  outside: (a, b) => a.x + a.width <= b.x + 1 || b.x + b.width <= a.x + 1 || a.y + a.height <= b.y + 1 || b.y + b.height <= a.y + 1,
  above: (a, b) => a.y + a.height <= b.y + 1,
  below: (a, b) => a.y >= b.y + b.height - 1,
  left: (a, b) => a.x + a.width <= b.x + 1,
  right: (a, b) => a.x >= b.x + b.width - 1,
};

async function expectAreaRelation(page: Page, target: ElementRef, a: string, relation: string, b: string): Promise<void> {
  if (!(relation in RELATION))
    throw new Error(`an area lies inside, outside, above or below another, or to the left or right of it — not "${relation}"`);
  const [ra, rb] = await v.areaRects(page, target, [a, b]);
  const where = relation === 'left' || relation === 'right' ? `to the ${relation} of` : relation;
  expect(RELATION[relation as Relation](ra, rb), `the "${a}" area of ${target.phrase} is at ${fmtBox(ra)} and does not lie ${where} ` +
    `the "${b}" area at ${fmtBox(rb)}`).toBe(true);
}

export const areaLies = Then('the {string} area of {widget} should lie {word} the {string} area',
  (page: Page, a: string, target: ElementRef, relation: string, b: string) => expectAreaRelation(page, target, a, relation, b),
  {description: 'inside, outside, above or below — a region title drawn in its region, or in the strip the layout reserved above the plot'});

export const areaLiesBeside = Then('the {string} area of {widget} should lie to the {word} of the {string} area',
  (page: Page, a: string, target: ElementRef, side: string, b: string) => expectAreaRelation(page, target, a, side, b),
  {description: 'left or right — a vertical band title in the strip the layout reserved to the right of the plot'});
