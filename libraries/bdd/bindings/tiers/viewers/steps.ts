/* The `viewers` tier: viewers on the current table view — adding them, their properties by
   caption, their context menus and hit areas, their canvas, their events, and the page's error
   floor. A project opts in with `{"tiers": ["viewers"]}` in its bdd.config.json.
   A step that reads or changes a viewer takes `{widget}` — a phrase naming a viewer or a widget —
   so the general shapes ("X of Y should be Z") stay free for what they say; the context menu and
   the pointer steps take any element.
   Every step here asks the platform (properties, hit areas, render and menu events) instead of
   guessing pixels or sleeping; viewers render immediately (no debounce) on a bdd page. */
import {expect, Page} from '@playwright/test';
import {Given, Then, When} from '../../../src/registry.js';
import type {ElementRef} from '../../../src/runtime/args.js';
import {normalizeKey} from '../../../src/runtime/gestures.js';
import {atFeatureEnd, takeErrors} from '../../../src/runtime/harness.js';
import * as v from '../../../src/runtime/viewers.js';

declare const grok: any;

const PATH = '"Group > Item" (groups separated by ">" or "|")';

export const openToolbox = When('user opens toolbox', async (page: Page) => {
  const shown = await page.evaluate(() => grok.shell.windows.showToolbox as boolean);
  if (!shown)
    await page.locator('[name="Toolbox"]').first().click();
  await page.locator('[name="div-section--Viewers"]').first().waitFor();
}, {tier: 'ui'});

export const viewerAdded = Then('{viewer} viewer should be added to the open tableview', async (page: Page, viewer: string) => {
  await expect(page.locator(`[name="viewer-${viewer.replace(/\s+/g, '-')}" i]`).first()).toBeVisible();
  const types: string[] = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map((x: any) => String(x.type)));
  expect(types.map((t) => t.toLowerCase())).toContain(viewer.toLowerCase());
}, {tier: 'api', description: 'both the DOM and grok.shell.tv.viewers'});

export const viewerCount = Then('the open tableview should have {int} {viewer} viewer(s)', async (page: Page, count: number, viewer: string) => {
  await expect.poll(() => page.evaluate((t) => Array.from(grok.shell.tv.viewers).filter((x: any) => String(x.type).toLowerCase() === t).length, viewer.toLowerCase()),
    {message: `${viewer} viewers of the open table view`}).toBe(count);
}, {tier: 'api', description: 'how many viewers of that type the view holds (grok.shell.tv.viewers) — what "added" and "closed" rest on'});

// --- adding and configuring -------------------------------------------------------------------------

export const addViewer = Given('user adds (a ){viewer} viewer', (page: Page, viewer: string) => v.addViewer(page, viewer),
  {tier: 'api', description: 'grok.shell.tv.addViewer — the UI path is the toolbox icon'});

export const addViewerWith = Given('user adds (a ){viewer} viewer with:', async (page: Page, viewer: string, table: string[][]) => {
  await v.addViewer(page, viewer);
  await v.setProperties(page, {phrase: `${viewer} viewer`}, table.map(([caption, value]) => [caption, value]));
}, {tier: 'api', description: '| property caption | value | rows applied right after adding'});

export const setProperty = When('user sets {string} property of {widget} to {string}', (page: Page, caption: string, target: ElementRef, value: string) =>
  v.setProperties(page, target, [[caption, value]]),
  {tier: 'api', description: 'a property by its caption ("Value", "Show Markers") or name; true/false, numbers, #rrggbb colors, "" for none, \\n for a line break'});

export const setProperties = When('user sets properties of {widget}:', (page: Page, target: ElementRef, table: string[][]) =>
  v.setProperties(page, target, table.map(([caption, value]) => [caption, value])),
  {tier: 'api', description: '| property caption | value | rows, one repaint for the group'});

export const propertyShouldBe = Then('{string} property of {widget} should be {string}', (page: Page, caption: string, target: ElementRef, value: string) =>
  v.expectProperty(page, target, caption, value), {description: 'read back through the property bag: "" for none, true/false, numbers as written'});

export const propertyShouldNotBe = Then('{string} property of {widget} should not be {string}', (page: Page, caption: string, target: ElementRef, value: string) =>
  v.expectProperty(page, target, caption, value, true));

export const propertyShouldContain = Then('{string} property of {widget} should contain {string}', async (page: Page, caption: string, target: ElementRef, text: string) => {
  await expect.poll(() => v.readProperty(page, target, caption), {message: `"${caption}" property of ${target.phrase}`}).toContain(text);
}, {description: 'a text property (a description a command writes its settings into) by substring'});

export const propertiesShouldBe = Then('properties of {widget} should be:', async (page: Page, target: ElementRef, table: string[][]) => {
  for (const [caption, value] of table)
    await v.expectProperty(page, target, caption, value);
}, {description: '| property caption | value | rows, each read back through the property bag'});

export const saveLayout = When('user saves the layout of the current table view', (page: Page) => v.saveLayout(page),
  {tier: 'api', description: 'tv.saveLayout(), kept for "loads the saved layout" later in the feature'});

export const saveLayoutToServer = When('user saves the layout of the current table view to the server', async (page: Page) => {
  const id = await v.saveLayoutToServer(page);
  atFeatureEnd(page, () => v.deleteLayout(page, id));
}, {tier: 'api', description: 'dapi.layouts.save; "loads the saved layout" then fetches what the server stored, so the round-trip covers its serialization; deleted when the feature ends'});

export const loadLayout = When('user loads the saved layout', (page: Page) => v.loadLayout(page), {tier: 'api'});

// --- context menus and hit areas ------------------------------------------------------------------

export const pickFromContextMenu = When('user picks {string} from the context menu of {element}', async (page: Page, path: string, target: ElementRef) => {
  await v.openContextMenuOf(page, target);
  await v.pickMenuPath(page, path);
}, {tier: 'ui', description: `right-clicks the element (its "view" area when it reports one) and clicks ${PATH}`});

export const pickFromAreaContextMenu = When('user picks {string} from the context menu of the {string} area of {widget}',
  async (page: Page, path: string, area: string, target: ElementRef) => {
    await v.openContextMenuOf(page, target, area);
    await v.pickMenuPath(page, path);
  }, {tier: 'ui', description: `right-clicks a named hit area of the viewer and clicks ${PATH}`});

export const openContextMenu = When('user opens the context menu of {element}', async (page: Page, target: ElementRef) => {
  await v.openContextMenuOf(page, target);
}, {tier: 'ui', description: 'leaves it open for menu item checks'});

export const rightClickArea = When('user right-clicks on the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  await v.openContextMenuOf(page, target, area);
}, {tier: 'ui', description: 'a named hit area the viewer reports (grok-bdd outline lists them); leaves the menu open'});

export const closeContextMenu = When('user closes the context menu', (page: Page) => v.closeContextMenu(page), {tier: 'ui'});

export const menuLists = Then('the open menu should list {string}', async (page: Page, path: string) => {
  const parts = path.split(/\s*[>|]\s*/).filter((s) => s.length > 0);
  const item = parts.pop();
  await expect.poll(() => v.menuLabels(page, parts.join(' > ')), {message: `the items under "${parts.join(' > ') || 'the menu'}"`}).toContain(item);
}, {tier: 'ui', description: `the last segment of ${PATH} is among the items, the ones before it are the groups opened to reach them — the menu must already be open`});

export const menuDoesNotList = Then('the open menu should not list {string}', async (page: Page, path: string) => {
  const parts = path.split(/\s*[>|]\s*/).filter((s) => s.length > 0);
  const item = parts.pop();
  // the group is opened first, so "not listed" cannot pass on a menu that never opened
  await expect.poll(() => v.menuLabels(page, parts.join(' > ')), {message: `the items under "${parts.join(' > ') || 'the menu'}"`}).not.toContain(item);
}, {tier: 'ui'});

export const clickArea = When('user clicks on the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  const c = v.centerOf(await v.hitArea(page, target, area, true));
  await page.mouse.click(c.x, c.y);
}, {tier: 'ui'});

export const doubleClickArea = When('user double-clicks on the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  const c = v.centerOf(await v.hitArea(page, target, area, true));
  await page.mouse.dblclick(c.x, c.y);
}, {tier: 'ui'});

export const hoverArea = When('user hovers over the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  const c = v.centerOf(await v.hitArea(page, target, area, true));
  await page.mouse.move(c.x - 3, c.y - 3);
  await page.mouse.move(c.x, c.y);
}, {tier: 'ui', description: 'snapshots the canvas first, so "should not have repainted" can follow'});

export const clickAreaHolding = When('user clicks on the {string} area of {widget} holding {key}', async (page: Page, area: string, target: ElementRef, key: string) => {
  const c = v.centerOf(await v.hitArea(page, target, area, true));
  // a chord is held key by key: Playwright's down() takes one key, never "Control+Shift"
  const keys = normalizeKey(key).split('+');
  for (const k of keys)
    await page.keyboard.down(k);
  try {
    await page.mouse.click(c.x, c.y);
  }
  finally {
    for (const k of keys.reverse())
      await page.keyboard.up(k);
  }
}, {tier: 'ui', description: 'a click with a key or a chord held: Control adds to the selection, Shift extends it, Control+Shift removes'});

export const dragSelectionOverArea = When('user drags a selection box over the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  const b = await v.hitArea(page, target, area, true);
  await page.keyboard.down('Shift');
  await page.mouse.move(b.x + b.width * 0.1, b.y + b.height * 0.1);
  await page.mouse.down();
  await page.mouse.move(b.x + b.width * 0.9, b.y + b.height * 0.9, {steps: 3});
  await page.mouse.up();
  await page.keyboard.up('Shift');
}, {tier: 'ui', description: 'a Shift-drag across the inner 80% of the area — the platform\'s rectangle selection'});

export const dragSelectionBetweenAreas = When('user drags a selection box from the {string} area to the {string} area of {widget}',
  async (page: Page, from: string, to: string, target: ElementRef) => {
    const a = v.centerOf(await v.hitArea(page, target, from, true));
    const b = v.centerOf(await v.hitArea(page, target, to));
    await page.keyboard.down('Shift');
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move(b.x, b.y, {steps: 3});
    await page.mouse.up();
    await page.keyboard.up('Shift');
  }, {tier: 'ui', description: 'a Shift-drag from the centre of one area to the centre of another — every area the rectangle touches is covered'});

export const dragAcrossArea = When('user drags across the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  const b = await v.hitArea(page, target, area, true);
  await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2);
  await page.mouse.down();
  await page.mouse.move(b.x + b.width * 0.75, b.y + b.height * 0.3, {steps: 2});
  await page.mouse.up();
}, {tier: 'ui', description: 'a plain drag from the centre of the area towards its upper right — a rotation on a 3D plot, a pan on a chart'});

export const dragAreaToArea = When('user drags the {string} area of {widget} to the {string} area', (page: Page, from: string, target: ElementRef, to: string) =>
  v.dragArea(page, target, from, to), {tier: 'ui', description: 'a plain drag from the centre of one hit area to the centre of another — a column header to a new place, a range handle onto a bin'});

export const dragAreaBy = When('user drags the {string} area of {widget} by {int} pixels to the {word}', (page: Page, area: string, target: ElementRef, px: number, direction: string) =>
  v.dragAreaBy(page, target, area, px, direction), {tier: 'ui', description: 'left, right, up or down — a column or row resizer, a splitter'});

export const dragDeselectionOverArea = When('user drags a deselection box over the {string} area of {widget}', (page: Page, area: string, target: ElementRef) =>
  v.dragBoxOverArea(page, target, area, ['Control', 'Shift']), {tier: 'ui', description: 'a Control+Shift-drag across the inner 80% of the area — the platform removes the rows inside from the selection'});

export const dragZoomOverArea = When('user drags a zoom box over the {string} area of {widget}', (page: Page, area: string, target: ElementRef) =>
  v.dragBoxOverArea(page, target, area, ['Alt']), {tier: 'ui', description: 'an Alt-drag across the inner 80% of the area — the platform zooms into the rectangle'});

export const enterIntoArea = When('user enters {string} into the {string} area of {widget}', (page: Page, text: string, area: string, target: ElementRef) =>
  v.typeIntoArea(page, target, area, text), {tier: 'ui', description: 'a hit area that holds an editor (a range input, a form field): a click on it, select all, the text, Enter'});

export const wheelOverArea = When('user scrolls the mouse wheel {word} over the {string} area of {widget}', async (page: Page, direction: string, area: string, target: ElementRef) => {
  if (direction !== 'up' && direction !== 'down')
    throw new Error(`the wheel scrolls up or down, not "${direction}"`);
  const c = v.centerOf(await v.hitArea(page, target, area, true));
  await page.mouse.move(c.x, c.y);
  await page.mouse.wheel(0, direction === 'up' ? -600 : 600);
}, {tier: 'ui', description: 'up or down, a few notches, with the pointer at the centre of the area'});

export const pointerAway = When('user moves the pointer away from {element}', async (page: Page, target: ElementRef) => {
  const box = await (await v.viewerLocator(page, target)).boundingBox();
  if (!box)
    throw new Error(`${target.phrase}: no box to leave`);
  await page.mouse.move(box.x + box.width / 2, Math.max(0, box.y - 40));
}, {tier: 'ui', description: 'above the element — tooltips close on leaving'});

// --- size ------------------------------------------------------------------------------------------

export const resizeTo = When('user resizes {widget} to {int} by {int}', (page: Page, target: ElementRef, width: number, height: number) =>
  v.resize(page, target, width, height), {tier: 'ui', description: 'width by height in pixels; "restores the size" undoes it'});

export const resizeWidth = When('user resizes {widget} to {int} wide', (page: Page, target: ElementRef, width: number) =>
  v.resize(page, target, width, null), {tier: 'ui'});

export const restoreSize = When('user restores the size of {widget}', (page: Page, target: ElementRef) => v.restoreSize(page, target), {tier: 'ui'});

// --- the canvas ------------------------------------------------------------------------------------

export const takeSnapshot = When('user takes a snapshot of {widget}', async (page: Page, target: ElementRef) => {
  await v.snapshot(page, target);
}, {tier: 'api', description: 'the baseline for "should have repainted" — every property set, menu pick and resize takes one by itself'});

export const repainted = Then('{widget} should have repainted', (page: Page, target: ElementRef) => v.expectRepainted(page, target),
  {description: 'the canvas differs from the snapshot taken before the last change (a property set, a menu pick, a gesture, a data step); checks never move the snapshot, so several can follow one change'});

export const repaintedBy = Then('{widget} should have repainted by at least {int} pixels', (page: Page, target: ElementRef, px: number) =>
  v.expectRepainted(page, target, px), {description: 'a change of at least that many pixels — for a toggle with no shape of its own (an axis, a selector)'});

export const lessInk = Then('{widget} should have less ink than before', (page: Page, target: ElementRef) => v.expectInk(page, target, 'less'),
  {description: 'fewer painted pixels than the snapshot before the last change'});

export const moreInk = Then('{widget} should have more ink than before', (page: Page, target: ElementRef) => v.expectInk(page, target, 'more'));

export const painted = Then('{widget} should be painted', (page: Page, target: ElementRef) => v.expectInk(page, target, 'some'),
  {description: 'the canvas has painted pixels'});

export const areaRepainted = Then('the {string} area of {widget} should have repainted', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaRepainted(page, target, area),
  {description: 'the pixels inside one hit area differ from the snapshot before the last change — that area repainted, not merely the canvas around it'});

export const areaLessInk = Then('the {string} area of {widget} should have less ink than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaInk(page, target, area, 'less'), {description: 'fewer painted pixels inside a hit area than the snapshot before the last change had there'});

export const areaMoreInk = Then('the {string} area of {widget} should have more ink than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaInk(page, target, area, 'more'));

export const areaColor = Then('the {string} area of {widget} should contain the color {string}', (page: Page, area: string, target: ElementRef, color: string) =>
  v.expectAreaColor(page, target, area, color), {description: 'pixels of a #rrggbb color (a shade of anti-aliasing allowed) inside a hit area'});

export const areasDiffer = Then('the {string} and {string} areas of {widget} should be painted in different colors',
  (page: Page, a: string, b: string, target: ElementRef) => v.expectAreasDiffer(page, target, a, b),
  {description: 'one area has a color the other does not — per-category coloring, not chrome'});

export const areasSame = Then('the {string} and {string} areas of {widget} should be painted in the same colors',
  (page: Page, a: string, b: string, target: ElementRef) => v.expectSameColors(page, target, a, b),
  {description: 'every significant color of either area has a match in the other — a linked color coding, a swatch and its category'});

export const areaNotColor = Then('the {string} area of {widget} should not contain the color {string}', (page: Page, area: string, target: ElementRef, color: string) =>
  v.expectAreaNotColor(page, target, area, color), {description: 'no pixel of the #rrggbb color (nor a shade of it) inside the hit area, read once'});

export const areaAtLeastTall = Then('the {string} area of {widget} should be at least {int} pixels tall', (page: Page, area: string, target: ElementRef, px: number) =>
  v.expectAreaSize(page, target, area, 'tall', px), {description: 'the rectangle the viewer reports for the area, in CSS pixels'});

export const areaAtLeastWide = Then('the {string} area of {widget} should be at least {int} pixels wide', (page: Page, area: string, target: ElementRef, px: number) =>
  v.expectAreaSize(page, target, area, 'wide', px));

export const areaTaller = Then('the {string} area of {widget} should be taller than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaGrew(page, target, area, 'taller'), {description: 'against the rectangle at the snapshot before the last change'});

export const areaWider = Then('the {string} area of {widget} should be wider than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaGrew(page, target, area, 'wider'));

export const areaShorter = Then('the {string} area of {widget} should be shorter than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaGrew(page, target, area, 'shorter'));

export const areaNarrower = Then('the {string} area of {widget} should be narrower than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaGrew(page, target, area, 'narrower'));

export const areaColors = Then('the {string} area of {widget} should be painted in at least {int} colors',
  (page: Page, area: string, target: ElementRef, count: number) => v.expectAreaColors(page, target, area, count),
  {description: 'distinct hues covering some pixels each, greys and white aside — a grid cell whose letters take their colors from the data, not a text cell'});

export const boundTable = Then('{widget} should be bound to table {string}', (page: Page, target: ElementRef, name: string) => v.expectBoundTable(page, target, name),
  {description: 'the table the viewer draws (viewer.dataFrame), not the Table property it was asked for'});

export const showsRows = Then('{widget} should show {int} rows', (page: Page, target: ElementRef, count: number) =>
  v.expectReading(page, target, 'rows shown', 'equal', count), {description: 'the rows the viewer draws after its own filter and the table\'s (the "rows shown" reading)'});

export const showsFewerRows = Then('{widget} should show fewer rows than before', (page: Page, target: ElementRef) =>
  v.expectReading(page, target, 'rows shown', 'lower'), {description: 'against the snapshot before the last change'});

export const showsMoreRows = Then('{widget} should show more rows than before', (page: Page, target: ElementRef) =>
  v.expectReading(page, target, 'rows shown', 'higher'));

export const readingIs = Then('the {string} reading of {widget} should be {float}', (page: Page, name: string, target: ElementRef, value: number) =>
  v.expectReading(page, target, name, 'equal', value),
{description: 'a reading the viewer reports (getWidgetStatus().values): "rows shown", the bar chart\'s "bars" / "stack segments" / "clipped bars", the 3D scatter plot\'s "camera distance"'});

export const readingReads = Then('the {string} reading of {widget} should be {string}', async (page: Page, name: string, target: ElementRef, value: string) => {
  await expect.poll(async () => String(await v.readValue(page, target, name)), {message: `"${name}" reading of ${target.phrase}`}).toBe(value);
}, {description: 'a reading as text (a source column, a signature, a true/false flag) by exact value'});

/** Two readings of one viewer, compared as text: 'same', or what each holds when they disagree. */
async function compareReadings(page: Page, target: ElementRef, a: string, b: string): Promise<string> {
  const x = String(await v.readValue(page, target, a));
  const y = String(await v.readValue(page, target, b));
  return x === y ? 'same' : `${a} is ${x}, ${b} is ${y}`;
}

export const readingsEqual = Then('the {string} and {string} readings of {widget} should be the same', async (page: Page, a: string, b: string, target: ElementRef) => {
  await expect.poll(() => compareReadings(page, target, a, b), {message: `${target.phrase}`}).toBe('same');
}, {description: 'two readings of the same viewer, as text — the colour of a linked column\'s cell against its source\'s'});

export const readingsDiffer = Then('the {string} and {string} readings of {widget} should differ', async (page: Page, a: string, b: string, target: ElementRef) => {
  await expect.poll(() => compareReadings(page, target, a, b), {message: `${target.phrase}`}).not.toBe('same');
});

export const readingDoesNotRead = Then('the {string} reading of {widget} should not be {string}', async (page: Page, name: string, target: ElementRef, value: string) => {
  await expect.poll(async () => String(await v.readValue(page, target, name)), {message: `"${name}" reading of ${target.phrase}`}).not.toBe(value);
});

export const readingFinite = Then('the {string} reading of {widget} should be a finite number', async (page: Page, name: string, target: ElementRef) => {
  await expect.poll(async () => {
    const value = await v.readValue(page, target, name);
    return typeof value === 'number' && isFinite(value) ? 'finite' : String(value);
  }, {message: `"${name}" reading of ${target.phrase}`}).toBe('finite');
}, {description: 'not null, NaN or Infinity — a viewer that let a NaN or an Infinity of the data into its axis reports no number there'});

export const readingBetween = Then('the {string} reading of {widget} should be between {float} and {float}',
  async (page: Page, name: string, target: ElementRef, lo: number, hi: number) => {
    await expect.poll(async () => Number(await v.readValue(page, target, name)),
      {message: `"${name}" reading of ${target.phrase}`}).toBeGreaterThanOrEqual(lo);
    await expect.poll(async () => Number(await v.readValue(page, target, name)),
      {message: `"${name}" reading of ${target.phrase}`}).toBeLessThanOrEqual(hi);
  }, {description: 'a reading that carries float noise or depends on the layout, bounded on both sides'});

export const pickColorSwatch = When('user picks the color {string} in the color picker dialog', async (page: Page, hex: string) => {
  const swatch = page.locator(`.d4-dialog [name="color-${hex.replace('#', '')}" i]`).filter({visible: true}).first();
  await expect(swatch, `a "${hex}" swatch in the open colour dialog`).toBeVisible();
  await swatch.click();
}, {tier: 'ui', description: 'a swatch of the open colour dialog by its #rrggbb — the dialog every categorical legend opens'});

export const readingAtLeast = Then('the {string} reading of {widget} should be at least {float}', async (page: Page, name: string, target: ElementRef, value: number) => {
  await expect.poll(() => v.readValue(page, target, name), {message: `"${name}" reading of ${target.phrase}`}).toBeGreaterThanOrEqual(value);
});

export const readingLower = Then('the {string} reading of {widget} should be lower than before', (page: Page, name: string, target: ElementRef) =>
  v.expectReading(page, target, name, 'lower'), {description: 'against the snapshot before the last change'});

export const readingHigher = Then('the {string} reading of {widget} should be higher than before', (page: Page, name: string, target: ElementRef) =>
  v.expectReading(page, target, name, 'higher'));

export const readingDiffers = Then('the {string} reading of {widget} should differ from before', (page: Page, name: string, target: ElementRef) =>
  v.expectReading(page, target, name, 'differ'), {description: 'not what the snapshot before the last change held — the 3D scatter plot\'s "scene signature" is what "repainted" means on a WebGL canvas'});

export const readingSame = Then('the {string} reading of {widget} should be the same as before', (page: Page, name: string, target: ElementRef) =>
  v.expectReading(page, target, name, 'same'), {description: 'read once the viewer is quiet, and read once'});

export const rememberReading = When('user remembers the {string} reading of {widget}', (page: Page, name: string, target: ElementRef) =>
  v.rememberReading(page, target, name), {tier: 'api', description: 'kept by viewer type, so "as remembered" holds across a close and a reopen (a layout or project round-trip)'});

export const readingAsRemembered = Then('the {string} reading of {widget} should be as remembered', (page: Page, name: string, target: ElementRef) =>
  v.expectRememberedReading(page, target, name));

export const readingNotAsRemembered = Then('the {string} reading of {widget} should not be as remembered', (page: Page, name: string, target: ElementRef) =>
  v.expectRememberedReading(page, target, name, true), {description: 'the change the step in between was supposed to make actually reached the reading'});

// --- the legend ------------------------------------------------------------------------------------

export const legendSide = Then('the legend of {widget} should be on the {word}', (page: Page, target: ElementRef, side: string) =>
  v.expectLegendSide(page, target, side), {description: 'left, right, top or bottom — the side the viewer laid its legend out on'});

export const legendLists = Then('the legend of {widget} should list {int} item(s)', (page: Page, target: ElementRef, count: number) =>
  v.expectLegendItems(page, target, count), {description: 'the item total the legend publishes (data-legend-items), every section, rendered or scrolled out — "should have N items" counts the rendered rows only'});

export const legendFewer = Then('the legend of {widget} should list fewer items than before', (page: Page, target: ElementRef) =>
  v.expectLegendItemsChange(page, target, 'fewer'), {description: 'against the snapshot before the last change'});

export const legendSameItems = Then('the legend of {widget} should list the same items as before', (page: Page, target: ElementRef) =>
  v.expectLegendItemsChange(page, target, 'same'), {description: 'the same total and the same rendered keys as at the snapshot before the last change'});

export const legendDocked = Then('the legend of {widget} should be docked', (page: Page, target: ElementRef) => v.expectLegendMode(page, target, 'docked'),
  {description: 'the mode the legend publishes: docked at a side, in a corner over the plot, collapsed to the mini icon, shown in the tooltip'});

export const legendCorner = Then('the legend of {widget} should be in a corner', (page: Page, target: ElementRef) => v.expectLegendMode(page, target, 'corner'));

export const legendMini = Then('the legend of {widget} should be collapsed to the mini icon', (page: Page, target: ElementRef) => v.expectLegendMode(page, target, 'mini icon'));

export const legendTooltip = Then('the legend of {widget} should be shown in the tooltip', (page: Page, target: ElementRef) => v.expectLegendMode(page, target, 'tooltip'));

export const legendSlot = Then('the legend of {widget} should be in the {string} slot', (page: Page, target: ElementRef, slot: string) =>
  v.expectLegendSlot(page, target, slot), {description: 'left, right, top, bottom, leftTop, leftBottom, rightTop, rightBottom'});

export const legendNotSlot = Then('the legend of {widget} should not be in the {string} slot', (page: Page, target: ElementRef, slot: string) =>
  v.expectLegendSlot(page, target, slot, true));

export const legendPlacedAsBefore = Then('the legend of {widget} should be placed as before', (page: Page, target: ElementRef) =>
  v.expectLegendPlacedAsBefore(page, target), {description: 'the same mode and slot as at the snapshot before the last change, read once the viewer is quiet'});

export const clickLegendItem = When('user clicks on {string} item in the legend of {widget}', (page: Page, label: string, target: ElementRef) =>
  v.clickLegendItem(page, target, label), {tier: 'ui', description: 'the item by its label; the category then filters the viewer — the baseline is taken before and the viewer settles after'});

export const clickLegendItemHolding = When('user clicks on {string} item in the legend of {widget} holding {key}', (page: Page, label: string, target: ElementRef, key: string) =>
  v.clickLegendItem(page, target, label, {key: normalizeKey(key)}), {tier: 'ui', description: 'Control adds the category to the legend selection'});

export const clickLegendCross = When('user clicks on the cross of {string} item in the legend of {widget}', (page: Page, label: string, target: ElementRef) =>
  v.clickLegendItem(page, target, label, {cross: true}), {tier: 'ui', description: 'the cross a selected item shows on hover — removes it from the selection'});

export const legendItemColor = Then('the {string} item in the legend of {widget} should be colored {string}', (page: Page, label: string, target: ElementRef, color: string) =>
  v.expectLegendItemColor(page, target, label, color), {description: 'the color the item is drawn in (a shade of anti-aliasing allowed), #rrggbb'});

export const legendItemNotColor = Then('the {string} item in the legend of {widget} should not be colored {string}', (page: Page, label: string, target: ElementRef, color: string) =>
  v.expectLegendItemColor(page, target, label, color, true));

export const legendItemsDiffer = Then('the {string} and {string} items in the legend of {widget} should be colored differently',
  (page: Page, a: string, b: string, target: ElementRef) => v.expectLegendItemsDiffer(page, target, a, b));

export const dragLegendSplitter = When('user drags the legend splitter of {widget} by {int} pixels', (page: Page, target: ElementRef, px: number) =>
  v.dragLegendSplitter(page, target, px), {tier: 'ui', description: 'along the axis the bar resizes, positive away from the plot; the viewer settles after'});

export const notRepainted = Then('{widget} should not have repainted', (page: Page, target: ElementRef) => v.expectNotRepainted(page, target),
  {description: 'no canvas change since the snapshot before the last gesture, read once the viewer is quiet (a render pass that draws the same picture — the mouse-over row — is not a repaint)'});

export const paintedInColors = Then('{widget} should be painted in at least {int} colors', (page: Page, target: ElementRef, count: number) =>
  v.expectPalette(page, target, count), {description: 'distinct colors covering 500 pixels or more, white aside'});

export const areaPainted = Then('the {string} area of {widget} should be painted', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaPainted(page, target, area), {description: 'painted pixels inside a hit area the viewer reports ("M values", "bar Asian")'});

export const hasArea = Then('{widget} should have a(n) {string} area', (page: Page, target: ElementRef, area: string) => v.expectHasArea(page, target, area),
  {description: 'the viewer reports the hit area right now — the statistics strip, a category, a comparison row'});

export const hasNoArea = Then('{widget} should not have a(n) {string} area', (page: Page, target: ElementRef, area: string) =>
  v.expectHasArea(page, target, area, true));

export const moreHighlight = Then('{widget} should show more selection highlight than before', (page: Page, target: ElementRef) =>
  v.expectHighlight(page, target, 'more'), {description: 'more pixels in the selected-rows color than the snapshot before the last change'});

export const lessHighlight = Then('{widget} should show less selection highlight than before', (page: Page, target: ElementRef) =>
  v.expectHighlight(page, target, 'less'));

export const someHighlight = Then('{widget} should show a selection highlight', (page: Page, target: ElementRef) => v.expectHighlight(page, target, 'some'),
  {description: 'pixels in the selected-rows color'});

export const noHighlight = Then('{widget} should show no selection highlight', (page: Page, target: ElementRef) => v.expectHighlight(page, target, 'none'),
  {description: 'not one pixel in the selected-rows color'});

// --- the value range ---------------------------------------------------------------------------------

export const narrowerRange = Then('{widget} should show a narrower value range than before', (page: Page, target: ElementRef) =>
  v.expectValueRange(page, target, 'narrower'), {description: 'the viewport the viewer reports, against the snapshot before the last change'});

export const sameRange = Then('{widget} should show the same value range as before', (page: Page, target: ElementRef) =>
  v.expectValueRange(page, target, 'same'));

export const widerRange = Then('{widget} should show a wider value range than before', (page: Page, target: ElementRef) =>
  v.expectValueRange(page, target, 'wider'));

export const rangeWithinColumn = Then('the value range of {widget} should lie within {string} column', (page: Page, target: ElementRef, column: string) =>
  v.expectValueRangeWithin(page, target, column), {description: 'no empty space beyond the column\'s min and max (a tenth of slack)'});

export const scaleNarrower = Then('the color scale of {widget} should cover a narrower range than before', (page: Page, target: ElementRef) =>
  v.expectScaleRange(page, target, 'narrower'), {description: 'the range the color scale labels (the filtered rows\' when they narrow it), against the snapshot before the last change'});

export const scaleWider = Then('the color scale of {widget} should cover a wider range than before', (page: Page, target: ElementRef) =>
  v.expectScaleRange(page, target, 'wider'));

export const scaleSame = Then('the color scale of {widget} should cover the same range as before', (page: Page, target: ElementRef) =>
  v.expectScaleRange(page, target, 'same'));

export const rememberRange = When('user remembers the value range of {widget}', (page: Page, target: ElementRef) => v.rememberRange(page, target),
  {tier: 'api', description: 'kept by viewer type, so "the remembered value range" holds across a close and a reopen (a project round-trip)'});

export const rememberedRange = Then('{widget} should show the remembered value range', (page: Page, target: ElementRef) => v.expectRememberedRange(page, target));

// --- events and errors -----------------------------------------------------------------------------

export const listenFor = Given('user listens for {string} event on {widget}', (page: Page, event: string, target: ElementRef) =>
  v.listenFor(page, target, event), {tier: 'api', description: 'a viewer event by name (d4-boxplot-reset-view), subscribed until "should have fired" reads it or the viewer closes'});

export const eventFired = Then('{string} event should have fired on {widget}', (page: Page, event: string, target: ElementRef) =>
  v.expectFired(page, target, event), {description: 'at least once since "listens for"; reading it ends the subscription'});

export const eventNotFired = Then('{string} event should not have fired on {widget}', (page: Page, event: string, target: ElementRef) =>
  v.expectNotFired(page, target, event), {description: 'not once since "listens for"; the subscription stays'});

export const noErrors = Then('no errors should have been logged', (page: Page) => {
  expect(takeErrors(page), 'console errors and page errors since the last check').toEqual([]);
}, {description: 'console errors and uncaught exceptions since the previous check, the scenario start or the login; checking clears them'});

export const noBalloons = Then('no error or warning balloon should have been shown', async (page: Page) => {
  const shown = (await v.takeBalloons(page)).filter((b) => b.type === 'error' || b.type === 'warning');
  expect(shown.map((b) => `${b.type}: ${b.message}`), 'error and warning balloons since the last check').toEqual([]);
}, {description: 'the platform\'s balloons (d4-balloon-shown) since the previous check, the scenario start or the login; checking clears them'});

/** The balloons of a type since the last read, polled: a balloon a command raises lands a task
 * after the gesture. */
async function expectBalloon(page: Page, type: string, text?: string): Promise<void> {
  let shown: string[] = [];
  await expect.poll(async () => {
    shown = shown.concat((await v.takeBalloons(page)).map((b) => `${b.type}: ${b.message}`));
    return shown.some((s) => s.startsWith(`${type}: `) && (text === undefined || s.includes(text)));
  }, {timeout: 5000, message: `${text === undefined ? `an ${type} balloon` : `an ${type} balloon containing "${text}"`}; balloons since the last check: ${shown.join(' | ') || 'none'}`}).toBe(true);
}

export const errorBalloon = Then('an error balloon should have been shown', (page: Page) => expectBalloon(page, 'error'),
  {description: 'since the previous balloon check; reading clears the balloons'});

export const warningBalloon = Then('a warning balloon should have been shown', (page: Page) => expectBalloon(page, 'warning'));

export const errorBalloonText = Then('an error balloon containing {string} should have been shown', (page: Page, text: string) => expectBalloon(page, 'error', text));

export const warningBalloonText = Then('a warning balloon containing {string} should have been shown', (page: Page, text: string) => expectBalloon(page, 'warning', text));

// --- tooltips --------------------------------------------------------------------------------------

export const tooltipColumns = Then('the tooltip should show columns {string}', (page: Page, list: string) => v.expectTooltipColumns(page, list),
  {description: 'the row tooltip lists exactly these columns (comma-separated, any order)'});

export const tooltipNotColumns = Then('the tooltip should not show columns {string}', (page: Page, list: string) => v.expectTooltipColumns(page, list, true),
  {description: 'the row tooltip lists a different set of columns'});

export const tooltipValue = Then('the tooltip should show {string} as {string}', (page: Page, column: string, value: string) => v.expectTooltipValue(page, column, value),
  {description: 'the row tooltip\'s value for a column, as text'});

export const tooltipSomeColumns = Then('the tooltip should show some columns', async (page: Page) => {
  await expect.poll(() => v.tooltipColumns(page), {timeout: 5000, message: 'the row tooltip lists no column'}).not.toEqual([]);
}, {description: 'the row tooltip lists at least one column — the table\'s default tooltip, whatever it holds'});
