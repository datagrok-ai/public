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
import {takeErrors} from '../../../src/runtime/harness.js';
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

export const propertiesShouldBe = Then('properties of {widget} should be:', async (page: Page, target: ElementRef, table: string[][]) => {
  for (const [caption, value] of table)
    await v.expectProperty(page, target, caption, value);
}, {description: '| property caption | value | rows, each read back through the property bag'});

export const saveLayout = When('user saves the layout of the current table view', (page: Page) => v.saveLayout(page),
  {tier: 'api', description: 'tv.saveLayout(), kept for "loads the saved layout" later in the feature'});

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
  await page.keyboard.down(normalizeKey(key));
  await page.mouse.click(c.x, c.y);
  await page.keyboard.up(normalizeKey(key));
}, {tier: 'ui', description: 'a click with a modifier held: Control adds to the selection, Shift extends it'});

export const dragSelectionOverArea = When('user drags a selection box over the {string} area of {widget}', async (page: Page, area: string, target: ElementRef) => {
  const b = await v.hitArea(page, target, area, true);
  await page.keyboard.down('Shift');
  await page.mouse.move(b.x + b.width * 0.1, b.y + b.height * 0.1);
  await page.mouse.down();
  await page.mouse.move(b.x + b.width * 0.9, b.y + b.height * 0.9, {steps: 3});
  await page.mouse.up();
  await page.keyboard.up('Shift');
}, {tier: 'ui', description: 'a Shift-drag across the inner 80% of the area — the platform\'s rectangle selection'});

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

export const areaLessInk = Then('the {string} area of {widget} should have less ink than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaInk(page, target, area, 'less'), {description: 'fewer painted pixels inside a hit area than the snapshot before the last change had there'});

export const areaMoreInk = Then('the {string} area of {widget} should have more ink than before', (page: Page, area: string, target: ElementRef) =>
  v.expectAreaInk(page, target, area, 'more'));

export const areaColor = Then('the {string} area of {widget} should contain the color {string}', (page: Page, area: string, target: ElementRef, color: string) =>
  v.expectAreaColor(page, target, area, color), {description: 'pixels of a #rrggbb color (a shade of anti-aliasing allowed) inside a hit area'});

export const areasDiffer = Then('the {string} and {string} areas of {widget} should be painted in different colors',
  (page: Page, a: string, b: string, target: ElementRef) => v.expectAreasDiffer(page, target, a, b),
  {description: 'one area has a color the other does not — per-category coloring, not chrome'});

export const boundTable = Then('{widget} should be bound to table {string}', (page: Page, target: ElementRef, name: string) => v.expectBoundTable(page, target, name),
  {description: 'the table the viewer draws (viewer.dataFrame), not the Table property it was asked for'});

export const showsRows = Then('{widget} should show {int} rows', (page: Page, target: ElementRef, count: number) =>
  v.expectReading(page, target, 'rows shown', 'equal', count), {description: 'the rows the viewer draws after its own filter and the table\'s (the "rows shown" reading)'});

export const showsFewerRows = Then('{widget} should show fewer rows than before', (page: Page, target: ElementRef) =>
  v.expectReading(page, target, 'rows shown', 'lower'), {description: 'against the snapshot before the last change'});

export const showsMoreRows = Then('{widget} should show more rows than before', (page: Page, target: ElementRef) =>
  v.expectReading(page, target, 'rows shown', 'higher'));

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

// --- tooltips --------------------------------------------------------------------------------------

export const tooltipColumns = Then('the tooltip should show columns {string}', (page: Page, list: string) => v.expectTooltipColumns(page, list),
  {description: 'the row tooltip lists exactly these columns (comma-separated, any order)'});

export const tooltipNotColumns = Then('the tooltip should not show columns {string}', (page: Page, list: string) => v.expectTooltipColumns(page, list, true),
  {description: 'the row tooltip lists a different set of columns'});

export const tooltipSomeColumns = Then('the tooltip should show some columns', async (page: Page) => {
  await expect.poll(() => v.tooltipColumns(page), {timeout: 5000, message: 'the row tooltip lists no column'}).not.toEqual([]);
}, {description: 'the row tooltip lists at least one column — the table\'s default tooltip, whatever it holds'});
