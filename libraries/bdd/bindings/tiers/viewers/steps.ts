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
  const c = v.centerOf(await v.hitArea(page, target, area));
  await page.mouse.move(c.x - 3, c.y - 3);
  await page.mouse.move(c.x, c.y);
}, {tier: 'ui'});

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
  {description: 'the canvas differs from the snapshot taken before the last change; the new state becomes the snapshot'});

export const lessInk = Then('{widget} should have less ink than before', (page: Page, target: ElementRef) => v.expectInk(page, target, 'less'),
  {description: 'fewer painted pixels than the snapshot before the last change'});

export const moreInk = Then('{widget} should have more ink than before', (page: Page, target: ElementRef) => v.expectInk(page, target, 'more'));

export const painted = Then('{widget} should be painted', (page: Page, target: ElementRef) => v.expectInk(page, target, 'some'),
  {description: 'the canvas has painted pixels'});

// --- events and errors -----------------------------------------------------------------------------

export const listenFor = Given('user listens for {string} event on {widget}', (page: Page, event: string, target: ElementRef) =>
  v.listenFor(page, target, event), {tier: 'api', description: 'a viewer event by name (d4-boxplot-reset-view), subscribed until "should have fired" reads it or the viewer closes'});

export const eventFired = Then('{string} event should have fired on {widget}', (page: Page, event: string, target: ElementRef) =>
  v.expectFired(page, target, event), {description: 'at least once since "listens for"; reading it ends the subscription'});

export const noErrors = Then('no errors should have been logged', (page: Page) => {
  expect(takeErrors(page), 'console errors and page errors since the last check').toEqual([]);
}, {description: 'console errors and uncaught exceptions since the previous check (or the scenario start); checking clears them'});

// --- tooltips --------------------------------------------------------------------------------------

export const tooltipColumns = Then('the tooltip should show columns {string}', (page: Page, list: string) => v.expectTooltipColumns(page, list),
  {description: 'the row tooltip lists exactly these columns (comma-separated, any order)'});

export const tooltipNotColumns = Then('the tooltip should not show columns {string}', (page: Page, list: string) => v.expectTooltipColumns(page, list, true),
  {description: 'the row tooltip lists a different set of columns'});
