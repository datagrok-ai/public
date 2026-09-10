/* The steps only the tile viewer needs. Everything else its features use is the library's
   `viewers` tier and the platform's data steps, read from the hit areas and readings the viewer
   reports (`view`, `viewer menu`, `lane "<name>"`, `lane header "<name>"`, `lane content "<name>"`,
   `tile of row <r>`, `field "<COL>" of row <r>`, `label "<COL>" of row <r>`; `lanes`, `lane names`,
   `lanes list`, `single lane`, `tiles`, `tiles in lane "<name>"`, `lane of row <r>`, `fields`,
   `fields shown`, `<COL> of row <r>`, `current row`, `rows selected`, `selected rows shown`,
   `auto generate`, `form designed`, `scroll of lane "<name>"` — see
   `core/client/d4/lib/src/viewers/tile_viewer/CLAUDE.md`).

   What is left here is this viewer's own: the sketch form designer, whose hosts are a view portaled
   out of the viewer, and the lane drag. Its card-content claims, its field-composition contrast,
   its viewer menu, its list readings and its description placement were every card viewer's and
   are in the library now (`bindings/tiers/viewers/widgets.ts`) — the forms viewer joins them by
   reporting the same `fields` and `<COL> of row <r>` readings. */
import {expect, Locator, Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {ElementRef, el, viewers} from '@datagrok-libraries/bdd/runtime';

/* The Dart name of the viewer is `viewer-Tile-Viewer`, which the generic `viewer` kind would only
   reach through the phrase "tile viewer viewer"; the designer is a view of its own, portaled out
   of the viewer (its EDIT / RESET / CLOSE AND APPLY buttons live in the ribbon, not inside it). */
const PANEL = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]';
element('tile viewer', {selector: '[name="viewer-Tile-Viewer"]', parts: {
  title: `${PANEL}//*[contains(@class, "panel-titlebar-text")]`,
  'settings icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="icon-font-icon-settings"]`,
  'menu icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="icon-font-icon-menu"]`,
  'close icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="Close" or @name="icon-font-icon-close"]`,
  description: '.d4-viewer-description',
}});
element('form designer', {selector: '.grok-view-sketch'});

const TILE_VIEWER = 'Tile Viewer';

// --- the sketch form designer -----------------------------------------------------------------------

/* The designer names a host after the column with every non-alphanumeric run turned into a dash
   (`DIS_POP` is `div-DIS-POP`), and holds two hosts per column: the value (an `input-<slug>`) and
   its caption (a `.d4-sketch-column-name`). The two sets are read apart because RESET restores
   both and a summed host count would hide one channel changing while the other did not. */
const slug = (name: string): string => name.replace(/[^A-Za-z0-9]+/g, '-');

function designerHost(page: Page, column: string, kind: 'value' | 'label'): Locator {
  const host = page.locator(`.grok-view-sketch .d4-host[name="div-${slug(column)}"]`);
  const value = page.locator(`input[name="input-${slug(column)}"]`);
  return kind === 'value'
    ? host.filter({has: value}).first()
    : host.filter({has: page.locator('input.d4-sketch-column-name')}).filter({hasNot: value}).first();
}

async function designerFields(page: Page): Promise<{values: string[]; labels: string[]}> {
  return page.evaluate(() => {
    const view = document.querySelector('.grok-view-sketch');
    if (view === null)
      return {values: [], labels: []};
    const hosts = Array.from(view.querySelectorAll('.d4-host[name^="div-"]'));
    const nameOf = (h: Element): string => (h.getAttribute('name') ?? '').replace('div-', '');
    const values = hosts.filter((h) => h.querySelector('input[name^="input-"]') !== null);
    const labels = hosts.filter((h) => !values.includes(h) && h.querySelector('input.d4-sketch-column-name') !== null);
    return {values: values.map(nameOf).sort(), labels: labels.map(nameOf).sort()};
  });
}

async function expectDesignerFields(page: Page, kind: 'values' | 'labels', expected: string): Promise<void> {
  const want = (expected === '' ? [] : expected.split(/\s*,\s*/)).map(slug).sort();
  let now: {values: string[]; labels: string[]} = {values: [], labels: []};
  await expect.poll(async () => {
    now = await designerFields(page);
    return now[kind].join(', ');
  }, {timeout: 5000, message: `the ${kind === 'values' ? 'value' : 'label'} fields of the form designer ` +
    `(its value fields are ${now.values.join(', ') || 'none'}, its label fields ${now.labels.join(', ') || 'none'})`})
    .toBe(want.join(', '));
}

export const designerValueFields = Then('the form designer should show value fields {string}',
  (page: Page, expected: string) => expectDesignerFields(page, 'values', expected),
{description: 'the whole set of value hosts, by column name and regardless of order'});

export const designerLabelFields = Then('the form designer should show label fields {string}',
  (page: Page, expected: string) => expectDesignerFields(page, 'labels', expected),
{description: 'the caption hosts, read apart from the value hosts — a summed count would hide one channel changing'});

async function deleteHost(page: Page, column: string, kind: 'value' | 'label'): Promise<void> {
  const host = designerHost(page, column, kind);
  await host.waitFor({state: 'visible', timeout: 5000}).catch(() => {
    throw new Error(`the form designer shows no "${column}" ${kind} field`);
  });
  await host.click();
  await page.keyboard.press('Delete');
  await host.waitFor({state: 'detached', timeout: 5000}).catch(() => undefined);
}

export const deleteValueField = When('user deletes the {string} value field in the form designer',
  (page: Page, column: string) => deleteHost(page, column, 'value'),
{tier: 'ui', description: 'a real click then Delete — a synthetic event does not select a sketch host'});

export const deleteLabelField = When('user deletes the {string} label field in the form designer',
  (page: Page, column: string) => deleteHost(page, column, 'label'),
{tier: 'ui', description: 'the caption host of that column, the one that carries no value input'});

// --- the drag between lanes ---------------------------------------------------------------------

/* The library's `drags the {string} area … to the {string} area` moves down, three steps across
   and up. The platform's drop zone is a body-level overlay the FIRST move past five pixels
   creates (`ui.makeDroppable`), and the browser coalesces moves queued while the main thread is
   busy — so a run where the crossing move and the release collapse into one lands the release
   before the zone exists and the drop is silently lost. This gesture separates the three moves the
   drag needs: one to start it, the crossing, and one inside the target while the zone is there. */
export const dragCardIntoLane = When('user drags the card of row {int} of {widget} into lane {string}',
  async (page: Page, row: number, target: ElementRef, lane: string) => {
    const from = viewers.centerOf(await viewers.hitArea(page, target, `tile of row ${row}`, true));
    const to = viewers.centerOf(await viewers.hitArea(page, target, `lane content "${lane}"`));
    await page.mouse.move(from.x, from.y);
    await page.mouse.down();
    await page.mouse.move(from.x + 12, from.y + 12);
    await page.mouse.move(to.x, to.y, {steps: 6});
    await page.mouse.move(to.x + 1, to.y + 1);
    await page.mouse.up();
    await viewers.settle(page, target, 1000);
  }, {tier: 'ui', description: 'a card picked up and dropped on another lane, the way the Kanban move is made'});

// --- adding the viewer ------------------------------------------------------------------------------

/* `user adds a {viewer} viewer` would pass "tile" to the platform, whose type is "Tile Viewer". */
export const addTileViewer = Given('user adds a tile viewer', (page: Page) => viewers.addViewer(page, TILE_VIEWER),
  {tier: 'api', description: 'grok.shell.tv.addViewer("Tile Viewer")'});

export const addTileViewerWith = Given('user adds a tile viewer with:', async (page: Page, table: string[][]) => {
  await viewers.addViewer(page, TILE_VIEWER);
  await viewers.setProperties(page, el('tile viewer'), table.map(([caption, value]) => [caption, value] as [string, string]));
}, {tier: 'api', description: '| property caption | value | rows applied right after adding'});
