/* The steps only the tile viewer needs. Everything else its features use is the library's
   `viewers` tier and the platform's data steps, read from the hit areas and readings the viewer
   reports (`view`, `viewer menu`, `lane "<name>"`, `lane header "<name>"`, `lane content "<name>"`,
   `tile of row <r>`, `field "<COL>" of row <r>`, `label "<COL>" of row <r>`; `lanes`, `lane names`,
   `lanes list`, `single lane`, `tiles`, `tiles in lane "<name>"`, `lane of row <r>`, `fields`,
   `fields shown`, `<COL> of row <r>`, `current row`, `rows selected`, `selected rows shown`,
   `auto generate`, `form designed`, `scroll of lane "<name>"` — see
   `core/client/d4/lib/src/viewers/tile_viewer/CLAUDE.md`).

   Four of these are not tile viewer business and should be promoted to the library:
   `the {string} reading of {widget} should (not )contain {string}` (a reading that is a list —
   the tile viewer's `fields` and `lane names`, the pc plot's axes), `user opens the viewer menu of
   {widget}` and `user picks {string} from the viewer menu of {widget}` (they honour the
   `ContextMenu` shortcut a widget declares in `getWidgetStatus`, which is the generic answer to a
   viewer whose centre is covered by something with a menu of its own), and
   `the description of {widget} should be above/below its content` (any viewer with a Description
   Position). The rest — the per-tile column checks, the auto-vs-designed refill contrast and the
   sketch form designer — are this viewer's own. */
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

/** Every reading the viewer reports right now, keyed as it names them. */
async function readings(page: Page, target: ElementRef): Promise<Record<string, unknown>> {
  await viewers.installViewerRuntime(page);
  const loc = await viewers.viewerLocator(page, target);
  return loc.evaluate((e) => (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {});
}

/** The rows the lanes have laid out and what their cards show for one column: the `<COL> of row
 * <r>` readings, which is the grid's display string for that cell. */
async function shownColumn(page: Page, target: ElementRef, column: string): Promise<{row: number; text: string}[]> {
  const values = await readings(page, target);
  const out: {row: number; text: string}[] = [];
  for (const key of Object.keys(values)) {
    const m = new RegExp(`^${column.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')} of row (\\d+)$`).exec(key);
    if (m)
      out.push({row: Number(m[1]), text: String(values[key])});
  }
  return out;
}

async function expectEveryTile(page: Page, target: ElementRef, column: string, holds: (text: string) => boolean, what: string): Promise<void> {
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
  }, {timeout: 5000, message: seen === 0
    ? `${target.phrase} has laid out no card with a "${column}" field`
    : `not every card of ${target.phrase} shows ${what} in "${column}": ${bad}`}).toBe(true);
}

export const everyTileShows = Then('every tile of {widget} should show {string} in {string}',
  (page: Page, target: ElementRef, value: string, column: string) =>
    expectEveryTile(page, target, column, (t) => t === value, `"${value}"`),
{description: 'every card the lanes have laid out shows that text for the column — the "all tiles are M" claim'});

export const everyTileBetween = Then('every tile of {widget} should show a value between {float} and {float} in {string}',
  (page: Page, target: ElementRef, lo: number, hi: number, column: string) =>
    expectEveryTile(page, target, column, (t) => Number(t) >= lo && Number(t) <= hi, `a value between ${lo} and ${hi}`),
{description: 'the "every age is over 50" claim, read off the cards and not off the table'});

// --- readings that are lists ----------------------------------------------------------------------

async function expectReadingContains(page: Page, target: ElementRef, name: string, item: string, negate: boolean): Promise<void> {
  const has = async (): Promise<boolean> => {
    const text = String(await viewers.readValue(page, target, name));
    return text.split(/\s*,\s*/).includes(item);
  };
  let last = '';
  const poll = expect.poll(async () => {
    last = String(await viewers.readValue(page, target, name));
    return await has();
  }, {timeout: 5000, message: `"${name}" of ${target.phrase} is "${last}"`});
  await (negate ? poll.not : poll).toBe(true);
}

export const readingContains = Then('the {string} reading of {widget} should contain {string}',
  (page: Page, name: string, target: ElementRef, item: string) => expectReadingContains(page, target, name, item, false),
{description: 'a reading the viewer reports as a comma-separated list ("fields", "lane names", "lanes list") holds that member — should be promoted to the library'});

export const readingNotContains = Then('the {string} reading of {widget} should not contain {string}',
  (page: Page, name: string, target: ElementRef, item: string) => expectReadingContains(page, target, name, item, true),
{description: 'the negative of the above; a feature pairs it with a positive one that proves the member was there to lose'});

// --- the viewer menu ------------------------------------------------------------------------------

/** The region the widget declares as its context-menu target (`getWidgetStatus().shortcuts`
 * `ContextMenu` → a hit-area name), or nothing when it declares none. The tile viewer needs it:
 * a right-click on a card opens the *column* menu of the field under the pointer, so the viewer's
 * own menu is only reachable where no card is — a lane header, or the free strip under a lane. */
async function menuArea(page: Page, target: ElementRef): Promise<string | undefined> {
  await viewers.installViewerRuntime(page);
  const loc = await viewers.viewerLocator(page, target);
  const name: string | null = await loc.evaluate((e) =>
    (window as any).__bdd.viewerOf(e).getWidgetStatus()?.shortcuts?.['ContextMenu'] ?? null);
  return name ?? undefined;
}

export const openViewerMenu = When('user opens the viewer menu of {widget}', async (page: Page, target: ElementRef) => {
  await viewers.openContextMenuOf(page, target, await menuArea(page, target));
}, {tier: 'ui', description: 'right-clicks the region the widget declares as its ContextMenu shortcut — should be promoted to the library'});

export const pickFromViewerMenu = When('user picks {string} from the viewer menu of {widget}', async (page: Page, path: string, target: ElementRef) => {
  await viewers.snapshot(page, target);
  await viewers.openContextMenuOf(page, target, await menuArea(page, target));
  await viewers.pickMenuPath(page, path);
}, {tier: 'ui', description: 'the same, then the path picked — should be promoted to the library'});

// --- the card composition -------------------------------------------------------------------------

let rememberedFields: string[] | null = null;

async function fieldsOf(page: Page, target: ElementRef): Promise<string[]> {
  const text = String(await viewers.readValue(page, target, 'fields'));
  return text === '' ? [] : text.split(/\s*,\s*/);
}

export const rememberFields = When('user remembers the fields of {widget}', async (page: Page, target: ElementRef) => {
  rememberedFields = await fieldsOf(page, target);
  if (rememberedFields.length === 0)
    throw new Error(`${target.phrase} shows no fields to remember`);
}, {tier: 'api', description: 'the card\'s composition, for the auto-versus-designed refill contrast that follows'});

/** Both halves of the contrast compare the composition with the remembered one; neither names the
 * column that left or the one that took its place, because the form orders its fields by relevance
 * and the excluded one is whatever the score left over. */
async function expectComposition(page: Page, target: ElementRef, refilled: boolean): Promise<void> {
  if (rememberedFields === null)
    throw new Error('nothing was remembered: put "user remembers the fields of tile viewer" before the change');
  const before = rememberedFields;
  let now: string[] = [];
  const holds = async (): Promise<boolean> => {
    now = await fieldsOf(page, target);
    const gone = before.filter((f) => !now.includes(f));
    const gained = now.filter((f) => !before.includes(f));
    return refilled
      ? now.length === before.length && gone.length === 1 && gained.length === 1
      : now.length === before.length && gone.length === 0 && gained.length === 0;
  };
  await expect.poll(holds, {timeout: 8000, message: refilled
    ? `${target.phrase} did not refill the freed slot: it showed ${before.join(', ')} and now shows ${now.join(', ')}`
    : `${target.phrase} did not keep its composition: it showed ${before.join(', ')} and now shows ${now.join(', ')}`}).toBe(true);
}

export const fieldsRefilled = Then('the fields of {widget} should have refilled the freed slot',
  (page: Page, target: ElementRef) => expectComposition(page, target, true),
{description: 'the auto-generated card: the column that left took a field with it and a column that had none took its place, so the count is unchanged'});

export const fieldsAsRemembered = Then('the fields of {widget} should be as remembered',
  (page: Page, target: ElementRef) => expectComposition(page, target, false),
{description: 'the designed card: the same fields as before, so no column that had none gains one — the field of a column that left stays, empty'});

// --- the description's place ------------------------------------------------------------------------

async function expectDescriptionPlace(page: Page, target: ElementRef, above: boolean): Promise<void> {
  const loc = await viewers.viewerLocator(page, target);
  const view = await viewers.hitArea(page, target, 'view');
  const box = await loc.locator('.d4-viewer-description').filter({visible: true}).first().boundingBox();
  if (box === null)
    throw new Error(`${target.phrase} shows no description`);
  const what = above ? 'above' : 'below';
  const ok = above ? box.y + box.height <= view.y + 1 : box.y + 1 >= view.y + view.height;
  expect(ok, `the description of ${target.phrase} is not ${what} its content: the description spans ` +
    `${Math.round(box.y)}..${Math.round(box.y + box.height)} and the tiles ${Math.round(view.y)}..${Math.round(view.y + view.height)}`).toBe(true);
}

export const descriptionAbove = Then('the description of {widget} should be above its content',
  (page: Page, target: ElementRef) => expectDescriptionPlace(page, target, true),
{description: 'Description Position Top: the description box ends where the content begins — should be promoted to the library'});

export const descriptionBelow = Then('the description of {widget} should be below its content',
  (page: Page, target: ElementRef) => expectDescriptionPlace(page, target, false),
{description: 'Description Position Bottom: the description box starts where the content ends'});

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
    const loc = await viewers.viewerLocator(page, target);
    await loc.evaluate((e) => (window as any).__bdd.settle(e, 1000));
  }, {tier: 'ui', description: 'a card picked up and dropped on another lane, the way the Kanban move is made'});

// --- adding the viewer ------------------------------------------------------------------------------

/* `user adds a {viewer} viewer` would pass "tile" to the platform, whose type is "Tile Viewer". */
export const addTileViewer = Given('user adds a tile viewer', (page: Page) => viewers.addViewer(page, TILE_VIEWER),
  {tier: 'api', description: 'grok.shell.tv.addViewer("Tile Viewer")'});

export const addTileViewerWith = Given('user adds a tile viewer with:', async (page: Page, table: string[][]) => {
  await viewers.addViewer(page, TILE_VIEWER);
  await viewers.setProperties(page, el('tile viewer'), table.map(([caption, value]) => [caption, value] as [string, string]));
}, {tier: 'api', description: '| property caption | value | rows applied right after adding'});
