/* The catalog of the Forge view: its reloads (the models pane publishes aria-busy and data-forge-catalog-load,
   ARCHITECTURE.md "Automation surface"), the gestures on a model's row, what the grid lists and shows, the top menu's
   groups, and the Share dialog's second account. */
import type {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, type ElementRef, expect, gestures, locate, menuNames, openTopMenu, pollMs, viewers}
  from '@datagrok-libraries/bdd/runtime';
import {catalogRowOf, drawnTexts, modelArea} from './catalog-grid.js';
import {listOf, revokeSharingGrant, sharingAccount} from './grok-s.js';
import {slidePointer} from './pointer.js';

declare const DG: any;
declare const grok: any;

interface CatalogLoad {load: number; busy: string | null; isShown: boolean}

/** The models pane of every open Forge view, in the order of the shell's views: the load it shows, aria-busy, whether
 * it is shown. A hidden view's root is detached from the page, so the panes are read through the views. */
function catalogLoads(page: Page): Promise<CatalogLoad[]> {
  return page.evaluate(() => Array.from(grok.shell.views as any[])
    .map((v) => v.root?.querySelector('.forge-catalog') as Element | null)
    .filter((e): e is Element => e !== null && e !== undefined)
    .map((e) => ({load: Number(e.getAttribute('data-forge-catalog-load')), busy: e.getAttribute('aria-busy'),
      isShown: e.getClientRects().length > 0})));
}

const rememberedLoads = new WeakMap<Page, CatalogLoad[]>();

export const rememberCatalogLoads = When('user remembers the catalog loads', async (page: Page) => {
  let loads: CatalogLoad[] = [];
  await expect.poll(async () => {
    loads = await catalogLoads(page);
    return loads.length > 0 && loads.every((l) => l.busy !== 'true');
  }, {message: 'a Forge catalog on the page, none of them reloading'}).toBe(true);
  rememberedLoads.set(page, loads);
}, {tier: 'api', description: 'once no catalog is reloading (aria-busy not true), the load every catalog of the page shows'});

async function expectReloaded(page: Page, isShownOnly: boolean): Promise<void> {
  const before = rememberedLoads.get(page);
  if (before === undefined)
    throw new Error('the catalog loads were not remembered');
  await expect.poll(async () => {
    const now = await catalogLoads(page);
    if (now.length !== before.length)
      return `${now.length} catalogs, ${before.length} remembered`;
    if (isShownOnly && !now.some((l) => l.isShown))
      return 'no catalog is shown';
    const late = now.map((l, i) => ({l, i})).filter(({l, i}) => (!isShownOnly || l.isShown) &&
      (l.load <= before[i].load || l.busy === 'true'));
    return late.length === 0 ? 'reloaded' :
      late.map(({l, i}) => `catalog ${i + 1}: load ${l.load} (remembered ${before[i].load}), busy ${l.busy}`).join('; ');
  }, {message: 'the catalogs\' loads against the remembered ones'}).toBe('reloaded');
}

export const everyCatalogReloaded = Then('every catalog should have reloaded', (page: Page) => expectReloaded(page, false),
  {description: 'the models pane of every open Forge view (a hidden one too) shows a later load than remembered and is ' +
  'not reloading: the reload a change of the models starts in each open catalog has ended'});

export const shownCatalogReloaded = Then('the shown catalog should have reloaded', (page: Page) => expectReloaded(page, true),
  {description: 'the shown models pane shows a later load than remembered and is not reloading (Refresh reloads one)'});

export const firstCatalogShown = Then('the first Forge view should be shown', async (page: Page) => {
  await expect.poll(async () => (await catalogLoads(page)).map((l) => l.isShown ? 'shown' : 'hidden').join(', '),
    {message: 'the open Forge views, in the order of the shell\'s views'}).toMatch(/^shown(, hidden)*$/);
}, {description: 'of the open Forge views, in the order the shell lists them, the first one is shown and no other'});

/** Clicks the row header of the row that shows [name] holding Control: the grid adds the row to its selection or
 * takes it out. */
async function toggleRow(page: Page, target: ElementRef, name: string): Promise<void> {
  const c = await modelArea(page, target, name, (row) => `row header ${row}`);
  await gestures.withKeys(page, gestures.keysOf('Control'), () => page.mouse.click(c.x, c.y));
}

export const selectModels = When('user selects the Forge models {string} in {widget}', async (page: Page, list: string,
  target: ElementRef) => {
  for (const name of listOf(list))
    await toggleRow(page, target, name);
}, {tier: 'ui', description: 'a Control-click on the row header of each model\'s row, in the order given'});

export const toggleModel = When('user toggles the selection of the Forge model {string} in {widget}', (page: Page,
  name: string, target: ElementRef) => toggleRow(page, target, name),
{tier: 'ui', description: 'a Control-click on the row header of the model\'s row'});

export const rightClickModel = When('user right-clicks on the Forge model {string} in {widget}', async (page: Page,
  name: string, target: ElementRef) => {
  const row = await catalogRowOf(page, target, name);
  await viewers.openContextMenuOf(page, target, `cell ${row} of name`);
}, {tier: 'ui', description: 'opens the context menu on the name cell of the model\'s row and leaves it open'});

export const hoverModel = When('user hovers over the Forge model {string} in {widget}', async (page: Page, name: string,
  target: ElementRef) => {
  const row = await catalogRowOf(page, target, name);
  const box = await viewers.hitArea(page, target, `cell ${row} of name`, true);
  const c = viewers.centerOf(box);
  // A pointer that rests in the cell (the click that chose the row) asks for no tooltip: it enters from the left.
  await slidePointer(page, c, {x: box.x - 8, y: c.y});
  await viewers.settle(page, target);
}, {tier: 'ui', description: 'the pointer enters the name cell of the model\'s row from its left in five moves, the grid ' +
  'settled'});

export const modelCellReads = Then('the {string} cell of the Forge model {string} in {widget} should read {string}',
  async (page: Page, column: string, name: string, target: ElementRef, text: string) => {
    await expect.poll(() => viewers.onViewer(page, target, (el, [n, c]) => {
      const values: Record<string, unknown> = DG.Widget.find(el).getWidgetStatus().values ?? {};
      const key = Object.keys(values).find((k) => /^text of cell \d+ of name$/.test(k) && values[k] === n);
      if (key === undefined)
        return 'no drawn row of the model';
      const cell = `${key.replace(/ of name$/, '')} of ${c}`;
      return cell in values ? String(values[cell]) : `(no drawn ${cell})`;
    }, [name, column] as [string, string]), {message: `the "${column}" cell of "${name}" in ${target.phrase}`}).toBe(text);
  }, {description: 'the drawn text of the cell of the frame\'s column (tags, engine_name, ...) in the model\'s row'});

/* The grid reports only the rows it draws: a negative read of them passes for a row scrolled out of sight. */
export const catalogListsNot = Then('{widget} should not list the Forge model {string}', async (page: Page,
  target: ElementRef, name: string) => {
  await viewers.settle(page, target);
  const rows = await viewers.onViewer(page, target, (el, n) => {
    const df = DG.Widget.find(el).dataFrame;
    const names = df.col('name');
    let count = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (df.filter.get(i) && names.get(i) === n)
        count++;
    }
    return count;
  }, name);
  expect(rows, `the rows of ${target.phrase} that show "${name}"`).toBe(0);
}, {description: 'once the grid has settled, no row of its table that passes its filter is the model: read once, after ' +
  'the end signal of the change (a reload, a filter)'});

/* LIB-PROMOTE: the library scrolls the current table view's grid only; a grid in a narrow pane (the model's Activity)
   draws, and so reports, only the columns in sight. */
export const scrollGridToColumn = When('user scrolls {widget} to the {string} column', async (page: Page,
  target: ElementRef, column: string) => {
  await viewers.onViewer(page, target, (el, c) => {
    const grid = DG.Widget.find(el);
    if (grid.dataFrame.col(c) === null)
      throw new Error(`the grid has no "${c}" column`);
    grid.scrollToCell(c, 0);
  }, column);
  await viewers.settle(page, target);
}, {tier: 'api', description: 'grid.scrollToCell on the first row of the column, then the grid settles'});

/* LIB-PROMOTE: the captions a grid draws over its columns, in order: what a reader sees as headers. */
export const gridShowsColumns = Then('{widget} should show the columns {string}', async (page: Page, target: ElementRef,
  list: string) => {
  await expect.poll(async () => {
    const drawn = await drawnTexts(page, target);
    return drawn.isFriendlyNameShown ? 'showFriendlyName is on' :
      drawn.headers.map((h) => h ?? '(no header drawn)').join(', ');
  }, {message: `the column captions ${target.phrase} draws`}).toBe(list);
}, {description: 'the visible columns of the grid in order (no row header), each by the text of the header cell the ' +
  'grid finds at its header\'s area; showFriendlyName off'});

/* LIB-PROMOTE: a word that must appear nowhere a grid draws text: its drawn cells and its headers. Rows below the
   drawn ones and tooltips are not read. */
export const gridShowsNoText = Then('no drawn cell or header of {widget} should contain {string}', async (page: Page,
  target: ElementRef, text: string) => {
  await viewers.settle(page, target);
  const drawn = await drawnTexts(page, target);
  const texts = [...drawn.cells, ...drawn.headers.filter((h): h is string => h !== null)];
  expect({isFriendlyNameShown: drawn.isFriendlyNameShown, hasCells: drawn.cells.length > 0},
    `${target.phrase}: showFriendlyName and drawn cells`).toEqual({isFriendlyNameShown: false, hasCells: true});
  expect(texts.filter((s) => s.toLowerCase().includes(text.toLowerCase())),
    `the drawn cells and headers of ${target.phrase} holding "${text}" (ignoring case)`).toEqual([]);
}, {description: 'once the grid has settled, no drawn cell text and no drawn header caption contains the text, ' +
  'ignoring case; showFriendlyName off'});

/* LIB-PROMOTE: a selection the steps made and took back, read from the grid's table. */
export const gridSelectedRows = Then('{widget} should have {int} selected row(s)', async (page: Page, target: ElementRef,
  n: number) => {
  await expect.poll(() => viewers.onViewer(page, target, (el) => DG.Widget.find(el).dataFrame.selection.trueCount as number),
    {message: `the selected rows of ${target.phrase}`}).toBe(n);
}, {description: 'the selection of the grid\'s table'});

/* LIB-PROMOTE: "should be current" names the view in front; a command that may open a view again is claimed by the
   count of the open views of that name. */
export const openViewsNamed = Then('there should be {int} open view(s) named {string}', async (page: Page, n: number,
  name: string) => {
  await expect.poll(() => page.evaluate((v) => Array.from(grok.shell.views as any[]).filter((x) => x.name === v).length,
    name), {message: `the open views named "${name}"`}).toBe(n);
}, {description: 'grok.shell.views with that name, hidden ones included'});

/* LIB-PROMOTE: the library's top-menu claim finds the items it names; "only these" needs the group's items read. */
export const topMenuGroupExactly = Then('the {string} top menu group should list exactly {string}', async (page: Page,
  group: string, list: string) => {
  await openTopMenu(page, group, false);
  const {names} = menuNames(group);
  // A vertical group holds its items in a fixed wrapper; a group of the bar holds them directly.
  const parent = `[name="${names[names.length - 1]}"]`;
  const items = '.d4-menu-item-container > .d4-menu-item > .d4-menu-item-label';
  const labels = page.locator(`${parent} > .d4-menu-item-container-fixed > ${items}, ${parent} > ${items}`)
    .filter({visible: true});
  try {
    await expect.poll(async () => (await labels.allTextContents()).map((s) => s.trim()).filter((s) => s !== '').join(', '),
      {message: `the items of the ${group} menu`}).toBe(list);
  }
  finally {
    const size = page.viewportSize() ?? {width: 1280, height: 800};
    await page.mouse.move(size.width / 2, size.height / 2);
    await expect(page.locator('[role="menubar"] .d4-vert-menu').filter({visible: true}), 'the top menu closed').toHaveCount(0);
  }
}, {tier: 'ui', description: 'the group opened by the pointer: its visible items, in order, are exactly the list; the menu ' +
  'closes after'});

/* LIB-PROMOTE: the library's sharing steps find the second account by its login with the punctuation stripped
   ("forgetester"); the platform shows the account by its friendly name ("Forge Tester"). The model shown in the
   context panel is the one shared: its grant is taken back first at feature end, should the scenario fail before it
   removes the grant itself; a deleted model's entity keeps its grants but no longer reaches the page, so after the
   delete a grant left fails the cleanup instead. */
export const pickSharingUserByName = When('user picks the sharing user by name in {element}', async (page: Page,
  target: ElementRef) => {
  const model = await page.evaluate(() => {
    const o: any = grok.shell.o;
    return o instanceof DG.DomainRow && o.typeName === 'forge.model' ? {id: String(o.id), name: String(o.values.name)} :
      null;
  });
  if (model === null)
    throw new Error('the context panel shows no Forge model to share');
  atFeatureEnd(page, () => revokeSharingGrant(page, model.id, model.name), true);
  const {login, name} = await sharingAccount();
  const editor = await gestures.editorOf(page, target);
  await editor.click();
  // The typeahead searches by the start of a word: the login's first word finds the account.
  await editor.pressSequentially(login.split(/[^a-z0-9]/i)[0]);
  const row = page.locator('.d4-user-selector-drop-down tr, .d4-tags-selector-drop-down tr').filter({hasText: name}).first();
  await expect(row, `the "${name}" row of the user typeahead`).toBeVisible({timeout: pollMs(15000)});
  await row.click();
  await expect((await locate(page, target)).first(), `${target.phrase} after the pick`).toContainText(name);
}, {tier: 'ui', description: 'types the first word of DATAGROK_SHARING_LOGIN, takes the typeahead row that shows the ' +
  'account\'s name and sees its chip'});

export const removeSharingUserByName = When('user removes the sharing user by name from {element}', async (page: Page,
  target: ElementRef) => {
  const {name} = await sharingAccount();
  const row = (await locate(page, target)).locator('[name^="div-permissions-row-"]').filter({hasText: name}).first();
  await expect(row, `the grant row of "${name}"`).toBeVisible({timeout: pollMs(15000)});
  await row.hover();
  await row.locator('[name="button-Remove"]').click();
  await expect(row, `the grant row of "${name}" after Remove`).toBeHidden();
}, {tier: 'ui', description: 'hovers the grant row of DATAGROK_SHARING_LOGIN\'s shown name and clicks its Remove button'});

export const namesSharingUser = Then('{element} should name the sharing user', async (page: Page, target: ElementRef) => {
  const {name} = await sharingAccount();
  await expect((await locate(page, target)).first(), `${target.phrase} and "${name}"`).toContainText(name);
}, {description: 'the element\'s text holds the shown name of DATAGROK_SHARING_LOGIN (its friendly name)'});
