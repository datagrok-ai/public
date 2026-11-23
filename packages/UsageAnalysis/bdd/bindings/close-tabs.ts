/* Closing views and tables safely (core#456, core#531): the fixtures the close features set up through
   the API, and the close gestures no platform kind names: the view selector the simple mode closes
   views from, the table's menu on its name in the status bar, and an open project's menu in the
   Dashboards pane. A table link and a second view of a table are fixtures: the Link Tables dialog and
   the view menu are other features' subjects. */
import {Locator, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, escapeRegExp, expect} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

export const linkTables = Given('tables {string} and {string} are linked on {string}', (page: Page, t1: string, t2: string, key: string) =>
  page.evaluate(([a, b, k]) => grok.data.linkTables(grok.shell.table(a), grok.shell.table(b), [k], [k],
    [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER]), [t1, t2, key]),
{tier: 'api', description: 'the current row of the first table filters the second on the key column both tables have'});

export const openSecondView = Given('user opens a second view of table {string} named {string}', async (page: Page, table: string, name: string) => {
  await page.evaluate(([t, n]) => { const v = DG.TableView.create(grok.shell.table(t), true); grok.shell.addView(v); v.name = n; }, [table, name]);
  await expect.poll(() => page.evaluate((t) => Array.from(grok.shell.tableViews).filter((v: any) => v.dataFrame?.name === t).length, table),
    {message: `two views of table "${table}"`}).toBe(2);
}, {tier: 'api'});

const restoring = new WeakSet<Page>();

export const askToSaveOnClose = Given('the "Ask to save on close" setting is {word} until the feature ends', async (page: Page, state: string) => {
  if (state !== 'on' && state !== 'off')
    throw new Error(`"${state}": the setting is on or off`);
  const was = await page.evaluate((on) => {
    const w = grok.shell.settings.askToSaveOnClose;
    grok.shell.settings.askToSaveOnClose = on;
    return w;
  }, state === 'on');
  // a Background step runs before every scenario: the first one that changes the setting puts it back
  if (was !== (state === 'on') && !restoring.has(page)) {
    restoring.add(page);
    atFeatureEnd(page, async () => {
      restoring.delete(page);
      await page.evaluate((w) => { grok.shell.settings.askToSaveOnClose = w; }, was);
    });
  }
}, {tier: 'api', description: 'whether a close also asks about unsaved changes (off by default); what depends on a closing table is asked regardless'});

export const clickDialogLink = When('user clicks on the {string} link in {string} dialog', async (page: Page, text: string, title: string) => {
  const dialog = page.locator('.d4-dialog').filter({has: page.locator('.d4-dialog-title', {hasText: title})}).filter({visible: true}).first();
  const link = dialog.locator('a', {hasText: new RegExp(`^\\s*${escapeRegExp(text)}\\s*$`)}).first();
  await expect(link, `the "${text}" link in the "${title}" dialog`).toBeVisible();
  await link.click();
}, {tier: 'ui', description: 'a link in a dialog\'s text, such as "save the project" in the close dialog, which saves first and then closes'});

export const simpleModeOn = Given('simple mode is on', async (page: Page) => {
  const was = await page.evaluate(() => { const s = grok.shell.windows.simpleMode; grok.shell.windows.simpleMode = true; return s; });
  await expect(page.locator('[name="view selector"]').filter({visible: true}).first(), 'the view selector').toBeVisible({timeout: 15000});
  if (!was)
    atFeatureEnd(page, () => page.evaluate(() => { grok.shell.windows.simpleMode = false; }));
}, {tier: 'api', description: 'the simple shell, views switched and closed from the ribbon\'s view selector; a filmed run otherwise shows the full shell'});

export const setColumnTag =When('code sets the {string} tag of column {string} of table {string} to {string}', (page: Page, tag: string, column: string, table: string, value: string) =>
  page.evaluate(([k, c, t, v]) => { grok.shell.table(t).col(c).setTag(k, v); }, [tag, column, table, value]),
{tier: 'api', description: 'a column tag set the way a package or the column\'s context panel sets it; which tags count as the user\'s work is by the key'});

export const closeFromViewSelector = When('user closes the {string} view from the view selector', async (page: Page, view: string) => {
  await page.locator('[name="view selector"]').filter({visible: true}).first().click();
  const close = page.locator(`.d4-list-item[name="${view}"] [aria-label="Close view"]`).filter({visible: true}).first();
  await expect(close, `the close icon of "${view}" in the view selector`).toBeVisible();
  await close.click();
}, {tier: 'ui', description: 'the x beside a view in the ribbon view selector, the way the simple mode closes a view'});

export const closeCurrentTab = When('user closes the tab of the current view', async (page: Page) => {
  const name = await page.evaluate(() => grok.shell.v.name);
  await page.locator(`.tab-handle[name="view-handle: ${name}"] > .tab-handle-close-button`).filter({visible: true}).first().click();
}, {tier: 'ui', description: 'the close button on the tab of whatever view is current, for a view whose name the scenario does not know'});

export const noViewOfType = Then('no {word} view should be open', (page: Page, type: string) =>
  expect.poll(() => page.evaluate((t) => [...grok.shell.views].filter((v: any) => v.type === t).map((v: any) => v.name), type),
    {message: `open ${type} views`}).toEqual([]),
{description: 'by the view type, as "the current view should be a {word} view" names it'});

async function pickFromMenuOf(page: Page, target: Locator, what: string, item: string, at?: (box: {width: number; height: number}) => {x: number; y: number}): Promise<void> {
  await expect(target, what).toBeVisible();
  const box = at ? await target.boundingBox() : null;
  await target.click({button: 'right', position: box ? at!(box) : undefined});
  const entry = page.locator('.d4-menu-popup .d4-menu-item')
    .filter({has: page.locator(':scope > .d4-menu-item-label', {hasText: new RegExp(`^\\s*${escapeRegExp(item)}\\s*$`)})})
    .filter({visible: true}).first();
  await expect(entry, `"${item}" in the context menu of ${what}`).toBeVisible();
  await entry.click();
}

export const pickFromSidebarMenu = When('user picks {string} from the context menu of the sidebar', (page: Page, item: string) =>
  pickFromMenuOf(page, page.locator('.d4-tab-header-stripe.vertical').filter({visible: true}).first(), 'the sidebar', item,
    (box) => ({x: box.width / 2, y: box.height * 0.6})),
{tier: 'ui', description: 'the menu of the left sidebar\'s empty strip between its icons, where Close All is'});

export const pickFromTableMenu = When('user picks {string} from the context menu of the table name in the status bar', (page: Page, item: string) =>
  pickFromMenuOf(page, page.locator('.d4-view-status-panel [name^="div-"]').filter({visible: true}).first(), 'the table name in the status bar', item),
{tier: 'ui', description: 'the current table\'s own menu (Add View, Close, Rename...), opened on its name in the status bar'});

export const pickFromProjectMenu = When('user picks {string} from the context menu of the {string} project in the Dashboards pane', async (page: Page, item: string, project: string) => {
  await page.locator('.layout-sidebar [name="Dashboards"], .grok-sidebar [name="Dashboards"]').filter({visible: true}).first().click();
  await pickFromMenuOf(page, page.locator(`[name="tree-${project}"]`).filter({visible: true}).first(), `the "${project}" project in the Dashboards pane`, item);
}, {tier: 'ui', description: 'an open project\'s menu in the Dashboards pane of the sidebar, which lists the open projects and their views'});
