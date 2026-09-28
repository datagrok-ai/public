/* What only the Projects regression features (fixed Jira tickets) need. */
import type {Page} from '@playwright/test';
import {element, Then, When} from '@datagrok-libraries/bdd';
import {expect, pollMs, takeErrors} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

const listOf = (text: string): string[] => text.split(',').map((s) => s.trim()).filter(Boolean);

/* The tabs of the document area left to right, as the dock manager shows them — a drag reorders
   them without touching grok.shell.views, so the order is read off the tab strip; only the tabs of
   table views count (Home and the Dashboards gallery keep their own). */
export const tableViewTabs = Then('the table view tabs should read {string}', async (page: Page, names: string) => {
  await expect.poll(() => page.evaluate(() => {
    const tableViews = new Set((Array.from(grok.shell.views) as any[]).filter((v) => v.type === 'TableView').map((v) => v.name));
    return (Array.from(document.querySelectorAll('.grok-tab-host .tab-handle[name^="view-handle: "]')) as HTMLElement[])
      .filter((e) => e.offsetParent !== null)
      .map((e) => e.getAttribute('name')!.slice('view-handle: '.length))
      .filter((n) => tableViews.has(n))
      .join(', ');
  }), {message: 'the table view tabs, left to right', timeout: pollMs(30000)}).toBe(listOf(names).join(', '));
}, {description: 'the visible tab handles of the document area that belong to table views, in their on-screen order'});

/* The Save project dialog renders a picture of every view by cloning it into an iframe, which logs
   "Unable to find element in cloned iframe" for some views (GROK-18606, won't fix; the operator's
   ruling of 2026-07-24: known noise that breaks nothing). Only that message is let through, only
   here — the checks after a reopen stay strict. */
const PREVIEW_NOISE = 'Unable to find element in cloned iframe';

export const noErrorsButPreviewNoise = Then('no errors but the project preview\'s should have been logged', async (page: Page) => {
  expect(takeErrors(page).filter((e) => !e.includes(PREVIEW_NOISE)), 'console errors and page errors since the last check, the Save dialog\'s preview noise aside').toEqual([]);
}, {description: 'the error floor since the last check, less the Save dialog\'s "cloned iframe" message (GROK-18606); checking clears it'});

export const closeViewByTab = When('user closes the {string} view by the cross on its tab', async (page: Page, name: string) => {
  const tab = page.locator(`.grok-tab-host .tab-handle[name="view-handle: ${name}"]`).filter({visible: true}).first();
  await expect(tab, `the tab of the "${name}" view`).toBeVisible({timeout: pollMs(15000)});
  await tab.hover();
  await tab.locator('.tab-handle-close-button').click();
  await expect(tab, `the tab of the "${name}" view after its cross was clicked`).toBeHidden();
}, {tier: 'ui', description: 'the close cross of the view\'s tab handle in the document area, not Close All'});


/* Tabs off (the status bar's Tabs toggle) is the shell without its menu and view tabs — the simple
   mode; the switcher it shows instead is the ribbon's view selector. */
export const viewTabsHidden = Then('the view tabs should be hidden', async (page: Page) => {
  await expect(page.locator('.grok-tab-host .tab-handle[name^="view-handle: "]').filter({visible: true}), 'the visible view tabs').toHaveCount(0);
  expect(await page.evaluate(() => grok.shell.windows.simpleMode), 'the simple mode the Tabs toggle switches').toBe(true);
}, {description: 'no view tab is shown in the document area and the shell is in the simple mode'});

export const viewSelectorLists = Then('the view selector should list {string}', async (page: Page, names: string) => {
  const selector = page.locator('[name="view selector"]').filter({visible: true});
  await expect(selector, 'the view selector of the ribbon').toHaveCount(1, {timeout: pollMs(30000)});
  for (const name of listOf(names))
    await expect(selector.first().locator(`.d4-list-item[name="${name}"]`), `the "${name}" entry of the view selector`).toHaveCount(1);
}, {description: 'the ribbon\'s view switcher is shown and holds an entry for each view'});

/* A project card's thumbnail is the picture the Save dialog took of the project's view
   (`/api/entities/picture/<id>`): a card that shows another project's picture is the same URL. */
export const galleryPictures = Then('the {int} gallery cards should show {int} different pictures', async (page: Page, cards: number, pictures: number) => {
  await expect.poll(() => page.evaluate(() => {
    const urls = (Array.from(document.querySelectorAll('.grok-gallery-grid .d4-gallery-card')) as HTMLElement[])
      .filter((e) => e.offsetParent !== null)
      .map((e) => (e.querySelector('.grok-gallery-grid-item-thumbnail') as HTMLElement | null)?.style.backgroundImage ?? '');
    return `${urls.length} cards, ${new Set(urls.filter(Boolean)).size} different pictures`;
  }), {message: 'the cards of the gallery and their pictures', timeout: pollMs(30000)}).toBe(`${cards} cards, ${pictures} different pictures`);
}, {description: 'the visible cards of the gallery and how many distinct thumbnail URLs they carry'});

/* A table of one open project, when two open projects hold tables of the same names (the tables of
   a project opened twice from equal sources): found through grok.shell.projects, never by name
   among grok.shell.tables, which would give the first project's. */
function projectTable(page: Page, project: string, table: string, read: 'selected' | 'filtered' | 'select10'): Promise<string | number> {
  return page.evaluate(([p, t, r]) => {
    const projects = (Array.from(grok.shell.projects) as any[]).filter((x) => x.friendlyName === p || x.name === p);
    if (projects.length !== 1)
      return `${projects.length} open projects named "${p}"; open: ${(Array.from(grok.shell.projects) as any[]).map((x) => x.friendlyName).join(', ')}`;
    const infos = projects[0].children.filter((c: any) => c instanceof DG.TableInfo);
    const df = infos.find((c: any) => c.friendlyName === t)?.dataFrame;
    if (!df)
      return `the "${p}" project has no open "${t}" table; it holds: ${infos.map((c: any) => c.friendlyName).join(', ')}`;
    if (r === 'select10') {
      df.selection.init((i: number) => i < 10);
      return 10;
    }
    return r === 'selected' ? df.selection.trueCount : df.filter.trueCount;
  }, [project, table, read] as const);
}

export const selectInProjectTable = When('user selects the first 10 rows of the {string} table of the open {string} project',
  async (page: Page, table: string, project: string) => {
    expect(await projectTable(page, project, table, 'select10'), `the "${table}" table of the open "${project}" project`).toBe(10);
  }, {tier: 'api', description: 'the selection bitset of that project\'s table — the table is named the same in both open projects, so no view or grid of it can be told apart by name'});

export const projectTableSelected = Then('the {string} table of the open {string} project should have {int} selected rows',
  async (page: Page, table: string, project: string, rows: number) => {
    await expect.poll(() => projectTable(page, project, table, 'selected'),
      {message: `selected rows of the "${table}" table of the open "${project}" project`, timeout: pollMs(15000)}).toBe(rows);
  });

export const projectTableFiltered = Then('the {string} table of the open {string} project should have {int} rows passing the filter',
  async (page: Page, table: string, project: string, rows: number) => {
    await expect.poll(() => projectTable(page, project, table, 'filtered'),
      {message: `rows passing the filter in the "${table}" table of the open "${project}" project`, timeout: pollMs(15000)}).toBe(rows);
  });

/* A query's table is refilled in place by REFRESH and by the data sync of a reopen, a moment after the
   gesture: the value is polled, and the table has to hold exactly one row. */
export const onlyRowReads = Then('the only row of the table should read {string} in the {string} column', async (page: Page, value: string, column: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const t = grok.shell.t;
    if (!t)
      return 'no table is open';
    const col = t.col(c);
    if (!col)
      return `no "${c}" column; it has: ${t.columns.names().join(', ')}`;
    return t.rowCount === 1 ? String(col.get(0)) : `${t.rowCount} rows`;
  }, column), {message: `the "${column}" value of the current table's only row`, timeout: pollMs(30000)}).toBe(value);
}, {description: 'the current table has one row, and its value in that column is the one named — polled while a refresh lands'});

/* The Save project dialog's "Share link:" line and Toolbox > Source's "Choose which parameters the
   dashboard link carries" icon (project_entity_move.dart, table_view.dart). */
element('save dialog share link', {selector: '[name="dialog-Save-project"] .d4-url-params-share-url',
  description: 'the URL of the "Share link:" line of the Save project dialog, shown for a table of a parameterized query'});
element('url parameters icon', {selector: '.d4-toolbox[caption] .d4-url-params-toolbar-icon',
  description: 'the sliders icon next to REFRESH in Toolbox > Source of a saved dashboard'});

/* The Where row of the visual query builder: a column's tag holds a condition field and the "Expose
   as function parameter" checkbox, named after the column (visual query editor). */
element('visual query name condition', {selector: '[name="input-where-condition-name"]',
  description: 'the condition of the "name" tag in the Where row of the visual query'});
element('visual query name parameter checkbox', {selector: '[name="input-where-param-name"]',
  description: 'the checkbox of the "name" tag in the Where row that exposes the condition as a parameter of the query'});

