/* Domain schemas over an external database: Create Domain Schema (PowerPack:createDomainBinding) and the
   generic /domains/<schema>/<table> pages PowerPack hosts (PowerPack:domainRouteView). The setting that
   turns them on, the schema cleanup, the server-side counts and the rows of a list or a table page. */
import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, expect, locate, pollMs, silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

element('table in schema design', {selector: '[data-u2-item^="t:"] .u2-tree-label',
  description: 'a table\'s row in the design tree of the Create domain schema wizard'});
element('domain rows', {selector: '.grok-view-container [data-u2="domain-list"] [data-u2="list"], ' +
  '.grok-view-container .u2-domain-data-table',
  description: 'the rows of a domain table\'s list page: its cards, or its table of columns'});
element('domain search box', {selector: '.grok-view-container > .d4-ribbon input[type="text"]',
  description: 'the ribbon\'s search box over the table\'s searchable columns'});
element('domain status count', {selector: '.u2-domain-status-count',
  description: 'the count slot of the status bar: "50 of N" while the list has fetched one page'});
element('domain form', {selector: '[data-u2="domain-form"]', description: 'the entity page\'s fields'});
element('domain children', {selector: '[data-u2="domain-children"]',
  description: 'the entity page\'s tabs of the rows that refer to it'});
element('domain breadcrumbs', {selector: '[data-u2="breadcrumbs"]',
  description: 'the entity page\'s trail: the table, then the row'});
element('domain filters', {selector: '.grok-view-container .u2-domain-filters',
  description: 'the condition builder beside the rows'});

async function dropSchema(page: Page, name: string): Promise<void> {
  const left = await page.evaluate(async (n) => {
    const all = async () => (await grok.dapi.domains.schemas.list()).map((s: any) => s.name);
    if ((await all()).includes(n))
      await grok.dapi.domains.schema(n).delete();
    return (await all()).includes(n);
  }, name);
  if (left)
    throw new Error(`the domain schema "${name}" is still registered after the delete`);
}

export const settingOn = Given('the {string} shell setting is on', async (page: Page, name: string) => {
  silent(page);
  await page.evaluate((n) => { (grok.shell.settings as any)[n] = true; }, name);
  await expect.poll(() => page.evaluate((n) => (grok.shell.settings as any)[n], name),
    {message: `the "${name}" shell setting`, timeout: pollMs(15000)}).toBe(true);
}, {tier: 'api', description: 'turns an account setting on (Settings › Beta), read back; pair it with the put-back step'});

export const schemaGone = Given('the domain schema {string} is removed now and at feature end', async (page: Page,
  name: string) => {
  silent(page);
  await dropSchema(page, name);
  // before the feature's connections go: a binding must not outlive the connection it reads through
  atFeatureEnd(page, () => dropSchema(page, name), true);
}, {tier: 'api', description: 'the registered schema (an external binding: its registry entry, never the database) is ' +
  'deleted now and at feature end, before any connection cleanup, and verified gone'});

export const tableCount = Then('the {string} domain table should hold {int} rows on the server', async (page: Page,
  address: string, n: number) => {
  await expect.poll(() => page.evaluate((a) => grok.dapi.domains.table(a).count(), address),
    {message: `rows of ${address}`, timeout: pollMs(60000)}).toBe(n);
}, {tier: 'api', description: 'counts the table through the domains API, which reads the bound database'});

export const openRow = When('user opens the {string} row of domain rows', async (page: Page, text: string) => {
  const rows = await locate(page, el('domain rows'));
  // the list and its table both stay in the page; only the one on show has rows to open
  const row = rows.locator('.u2-list-row, .u2-data-table-row').filter({visible: true})
    .filter({has: page.getByText(text, {exact: true})}).first();
  await expect(row).toBeVisible({timeout: pollMs(30000)});
  await row.dblclick();
}, {tier: 'ui', description: 'double-clicks the card or table row that shows exactly that text: its entity page opens'});

/** The design tree's row of [table], rendered: the tree is virtual, and a table further down exists only once the
 * list scrolls to it. Centred from inside the page, so an action on it does not scroll the list again. */
async function treeRow(page: Page, table: string, part: string) {
  const target = page.locator(`[data-u2-item="t:${table}"] ${part}`).first();
  const list = page.locator('[data-u2="dialog"] .u2-tree .u2-list').first();
  await list.evaluate((el) => el.scrollTop = 0);
  for (let i = 0; i < 100 && await target.count() === 0; i++) {
    await list.evaluate((el) => {
      el.scrollTop += el.clientHeight / 2;
      return new Promise((frame) => requestAnimationFrame(frame));
    });
  }
  return target;
}

async function clickInTree(page: Page, table: string, part: string, done: () => Promise<void>): Promise<void> {
  // the virtual list renders its rows again whenever it scrolls or the selection moves: each try finds the row anew,
  // centres it from inside the page and clicks where it is
  await expect(async () => {
    const target = await treeRow(page, table, part);
    const at = await target.evaluate(async (el) => {
      el.scrollIntoView({block: 'center'});
      await new Promise((frame) => requestAnimationFrame(frame));
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }, undefined, {timeout: 2000});
    await page.mouse.click(at.x, at.y);
    await done();
  }).toPass({timeout: pollMs(20000)});
}

export const leaveOut = When('user leaves the {string} table out of schema design', async (page: Page, table: string) => {
  await clickInTree(page, table, 'input[type="checkbox"]',
    async () => expect(await treeRow(page, table, 'input[type="checkbox"]')).not.toBeChecked({timeout: 2000}));
}, {tier: 'ui', description: 'unticks a table in the wizard\'s design tree: it is not bound, and links into it become plain values'});

export const pickTable = When('user clicks on the {string} table in schema design', async (page: Page, table: string) => {
  const dialog = page.locator('[data-u2="dialog"]').first();
  await clickInTree(page, table, '.u2-tree-label',
    () => expect(dialog).toContainText(new RegExp(`Table\\s*${table}(?![a-z0-9_])`), {timeout: 2000}));
}, {tier: 'ui', description: 'scrolls the wizard\'s design tree to a table and selects it, so the panel beside it edits that table'});

export const followLink = When('user follows the {string} link in domain form', async (page: Page, text: string) => {
  const form = await locate(page, el('domain form'));
  const link = form.locator('a, .u2-link, [role="link"]').filter({hasText: text}).first();
  await expect(link).toBeVisible({timeout: pollMs(30000)});
  await link.click();
}, {tier: 'ui', description: 'clicks a reference field\'s link on the entity page: the referenced row\'s page opens'});

export const columnSemType = Then('the {string} column of the {string} domain table should be recognized as {string}',
  async (page: Page, column: string, address: string, semType: string) => {
    await expect.poll(() => page.evaluate(async ([a, c]) => {
      const properties = await grok.dapi.domains.registry.rowProperties(a);
      return properties.find((p: any) => p.name === c)?.semType ?? null;
    }, [address, column] as const), {message: `the semantic type of ${address}.${column}`, timeout: pollMs(30000)})
      .toBe(semType);
  }, {tier: 'api', description: 'reads the bound column\'s semantic type from the registry: what Create Domain Schema ' +
    'recognized in a sample of its values'});

export const goBack = When('user goes back', async (page: Page) => {
  await page.goBack({waitUntil: 'commit'});
}, {tier: 'ui', description: 'the browser\'s Back: the domain app follows the address to the page before'});

export const childRows = Then('domain children should show their rows', async (page: Page) => {
  const grid = page.locator('[data-u2="domain-children"] .u2-domain-grid').filter({visible: true}).first();
  await expect.poll(async () => (await grid.boundingBox())?.height ?? 0,
    {message: 'the height of the open tab\'s grid', timeout: pollMs(60000)}).toBeGreaterThan(0);
}, {tier: 'ui', description: 'the open tab\'s grid has its first page: a tab over millions of rows counts them before it fetches them'});
