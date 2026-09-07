import { Page, expect } from '@playwright/test';

export const AUTH_STATE = 'e2e/.auth.json';

export const PG_SERVER = process.env.DG_PG_SERVER ?? 'db.datagrok.ai';
export const PG_PORT = process.env.DG_PG_PORT ?? '54322';
export const PG_DB = process.env.DG_PG_DB ?? 'northwind';
export const PG_LOGIN = process.env.DG_PG_LOGIN ?? 'datagrok';
export const PG_PASSWORD = process.env.DG_PG_PASSWORD ?? '';

export const PG_EXT_SERVER = process.env.DG_PG_EXT_SERVER ?? 'db.datagrok.ai';
export const PG_EXT_PORT = process.env.DG_PG_EXT_PORT ?? '54327';
export const PG_EXT_DB = process.env.DG_PG_EXT_DB ?? 'test';
export const PG_EXT_LOGIN = process.env.DG_PG_EXT_LOGIN ?? 'superuser';
export const PG_EXT_PASSWORD = process.env.DG_PG_EXT_PASSWORD ?? '';

export const SPARQL_ENDPOINT =
  process.env.DG_SPARQL_ENDPOINT ?? 'http://data.ontotext.com/repositories/data-last';

export async function goHome(page: Page): Promise<void> {
  await page.goto('/', { waitUntil: 'domcontentloaded', timeout: 30_000 });
  await page.locator('[name="Browse"]').first().waitFor({ state: 'visible', timeout: 60_000 });
  await page.waitForTimeout(500);

  await page.addStyleTag({ content: '.d4-tooltip { display: none !important; }' });

  const databasesRoot = treeNodeLocator(page, 'tree-Databases');
  if (!(await databasesRoot.isVisible().catch(() => false))) {
    await page.locator('[name="Browse"]').first().click();
    await databasesRoot.waitFor({ state: 'visible', timeout: 10_000 });
  }
}

function treeNodeLocator(page: Page, nodeName: string) {
  return page.locator(`[name="${nodeName}"]:not(.d4-tree-view-list-more)`).first();
}

export async function expandTreeNode(page: Page, nodeName: string): Promise<void> {
  const node = treeNodeLocator(page, nodeName);
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();
  const tri = node.locator('.d4-tree-view-tri').first();
  const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded'))
    .catch(() => false);
  if (!expanded)
    await tri.click();
  await page.waitForTimeout(800);
}

export async function expandDbProvider(page: Page, provider: string): Promise<void> {

  await expandTreeNode(page, 'tree-Databases');
  await expandTreeNode(page, providerNodeName(provider));
}

export async function expandDbConnection(
  page: Page, provider: string, connection: string,
): Promise<void> {
  await expandTreeNode(page, connectionNodeName(provider, connection));
}

export function providerNodeName(provider: string): string {
  return `tree-Databases---${provider.replace(/ /g, '-')}`;
}

export function connectionNodeName(provider: string, connection: string): string {
  return `${providerNodeName(provider)}---${connection.replace(/_/g, '-').replace(/ /g, '-')}`;
}

export async function showAllDatabaseProviders(page: Page): Promise<void> {
  await expandTreeNode(page, 'tree-Databases');

  let clicked = false;
  const moreByClass = page.locator('[name="tree-Databases"] .d4-tree-view-list-more').first();
  if (await moreByClass.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await moreByClass.scrollIntoViewIfNeeded();
    await moreByClass.click();
    clicked = true;
  }
  else {
    clicked = await page.evaluate(() => {
      const root = document.querySelector('[name="tree-Databases"]');
      if (!root) return false;
      const candidates = Array.from(root.querySelectorAll('div, span, label')) as HTMLElement[];
      const more = candidates.find((el) => {
        const t = el.textContent?.trim() ?? '';
        return (t === '...' || t === 'Show more' || t === 'more') && el.offsetParent !== null;
      });
      if (!more) return false;
      more.scrollIntoView({ block: 'center' });
      more.click();
      return true;
    });
  }
  if (clicked) await page.waitForTimeout(1500);
}

export async function clickTreeExpander(page: Page, expanderName: string): Promise<void> {
  await page.waitForSelector(`[name="${expanderName}"]`, { state: 'attached', timeout: 15_000 });
  await page.evaluate((name) => {
    const el = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
    if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
  }, expanderName);
  await page.waitForTimeout(1200);
}

export async function expandDbGroupWrapper(
  page: Page, provider: string, connServerName: string, group: 'Catalogs' | 'Schemas',
): Promise<string> {
  const wrapper = `[name="div-${provider.replace(/ /g, '-')}-${connServerName}-${group}"]`;

  await page.waitForSelector(wrapper, { state: 'attached', timeout: 45_000 });
  await page.evaluate((sel) => {
    const w = document.querySelector(sel);
    const innerNode = w?.querySelector(':scope > .d4-tree-view-node');
    (innerNode as HTMLElement | null)?.scrollIntoView({ block: 'center' });
    const tri = innerNode?.querySelector(':scope > .d4-tree-view-tri');
    if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) (tri as HTMLElement).click();
  }, wrapper);
  await page.waitForTimeout(1500);
  return wrapper;
}

export async function expandChildByLabel(
  page: Page, parentName: string, label: string,
): Promise<void> {
  await page.waitForSelector(`[name="${parentName}"]`, { state: 'attached', timeout: 15_000 });
  await page.evaluate(({ p, l }) => {
    const parent = document.querySelector(`[name="${p}"]`);
    if (!parent) throw new Error(`parent tree node not found: ${p}`);

    const candidates = Array.from(parent.querySelectorAll(
      '.d4-tree-view-group-label, .d4-tree-view-item-label',
    )) as HTMLElement[];
    const match = candidates.find((el) => el.textContent?.trim() === l);
    if (!match) throw new Error(`tree node "${l}" not found under ${p}`);
    const node = match.closest('.d4-tree-view-node') as HTMLElement | null;
    if (!node) throw new Error(`tree-view-node ancestor missing for "${l}"`);
    node.scrollIntoView({ block: 'center' });
    const tri = node.querySelector(':scope > .d4-tree-view-tri') as HTMLElement | null;
    if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) tri.click();
  }, { p: parentName, l: label });
  await page.waitForTimeout(1500);
}

export async function rightClickChildByLabel(
  page: Page, parentName: string, label: string,
): Promise<void> {
  await page.evaluate(({ p, l }) => {
    const parent = document.querySelector(`[name="${p}"]`);
    if (!parent) throw new Error(`parent tree node not found: ${p}`);
    const candidates = Array.from(parent.querySelectorAll(
      '.d4-tree-view-group-label, .d4-tree-view-item-label',
    )) as HTMLElement[];
    const match = candidates.find((el) => el.textContent?.trim() === l);
    if (!match) throw new Error(`tree node "${l}" not found under ${p}`);
    const node = match.closest('.d4-tree-view-node') as HTMLElement | null;
    if (!node) throw new Error(`tree-view-node ancestor missing for "${l}"`);
    node.scrollIntoView({ block: 'center' });
    node.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
    }));
  }, { p: parentName, l: label });
  await page.waitForSelector('.d4-menu-popup', { timeout: 5_000 });
}

export async function openDbConnectionView(
  page: Page, provider: string, connection: string,
): Promise<void> {
  const label = page
    .locator(
      `[name="${connectionNodeName(provider, connection)}"]:not(.d4-tree-view-list-more)`
      + ' .d4-tree-view-group-label',
    )
    .first();
  await label.waitFor({ state: 'visible', timeout: 15_000 });
  await label.scrollIntoViewIfNeeded();
  await label.click();
}

export async function rightClickTreeNode(page: Page, nodeName: string): Promise<void> {
  const node = treeNodeLocator(page, nodeName);
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();

  try {
    await node.click({ button: 'right', position: { x: 80, y: 8 }, timeout: 5_000 });
  }
  catch {
    await node.evaluate((el) => el.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
    })));
  }
  await page.waitForSelector('.d4-menu-popup', { timeout: 5_000 });
}

export async function clickMenuItemExact(page: Page, text: string): Promise<void> {
  const escaped = text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  await page.locator('.d4-menu-popup .d4-menu-item-label')
    .filter({ hasText: new RegExp(`^${escaped}$`) })
    .first()
    .click();
}

export async function readMenuItems(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
    .map((el) => (el as HTMLElement).textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

export async function closeMenuPopup(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.waitForTimeout(150);
}

export async function openAddConnectionDialog(page: Page, provider: string): Promise<void> {
  await expandDbProvider(page, provider);
  await rightClickTreeNode(page, providerNodeName(provider));

  const items = await readMenuItems(page);
  const label = items.find((t) => /^Add connection\.\.\.$|^New connection\.\.\.$/.test(t));
  if (!label) throw new Error(`Add/New connection menu item not found for ${provider}: got ${JSON.stringify(items)}`);
  await clickMenuItemExact(page, label);
  await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });
}

export async function fillConnectionField(
  page: Page, host: string, value: string,
): Promise<void> {
  const input = page.locator(`.d4-dialog [name="input-host-${host}"] input`).first();
  await input.waitFor({ state: 'visible', timeout: 10_000 });
  await input.click({ clickCount: 3 });
  await page.keyboard.type(value);

  await expect(input).toHaveValue(value, { timeout: 5_000 });
}

export async function readConnectionField(page: Page, host: string): Promise<string> {
  const input = page.locator(`.d4-dialog [name="input-host-${host}"] input`).first();
  return input.inputValue();
}

export async function selectConnectionField(
  page: Page, host: string, value: string,
): Promise<void> {
  const scoped = page.locator(`.d4-dialog [name="input-host-${host}"] select`).first();
  const fallback = page.locator('.d4-dialog select').first();
  const target = (await scoped.count()) ? scoped : fallback;
  await target.waitFor({ state: 'visible', timeout: 10_000 });
  await target.selectOption(value);
}

export async function clickConnectionTest(page: Page, timeout = 60_000): Promise<void> {

  await page.evaluate(() => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
    .forEach((b) => (b as HTMLElement).remove()));
  await page.locator('.d4-dialog [name="button-TEST"]').click();

  await page.waitForFunction(
    () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
      .some((b) => (b as HTMLElement).textContent && (b as HTMLElement).textContent!.trim().length > 0),
    undefined,
    { timeout },
  );
}

export async function clickConnectionOk(page: Page): Promise<void> {
  await page.locator('.d4-dialog [name="button-OK"]').click();
  await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 15_000 });
}

export async function clickConnectionSave(page: Page): Promise<void> {

  const ok = page.locator('.d4-dialog [name="button-OK"]');
  const save = page.locator('.d4-dialog [name="button-SAVE"]');
  if (await save.isVisible().catch(() => false))
    await save.click();
  else
    await ok.click();
  await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 15_000 });
}

export async function clickConnectionCancel(page: Page): Promise<void> {
  const cancel = page.locator('.d4-dialog [name="button-CANCEL"]').first();
  if (await cancel.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await cancel.click();
    await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 5_000 }).catch(() => null);
  }
}

export async function deleteTreeNodeViaContext(page: Page, nodeName: string): Promise<void> {
  await rightClickTreeNode(page, nodeName);
  await clickMenuItemExact(page, 'Delete...');

  await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });
  const del = page.locator('.d4-dialog [name="button-DELETE"]');
  const yes = page.locator('.d4-dialog [name="button-YES"]');
  if (await del.isVisible().catch(() => false))
    await del.click();
  else
    await yes.click();
  await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 10_000 });
}

export async function testConnectionViaContext(
  page: Page, nodeName: string, timeout = 60_000,
): Promise<void> {

  await page.evaluate(() => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
    .forEach((b) => (b as HTMLElement).remove()));
  await rightClickTreeNode(page, nodeName);
  await clickMenuItemExact(page, 'Test connection');
  await page.waitForFunction(
    () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
      .some((b) => (b as HTMLElement).textContent && (b as HTMLElement).textContent!.trim().length > 0),
    undefined,
    { timeout },
  );
}

export async function getAllBalloonsText(page: Page): Promise<string> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
    .map((b) => (b as HTMLElement).textContent ?? '')
    .join(' | '));
}

export async function findConnectionByFriendlyName(
  page: Page, friendlyName: string,
): Promise<{ id: string; name: string; friendlyName: string; dataSource: string } | null> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const cs = await g.dapi.connections.filter(`friendlyName = "${fn}"`).list();
    if (!cs.length) return null;
    const c = cs[0];
    return { id: c.id, name: c.name, friendlyName: c.friendlyName, dataSource: c.dataSource };
  }, friendlyName);
}

export async function deleteConnectionByFriendlyName(
  page: Page, friendlyName: string,
): Promise<number> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const cs = await g.dapi.connections.filter(`friendlyName = "${fn}"`).list();
    for (const c of cs) await g.dapi.connections.delete(c);
    return cs.length;
  }, friendlyName);
}

export async function refreshBrowseTree(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;

    g.shell.windows.simpleMode = true;
  });
  await page.locator('[name="Browse"]').first().click();
  await page.waitForTimeout(500);
}

export async function applyAutomationSetup(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
  });
}

export async function showContextPanel(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    g.shell.windows.showProperties = true;
  });
}

export async function clickContextPanelSection(
  page: Page, section: string,
): Promise<{ paneTextContent: () => Promise<string> }> {
  const escaped = section.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

  const header = page.locator('.d4-accordion-pane-header')
    .filter({ hasText: new RegExp(`^${escaped}`) })
    .first();
  await header.waitFor({ state: 'visible', timeout: 10_000 });
  await header.scrollIntoViewIfNeeded();
  await header.click();
  await page.waitForTimeout(400);
  return {
    paneTextContent: async () => page.evaluate((t) => {
      const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'));
      const h = headers.find((el) => (el.textContent ?? '').trim().startsWith(t)) as HTMLElement | undefined;
      return h?.parentElement?.textContent ?? '';
    }, section),
  };
}
