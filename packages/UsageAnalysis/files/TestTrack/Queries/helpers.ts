import { Page, expect } from '@playwright/test';

export const POSTGRES_CONNECTION = 'Northwind';
export const MS_SQL_CONNECTION = 'Northwind';

export const AUTH_STATE = 'e2e/.auth.public.json';

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
  await page.waitForTimeout(1000);
}

export async function expandDbProvider(page: Page, provider: string): Promise<void> {

  await expandTreeNode(page, 'tree-Databases');
  await expandTreeNode(page, `tree-Databases---${provider.replace(/ /g, '-')}`);
}

export async function expandDbConnection(page: Page, provider: string, connection: string): Promise<void> {
  await expandTreeNode(page, `tree-Databases---${provider.replace(/ /g, '-')}---${connection.replace(/ /g, '-')}`);
}

export async function openDbConnectionView(page: Page, provider: string, connection: string): Promise<void> {
  const nodeName = `tree-Databases---${provider.replace(/ /g, '-')}---${connection.replace(/ /g, '-')}`;
  const label = page
    .locator(`[name="${nodeName}"]:not(.d4-tree-view-list-more) .d4-tree-view-group-label`)
    .first();
  await label.waitFor({ state: 'visible', timeout: 15_000 });
  await label.scrollIntoViewIfNeeded();
  await label.click();
}

export async function showContextPanel(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: { shell: { windows: { showProperties: boolean } } } }).grok;
    g.shell.windows.showProperties = true;
  });
}

export async function openTransformationsTab(page: Page): Promise<void> {
  await page.locator('[name="Transformations"]').first().click();
  await page.waitForSelector('[data-source="tab-content-Transformations"]', { state: 'attached', timeout: 10_000 });

  await page.waitForFunction(() => {
    const firstStep = document.querySelector('.grok-action-editor');

    return firstStep !== null && !firstStep.classList.contains('running');
  }, undefined, { timeout: 15_000 });
}

export async function clickTransformationAction(page: Page, actionName: string): Promise<void> {

  await page.evaluate((name) => {
    const label = Array.from(document.querySelectorAll('.grok-actions-browser label'))
      .find((el) => el.textContent?.trim() === name) as HTMLElement | undefined;
    label?.click();
  }, actionName);
  await page.waitForSelector('.d4-dialog', { timeout: 10_000 });
}

export async function addNewColumnTransformation(page: Page, expression: string): Promise<void> {
  await clickTransformationAction(page, 'Add New Column');

  const cm = page.locator('.add-new-column-dialog-cm-div .cm-content');
  await cm.focus();
  await page.keyboard.type(expression);
  await page.locator('[name="button-Add-New-Column---OK"]').click();
  await expect(page.locator('.d4-dialog')).toHaveCount(0, { timeout: 10_000 });
}

export async function transformationStepNames(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.grok-action-editor'))
    .map((el) => el.querySelector('#name')?.textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

export async function deleteTransformationStep(page: Page, index: number): Promise<void> {
  const step = page.locator('.grok-action-editor').nth(index);
  await step.scrollIntoViewIfNeeded();
  await step.hover();
  await step.locator('[name="icon-delete"]').click({ force: true });
  await page.waitForTimeout(500);
}

export async function findProjectByFriendlyName(page: Page, friendlyName: string): Promise<{
  id: string; name: string; friendlyName: string;
} | null> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const ps = await g.dapi.projects.filter(`friendlyName = "${fn}"`).list();
    if (!ps.length) return null;
    return { id: ps[0].id, name: ps[0].name, friendlyName: ps[0].friendlyName };
  }, friendlyName);
}

export async function deleteProjectByFriendlyName(page: Page, friendlyName: string): Promise<number> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const ps = await g.dapi.projects.filter(`friendlyName = "${fn}"`).list();
    for (const p of ps) await g.dapi.projects.delete(p);
    return ps.length;
  }, friendlyName);
}

export async function getConnectionServerName(page: Page, provider: string, friendlyName: string): Promise<string> {
  return page.evaluate(async ({ p, f }) => {
    const g = (window as unknown as { grok: any }).grok;
    const c = (await g.dapi.connections.filter(`friendlyName = "${f}" and dataSource = "${p}"`).list())[0];
    return c?.name as string;
  }, { p: provider, f: friendlyName });
}

export async function expandDbSchemas(page: Page, provider: string, connServerName: string): Promise<void> {
  const groupSelector = `[name="div-${provider}-${connServerName}-Schemas"]`;
  await page.waitForSelector(groupSelector, { state: 'attached', timeout: 15_000 });
  await page.evaluate((sel) => {
    const group = document.querySelector(sel);
    const innerNode = group?.querySelector(':scope > .d4-tree-view-node');
    (innerNode as HTMLElement | null)?.scrollIntoView({ block: 'center' });
    const tri = innerNode?.querySelector(':scope > .d4-tree-view-tri');
    if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) (tri as HTMLElement).click();
  }, groupSelector);
  await page.waitForTimeout(1500);
}

export async function selectTreeNodeAsCurrentObject(page: Page, nodeName: string): Promise<void> {
  const node = page.locator(`[name="${nodeName}"]:not(.d4-tree-view-list-more)`).first();
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();

  await node.click({ position: { x: 80, y: 8 } });
  await page.waitForTimeout(250);
}

export async function getCurrentObjectName(page: Page): Promise<string | null> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    const o = g.shell.o;
    if (!o) return null;
    return o.name ?? String(o);
  });
}

export async function getVisibleErrorBalloons(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-balloon-error, .grok-balloon-error, [class*="balloon"][class*="error"]'))
    .filter((b) => (b as HTMLElement).getBoundingClientRect().width > 0)
    .map((b) => (b as HTMLElement).textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

export async function listDbTableColumnNodeNames(page: Page, tableNodeName: string): Promise<string[]> {
  return page.evaluate((t) => Array.from(document.querySelectorAll(`[name^="${t}---"]`))
    .map((n) => n.getAttribute('name')!)
    .filter((s) => !!s), tableNodeName);
}

export async function focusBrowseSidebar(page: Page): Promise<void> {
  await page.locator('[name="Browse"]').first().click();
  await page.waitForTimeout(500);
}

export function queryTreeNodeSuffix(friendlyName: string): string {
  return friendlyName.replace(/_/g, '-');
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

export async function setQueryName(page: Page, name: string): Promise<void> {
  const input = page.locator('[name="input-Name"]');
  await input.waitFor({ timeout: 10_000 });

  await input.click({ clickCount: 3 });
  await page.keyboard.type(name);
  await expect(input).toHaveValue(name);
}

export async function typeQuerySql(page: Page, sql: string): Promise<void> {
  await page.waitForSelector('.CodeMirror', { state: 'visible', timeout: 20_000 });

  await page.waitForFunction(() => {
    const el = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: unknown } | null;
    return !!el?.CodeMirror;
  }, undefined, { timeout: 10_000 });

  const surface = page.locator('.CodeMirror .CodeMirror-code').first();
  await surface.click();

  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.keyboard.type(sql, { delay: 15 });

  await expect.poll(async () => page.evaluate(() => {
    const el = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: { getValue: () => string } } | null;
    return el?.CodeMirror?.getValue() ?? '';
  }), { timeout: 5_000 }).toBe(sql);
}

export async function setQuerySql(page: Page, sql: string): Promise<void> {

  await page.waitForSelector('.CodeMirror', { state: 'attached', timeout: 20_000 });
  await page.waitForFunction(() => {
    const cm = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: unknown } | null;
    return !!cm?.CodeMirror;
  }, undefined, { timeout: 10_000 });
  await page.evaluate((value) => {
    const cm = document.querySelector('.CodeMirror') as unknown as { CodeMirror: { setValue: (v: string) => void } };
    cm.CodeMirror.setValue(value);
  }, sql);
}

export async function setQueryPostProcessScript(page: Page, friendlyName: string, script: string): Promise<void> {
  await page.evaluate(async ({ fn, s }) => {
    const grok = (window as unknown as { grok: any }).grok;
    const queries = await grok.dapi.queries.filter(`friendlyName = "${fn}"`).list();
    if (!queries.length) throw new Error(`query "${fn}" not found on server`);
    const q = queries[0];
    q.postProcessScript = s;
    await grok.dapi.queries.save(q, false);
  }, { fn: friendlyName, s: script });
}

export async function waitForQuerySql(page: Page, expected: string): Promise<void> {
  await page.waitForSelector('.CodeMirror', { state: 'attached', timeout: 20_000 });
  await expect.poll(async () => page.evaluate(() => {
    const cm = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: { getValue: () => string } } | null;
    return cm?.CodeMirror?.getValue() ?? '';
  }), { timeout: 20_000 }).toBe(expected);
}

export async function runQueryViaPlay(page: Page): Promise<void> {
  await page.locator('[name="icon-play"]').first().click();
  try {
    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 30_000 });
  }
  catch {

    await page.locator('[name="icon-play"]').first().click();
    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 30_000 });
  }
}

export async function runQueryViaActions(page: Page, queryName: string): Promise<void> {

  const selector = `[name="view-handle: ${queryName}"]`;
  const before = await page.locator(selector).count();
  await page.locator('label.d4-link-action')
    .filter({ hasText: /^Run query\.\.\.$/ })
    .first()
    .click();
  await page.waitForFunction(({ sel, prev }) =>
    document.querySelectorAll(sel).length > prev, { sel: selector, prev: before }, { timeout: 30_000 });
}

export async function focusQueryEditorTab(page: Page, queryName: string): Promise<void> {
  await page.evaluate((name) => {
    const handles = Array.from(document.querySelectorAll(`[name="view-handle: ${name}"]`));
    const editor = handles.find((h) => h.querySelector('[name="icon-data-query"]')) as HTMLElement | undefined;
    editor?.click();
  }, queryName);
  await page.waitForTimeout(400);
}

export async function saveQuery(page: Page, friendlyName: string): Promise<void> {
  await page.locator('[name="button-Save"]').first().click();

  await expect.poll(async () =>
    (await findQueryByFriendlyName(page, friendlyName)) !== null,
  { timeout: 30_000 }).toBe(true);
}

export async function findQueryByFriendlyName(page: Page, friendlyName: string): Promise<{
  name: string; friendlyName: string; connection: string | undefined;
} | null> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const qs = await g.dapi.queries.filter(`friendlyName = "${fn}"`).list();
    if (!qs.length) return null;
    const q = qs[0];
    return { name: q.name, friendlyName: q.friendlyName, connection: q.connection?.name };
  }, friendlyName);
}

export async function deleteQueryByFriendlyName(page: Page, friendlyName: string): Promise<number> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const qs = await g.dapi.queries.filter(`friendlyName = "${fn}"`).list();
    for (const q of qs) await g.dapi.queries.delete(q);
    return qs.length;
  }, friendlyName);
}
