import { test, expect } from '@playwright/test';
import {
  applyAutomationSetup,
  clickContextPanelSection,
  clickMenuItemExact,
  closeMenuPopup,
  expandDbConnection,
  expandDbGroupWrapper,
  expandDbProvider,
  expandTreeNode,
  goHome,
  readMenuItems,
  rightClickTreeNode,
  showContextPanel,
} from './helpers';

const PROVIDER = 'Postgres';
const CONNECTION = 'Datagrok';

const CATALOG = 'datagrok';

const CATALOGS_ROOT = `tree-Databases---${PROVIDER}---${CONNECTION}---Catalogs`;
const CATALOG_NODE = `${CATALOGS_ROOT}---${CATALOG}`;

test.describe.serial('Connections / Catalogs (Postgres / Datagrok substitution)', () => {
  test('1. Catalogs node visible; expand catalog → schema → tables', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONNECTION, 'Catalogs');

    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOGS_ROOT,
      { timeout: 30_000 },
    );

    await page.evaluate((n) => {
      const el = document.querySelector(`[name="tree-expander-${n}"]`) as HTMLElement | null;
      if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
    }, CATALOG_NODE.replace('tree-', ''));
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOG_NODE,
      { timeout: 30_000 },
    );

    const firstSchema = page.locator(`[name^="${CATALOG_NODE}---"]:not([name*="tree-expander-"])`).first();
    const firstSchemaName = (await firstSchema.getAttribute('name'))!;
    await page.evaluate((n) => {
      const el = document.querySelector(`[name="tree-expander-${n.replace(/^tree-/, '')}"]`) as HTMLElement | null;
      if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
    }, firstSchemaName);

    const firstTable = page.locator(`[name^="${firstSchemaName}---"]:not([name*="tree-expander-"])`).first();
    await firstTable.waitFor({ state: 'visible', timeout: 30_000 });
  });

  test('2. Click a catalog → Context Panel shows preview with the catalog name', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONNECTION, 'Catalogs');
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOGS_ROOT,
      { timeout: 30_000 },
    );

    const label = page.locator(
      `[name="${CATALOG_NODE}"]:not(.d4-tree-view-list-more) .d4-tree-view-group-label`,
    ).first();
    await label.click({ position: { x: 80, y: 8 }, force: true });
    await showContextPanel(page);
    await page.waitForTimeout(1000);

    const panelText = await page.evaluate(() => {
      const panel = document.querySelector('[name="div-grok-prop-panel"]')
        ?? document.querySelector('.grok-prop-panel')
        ?? document.querySelector('.d4-accordion');
      return panel?.textContent ?? '';
    });
    expect(panelText).toContain(CATALOG);
  });

  test('3. Set Comment + LLM-Comment; verify persistence after re-select', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONNECTION, 'Catalogs');
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOGS_ROOT,
      { timeout: 30_000 },
    );

    const label = page.locator(
      `[name="${CATALOG_NODE}"]:not(.d4-tree-view-list-more) .d4-tree-view-group-label`,
    ).first();
    await label.click({ position: { x: 80, y: 8 }, force: true });
    await showContextPanel(page);

    const stamp = Date.now();
    const comment = `pw-comment-${stamp}`;
    const llm = `pw-llm-${stamp}`;

    await clickContextPanelSection(page, 'Database meta');
    const commentInput = page.locator(
      '[name="input-host-Comment"] input, [name="input-host-Comment"] textarea',
    ).first();
    await commentInput.waitFor({ state: 'visible', timeout: 10_000 });
    await commentInput.click({ clickCount: 3 });
    await page.keyboard.type(comment);
    await page.keyboard.press('Tab');

    const llmInput = page.locator(
      '[name="input-host-LLM-Comment"] input, [name="input-host-LLM-Comment"] textarea',
    ).first();
    await llmInput.waitFor({ state: 'visible', timeout: 10_000 });
    await llmInput.click({ clickCount: 3 });
    await page.keyboard.type(llm);
    await page.keyboard.press('Tab');

    const siblings = page.locator(`[name^="${CATALOGS_ROOT}---"]:not([name*="tree-expander-"])`);
    const siblingCount = await siblings.count();
    if (siblingCount >= 2) {
      const second = siblings.nth(1);
      await second.click({ position: { x: 80, y: 8 }, force: true });
      await page.waitForTimeout(800);
    }
    await label.click({ position: { x: 80, y: 8 }, force: true });
    await page.waitForTimeout(1000);
    await clickContextPanelSection(page, 'Database meta');

    const persisted = await page.evaluate(() => {
      const ci = document.querySelector(
        '[name="input-host-Comment"] input, [name="input-host-Comment"] textarea',
      ) as HTMLInputElement | HTMLTextAreaElement | null;
      const li = document.querySelector(
        '[name="input-host-LLM-Comment"] input, [name="input-host-LLM-Comment"] textarea',
      ) as HTMLInputElement | HTMLTextAreaElement | null;
      return { c: ci?.value ?? '', l: li?.value ?? '' };
    });
    expect(persisted.c, 'Comment value persists across re-select').toBe(comment);
    expect(persisted.l, 'LLM-Comment value persists across re-select').toBe(llm);
  });

  test('4. Right-click catalog: Browse and Open as table are present and openable', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONNECTION, 'Catalogs');
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOGS_ROOT,
      { timeout: 30_000 },
    );

    await rightClickTreeNode(page, CATALOG_NODE);
    const items = await readMenuItems(page);
    expect(items).toEqual(expect.arrayContaining(['Browse']));
    expect(items.some((t) => /^Open as table$|^Open schema as table$/.test(t))).toBe(true);

    await clickMenuItemExact(page, 'Browse');

    await page.waitForTimeout(2000);

    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONNECTION, 'Catalogs');
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOGS_ROOT,
      { timeout: 30_000 },
    );
    await rightClickTreeNode(page, CATALOG_NODE);
    const items2 = await readMenuItems(page);
    const openAsTable = items2.find((t) => /^Open as table$|^Open schema as table$/.test(t));
    if (openAsTable) {
      await clickMenuItemExact(page, openAsTable);

      await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 60_000 });
    }
    else {
      await closeMenuPopup(page);
    }
  });
});
