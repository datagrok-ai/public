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

// Manual scenario `catalogs.md` (order 10).
//
// Catalog browsing: tree shows `Catalogs` node, expand a catalog → schemas → tables
// Catalog preview: click catalog → Context Panel shows preview with name
// Catalog meta: Comment + LLM-Comment can be set and persist
// Catalog schema view: context menu has Browse + Open as table
//
// SCOPE NOTE: the original Test Track scenario uses `MS SQL > NorthwindTest`.
// `opavlenko+playwright@datagrok.ai` has no MS SQL connection accessible (0 of
// 18 connections this user sees are MS SQL), so we substitute with the
// `Postgres > Datagrok` connection — a catalog-capable Postgres on dev that
// 07-schema test 2 already drills into. Its catalogs are `datagrok / postgres /
// template1`; the Context Panel exposes the same `General` + `Database meta`
// accordion with Comment / LLM-Comment input hosts. The structural assertions
// (catalog tree → schema → tables, catalog preview, meta-prop persistence,
// catalog right-click menu) are independent of provider.
//
// MS SQL fidelity: the catalog UX is identical between Postgres and MS SQL on
// the platform side, so this substitution covers the same user-facing
// behaviour. Switching back to MS SQL is a one-line PROVIDER/CONNECTION change
// once an MS SQL fixture is shared with the playwright user.

const PROVIDER = 'Postgres';
const CONNECTION = 'Datagrok';
// Postgres `Datagrok` connection has 3 catalogs on dev — pick `datagrok` for
// the drill-down; siblings (`postgres`, `template1`) exercise the re-select
// path in test 3.
const CATALOG = 'datagrok';

const CATALOGS_ROOT = `tree-Databases---${PROVIDER}---${CONNECTION}---Catalogs`;
const CATALOG_NODE = `${CATALOGS_ROOT}---${CATALOG}`;

// CI: the System:Datagrok connection on the ephemeral test Datlas points to
// the single per-build platform DB (`datagrok_<network>`), so the platform
// doesn't render a "Catalogs" section under the connection at all (it
// appears only when the connection enumerates multiple catalog DBs as it
// does on dev). Skip the whole describe on the CI server; the dev-targeted
// playwright-tests/ copy still runs there against dev's multi-catalog
// Postgres > Datagrok connection.
test.describe.skip('Connections / Catalogs (Postgres / Datagrok substitution)', () => {
  test.describe.configure({ mode: 'serial' });
  test('1. Catalogs node visible; expand catalog → schema → tables', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONNECTION, 'Catalogs');

    // Catalogs root should now have catalog children.
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOGS_ROOT,
      { timeout: 30_000 },
    );

    // Expand `datagrok` catalog (clicking the dedicated tree-expander avoids
    // hover-flake on the inner `.d4-tree-view-tri`).
    await page.evaluate((n) => {
      const el = document.querySelector(`[name="tree-expander-${n}"]`) as HTMLElement | null;
      if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
    }, CATALOG_NODE.replace('tree-', ''));
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      CATALOG_NODE,
      { timeout: 30_000 },
    );

    // Pick the first schema under the catalog and drill in.
    const firstSchema = page.locator(`[name^="${CATALOG_NODE}---"]:not([name*="tree-expander-"])`).first();
    const firstSchemaName = (await firstSchema.getAttribute('name'))!;
    await page.evaluate((n) => {
      const el = document.querySelector(`[name="tree-expander-${n.replace(/^tree-/, '')}"]`) as HTMLElement | null;
      if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
    }, firstSchemaName);

    // Tables should appear under the schema.
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

    // Click the catalog label so it becomes the current object.
    const label = page.locator(
      `[name="${CATALOG_NODE}"]:not(.d4-tree-view-list-more) .d4-tree-view-group-label`,
    ).first();
    await label.click({ position: { x: 80, y: 8 }, force: true });
    await showContextPanel(page);
    await page.waitForTimeout(1000);

    // The catalog name should appear somewhere on the right-side property panel.
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

    // Comment + LLM-Comment inputs live inside the "Database meta" accordion pane.
    // Input host names: `input-host-Comment` and `input-host-LLM-Comment` (capital C).
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

    // Click away to a sibling catalog, then re-select — values should remain.
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

    // Read the current input values directly — panel.textContent may not include
    // input values, only labels.
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
    // A schema/catalog browse view opens — wait briefly for the view handle update.
    await page.waitForTimeout(2000);

    // Reopen the menu and pick "Open as table" / "Open schema as table".
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
      // Result table view should render with a grid canvas.
      await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 60_000 });
    }
    else {
      await closeMenuPopup(page);
    }
  });
});
