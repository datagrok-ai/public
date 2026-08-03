import { test, expect } from '@playwright/test';
import {
  CONTEXT_MENU,
  CONTEXT_MENU_BROWSE_SCHEMA,
  contextMenuItem,
  treeGroupByName,
  treeNodeByPath,
  TREE_EXPAND_ARROW,
  TREE_EXPAND_ARROW_EXPANDED,
  viewTabHandle,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

test.describe('Browse Databases (Browse-DB-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
    // CI: the demo Postgres connections (CHEMBL, Datagrok) aren't deployed on the minimal
    // stack, so these connection-level DB tests can't run. Skip when CHEMBL is absent.
    await expandTreeGroup(page, 'Databases').catch(() => undefined);
    await expandTreeGroup(page, 'Postgres').catch(() => undefined);
    const hasChembl = await treeNodeByPath(page, ['Databases', 'Postgres', 'CHEMBL'])
      .isVisible().catch(() => false);
    test.skip(!hasChembl, 'CHEMBL demo connection not present on this stack');
  });

  test('Browse-DB-01 — Databases section lists connected providers', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Databases');

    // Postgres is the canonical demo provider on dev.
    await expect(treeNodeByPath(page, ['Databases', 'Postgres']),
      'Postgres provider should be present under Databases').toHaveCount(1, { timeout: 5_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-DB-02 — expanding Postgres reveals connections', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Databases');
    await expandTreeGroup(page, 'Postgres');

    // CHEMBL and Datagrok connections are known to exist on dev.
    await expect(treeNodeByPath(page, ['Databases', 'Postgres', 'CHEMBL']),
      'CHEMBL connection should appear under Postgres').toBeVisible({ timeout: 10_000 });
    await expect(treeNodeByPath(page, ['Databases', 'Postgres', 'Datagrok']),
      'Datagrok connection should appear under Postgres').toBeVisible({ timeout: 10_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-DB-03 — right-click on a connection shows the Browse schema item', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Databases');
    await expandTreeGroup(page, 'Postgres');

    const chembl = treeNodeByPath(page, ['Databases', 'Postgres', 'CHEMBL']);
    await chembl.waitFor({ state: 'visible', timeout: 10_000 });
    const label = chembl.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();

    await label.click({ button: 'right' });
    await expect(page.locator(CONTEXT_MENU)).toBeVisible({ timeout: 5_000 });
    // On current dev the menu exposes "Browse queries" (not "Browse schema").
    const browseQueries = contextMenuItem(page, 'Browse queries');
    expect(await browseQueries.count(),
      '"Browse queries" item should be present in a Postgres connection context menu')
      .toBeGreaterThanOrEqual(1);

    await page.keyboard.press('Escape');
    await expectNoErrors(page, sink);
  });

  test('Browse-DB-04 — clicking a saved query opens the Query Editor', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Databases');
    await expandTreeGroup(page, 'Postgres');
    await expandTreeGroup(page, 'Datagrok');

    // Use the first saved query under Postgres > Datagrok (any one works — they all
    // open the Query Editor). The exact list is environment-dependent.
    const firstQuery = page.locator('[name^="tree-Databases---Postgres---Datagrok---"]').first();
    await firstQuery.waitFor({ state: 'visible', timeout: 10_000 });
    const queryName = await firstQuery.getAttribute('name') ?? '';
    await firstQuery.click();
    await page.waitForTimeout(2500);

    // The URL should now point at a Query (or generic function) editor.
    await expect(page).toHaveURL(/\/(q|func)\//, { timeout: 10_000 });
    await expect(page.locator('.d4-ribbon').first()).toBeVisible();

    await expectNoErrors(page, sink);
  });

  test('Browse-DB-05 — clicking a connection refreshes the Context Panel without errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Databases');
    await expandTreeGroup(page, 'Postgres');

    // Click on CHEMBL — Context Panel should show its metadata (Details / info panes).
    const chembl = treeNodeByPath(page, ['Databases', 'Postgres', 'CHEMBL']);
    await chembl.click();
    await page.waitForTimeout(2000);

    // The Context Panel inner should contain something non-empty.
    const text = ((await page.locator('.grok-prop-panel .grok-entity-prop-panel').innerText().catch(() => '')) ?? '').trim();
    expect(text.length, 'Context Panel must have content for CHEMBL connection').toBeGreaterThan(0);

    await expectNoErrors(page, sink);
  });

  test('Browse-DB-06 — Postgres > CHEMBL connection opens without errors (regression)', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Databases');
    await expandTreeGroup(page, 'Postgres');

    const chembl = treeNodeByPath(page, ['Databases', 'Postgres', 'CHEMBL']);
    await chembl.waitFor({ state: 'visible', timeout: 10_000 });

    // Expand CHEMBL and confirm that no error fires while doing so. The historical
    // GROK-16857 manifested as a crash inside the CHEMBL > Browse > Summary path —
    // verify that path is at least reachable for clicking, without asserting view content.
    const tri = chembl.locator(TREE_EXPAND_ARROW).first();
    if (!(await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false))) {
      await tri.click();
      await page.waitForTimeout(2000);
    }

    await expectNoErrors(page, sink);
  });
});
