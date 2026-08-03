import { test, expect } from '@playwright/test';
import {
  CONTEXT_PANEL,
  CONTEXT_PANEL_INNER,
  LIST_SEARCH_INPUT,
  treeItemByName,
  viewTabHandle,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
} from './helpers';

test.describe('Browse Dashboards (Browse-Dash-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Dash-01 — Dashboards node opens the Dashboards view', async ({ page }) => {
    const sink = watchErrors(page);

    // Dashboards is a top-level item (not a group): clicking opens the Dashboards view.
    const dashboards = treeItemByName(page, 'Dashboards');
    await expect(dashboards, 'Dashboards should be visible at the top level')
      .toBeVisible({ timeout: 10_000 });
    await dashboards.click();
    await page.waitForTimeout(1500);

    // Clicking the "Dashboards" tree item navigates to the Projects view (/projects).
    await expect(viewTabHandle(page, 'Projects'),
      'Projects view (opened from the Dashboards tree node) should be selected')
      .toHaveClass(/tab-handle-selected/, { timeout: 10_000 });
    await expect(page).toHaveURL(/\/projects\b/);

    await expectNoErrors(page, sink);
  });

  test('Browse-Dash-02 — opening chemical_space_demo and demo-datagrok-api dashboards', async ({ page }) => {
    const sink = watchErrors(page);

    // Both dashboards confirmed present on dev (friendlyName ↔ name mapping):
    //   chemical_space_demo → ChemicalSpaceDemo
    //   demo-datagrok-api   → DemoDatagrokApi
    for (const fullName of ['Datagrok.ChemicalSpaceDemo', 'Datagrok.DemoDatagrokApi']) {
      await page.goto(`${process.env.DATAGROK_URL!}/p/${fullName}`);
      await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
      await page.waitForTimeout(2500);

      // The dashboard opens as a project view — assert there's no error balloon
      // and the ribbon is still present (= app did not crash).
      await expect(page.locator('.d4-ribbon').first(), `Ribbon visible for ${fullName}`).toBeVisible();
      await expectNoErrors(page, sink);
    }
  });

  test('Browse-Dash-03 — Context Panel content differs when switching between dashboards', async ({ page }) => {
    const sink = watchErrors(page);

    // Open the Projects view. The gallery only shows featured projects by default —
    // use the gallery search to surface specific dashboards.
    await page.goto(`${process.env.DATAGROK_URL!}/projects`);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
    await page.waitForTimeout(1500);

    /**
     * Search for a project and click it, then poll the Context Panel until it
     * actually reflects the selected entity (substring match) — avoids stale
     * content from a previous selection.
     */
    async function selectDashboard(searchTerm: string, displayName: RegExp, expectedMarker: RegExp): Promise<string> {
      const searchInput = page.locator('input[placeholder^="Search projects"]').first();
      await searchInput.fill('');
      await searchInput.fill(searchTerm);
      await page.waitForTimeout(2000);
      const proj = page.locator('.grok-gallery-grid-item-title', { hasText: displayName }).first();
      await expect(proj, `Project matching "${searchTerm}" must be visible`).toBeVisible({ timeout: 10_000 });
      await proj.click();

      // Wait for Context Panel to update — it should contain the selected entity's marker.
      await expect.poll(
        async () => ((await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '').toLowerCase(),
        {
          message: `Context Panel must reflect the selected project "${searchTerm}"`,
          timeout: 15_000,
          intervals: [400, 700, 1000, 1500],
        },
      ).toMatch(expectedMarker);

      return ((await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '').toLowerCase();
    }

    const contentA = await selectDashboard(
      'chemical_space_demo', /^chemical_space_demo$/i, /chemical[\s_-]?space/i,
    );
    const contentB = await selectDashboard(
      'demo-datagrok-api', /^demo-datagrok-api$/i, /demo[\s_-]?datagrok[\s_-]?api/i,
    );

    expect(contentA, 'Context Panel for A should be non-empty').not.toBe('');
    expect(contentB, 'Context Panel for B should be non-empty').not.toBe('');
    expect(contentA, 'Context Panel content must differ for two distinct dashboards (ref: GROK-19934)')
      .not.toEqual(contentB);

    await expectNoErrors(page, sink);
  });
});
