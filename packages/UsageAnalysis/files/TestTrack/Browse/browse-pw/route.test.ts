import { test, expect } from '@playwright/test';
import { BALLOON_CONTAINER, RIBBON } from './selectors';
import { goHome, watchErrors, expectNoErrors } from './helpers';

const BASE: string = process.env.DATAGROK_URL!;

test.describe('Browse routing (Browse-Route-*)', () => {
  test('Browse-Route-01 — opening a saved dashboard URL restores the view', async ({ page }) => {
    const sink = watchErrors(page);

    // Direct URL to a known dashboard.
    await page.goto(`${BASE}/p/Datagrok.ChemicalSpaceDemo`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });
    await page.waitForTimeout(2500);

    await expect(page).toHaveURL(/\/p\/Datagrok\.ChemicalSpaceDemo/);
    await expect(page.locator(RIBBON).first(),
      'Ribbon must be visible after restoring a dashboard URL').toBeVisible();

    await expectNoErrors(page, sink);
  });

  test('Browse-Route-02 — navigating to a Files folder URL opens the folder view', async ({ page }) => {
    const sink = watchErrors(page);

    await page.goto(`${BASE}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });

    // We should land on a Files folder view; the URL stays under /files/.
    await expect(page).toHaveURL(/\/files\//);

    // The folder content (e.g., demog.csv) should be visible.
    const demog = page.locator('label, .d4-tree-view-item-label', { hasText: /^demog\.csv$/i }).first();
    await expect(demog, 'demog.csv should be visible inside the opened Demo folder')
      .toBeVisible({ timeout: 15_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Route-03 — opening a saved query URL navigates to the query editor', async ({ page }) => {
    const sink = watchErrors(page);

    // Navigate to a saved query that we know exists on dev (Postgres > Datagrok > World).
    await page.goto(`${BASE}/q/Datagrok.Datagrok.World`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });
    await page.waitForTimeout(2500);

    await expect(page).toHaveURL(/\/q\/Datagrok\.Datagrok\.World/);
    await expect(page.locator(RIBBON).first(),
      'Ribbon must be visible after opening a saved query URL').toBeVisible();

    await expectNoErrors(page, sink);
  });

  test('Browse-Route-04 — invalid URL does not crash the app', async ({ page }) => {
    const sink = watchErrors(page);

    await page.goto(`${BASE}/p/no.such_project/none`);
    // The app should still be loaded (ribbon visible), even if the route is invalid.
    await expect(page.locator(RIBBON).first(), 'App ribbon should remain visible after invalid URL')
      .toBeVisible({ timeout: 30_000 });

    // We allow info / warning balloons; specifically we forbid hard JS errors.
    await expectNoErrors(page, sink);
  });
});
