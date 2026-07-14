import { test, expect } from '@playwright/test';
import { treeGroupByName, treeNodeByPath, BALLOON_CONTAINER, BROWSE_HEADER } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

/**
 * Permission / role-based tests — run under DATAGROK_SHARING_LOGIN (non-admin).
 * Gated by the `chromium-sharing` Playwright project (see playwright.config.ts).
 */
test.describe('Browse permissions (sharing user, non-admin)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Platform-03 — non-admin user does not see Platform (or no destructive actions)', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page.locator(BROWSE_HEADER), 'Browse panel should be visible').toBeVisible();

    // If Platform is visible, clicking it must not produce errors / error balloons.
    const platform = treeGroupByName(page, 'Platform');
    if (await platform.isVisible().catch(() => false)) {
      await platform.click();
      await page.waitForTimeout(1500);
    }

    await expect(
      page.locator(`${BALLOON_CONTAINER} .d4-balloon-error, ${BALLOON_CONTAINER} .grok-balloon-error`),
      'Non-admin must not see an error balloon when interacting with Platform',
    ).toHaveCount(0);

    await expectNoErrors(page, sink);
  });

  test('Browse-Tree-07 — clicking an object the user does not have access to does not crash', async ({ page }) => {
    const sink = watchErrors(page);

    // Open Databases > Postgres > CHEMBL. Non-admin may have limited access on it
    // depending on share settings; the test asserts only that no crash / error balloon
    // appears when clicking it. If the entity isn't visible at all (also a valid
    // permission outcome) — that's fine too.
    await expandTreeGroup(page, 'Databases');
    const chembl = treeNodeByPath(page, ['Databases', 'Postgres', 'CHEMBL']);
    if (await chembl.isVisible().catch(() => false)) {
      await chembl.click();
      await page.waitForTimeout(1500);
    }

    // No error balloons.
    await expect(
      page.locator(`${BALLOON_CONTAINER} .d4-balloon-error, ${BALLOON_CONTAINER} .grok-balloon-error`),
      'A restricted-access click must not produce an error balloon',
    ).toHaveCount(0);

    await expectNoErrors(page, sink);
  });

  test('Browse-MyStuff-04 — Shared with me node opens without errors', async ({ page }) => {
    const sink = watchErrors(page);

    // Without a controlled sharing fixture we can't assert exact grouping by user;
    // verify that the Shared-with-me subnode is reachable and clicking it does not
    // produce console / pageerror / balloon errors under the sharing user.
    await expandTreeGroup(page, 'My stuff');
    const shared = treeNodeByPath(page, ['My-stuff', 'Shared-with-me']);
    await expect(shared, 'Shared with me node must exist').toHaveCount(1, { timeout: 5_000 });
    await shared.scrollIntoViewIfNeeded();
    await shared.click();
    await page.waitForTimeout(1500);

    await expectNoErrors(page, sink);
  });
});
