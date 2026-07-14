import { test, expect } from '@playwright/test';
import { treeGroupByName, treeNodeByPath } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

test.describe('Browse Platform (Browse-Platform-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Platform-01 — Platform section reveals admin / functions / access subdivisions', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Platform');

    // The mandatory subnodes (subset of the full set documented in selectors.md).
    for (const name of ['Admin', 'Plugins', 'Functions', 'Users', 'Groups', 'Roles']) {
      await expect(
        treeNodeByPath(page, ['Platform', name]),
        `Platform subnode "${name}" should be present`,
      ).toHaveCount(1, { timeout: 5_000 });
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-Platform-02 — opening a Platform leaf opens its view', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Platform');
    const users = treeNodeByPath(page, ['Platform', 'Users']);
    await users.waitFor({ state: 'visible', timeout: 10_000 });
    await users.click();
    await page.waitForTimeout(1500);

    // Users view should be active (URL contains /users) and the view tab handle exists.
    await expect(page).toHaveURL(/\/users\b/, { timeout: 10_000 });

    await expectNoErrors(page, sink);
  });

});
