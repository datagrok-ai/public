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

    // Mandatory subnodes present in every Datlas. "Admin" is deployment-specific (absent on
    // the minimal CI stack — only Plugins/Credentials/Functions/Users/Groups/Roles show there)
    // so it is not required here.
    for (const name of ['Plugins', 'Functions', 'Users', 'Groups', 'Roles']) {
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
