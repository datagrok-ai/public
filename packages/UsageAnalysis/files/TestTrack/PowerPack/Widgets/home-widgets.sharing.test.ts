import { test, expect } from './fixtures';
import { widgetByTitle } from './selectors';
import { resetHome, watchErrors, expectNoErrors } from './helpers';

/**
 * canView gating — runs under DATAGROK_SHARING_LOGIN (non-admin) via the `chromium-sharing`
 * Playwright project (see playwright.config.ts). Usage / Reports are restricted to
 * Developers / Administrators, so a non-admin user must not see them.
 */
test.describe('Home page Widgets — permissions (non-admin)', () => {
  test.beforeEach(async ({ homePage: page }) => {
    test.setTimeout(180_000);
    await resetHome(page);
  });

  test('Widgets-Perm-01 — Usage and Reports are hidden for a non-admin user', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    // Gated widgets must be absent.
    await expect(widgetByTitle(page, 'Usage'), 'non-admin must not see Usage').toHaveCount(0);
    await expect(widgetByTitle(page, 'Reports'), 'non-admin must not see Reports').toHaveCount(0);

    // Ungated widgets are still shown.
    await expect(widgetByTitle(page, 'Spotlight')).toBeVisible({ timeout: 15_000 });
    await expect(widgetByTitle(page, 'Community')).toBeVisible({ timeout: 15_000 });

    await expectNoErrors(page, sink);
  });
});
