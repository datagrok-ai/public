import { test, expect } from '@playwright/test';
import { contextMenuItemByName, BALLOON_CONTAINER } from './selectors';
import {
  openPlatformView,
  openCardContextMenu,
  searchAndWaitCard,
  watchErrors,
  expectNoErrors,
  sweepGroupsByPrefix,
  createGroupAsAdmin,
} from './helpers';

const TARGET_PREFIX = 'qa_autotest_sharing_target_';
let TARGET = '';

test.describe('Groups View — membership requests (Groups-16, non-admin)', () => {
  test.beforeAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, TARGET_PREFIX);            
    TARGET = await createGroupAsAdmin(browser, `${TARGET_PREFIX}${Date.now()}`);
  });
  test.afterAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, TARGET_PREFIX);            
  });

  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Groups');
  });

  test('Groups-16 — request membership in a group as a non-admin', async ({ page }) => {
    const sink = watchErrors(page);

    await searchAndWaitCard(page, 'groups', TARGET);
    await openCardContextMenu(page, TARGET);

    const request = contextMenuItemByName(page, 'Request membership');
    await expect(request, '"Request membership" should be available to a non-member').toBeVisible({ timeout: 5_000 });
    await request.click();
    await page.waitForTimeout(1500);

    const balloon = page.locator(BALLOON_CONTAINER);
    await expect(balloon, 'a notification should confirm the request was sent').toBeVisible({ timeout: 5_000 });

    await expectNoErrors(page, sink);
  });
});
