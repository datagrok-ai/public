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

// Runs ONLY under the `chromium-sharing` project (secondary, non-admin account
// `opavlenko+pwsharing`) — see playwright.config.ts. A non-admin can request to join a group it
// is not a member of; that is the cross-account action this file covers.
//
// The target is a dedicated group created fresh per run (in an ADMIN context) rather than a shared
// fixture group. This is deliberate: the server's "notify admins of the request" step crashes with
// `NoSuchMethodError: getter 'id' on null` when any admin of the target group is itself a GROUP
// (not a user) — `child.user` is null for a group-admin (datlas groups_service requestMembership).
// A freshly created group's only admin is its creator (a real user), so the request succeeds. The
// underlying platform bug (group-admins break membership requests) is tracked separately.
const TARGET_PREFIX = 'qa_autotest_sharing_target_';
let TARGET = '';

test.describe('Groups View — membership requests (Groups-16, non-admin)', () => {
  test.beforeAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, TARGET_PREFIX);            // clear leftover from prior runs
    TARGET = await createGroupAsAdmin(browser, `${TARGET_PREFIX}${Date.now()}`);
  });
  test.afterAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, TARGET_PREFIX);            // delete target (+ its requests)
  });

  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Groups');
  });

  test('Groups-16 — request membership in a group as a non-admin', async ({ page }) => {
    const sink = watchErrors(page);

    // Find the freshly seeded target group (the non-admin is not a member of it).
    await searchAndWaitCard(page, 'groups', TARGET);
    await openCardContextMenu(page, TARGET);

    // A non-admin, non-member sees "Request membership".
    const request = contextMenuItemByName(page, 'Request membership');
    await expect(request, '"Request membership" should be available to a non-member').toBeVisible({ timeout: 5_000 });
    await request.click();
    await page.waitForTimeout(1500);

    // A confirmation balloon ("Request sent...") should appear, with no server error.
    const balloon = page.locator(BALLOON_CONTAINER);
    await expect(balloon, 'a notification should confirm the request was sent').toBeVisible({ timeout: 5_000 });

    await expectNoErrors(page, sink);
  });
});
