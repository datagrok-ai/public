import { test, expect } from '@playwright/test';
import {
  GALLERY_GRID,
  GALLERY_COUNTS,
  DIALOG,
  DIALOG_TITLE,
  dialogInput,
  dialogButton,
  contextMenuItemByName,
  galleryCardByName,
  infoPaneByName,
  ribbonButtonByText,
  VIEW_TOGGLE_BRIEF,
  VIEW_TOGGLE_CARD,
  VIEW_TOGGLE_GRID,
} from './selectors';
import {
  openPlatformView,
  readGalleryCount,
  searchGallery,
  searchAndWaitCard,
  clearGallerySearch,
  openCardContextMenu,
  closeMenu,
  createGroup,
  deleteEntityViaContextMenu,
  selectCard,
  openManageFromPane,
  addMembershipBySearch,
  removeMembershipRow,
  setMemberRowToggle,
  saveDialog,
  watchErrors,
  expectNoErrors,
  apiDeleteGroupsByPrefix,
  sweepGroupsByPrefix,
} from './helpers';

// All groups created by this suite share this prefix; used for the pre/post-run cleanup sweeps.
const GROUP_PREFIX = 'qa_autotest_g_';

// Groups can be created and deleted freely on dev (self-cleaning), so this set exercises the
// full lifecycle. Unique names per run; everything is deleted in cleanup.
const STAMP = Date.now();
const GROUP = `qa_autotest_g_${STAMP}`;
const GROUP_RENAMED = `${GROUP}_renamed`;
// Partitioned per file (see users.test.ts) so files can run on separate workers in parallel:
// this file owns opavlenko454 + nests "UsersTest".
const MEMBER_USER = 'opavlenko454';    // existing user to add as a member
const CHILD_GROUP = 'UsersTest';        // existing group to nest

test.describe.configure({ mode: 'serial' });

test.describe('Groups View (Groups-*)', () => {
  // Remove any leftover autotest groups from earlier (possibly crashed) runs before/after this file.
  test.beforeAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, GROUP_PREFIX);
  });
  test.afterAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, GROUP_PREFIX);
  });

  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Groups');
  });

  // Groups-01 + Groups-02 + Groups-07: open from tree, layout, view modes.
  test('Groups-01/02/07 — open view: gallery, count, toolbar, view modes', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page).toHaveURL(/\/groups\b/);
    await expect(page.locator(GALLERY_GRID).first()).toBeVisible();
    expect((await readGalleryCount(page)).total).toBeGreaterThan(0);

    await expect(ribbonButtonByText(page, 'NEW GROUP')).toBeVisible();
    await expect(page.locator('input[placeholder^="Search groups"]')).toBeVisible();
    await expect(page.locator(GALLERY_COUNTS)).toBeVisible();

    await page.locator(VIEW_TOGGLE_CARD).first().click();
    await expect(page.locator(VIEW_TOGGLE_CARD).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_GRID).first().click();
    await expect(page.locator(VIEW_TOGGLE_GRID).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_BRIEF).first().click();
    await expect(page.locator(VIEW_TOGGLE_BRIEF).first()).toHaveClass(/d4-current/);

    await expectNoErrors(page, sink);
  });

  // Groups-03 + Groups-04 + Groups-05 + Groups-09 + Groups-14: full CRUD lifecycle.
  test('Groups-03/04/05/09/14 — group lifecycle: create dialog, validation, create, rename, delete', async ({ page }) => {
    try {
      // Dialog fields + cancel (no creation). NOTE: unlike the New User dialog, the Create New Group
      // dialog does NOT disable OK for an empty name (no client-side validation), so we only verify
      // the fields render and CANCEL is non-destructive.
      await ribbonButtonByText(page, 'NEW GROUP').click();
      await expect(page.locator(DIALOG_TITLE)).toHaveText(/Create New Group/i);
      await expect(dialogInput(page, 'Name')).toBeVisible();
      await expect(dialogInput(page, 'Description')).toBeVisible();
      await dialogButton(page, 'CANCEL').click();
      await expect(page.locator(DIALOG)).toHaveCount(0);

      // Create.
      await createGroup(page, GROUP, 'created by autotest');
      await openPlatformView(page, 'Groups');
      await searchAndWaitCard(page, 'groups', GROUP);

      // Rename via Properties...
      await openCardContextMenu(page, GROUP);
      await contextMenuItemByName(page, 'Properties...').click();
      await expect(page.locator(DIALOG_TITLE)).toContainText(/Properties/i);
      await dialogInput(page, 'Name').fill(GROUP_RENAMED);
      await dialogButton(page, 'OK').click();
      await page.waitForTimeout(1000);
      await openPlatformView(page, 'Groups');
      await searchAndWaitCard(page, 'groups', GROUP_RENAMED);

      // Delete.
      await deleteEntityViaContextMenu(page, GROUP_RENAMED);
      await openPlatformView(page, 'Groups');
      await searchGallery(page, 'groups', GROUP_RENAMED);
      await expect(galleryCardByName(page, GROUP_RENAMED), 'deleted group should be gone').toHaveCount(0);
    } finally {
      // Safety net: delete the group under either name (`GROUP` is a prefix of `GROUP_RENAMED`)
      // if the test failed before its own UI delete ran.
      await apiDeleteGroupsByPrefix(page, GROUP);
    }
  });

  test('Groups-06 — search by name filters list and updates URL', async ({ page }) => {
    const before = await readGalleryCount(page);
    await searchGallery(page, 'groups', 'QA');
    await expect(page).toHaveURL(/\/groups\?q=QA/);
    expect((await readGalleryCount(page)).shown).toBeLessThan(before.total);
    await clearGallerySearch(page, 'groups');
    expect((await readGalleryCount(page)).shown).toBe(before.total);
  });

  // Groups-08 + Groups-10: context menu items and context-panel info panes.
  test('Groups-08/10 — group context menu items and Context Panel info panes', async ({ page }) => {
    await openCardContextMenu(page, CHILD_GROUP);
    for (const name of ['Properties...', 'Delete'])
      await expect(contextMenuItemByName(page, name), `menu item "${name}" should be present`)
        .toBeVisible({ timeout: 5_000 });
    await closeMenu(page);

    await selectCard(page, CHILD_GROUP);
    for (const pane of ['Members', 'Permissions'])
      await expect(infoPaneByName(page, pane), `pane "${pane}" should render`).toBeVisible({ timeout: 10_000 });
  });

  // Groups-11 + Groups-12 + Groups-13 + Groups-15: member management via MANAGE.
  test('Groups-11/12/13/15 — manage members: add user, admin toggle, nest group, remove', async ({ page }) => {
    await createGroup(page, GROUP, 'members test');
    try {
      await openPlatformView(page, 'Groups');
      await searchAndWaitCard(page, 'groups', GROUP);
      await selectCard(page, GROUP);

      // Add a user member, make them admin, and nest an existing group.
      await openManageFromPane(page, 'Members');
      await addMembershipBySearch(page, MEMBER_USER);
      await setMemberRowToggle(page, MEMBER_USER, true);
      await addMembershipBySearch(page, CHILD_GROUP);
      await saveDialog(page);

      // Verify both are now members.
      await selectCard(page, GROUP);
      await openManageFromPane(page, 'Members');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: MEMBER_USER })).toBeVisible({ timeout: 5_000 });
      await expect(page.locator('.d4-dialog .membership-row', { hasText: CHILD_GROUP })).toBeVisible({ timeout: 5_000 });

      // Remove the user member, keep the group; save.
      await removeMembershipRow(page, MEMBER_USER);
      await saveDialog(page);
      await selectCard(page, GROUP);
      await openManageFromPane(page, 'Members');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: MEMBER_USER })).toHaveCount(0);
      await page.locator('.d4-dialog .ui-btn', { hasText: /^CANCEL$/i }).first().click();
    } finally {
      // Cleanup: delete the group (removes all memberships).
      await openPlatformView(page, 'Groups');
      await searchGallery(page, 'groups', GROUP);
      if (await galleryCardByName(page, GROUP).isVisible().catch(() => false))
        await deleteEntityViaContextMenu(page, GROUP);
    }
  });

  // Groups-18 (personal group exists) is intentionally NOT automated here: personal groups DO
  // exist (verifiable via the API — every user's `group` is their personal group), but the Groups
  // View gallery does NOT list them (searching a user's name returns 0). So there is no UI to assert
  // against. Covered by API-level tests instead.
});
