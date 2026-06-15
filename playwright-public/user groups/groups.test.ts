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
  searchAndWaitGone,
  selectCard,
  openManageFromPane,
  addMembershipBySearch,
  removeMembershipRow,
  setMemberRowToggle,
  saveDialog,
  ensureUserSeeded,
  ensureGroupSeeded,
  watchErrors,
  expectNoErrors,
} from './helpers';

// CI/CD variant of the dev `e2e/user groups/groups.test.ts`. Groups can be created and deleted
// freely (self-cleaning), so the lifecycle cases are unchanged. The member-user and the nested
// child group don't exist on a fresh CI DB, so they're seeded once through the UI in `beforeAll`.
const STAMP = Date.now();
const GROUP = `qa_autotest_g_${STAMP}`;
const GROUP_RENAMED = `${GROUP}_renamed`;
const MEMBER_USER = 'opavlenko454';    // seeded user added as a member
const CHILD_GROUP = 'UsersTest';        // seeded group nested into GROUP

test.describe.configure({ mode: 'serial' });

test.describe('Groups View (Groups-*)', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: 'e2e/.auth.json' });
    const page = await ctx.newPage();
    try {
      await ensureUserSeeded(page, MEMBER_USER);
      await ensureGroupSeeded(page, CHILD_GROUP);
    } finally {
      await ctx.close();
    }
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
    await searchAndWaitGone(page, 'groups', GROUP_RENAMED);
  });

  test('Groups-06 — search by name filters list and updates URL', async ({ page }) => {
    const before = await readGalleryCount(page);
    await searchGallery(page, 'groups', CHILD_GROUP);
    await expect(page).toHaveURL(new RegExp(`/groups\\?q=${CHILD_GROUP}`));
    const filtered = await readGalleryCount(page);
    expect(filtered.shown).toBeLessThan(before.total);
    expect(filtered.shown).toBeGreaterThan(0);
    await clearGallerySearch(page, 'groups');
    expect((await readGalleryCount(page)).shown).toBe(before.total);
  });

  // Groups-08 + Groups-10: context menu items and context-panel info panes.
  test('Groups-08/10 — group context menu items and Context Panel info panes', async ({ page }) => {
    await searchAndWaitCard(page, 'groups', CHILD_GROUP);
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
