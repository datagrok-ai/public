import { test, expect } from '@playwright/test';
import {
  CONTEXT_MENU,
  contextMenuItem,
  contextMenuItemByName,
  GALLERY_GRID,
  GALLERY_COUNTS,
  NEW_USER_BUTTON,
  VIEW_TOGGLE_BRIEF,
  VIEW_TOGGLE_CARD,
  VIEW_TOGGLE_GRID,
  SORT_LIST,
  TOGGLE_FILTERS,
  DIALOG,
  DIALOG_TITLE,
  dialogInput,
  dialogButton,
  galleryCardByName,
  infoPaneByName,
} from './selectors';
import {
  openPlatformView,
  readGalleryCount,
  searchGallery,
  clearGallerySearch,
  openCardContextMenu,
  openNewUserOption,
  openUserMembershipDialog,
  addMembershipBySearch,
  removeMembershipRow,
  saveDialog,
  getUserStatus,
  unblockUser,
  closeMenu,
  watchErrors,
  expectNoErrors,
} from './helpers';

// Existing entities used as edit targets (present on dev). The role/group are added to and
// removed from TARGET_USER within the test so its membership is unchanged afterwards.
// Resources are partitioned per file so the three files can run on separate workers in parallel
// without racing on shared server-side state (this file owns opavlenko45 + group "QA").
const TEST_GROUP = 'QA';
const TEST_ROLE = 'Super Chemist';

// An existing disposable test user to act on. Must NOT be the runner (opavlenko+playwright)
// or the sharing account (opavlenko+pwsharing). All changes are reverted in-test.
const TARGET_USER = 'opavlenko45';

// Several tests mutate the same TARGET_USER (favorites, memberships, block) — run serially so
// parallel workers don't race on that shared server-side state.
test.describe.configure({ mode: 'serial' });

test.describe('Users View (Users-*)', () => {
  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Users');
  });

  // Users-01 + Users-02 + Users-11: open from tree, layout, and view modes.
  test('Users-01/02/11 — open view: tree path, gallery, count, toolbar, view modes', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page).toHaveURL(/\/users\b/);
    await expect(page.locator(GALLERY_GRID).first()).toBeVisible();
    const { total } = await readGalleryCount(page);
    expect(total).toBeGreaterThan(0);

    // Toolbar controls.
    await expect(page.locator(NEW_USER_BUTTON)).toBeVisible();
    await expect(page.locator('input[placeholder^="Search users"]')).toBeVisible();
    await expect(page.locator(SORT_LIST).first()).toBeVisible();
    await expect(page.locator(TOGGLE_FILTERS).first()).toBeVisible();
    await expect(page.locator(GALLERY_COUNTS)).toBeVisible();

    // View modes: card -> grid -> brief, each becomes current.
    await page.locator(VIEW_TOGGLE_CARD).first().click();
    await expect(page.locator(VIEW_TOGGLE_CARD).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_GRID).first().click();
    await expect(page.locator(VIEW_TOGGLE_GRID).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_BRIEF).first().click();
    await expect(page.locator(VIEW_TOGGLE_BRIEF).first()).toHaveClass(/d4-current/);

    await expectNoErrors(page, sink);
  });

  // Users-03 + Users-04 + Users-06 + Users-08: every NEW create dialog, always cancelled
  // (this dev set never creates users — they cannot be hard-deleted).
  test('Users-03/04/06/08 — NEW create dialogs: options, fields, validation, Invite (no creation)', async ({ page }) => {
    const before = await readGalleryCount(page);

    // NEW menu lists all three options.
    await page.locator(NEW_USER_BUTTON).click();
    await page.waitForTimeout(400);
    for (const label of ['User...', 'Service User...', 'Invite a Friend...'])
      await expect(
        page.locator('.d4-menu-popup .d4-menu-item-label', { hasText: new RegExp(`^${label.replace(/\./g, '\\.')}$`) }),
      ).toBeVisible({ timeout: 5_000 });
    await closeMenu(page);

    // Create new user: fields + client-side validation; OK stays disabled for empty / invalid email.
    await openNewUserOption(page, 'User...');
    await expect(page.locator(DIALOG_TITLE)).toHaveText(/Create new user/i);
    for (const cap of ['Email', 'Login', 'First Name', 'Last Name'])
      await expect(dialogInput(page, cap), `field "${cap}" should be present`).toBeVisible();
    await expect(dialogButton(page, 'OK'), 'OK disabled with empty email').toHaveClass(/disabled/);
    await dialogInput(page, 'Email').fill('not-an-email');
    await page.waitForTimeout(300);
    await expect(dialogButton(page, 'OK'), 'OK disabled with invalid email').toHaveClass(/disabled/);
    await dialogButton(page, 'CANCEL').click();
    await expect(page.locator(DIALOG)).toHaveCount(0);

    // Invite a Friend: Email field, cancel (never send from CI).
    await openNewUserOption(page, 'Invite a Friend...');
    await expect(page.locator(DIALOG_TITLE)).toHaveText(/Invite a Friend/i);
    await expect(dialogInput(page, 'Email')).toBeVisible();
    await dialogButton(page, 'CANCEL').click();
    await expect(page.locator(DIALOG)).toHaveCount(0);

    // Nothing was created.
    const after = await readGalleryCount(page);
    expect(after.total).toBe(before.total);
  });

  test('Users-09 — search by name filters list and updates URL', async ({ page }) => {
    const before = await readGalleryCount(page);
    await searchGallery(page, 'users', 'opavlenko');
    await expect(page).toHaveURL(/\/users\?q=opavlenko/);
    const filtered = await readGalleryCount(page);
    expect(filtered.shown).toBeLessThan(before.total);
    expect(filtered.shown).toBeGreaterThan(0);
    await clearGallerySearch(page, 'users');
    expect((await readGalleryCount(page)).shown).toBe(before.total);
  });

  test('Users-14 — user context menu items', async ({ page }) => {
    await openCardContextMenu(page, TARGET_USER);
    for (const name of ['Details', 'Chat', 'Block', 'Groups...', 'Roles...', 'Add to favorites'])
      await expect(contextMenuItemByName(page, name), `menu item "${name}" should be present`)
        .toBeVisible({ timeout: 5_000 });
    await closeMenu(page);
  });

  // Users-15 + Users-17: single-click populates the Context Panel; double-click opens the profile view.
  test('Users-15/17 — user profile and Context Panel info panes', async ({ page }) => {
    const sink = watchErrors(page);
    const card = galleryCardByName(page, TARGET_USER);
    await card.waitFor({ state: 'visible', timeout: 15_000 });

    await card.click();
    await page.waitForTimeout(1500);
    for (const pane of ['Personal', 'Roles', 'Member-of'])
      await expect(infoPaneByName(page, pane), `pane "${pane}" should render in Context Panel`)
        .toBeVisible({ timeout: 10_000 });

    await card.dblclick();
    await page.waitForTimeout(2000);
    expect(await page.evaluate(() => (window as any).grok?.shell?.v?.type)).toBe('user_edit');

    await expectNoErrors(page, sink);
  });

  test('Users-21 — add / remove a user from favorites (roundtrip)', async ({ page }) => {
    const sink = watchErrors(page);

    // Pre-clean: if already a favorite, remove it first.
    await openCardContextMenu(page, TARGET_USER);
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
      await page.waitForTimeout(800);
    } else await closeMenu(page);

    // Add, then confirm the menu flips to "Remove from favorites".
    await openCardContextMenu(page, TARGET_USER);
    await contextMenuItemByName(page, 'Add to favorites').click();
    await page.waitForTimeout(1000);
    await openCardContextMenu(page, TARGET_USER);
    await expect(contextMenuItem(page, 'Remove from favorites'),
      'after adding, "Remove from favorites" must be present').toBeVisible({ timeout: 5_000 });

    // Cleanup.
    await contextMenuItem(page, 'Remove from favorites').click();
    await page.waitForTimeout(500);
    await expectNoErrors(page, sink);
  });

  // Users-18 + Users-19: edit an existing user's group memberships and roles via the
  // MembershipEditor, then revert — net membership unchanged.
  test('Users-18/19 — edit user group memberships and roles (existing user, reverted)', async ({ page }) => {
    // --- Groups ---
    await openUserMembershipDialog(page, TARGET_USER, 'Groups...');
    await addMembershipBySearch(page, TEST_GROUP);
    await saveDialog(page);
    try {
      await openUserMembershipDialog(page, TARGET_USER, 'Groups...');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: TEST_GROUP }),
        'added group should persist').toBeVisible({ timeout: 5_000 });
    } finally {
      // Revert: remove the group again.
      await removeMembershipRow(page, TEST_GROUP);
      await saveDialog(page);
    }

    // --- Roles ---
    await openUserMembershipDialog(page, TARGET_USER, 'Roles...');
    await addMembershipBySearch(page, TEST_ROLE);
    await saveDialog(page);
    try {
      await openUserMembershipDialog(page, TARGET_USER, 'Roles...');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: TEST_ROLE }),
        'added role should persist').toBeVisible({ timeout: 5_000 });
    } finally {
      await removeMembershipRow(page, TEST_ROLE);
      await saveDialog(page);
    }
  });

  // Users-20: block AND unblock a user fully through the UI. After Block the gallery card is stale
  // (it doesn't react to ENTITY_MODIFIED), so the menu only flips to "Unblock" after the view is
  // refreshed — the test reloads to exercise the real UI unblock. API is used only for verification
  // reads and as a finally safety-net.
  test('Users-20 — block then unblock a user via UI (refresh clears the stale card)', async ({ page }) => {
    const waitForStatus = async (want: string): Promise<string> => {
      let s = '';
      for (let i = 0; i < 15; i++) {
        s = await getUserStatus(page, TARGET_USER);
        if (s === want) break;
        await page.waitForTimeout(500);
      }
      return s;
    };
    try {
      // Block.
      await openCardContextMenu(page, TARGET_USER);
      await contextMenuItemByName(page, 'Block').click();
      await page.locator('.d4-dialog button[name="button-YES"]').click();
      expect(await waitForStatus('blocked'), 'user should be blocked').toBe('blocked');

      // Refresh the view so the card picks up the new state, then unblock via the UI.
      await page.reload();
      await page.locator(GALLERY_GRID).first().waitFor({ state: 'visible', timeout: 20_000 });
      await page.waitForTimeout(800);
      await openCardContextMenu(page, TARGET_USER);
      await contextMenuItemByName(page, 'Unblock').click();
      await page.locator('.d4-dialog button[name="button-YES"]').click();
      expect(await waitForStatus('active'), 'user should be unblocked via UI').toBe('active');
    } finally {
      // Safety net: guarantee the account is restored even if the UI path failed mid-way.
      if (await getUserStatus(page, TARGET_USER) !== 'active') await unblockUser(page, TARGET_USER);
    }
  });
});
