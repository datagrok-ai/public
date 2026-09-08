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

const TEST_GROUP = 'QA';
const TEST_ROLE = 'Super Chemist';

const TARGET_USER = 'opavlenko45';

test.describe.configure({ mode: 'serial' });

test.describe('Users View (Users-*)', () => {
  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Users');
  });

  test('Users-01/02/11 — open view: tree path, gallery, count, toolbar, view modes', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page).toHaveURL(/\/users\b/);
    await expect(page.locator(GALLERY_GRID).first()).toBeVisible();
    const { total } = await readGalleryCount(page);
    expect(total).toBeGreaterThan(0);

    await expect(page.locator(NEW_USER_BUTTON)).toBeVisible();
    await expect(page.locator('input[placeholder^="Search users"]')).toBeVisible();
    await expect(page.locator(SORT_LIST).first()).toBeVisible();
    await expect(page.locator(TOGGLE_FILTERS).first()).toBeVisible();
    await expect(page.locator(GALLERY_COUNTS)).toBeVisible();

    await page.locator(VIEW_TOGGLE_CARD).first().click();
    await expect(page.locator(VIEW_TOGGLE_CARD).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_GRID).first().click();
    await expect(page.locator(VIEW_TOGGLE_GRID).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_BRIEF).first().click();
    await expect(page.locator(VIEW_TOGGLE_BRIEF).first()).toHaveClass(/d4-current/);

    await expectNoErrors(page, sink);
  });

  test('Users-03/04/06/08 — NEW create dialogs: options, fields, validation, Invite (no creation)', async ({ page }) => {
    const before = await readGalleryCount(page);

    await page.locator(NEW_USER_BUTTON).click();
    await page.waitForTimeout(400);
    for (const label of ['User...', 'Service User...', 'Invite a Friend...'])
      await expect(
        page.locator('.d4-menu-popup .d4-menu-item-label', { hasText: new RegExp(`^${label.replace(/\./g, '\\.')}$`) }),
      ).toBeVisible({ timeout: 5_000 });
    await closeMenu(page);

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

    await openNewUserOption(page, 'Invite a Friend...');
    await expect(page.locator(DIALOG_TITLE)).toHaveText(/Invite a Friend/i);
    await expect(dialogInput(page, 'Email')).toBeVisible();
    await dialogButton(page, 'CANCEL').click();
    await expect(page.locator(DIALOG)).toHaveCount(0);

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

    await openCardContextMenu(page, TARGET_USER);
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
      await page.waitForTimeout(800);
    } else await closeMenu(page);

    await openCardContextMenu(page, TARGET_USER);
    await contextMenuItemByName(page, 'Add to favorites').click();
    await page.waitForTimeout(1000);
    await openCardContextMenu(page, TARGET_USER);
    await expect(contextMenuItem(page, 'Remove from favorites'),
      'after adding, "Remove from favorites" must be present').toBeVisible({ timeout: 5_000 });

    await contextMenuItem(page, 'Remove from favorites').click();
    await page.waitForTimeout(500);
    await expectNoErrors(page, sink);
  });

  test('Users-18/19 — edit user group memberships and roles (existing user, reverted)', async ({ page }) => {

    await openUserMembershipDialog(page, TARGET_USER, 'Groups...');
    await addMembershipBySearch(page, TEST_GROUP);
    await saveDialog(page);
    try {
      await openUserMembershipDialog(page, TARGET_USER, 'Groups...');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: TEST_GROUP }),
        'added group should persist').toBeVisible({ timeout: 5_000 });
    } finally {

      await removeMembershipRow(page, TEST_GROUP);
      await saveDialog(page);
    }

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

      await openCardContextMenu(page, TARGET_USER);
      await contextMenuItemByName(page, 'Block').click();
      await page.locator('.d4-dialog button[name="button-YES"]').click();
      expect(await waitForStatus('blocked'), 'user should be blocked').toBe('blocked');

      await page.reload();
      await page.locator(GALLERY_GRID).first().waitFor({ state: 'visible', timeout: 20_000 });
      await page.waitForTimeout(800);
      await openCardContextMenu(page, TARGET_USER);
      await contextMenuItemByName(page, 'Unblock').click();
      await page.locator('.d4-dialog button[name="button-YES"]').click();
      expect(await waitForStatus('active'), 'user should be unblocked via UI').toBe('active');
    } finally {

      if (await getUserStatus(page, TARGET_USER) !== 'active') await unblockUser(page, TARGET_USER);
    }
  });
});
