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

const GROUP_PREFIX = 'qa_autotest_g_';

const STAMP = Date.now();
const GROUP = `qa_autotest_g_${STAMP}`;
const GROUP_RENAMED = `${GROUP}_renamed`;

const MEMBER_USER = 'opavlenko454';    
const CHILD_GROUP = 'UsersTest';        

test.describe.configure({ mode: 'serial' });

test.describe('Groups View (Groups-*)', () => {

  test.beforeAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, GROUP_PREFIX);
  });
  test.afterAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, GROUP_PREFIX);
  });

  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Groups');
  });

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

  test('Groups-03/04/05/09/14 — group lifecycle: create dialog, validation, create, rename, delete', async ({ page }) => {
    try {

      await ribbonButtonByText(page, 'NEW GROUP').click();
      await expect(page.locator(DIALOG_TITLE)).toHaveText(/Create New Group/i);
      await expect(dialogInput(page, 'Name')).toBeVisible();
      await expect(dialogInput(page, 'Description')).toBeVisible();
      await dialogButton(page, 'CANCEL').click();
      await expect(page.locator(DIALOG)).toHaveCount(0);

      await createGroup(page, GROUP, 'created by autotest');
      await openPlatformView(page, 'Groups');
      await searchAndWaitCard(page, 'groups', GROUP);

      await openCardContextMenu(page, GROUP);
      await contextMenuItemByName(page, 'Properties...').click();
      await expect(page.locator(DIALOG_TITLE)).toContainText(/Properties/i);
      await dialogInput(page, 'Name').fill(GROUP_RENAMED);
      await dialogButton(page, 'OK').click();
      await page.waitForTimeout(1000);
      await openPlatformView(page, 'Groups');
      await searchAndWaitCard(page, 'groups', GROUP_RENAMED);

      await deleteEntityViaContextMenu(page, GROUP_RENAMED);
      await openPlatformView(page, 'Groups');
      await searchGallery(page, 'groups', GROUP_RENAMED);
      await expect(galleryCardByName(page, GROUP_RENAMED), 'deleted group should be gone').toHaveCount(0);
    } finally {

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

  test('Groups-11/12/13/15 — manage members: add user, admin toggle, nest group, remove', async ({ page }) => {
    await createGroup(page, GROUP, 'members test');
    try {
      await openPlatformView(page, 'Groups');
      await searchAndWaitCard(page, 'groups', GROUP);
      await selectCard(page, GROUP);

      await openManageFromPane(page, 'Members');
      await addMembershipBySearch(page, MEMBER_USER);
      await setMemberRowToggle(page, MEMBER_USER, true);
      await addMembershipBySearch(page, CHILD_GROUP);
      await saveDialog(page);

      await selectCard(page, GROUP);
      await openManageFromPane(page, 'Members');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: MEMBER_USER })).toBeVisible({ timeout: 5_000 });
      await expect(page.locator('.d4-dialog .membership-row', { hasText: CHILD_GROUP })).toBeVisible({ timeout: 5_000 });

      await removeMembershipRow(page, MEMBER_USER);
      await saveDialog(page);
      await selectCard(page, GROUP);
      await openManageFromPane(page, 'Members');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: MEMBER_USER })).toHaveCount(0);
      await page.locator('.d4-dialog .ui-btn', { hasText: /^CANCEL$/i }).first().click();
    } finally {

      await openPlatformView(page, 'Groups');
      await searchGallery(page, 'groups', GROUP);
      if (await galleryCardByName(page, GROUP).isVisible().catch(() => false))
        await deleteEntityViaContextMenu(page, GROUP);
    }
  });

});
