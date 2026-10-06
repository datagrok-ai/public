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
  createRole,
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

const ROLE_PREFIX = 'qa_autotest_r_';

const STAMP = Date.now();
const ROLE = `qa_autotest_r_${STAMP}`;
const ROLE_RENAMED = `${ROLE}_renamed`;

const ASSIGNEE_USER = 'opavlenko656';  

test.describe.configure({ mode: 'serial' });

test.describe('Roles View (Roles-*)', () => {

  test.beforeAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, ROLE_PREFIX);
  });
  test.afterAll(async ({ browser }) => {
    await sweepGroupsByPrefix(browser, ROLE_PREFIX);
  });

  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Roles');
  });

  test('Roles-01/02/07 — open view: gallery, count, toolbar, view modes', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page).toHaveURL(/\/roles\b/);
    await expect(page.locator(GALLERY_GRID).first()).toBeVisible();
    expect((await readGalleryCount(page)).total).toBeGreaterThan(0);

    await expect(ribbonButtonByText(page, 'New Role')).toBeVisible();
    await expect(page.locator('input[placeholder^="Search roles"]')).toBeVisible();
    await expect(page.locator(GALLERY_COUNTS)).toBeVisible();

    await page.locator(VIEW_TOGGLE_CARD).first().click();
    await expect(page.locator(VIEW_TOGGLE_CARD).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_GRID).first().click();
    await expect(page.locator(VIEW_TOGGLE_GRID).first()).toHaveClass(/d4-current/);
    await page.locator(VIEW_TOGGLE_BRIEF).first().click();
    await expect(page.locator(VIEW_TOGGLE_BRIEF).first()).toHaveClass(/d4-current/);

    await expectNoErrors(page, sink);
  });

  test('Roles-03/04/09/15 — role lifecycle: create dialog, create, rename, delete', async ({ page }) => {
    try {

      await ribbonButtonByText(page, 'New Role').click();
      await expect(page.locator(DIALOG_TITLE)).toHaveText(/Create New Role/i);
      await expect(dialogInput(page, 'Name')).toBeVisible();
      await expect(dialogInput(page, 'Description')).toBeVisible();
      await dialogButton(page, 'CANCEL').click();
      await expect(page.locator(DIALOG)).toHaveCount(0);

      await createRole(page, ROLE, 'created by autotest');
      await openPlatformView(page, 'Roles');
      await searchAndWaitCard(page, 'roles', ROLE);

      await openCardContextMenu(page, ROLE);
      await contextMenuItemByName(page, 'Properties...').click();
      await expect(page.locator(DIALOG_TITLE)).toContainText(/Properties/i);
      await dialogInput(page, 'Name').fill(ROLE_RENAMED);
      await dialogButton(page, 'OK').click();
      await page.waitForTimeout(1000);
      await openPlatformView(page, 'Roles');
      await searchAndWaitCard(page, 'roles', ROLE_RENAMED);

      await deleteEntityViaContextMenu(page, ROLE_RENAMED);
      await openPlatformView(page, 'Roles');
      await searchGallery(page, 'roles', ROLE_RENAMED);
      await expect(galleryCardByName(page, ROLE_RENAMED), 'deleted role should be gone').toHaveCount(0);
    } finally {

      await apiDeleteGroupsByPrefix(page, ROLE);
    }
  });

  test('Roles-06 — search by name filters list and updates URL', async ({ page }) => {
    const before = await readGalleryCount(page);
    await searchGallery(page, 'roles', 'Admin');
    await expect(page).toHaveURL(/\/roles\?q=Admin/);
    expect((await readGalleryCount(page)).shown).toBeLessThan(before.total);
    await clearGallerySearch(page, 'roles');
    expect((await readGalleryCount(page)).shown).toBe(before.total);
  });

  test('Roles-08/10 — role context menu items and Context Panel info panes', async ({ page }) => {

    await createRole(page, ROLE, 'menu test');
    try {
      await openPlatformView(page, 'Roles');
      await searchAndWaitCard(page, 'roles', ROLE);
      await openCardContextMenu(page, ROLE);
      for (const name of ['Properties...', 'Delete'])
        await expect(contextMenuItemByName(page, name), `menu item "${name}" should be present`)
          .toBeVisible({ timeout: 5_000 });
      await closeMenu(page);

      await selectCard(page, ROLE);
      for (const pane of ['Assigned to', 'Permissions'])
        await expect(infoPaneByName(page, pane.replace(/ /g, '-')), `pane "${pane}" should render`)
          .toBeVisible({ timeout: 10_000 });
    } finally {
      await openPlatformView(page, 'Roles');
      await searchGallery(page, 'roles', ROLE);
      if (await galleryCardByName(page, ROLE).isVisible().catch(() => false))
        await deleteEntityViaContextMenu(page, ROLE);
    }
  });

  test('Roles-11/12/13 — assign role to a user, Can-assign toggle, remove', async ({ page }) => {
    await createRole(page, ROLE, 'assignment test');
    try {
      await openPlatformView(page, 'Roles');
      await searchAndWaitCard(page, 'roles', ROLE);
      await selectCard(page, ROLE);

      await openManageFromPane(page, 'Assigned to');
      await addMembershipBySearch(page, ASSIGNEE_USER);
      await setMemberRowToggle(page, ASSIGNEE_USER, true);
      await saveDialog(page);

      await selectCard(page, ROLE);
      await openManageFromPane(page, 'Assigned to');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: ASSIGNEE_USER })).toBeVisible({ timeout: 5_000 });

      await removeMembershipRow(page, ASSIGNEE_USER);
      await saveDialog(page);
      await selectCard(page, ROLE);
      await openManageFromPane(page, 'Assigned to');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: ASSIGNEE_USER })).toHaveCount(0);
      await page.locator('.d4-dialog .ui-btn', { hasText: /^CANCEL$/i }).first().click();
    } finally {
      await openPlatformView(page, 'Roles');
      await searchGallery(page, 'roles', ROLE);
      if (await galleryCardByName(page, ROLE).isVisible().catch(() => false))
        await deleteEntityViaContextMenu(page, ROLE);
    }
  });

});
