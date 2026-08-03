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
  searchAndWaitGone,
  selectCard,
  openManageFromPane,
  addMembershipBySearch,
  removeMembershipRow,
  setMemberRowToggle,
  saveDialog,
  ensureUserSeeded,
  watchErrors,
  expectNoErrors,
} from './helpers';

// CI/CD variant of the dev `e2e/user groups/roles.test.ts`. Roles can be created and deleted
// freely. The assignee user does not exist on a fresh CI DB, so it is seeded once via the UI.
const STAMP = Date.now();
const ROLE = `qa_autotest_r_${STAMP}`;
const ROLE_RENAMED = `${ROLE}_renamed`;
const ASSIGNEE_USER = 'opavlenko656';  // seeded user the role is assigned to
// A built-in role present on a fresh stack, used by the search case.
const BUILTIN_ROLE_SEARCH = 'Admin';

test.describe.configure({ mode: 'serial' });

test.describe('Roles View (Roles-*)', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: 'e2e/.auth.json' });
    const page = await ctx.newPage();
    try {
      await ensureUserSeeded(page, ASSIGNEE_USER);
    } finally {
      await ctx.close();
    }
  });

  test.beforeEach(async ({ page }) => {
    await openPlatformView(page, 'Roles');
  });

  // Roles-01 + Roles-02 + Roles-07: open from tree, layout, view modes.
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

  // Roles-03 + Roles-04 + Roles-09 + Roles-15: full CRUD lifecycle.
  test('Roles-03/04/09/15 — role lifecycle: create dialog, create, rename, delete', async ({ page }) => {
    // Dialog fields + cancel (no creation).
    await ribbonButtonByText(page, 'New Role').click();
    await expect(page.locator(DIALOG_TITLE)).toHaveText(/Create New Role/i);
    await expect(dialogInput(page, 'Name')).toBeVisible();
    await expect(dialogInput(page, 'Description')).toBeVisible();
    await dialogButton(page, 'CANCEL').click();
    await expect(page.locator(DIALOG)).toHaveCount(0);

    // Create.
    await createRole(page, ROLE, 'created by autotest');
    await openPlatformView(page, 'Roles');
    await searchAndWaitCard(page, 'roles', ROLE);

    // Rename via Properties...
    await openCardContextMenu(page, ROLE);
    await contextMenuItemByName(page, 'Properties...').click();
    await expect(page.locator(DIALOG_TITLE)).toContainText(/Properties/i);
    await dialogInput(page, 'Name').fill(ROLE_RENAMED);
    await dialogButton(page, 'OK').click();
    await page.waitForTimeout(1000);
    await openPlatformView(page, 'Roles');
    await searchAndWaitCard(page, 'roles', ROLE_RENAMED);

    // Delete.
    await deleteEntityViaContextMenu(page, ROLE_RENAMED);
    await openPlatformView(page, 'Roles');
    await searchAndWaitGone(page, 'roles', ROLE_RENAMED);
  });

  test('Roles-06 — search by name filters list and updates URL', async ({ page }) => {
    const before = await readGalleryCount(page);
    await searchGallery(page, 'roles', BUILTIN_ROLE_SEARCH);
    await expect(page).toHaveURL(new RegExp(`/roles\\?q=${BUILTIN_ROLE_SEARCH}`));
    expect((await readGalleryCount(page)).shown).toBeLessThan(before.total);
    await clearGallerySearch(page, 'roles');
    expect((await readGalleryCount(page)).shown).toBe(before.total);
  });

  // Roles-08 + Roles-10: context menu items and context-panel info panes.
  test('Roles-08/10 — role context menu items and Context Panel info panes', async ({ page }) => {
    // Act on a freshly created custom role so Delete is available and we don't touch built-ins.
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

  // Roles-11 + Roles-12 + Roles-13: assign a role to a user, toggle "Can assign", remove.
  test('Roles-11/12/13 — assign role to a user, Can-assign toggle, remove', async ({ page }) => {
    await createRole(page, ROLE, 'assignment test');
    try {
      await openPlatformView(page, 'Roles');
      await searchAndWaitCard(page, 'roles', ROLE);
      await selectCard(page, ROLE);

      // Assign a user and grant "Can assign".
      await openManageFromPane(page, 'Assigned to');
      await addMembershipBySearch(page, ASSIGNEE_USER);
      await setMemberRowToggle(page, ASSIGNEE_USER, true);
      await saveDialog(page);

      // Verify the assignment persisted.
      await selectCard(page, ROLE);
      await openManageFromPane(page, 'Assigned to');
      await expect(page.locator('.d4-dialog .membership-row', { hasText: ASSIGNEE_USER })).toBeVisible({ timeout: 5_000 });

      // Remove the assignment.
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

  // Roles-14 (edit a role's permissions via the Global Permissions / Permissions panes) is NOT
  // automated yet: the editing control in those panes rendered lazily during recon and its exact
  // mechanism is unconfirmed. To be added once the control is verified on a first manual pass.
});
