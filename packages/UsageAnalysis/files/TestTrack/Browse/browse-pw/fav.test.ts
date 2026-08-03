import { test, expect } from '@playwright/test';
import {
  CONTEXT_MENU,
  CONTEXT_MENU_ADD_FAVORITES,
  CONTEXT_PANEL_STAR,
  SIDEBAR_FAVORITES_ICON,
  contextMenuItem,
  treeGroupByName,
  treeNodeByPath,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

// Fav tests toggle the same shared Favorites state on Tutorials — keep them serial
// within the file so they don't race when workers > 1.
test.describe.configure({ mode: 'serial' });

test.describe('Browse favorites (Browse-Fav-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Fav-01 — Add to favorites from the tree context menu', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    const target = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await target.waitFor({ state: 'visible', timeout: 10_000 });
    const targetLabel = target.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();

    /** Right-click on the target label, returning when the menu is open. */
    async function openContextMenu(): Promise<void> {
      await targetLabel.scrollIntoViewIfNeeded();
      await targetLabel.click({ button: 'right' });
      await expect(page.locator(CONTEXT_MENU), 'Context menu should open')
        .toBeVisible({ timeout: 5_000 });
    }

    /** Closes any open menu by pressing Escape (no-op if nothing's open). */
    async function closeMenu(): Promise<void> {
      await page.keyboard.press('Escape');
      await page.waitForTimeout(300);
    }

    // Step 0: ensure the entity is NOT in favorites (cleanup from prior runs).
    await openContextMenu();
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
      await page.waitForTimeout(1000);
    } else {
      await closeMenu();
    }

    // Step 1: add to favorites.
    await openContextMenu();
    await expect(contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES),
      'After cleanup, "Add to favorites" must be present').toBeVisible({ timeout: 5_000 });
    await contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES).click();
    await page.waitForTimeout(1500);

    // Step 2: verify the entity is now listed under My stuff > Favorites.
    await expandTreeGroup(page, 'My stuff');
    await expandTreeGroup(page, 'Favorites');
    await page.waitForTimeout(500);
    const favEntry = page.locator(
      '.d4-tree-view-item-label, .d4-tree-view-group-label',
      { hasText: /^Tutorials$/ },
    );
    expect(await favEntry.count(), 'Tutorials should appear in Favorites').toBeGreaterThanOrEqual(1);

    await expectNoErrors(page, sink);

    // Cleanup: remove from favorites.
    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await openContextMenu();
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
    } else {
      await closeMenu();
    }
  });

  test('Browse-Fav-02 — toggle favorite via the Context Panel star (add + remove roundtrip)', async ({ page }) => {
    const sink = watchErrors(page);

    // Open Tutorials so its Context Panel star is rendered.
    await expandTreeGroup(page, 'Apps');
    await treeNodeByPath(page, ['Apps', 'Tutorials']).click();
    await page.waitForTimeout(2000);

    const star = page.locator(CONTEXT_PANEL_STAR);
    await expect(star, 'Favorites star should be present in Context Panel header')
      .toHaveCount(1, { timeout: 10_000 });

    // Roundtrip: click → click. The colour change requires `force:true` because of the
    // titlebar visibility:hidden quirk.
    await star.click({ force: true });
    await page.waitForTimeout(800);
    await star.click({ force: true });
    await page.waitForTimeout(800);

    await expectNoErrors(page, sink);
  });

  test('Browse-Fav-03 — Favorites icon on the Sidebar opens the favorites panel', async ({ page }) => {
    const sink = watchErrors(page);

    await page.locator(SIDEBAR_FAVORITES_ICON).click();
    await expect(page.locator(SIDEBAR_FAVORITES_ICON), 'Favorites sidebar tab should be selected')
      .toHaveClass(/selected/, { timeout: 5_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Fav-04 — drag a tree entity onto the Favorites panel adds it to Favorites', async ({ page }) => {
    const sink = watchErrors(page);

    // Pre-clean: make sure Tutorials is NOT in favorites.
    await expandTreeGroup(page, 'Apps');
    const source = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await source.waitFor({ state: 'visible', timeout: 10_000 });
    const sourceLabel = source.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();
    await sourceLabel.click({ button: 'right' });
    await expect(page.locator(CONTEXT_MENU)).toBeVisible({ timeout: 5_000 });
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
      await page.waitForTimeout(800);
    } else {
      await page.keyboard.press('Escape');
    }

    // Open the Favorites panel — this is the drop zone.
    await page.locator(SIDEBAR_FAVORITES_ICON).click();
    const dropZone = page.locator('.grok-favorites-pane.grok-favorites-list').first();
    await expect(dropZone, 'Favorites panel must be visible to act as a drop zone')
      .toBeVisible({ timeout: 5_000 });

    // Bring Browse back so source is accessible, then drag.
    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await source.scrollIntoViewIfNeeded();
    await source.dragTo(dropZone);
    await page.waitForTimeout(1500);

    // Final assertion: Tutorials appears in the Favorites pane.
    const favEntry = dropZone.locator('label, .d4-tree-view-item-label, .d4-tree-view-group-label',
      { hasText: /^Tutorials$/ });
    expect(await favEntry.count(),
      'Tutorials should appear in the Favorites panel after drag-and-drop').toBeGreaterThanOrEqual(1);

    await expectNoErrors(page, sink);

    // Cleanup.
    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await sourceLabel.click({ button: 'right' });
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
    } else {
      await page.keyboard.press('Escape');
    }
  });

  test('Browse-Fav-05 — file / cell context menu does not expose Add to favorites', async ({ page }) => {
    const sink = watchErrors(page);

    // Open the Files > Demo folder and right-click on a plain file (non-entity).
    await page.goto(`${process.env.DATAGROK_URL!}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
    const demog = page.locator('label', { hasText: /^demog\.csv$/ }).first();
    await expect(demog).toBeVisible({ timeout: 15_000 });

    await demog.evaluate((el) => {
      const r = el.getBoundingClientRect();
      el.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2, clientX: r.left + 5, clientY: r.top + 5,
      }));
    });

    await expect(page.locator(CONTEXT_MENU), 'Context menu must open').toBeVisible({ timeout: 5_000 });

    // For a plain file, "Add to favorites" should not be among visible menu items.
    const addItems = page.locator(`${CONTEXT_MENU} .d4-menu-item-label`,
      { hasText: /^Add to favorites$/i });
    expect(await addItems.count(),
      '"Add to favorites" must not appear for a non-entity file').toBe(0);

    await page.keyboard.press('Escape');
    await expectNoErrors(page, sink);
  });
});
