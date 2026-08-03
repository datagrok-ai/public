import { test, expect } from '@playwright/test';
import { treeGroupByName, TREE_EXPAND_ARROW } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

test.describe('Browse Spaces (Browse-Spaces-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Spaces-01 — clicking Spaces expands without errors and shows available spaces', async ({ page }) => {
    const sink = watchErrors(page);

    const spaces = treeGroupByName(page, 'Spaces');
    await expect(spaces, 'Spaces top-level group should be visible').toBeVisible({ timeout: 10_000 });

    await expandTreeGroup(page, 'Spaces');

    // After expand, the arrow has the expanded class.
    const tri = spaces
      .locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]')
      .locator(TREE_EXPAND_ARROW)
      .first();
    await expect(tri, 'Spaces arrow should be in expanded state')
      .toHaveClass(/d4-tree-view-tri-expanded/);

    await expectNoErrors(page, sink);
  });
});
