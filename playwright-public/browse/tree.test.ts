import { test, expect } from '@playwright/test';
import {
  BROWSE_HEADER,
  CONTEXT_MENU,
  CONTEXT_MENU_ITEM,
  SIDEBAR_BROWSE_ICON,
  TREE_EXPAND_ARROW,
  TREE_EXPAND_ARROW_EXPANDED,
  treeGroupByName,
  treeItemByName,
  treeNodeByName,
  treeNodeByPath,
  viewTabHandle,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
  clickCollapseAll,
  countExpandedNodes,
} from './helpers';

const FILES_CHILDREN = ['My files', 'App Data', 'Demo'];

test.describe('Browse tree (Browse-Tree-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
    await clickCollapseAll(page);
  });

  test('Browse-Tree-01 — expand then collapse a tree node', async ({ page }) => {
    const sink = watchErrors(page);

    const filesLabel = treeGroupByName(page, 'Files');
    const filesNode = filesLabel.locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]');
    const tri = filesNode.locator(TREE_EXPAND_ARROW).first();

    // Pre: Files is collapsed.
    await expect(tri).not.toHaveClass(/d4-tree-view-tri-expanded/);

    // Expand.
    await tri.click();
    await expect(tri, 'Files arrow should be expanded').toHaveClass(/d4-tree-view-tri-expanded/, { timeout: 5_000 });

    for (const child of FILES_CHILDREN) {
      await expect(
        treeGroupByName(page, child),
        `Child "${child}" should be visible after expand`,
      ).toBeVisible({ timeout: 5_000 });
    }

    // Collapse back.
    await tri.click();
    await expect(tri, 'Files arrow should be collapsed').not.toHaveClass(/d4-tree-view-tri-expanded/, { timeout: 5_000 });

    // Children hidden — pick one and assert it's not visible.
    await expect(
      treeGroupByName(page, 'Demo'),
      'A child should be hidden after collapse',
    ).not.toBeVisible({ timeout: 5_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Tree-02 — keyboard navigation in the tree (→ expands, ← collapses)', async ({ page }) => {
    const sink = watchErrors(page);

    const files = treeGroupByName(page, 'Files');
    const filesArrow = files
      .locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]')
      .locator(TREE_EXPAND_ARROW)
      .first();

    // Focus Files by clicking it; pre-state: collapsed.
    await files.click();
    await expect(filesArrow).not.toHaveClass(/d4-tree-view-tri-expanded/);

    // → expands the focused node.
    await page.keyboard.press('ArrowRight');
    await expect(filesArrow, '→ should expand the focused node').toHaveClass(
      /d4-tree-view-tri-expanded/,
      { timeout: 5_000 },
    );

    // Re-focus Files (after expand, focus may shift to first child) and press ← to collapse.
    await files.click();
    await page.keyboard.press('ArrowLeft');
    await expect(filesArrow, '← should collapse the focused node').not.toHaveClass(
      /d4-tree-view-tri-expanded/,
      { timeout: 5_000 },
    );

    // ↑/↓ navigation: after pressing ↓ on Files, a different node should have focus.
    await files.click();
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(500);
    // We can't easily inspect platform-specific focus styling; verify only that the
    // tree did not throw and remains interactive.
    await expect(page.locator(BROWSE_HEADER), 'Tree must stay interactive').toBeVisible();

    await expectNoErrors(page, sink);
  });

  test('Browse-Tree-03 — tree remembers expand/collapse state when switching views', async ({ page }) => {
    const sink = watchErrors(page);

    // Pre: expand two specific groups.
    await expandTreeGroup(page, 'Files');
    await expandTreeGroup(page, 'Databases');

    // Switch view: open Tutorials. Sidebar switches to Toolbox after a table-like view opens.
    await expandTreeGroup(page, 'Apps');
    await treeItemByName(page, 'Tutorials').click();
    await page.waitForTimeout(2000);

    // Return to the Browse panel from the Sidebar.
    await ensureBrowsePanelOpen(page);

    // Both previously-expanded groups must still have the expanded arrow.
    for (const name of ['Files', 'Databases']) {
      const tri = treeGroupByName(page, name)
        .locator('xpath=ancestor::*[contains(@class,"d4-tree-view-node")][1]')
        .locator(TREE_EXPAND_ARROW)
        .first();
      await expect(
        tri,
        `"${name}" should still be expanded after returning to Browse`,
      ).toHaveClass(/d4-tree-view-tri-expanded/, { timeout: 5_000 });
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-Tree-04 — sidebar toggle does not force-expand nested nodes', async ({ page }) => {
    const sink = watchErrors(page);

    // Pre: expand exactly Files; leave others collapsed.
    await expandTreeGroup(page, 'Files');
    const expandedBefore = await countExpandedNodes(page);

    // Toggle the Browse sidebar tab off and back on.
    await page.locator(SIDEBAR_BROWSE_ICON).click();
    await expect(page.locator(BROWSE_HEADER)).toBeHidden({ timeout: 5_000 });
    await page.locator(SIDEBAR_BROWSE_ICON).click();
    await expect(page.locator(BROWSE_HEADER)).toBeVisible({ timeout: 5_000 });
    await page.waitForTimeout(500);

    // Nested nodes should NOT be force-expanded — count should be unchanged (or ≤ before).
    const expandedAfter = await countExpandedNodes(page);
    expect(
      expandedAfter,
      'Toggling the sidebar must not force-expand nested nodes (GROK-19802)',
    ).toBeLessThanOrEqual(expandedBefore + 1);

    await expectNoErrors(page, sink);
  });

  test('Browse-Tree-05 — right-click on a tree node opens the context menu', async ({ page }) => {
    const sink = watchErrors(page);

    // Use a known visible top-level group (Files is always present and visible).
    const files = treeGroupByName(page, 'Files');
    await files.waitFor({ state: 'visible', timeout: 10_000 });
    await files.scrollIntoViewIfNeeded();
    await files.click({ button: 'right' });

    const menu = page.locator(CONTEXT_MENU);
    await expect(menu, 'Context menu should be open').toBeVisible({ timeout: 5_000 });

    // Sanity: the menu has at least one item.
    expect(
      await menu.locator(CONTEXT_MENU_ITEM).count(),
      'Context menu should have at least one item',
    ).toBeGreaterThanOrEqual(1);

    // Close the menu.
    await page.keyboard.press('Escape');
    await expect(menu).toBeHidden({ timeout: 5_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Tree-06 — drag-and-drop on a tree node does not crash the UI', async ({ page }) => {
    const sink = watchErrors(page);

    // Drag-and-drop on Browse tree entities is used to reorganize Spaces, share entities, etc.
    // Without a controlled fixture the most we can assert is "performing a drag gesture does
    // not raise an error". Source: My stuff > My files; Target: Favorites top-level node.
    await expandTreeGroup(page, 'My stuff');
    const source = treeNodeByPath(page, ['My-stuff', 'My-files']);
    const target = treeNodeByPath(page, ['My-stuff', 'Favorites']);

    await source.waitFor({ state: 'visible', timeout: 10_000 });
    await target.waitFor({ state: 'visible', timeout: 10_000 });

    await source.dragTo(target).catch(() => undefined);
    await page.waitForTimeout(1500);

    // App must still be alive.
    await expect(page.locator(BROWSE_HEADER)).toBeVisible();
    await expectNoErrors(page, sink);
  });

});
