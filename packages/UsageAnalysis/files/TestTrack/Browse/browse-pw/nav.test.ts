import { test, expect } from '@playwright/test';
import {
  BROWSE_HEADER,
  BROWSE_HEADER_HOME,
  BROWSE_HEADER_IMPORT_FILE,
  BROWSE_HEADER_IMPORT_TEXT,
  BROWSE_HEADER_REFRESH,
  BROWSE_HEADER_COLLAPSE_ALL,
  BROWSE_HEADER_LOCATE,
  BROWSE_HEADER_CLOSE_PANEL,
  HOME_GLOBAL_SEARCH_INPUT,
  SIDEBAR_BROWSE_ICON,
  TREE_EXPAND_ARROW_EXPANDED,
  treeGroupByName,
  treeItemByName,
  treeNodeByName,
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

const REQUIRED_TOP_LEVEL = [
  'My stuff',
  'Spaces',
  'Apps',
  'Files',
  'Dashboards',
  'Databases',
  'Platform',
];

test.describe('Browse navigation (Browse-Nav-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Nav-01 — Browse panel opens and contains required top-level sections', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page.locator(BROWSE_HEADER), 'Browse header should be visible').toBeVisible();

    for (const name of REQUIRED_TOP_LEVEL) {
      const isGroup = name !== 'Dashboards';
      const locator = isGroup ? treeGroupByName(page, name) : treeItemByName(page, name);
      await expect(locator, `Top-level section "${name}" should be visible`).toBeVisible({ timeout: 10_000 });
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-02 — toggle Browse panel by clicking the Sidebar icon', async ({ page }) => {
    const sink = watchErrors(page);
    const header = page.locator(BROWSE_HEADER);
    const sidebarBrowse = page.locator(SIDEBAR_BROWSE_ICON);

    await expect(header).toBeVisible();
    await sidebarBrowse.click();
    await expect(header, 'Browse panel should hide after first click').toBeHidden({ timeout: 5_000 });
    await sidebarBrowse.click();
    await expect(header, 'Browse panel should reappear after second click').toBeVisible({ timeout: 5_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-03 — Home icon opens Home Page with global search', async ({ page }) => {
    const sink = watchErrors(page);

    await page.locator(BROWSE_HEADER_HOME).click();
    // The global search input is unique to Home Page — its visibility is the source of truth.
    await expect(page.locator(HOME_GLOBAL_SEARCH_INPUT), 'Global search should be present').toBeVisible({ timeout: 10_000 });
    // The sidebar Home handle should be marked as selected.
    await expect(viewTabHandle(page, 'Home'), 'Home handle should be selected').toHaveClass(/tab-handle-selected/);

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-04 — Collapse all collapses the tree', async ({ page }) => {
    const sink = watchErrors(page);

    // Pre-state: expand a couple of nodes so we have something to collapse.
    await expandTreeGroup(page, 'Files');
    await expandTreeGroup(page, 'Databases');
    const before = await countExpandedNodes(page);
    expect(before, 'Pre-condition: at least 2 nodes expanded').toBeGreaterThanOrEqual(2);

    await clickCollapseAll(page);
    const after = await countExpandedNodes(page);
    expect(after, 'Number of expanded nodes should drop sharply after Collapse all')
      .toBeLessThan(before);
    // Some root-level Browse panel containers may keep the "expanded" marker;
    // accept ≤ 1 as effective "everything collapsed".
    expect(after, 'No regular tree nodes should remain expanded').toBeLessThanOrEqual(1);

    // Required sections still visible after collapse.
    for (const name of REQUIRED_TOP_LEVEL) {
      const isGroup = name !== 'Dashboards';
      const locator = isGroup ? treeGroupByName(page, name) : treeItemByName(page, name);
      await expect(locator).toBeVisible();
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-05 — Locate current object highlights the node in the tree', async ({ page }) => {
    const sink = watchErrors(page);

    // Pre: open Tutorials app via UI (it's a top-level item under Apps).
    await expandTreeGroup(page, 'Apps');
    const tutorials = treeItemByName(page, 'Tutorials');
    await expect(tutorials, 'Tutorials node should be visible under Apps').toBeVisible({ timeout: 10_000 });
    await tutorials.click();
    await page.waitForTimeout(2000);

    // Toolbox may have replaced Browse — reopen and collapse for a clean Locate.
    await ensureBrowsePanelOpen(page);
    await clickCollapseAll(page);

    // Confirm Tutorials is NOT visible in the collapsed tree.
    await expect(
      treeItemByName(page, 'Tutorials'),
      'Tutorials must be hidden before Locate (tree is collapsed)',
    ).not.toBeVisible({ timeout: 3_000 });

    // Now click Locate — the tree should auto-expand to reveal Tutorials.
    await page.locator(BROWSE_HEADER_LOCATE).click();
    await page.waitForTimeout(1500);

    await expect(
      treeItemByName(page, 'Tutorials'),
      'Tutorials should be revealed in the tree after Locate',
    ).toBeVisible({ timeout: 10_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-06 — Refresh tree re-reads tree without page reload', async ({ page }) => {
    const sink = watchErrors(page);

    // Expand a couple of sections so we can verify the tree state survives.
    await expandTreeGroup(page, 'Files');
    await expandTreeGroup(page, 'Databases');
    const expandedBefore = await countExpandedNodes(page);

    const initialUrl = page.url();
    await page.locator(BROWSE_HEADER_REFRESH).click();
    await page.waitForTimeout(2500);

    // URL and session unchanged.
    expect(page.url(), 'URL must not change after Refresh tree').toBe(initialUrl);

    // Required sections still present.
    for (const name of REQUIRED_TOP_LEVEL) {
      const isGroup = name !== 'Dashboards';
      const locator = isGroup ? treeGroupByName(page, name) : treeItemByName(page, name);
      await expect(locator, `"${name}" must still be present after Refresh`).toBeVisible();
    }

    // Expanded state preserved (or close to it — within ±1 to allow for minor relayout).
    const expandedAfter = await countExpandedNodes(page);
    expect(Math.abs(expandedAfter - expandedBefore),
      'Expanded-state should be preserved across Refresh tree').toBeLessThanOrEqual(1);

    await expectNoErrors(page, sink);

    // NOTE: verifying that an externally-created object appears after Refresh
    // is tracked separately (requires a stable API-creation path on dev).
  });

  test('Browse-Nav-07 — Import file opens system file dialog', async ({ page }) => {
    const sink = watchErrors(page);

    // We can verify the trigger by listening for the filechooser event.
    const fileChooserPromise = page.waitForEvent('filechooser', { timeout: 5_000 });
    await page.locator(BROWSE_HEADER_IMPORT_FILE).click();
    const chooser = await fileChooserPromise;
    expect(chooser, 'File chooser dialog should appear').toBeTruthy();

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-08 — Import text opens a view with a text input', async ({ page }) => {
    const sink = watchErrors(page);

    await page.locator(BROWSE_HEADER_IMPORT_TEXT).click();
    await page.waitForTimeout(1500);

    // The "Import text" view opens; look for the editor textarea.
    const textArea = page.locator('textarea, [contenteditable="true"]').first();
    await expect(textArea, 'A text input/area should be available').toBeVisible({ timeout: 10_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Nav-09 — every icon in the Browse header is clickable without errors', async ({ page }) => {
    const sink = watchErrors(page);

    // (icon-name, optional verification before returning to baseline)
    const icons: Array<{ name: string; selector: string; afterClick?: () => Promise<void> }> = [
      { name: 'Home', selector: BROWSE_HEADER_HOME },
      {
        name: 'Refresh',
        selector: BROWSE_HEADER_REFRESH,
        afterClick: async () => {
          // Refresh doesn't open a view; just give the tree a moment.
          await page.waitForTimeout(800);
        },
      },
      {
        name: 'Collapse all',
        selector: BROWSE_HEADER_COLLAPSE_ALL,
        afterClick: async () => {
          await page.waitForTimeout(500);
        },
      },
      {
        name: 'Locate',
        selector: BROWSE_HEADER_LOCATE,
        afterClick: async () => {
          await page.waitForTimeout(500);
        },
      },
      {
        name: 'Import text',
        selector: BROWSE_HEADER_IMPORT_TEXT,
        afterClick: async () => {
          // Close the view it opens to keep state clean for the next icon.
          await page.waitForTimeout(1000);
          const close = page.locator(`[aria-label="Close view"]`).first();
          if (await close.isVisible().catch(() => false)) await close.click();
          await page.waitForTimeout(500);
        },
      },
    ];

    for (const icon of icons) {
      const locator = page.locator(icon.selector).first();
      await expect(locator, `Icon "${icon.name}" should be visible`).toBeVisible();
      await locator.click();
      if (icon.afterClick) await icon.afterClick();
      await expectNoErrors(page, sink);
    }

    // Import file is special — it opens a system dialog; close it without choosing.
    const fileChooserPromise = page.waitForEvent('filechooser', { timeout: 5_000 });
    await page.locator(BROWSE_HEADER_IMPORT_FILE).click();
    const chooser = await fileChooserPromise;
    await chooser.setFiles([]).catch(() => undefined); // dismiss
    await expectNoErrors(page, sink);

    // The close-panel icon hides the Browse panel — verify, then reopen for cleanup.
    if (await page.locator(BROWSE_HEADER_CLOSE_PANEL).isVisible().catch(() => false)) {
      await page.locator(BROWSE_HEADER_CLOSE_PANEL).click();
      await expect(page.locator(BROWSE_HEADER), 'Browse panel should be hidden').toBeHidden({ timeout: 5_000 });
      await ensureBrowsePanelOpen(page);
      await expectNoErrors(page, sink);
    }
  });
});
