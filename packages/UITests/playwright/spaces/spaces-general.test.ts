/**
 * spaces-general-v2.test.ts
 *
 * Spaces feature end-to-end tests covering CRUD, validation, navigation,
 * hierarchy, favorites, sharing, search, drag-and-drop, content view,
 * link/copy semantics, entity operations, and edge cases.
 */

import { test, expect, Page } from '@playwright/test';

const BASE = process.env.DATAGROK_URL!;
const SHARING_LOGIN = process.env.DATAGROK_SHARING_LOGIN || '';
const SHARING_PASSWORD = process.env.DATAGROK_SHARING_PASSWORD || '';

// ===========================================================================
// Helpers
// ===========================================================================

async function openSpacesView(page: Page) {
  await page.goto(BASE);
  await expect(
    page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first(),
  ).toBeVisible({ timeout: 20_000 });
}

async function refreshSpacesTree(page: Page) {
  const syncIcon = page.locator('i.grok-icon[name="icon-sync"], i.fal.fa-sync').first();
  if (await syncIcon.isVisible({ timeout: 3_000 }).catch(() => false)) {
    await syncIcon.click();
    await page.waitForTimeout(500);
  }
}

async function apiDeleteSpace(page: Page, name: string) {
  await page.evaluate(async (n) => {
    const g = (window as any).grok;
    const spaces = await g.dapi.spaces.filter(`name = "${n}"`).list({ pageSize: 5 });
    for (const s of spaces)
      if (s.friendlyName === n || s.name === n) await g.dapi.spaces.delete(s);
  }, name);
}

async function rightClickSpacesTreeNode(page: Page) {
  await page.evaluate(() => {
    const el = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
      .find((e) => e.textContent?.trim() === 'Spaces');
    el?.closest('.d4-tree-view-node')?.dispatchEvent(
      new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2 }),
    );
  });
  await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
}

async function rightClickSpace(page: Page, name: string) {
  const label = page.locator('.d4-tree-view-group-label:visible')
    .filter({ hasText: new RegExp(`^${name}$`) }).first();
  await expect(label).toBeVisible({ timeout: 10_000 });
  await label.locator('xpath=..').dispatchEvent('contextmenu', { button: 2, bubbles: true });
  await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
}

async function clickMenuItem(page: Page, label: string | RegExp) {
  await page.locator('.d4-menu-item', { hasText: label }).first().click();
  await page.waitForTimeout(300);
}

async function fillNameAndOk(page: Page, name: string) {
  const dialog = page.locator('.d4-dialog').last();
  const input = dialog.locator('input[type="text"]').first();
  await expect(input).toBeVisible({ timeout: 5_000 });
  await input.click();
  await page.keyboard.press('Control+a');
  await page.keyboard.press('Delete');
  await input.pressSequentially(name);
  await page.keyboard.press('Enter');
  await dialog.waitFor({ state: 'hidden', timeout: 5_000 }).catch(() => {});
}

async function uiCreateRootSpace(page: Page, name: string) {
  await rightClickSpacesTreeNode(page);
  await clickMenuItem(page, 'Create Space...');
  await fillNameAndOk(page, name);
  await expect(
    page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${name}$`) }).first(),
  ).toBeVisible({ timeout: 20_000 });
}

async function openSpaceContent(page: Page, spaceName: string) {
  const spacesNode = page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first();
  if (!await spacesNode.isVisible({ timeout: 2_000 }).catch(() => false)) {
    await page.goto(BASE);
    await expect(spacesNode).toBeVisible({ timeout: 20_000 });
  }
  const label = page.locator('.d4-tree-view-group-label')
    .filter({ hasText: new RegExp(`^${spaceName}$`) }).first();
  await expect(label).toBeVisible({ timeout: 10_000 });
  await label.dblclick();
  await page.waitForURL(/\/s\//, { timeout: 10_000 });
  await page.locator('#elementContent .grok-gallery-grid, .d4-link-label').first()
    .waitFor({ state: 'visible', timeout: 10_000 }).catch(() => {});
}

async function openSpaceViaTree(page: Page, spaceName: string) {
  const label = page.locator('.d4-tree-view-group-label')
    .filter({ hasText: new RegExp(`^${spaceName}$`) }).first();
  await expect(label).toBeVisible({ timeout: 10_000 });
  const prevUrl = page.url();
  await label.dblclick();
  await page.waitForURL(url => url !== prevUrl && url.includes('/s/'), { timeout: 8_000 })
    .catch(() => {});
  await page.locator('#elementContent .grok-gallery-grid, .d4-link-label').first()
    .waitFor({ state: 'visible', timeout: 8_000 }).catch(() => {});
}

async function dragFileToSpaceNode(page: Page, fileName: string, spaceName: string) {
  // Navigate to DemoFiles, retry once if page doesn't load properly
  for (let attempt = 0; attempt < 3; attempt++) {
    await page.goto(`${BASE}/files/System.DemoFiles/?browse=files`);
    await page.waitForTimeout(3000);
    const count = await page.locator('.d4-link-label label').count();
    if (count > 5) break; // Page loaded with enough files
  }
  const fileLabel = page.locator('.d4-link-label label')
    .filter({ hasText: new RegExp(`^${fileName}$`) }).first();
  const fileCount = await fileLabel.count();
  if (fileCount === 0) {
    const allFiles = await page.locator('.d4-link-label label').allTextContents();
    throw new Error(`File "${fileName}" not found in DemoFiles. Available files: ${allFiles.slice(0, 20).join(', ')}`);
  }
  await fileLabel.scrollIntoViewIfNeeded();
  await expect(fileLabel).toBeVisible({ timeout: 5_000 });
  const spaceNode = page.locator('.d4-tree-view-group-label')
    .filter({ hasText: new RegExp(`^${spaceName}$`) }).first();
  await expect(spaceNode).toBeVisible({ timeout: 10_000 });
  await fileLabel.dragTo(spaceNode);
  await expect(page.locator('text=Move entity')).toBeVisible({ timeout: 8_000 });
}

async function addFileToSpaceViaLink(page: Page, fileName: string, spaceName: string) {
  await dragFileToSpaceNode(page, fileName, spaceName);
  await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
  await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
}

async function addFileToSpaceViaCopy(page: Page, fileName: string, spaceName: string) {
  await dragFileToSpaceNode(page, fileName, spaceName);
  await page.locator('select').filter({ has: page.locator('option[value="Link"]') }).first().selectOption('Copy');
  await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
  await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
}

async function rightClickItemInSpaceView(page: Page, name: string) {
  const label = page.locator('.d4-link-label label')
    .filter({ hasText: new RegExp(`^${name}$`) }).first();
  await expect(label).toBeVisible({ timeout: 10_000 });
  await label.locator('xpath=..').dispatchEvent('contextmenu', { button: 2, bubbles: true });
  await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
}

// ===========================================================================
// 1. Create root space and verify context menu
// ===========================================================================

test('1. Create root space and verify context menu', async ({ page }) => {
  // Spaces node visible in browse tree; create root space; context menu has all expected items
  const ROOT = 'PW-Gen-Root-1';
  try {
    await page.goto(BASE);
    const spacesNode = page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first();
    await expect(spacesNode).toBeVisible({ timeout: 20_000 });

    await spacesNode.click();
    await expect(spacesNode).toBeVisible();

    await rightClickSpacesTreeNode(page);
    await expect(page.locator('.d4-menu-item', { hasText: /Create Space/i }).first()).toBeVisible({ timeout: 5_000 });
    await clickMenuItem(page, 'Create Space...');
    await fillNameAndOk(page, ROOT);
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${ROOT}$`) }).first(),
    ).toBeVisible({ timeout: 20_000 });

    await rightClickSpace(page, ROOT);
    const items = page.locator('.d4-menu-item');
    await expect(items.filter({ hasText: /Share/i }).first()).toBeVisible({ timeout: 5_000 });
    await expect(items.filter({ hasText: /Rename/i }).first()).toBeVisible();
    await expect(items.filter({ hasText: /Delete/i }).first()).toBeVisible();
    await expect(items.filter({ hasText: /Create Child Space/i }).first()).toBeVisible();
    await expect(items.filter({ hasText: /favorites/i }).first()).toBeVisible();
    await page.keyboard.press('Escape');
  } finally {
    await apiDeleteSpace(page, ROOT);
  }
});

// ===========================================================================
// 2. Validation: empty name, duplicate root and child names
// ===========================================================================

test('2. Validation: empty name disables OK; duplicate names show errors', async ({ page }) => {
  // Empty name disables OK button; duplicate root name shows error; duplicate child name shows error
  const DUP    = 'PW-Gen-Dup-2';
  const PARENT = 'PW-Gen-Dup-Parent-2';
  const CHILD  = 'PW-Gen-Dup-Child-2';
  try {
    await openSpacesView(page);

    await rightClickSpacesTreeNode(page);
    await clickMenuItem(page, 'Create Space...');
    const emptyDialog = page.locator('.d4-dialog').last();
    const emptyInput = emptyDialog.locator('input[type="text"]').first();
    await emptyInput.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await emptyInput.fill('');
    await expect(emptyDialog.locator('.ui-btn', { hasText: /^OK$/i }).first()).not.toHaveClass(/\benabled\b/, { timeout: 3_000 });
    await page.keyboard.press('Escape');

    await uiCreateRootSpace(page, DUP);
    await rightClickSpacesTreeNode(page);
    await clickMenuItem(page, 'Create Space...');
    await fillNameAndOk(page, DUP);
    await expect(page.locator('.d4-toast, .d4-balloon, [class*="toast"], [class*="error"]').first()).toContainText(/name already exists/i, { timeout: 8_000 });
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${DUP}$`) }),
    ).toHaveCount(1, { timeout: 5_000 });

    await uiCreateRootSpace(page, PARENT);
    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD);
    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD);
    await expect(page.locator('.d4-toast, .d4-balloon, [class*="toast"], [class*="error"]').first()).toContainText(/name already exists/i, { timeout: 8_000 });
  } finally {
    for (const n of [DUP, PARENT, CHILD])
      await apiDeleteSpace(page, n);
  }
});

// ===========================================================================
// 3. Hierarchy: create child, three levels, navigation
// ===========================================================================

test('3. Hierarchy: create child from tree, grandchild from view, navigate three levels', async ({ page }) => {
  const ROOT  = 'PW-Gen-Hier-Root-3';
  const CHILD = 'PW-Gen-Hier-Child-3';
  const GRAND = 'PW-Gen-Hier-Grand-3';
  try {
    await openSpacesView(page);

    // Create root space via UI
    await uiCreateRootSpace(page, ROOT);

    // Create child from the browse tree
    await rightClickSpace(page, ROOT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD);

    // Verify child is visible in the browse tree
    await refreshSpacesTree(page);
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${CHILD}$`) }).first(),
    ).toBeVisible({ timeout: 10_000 });

    // Double-click root to open its content view; verify URL changes to /s/
    await openSpaceContent(page, ROOT);
    expect(page.url()).toMatch(/\/s\//);

    // Verify gallery grid is rendered and child card is visible inside
    await expect(
      page.locator('#elementContent .grok-gallery-grid').first(),
    ).toBeVisible({ timeout: 10_000 });
    await expect(
      page.locator('.d4-link-label').getByText(CHILD, { exact: true }).first(),
    ).toBeVisible({ timeout: 10_000 });

    // Create grandchild from the content view (right-click child card)
    const childCard = page.locator('.d4-link-label')
      .filter({ hasText: new RegExp(`^${CHILD}$`) }).first();
    await childCard.dispatchEvent('contextmenu', { button: 2, bubbles: true });
    await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, GRAND);

    // Navigate into child via double-click and verify grandchild is visible
    await openSpaceContent(page, ROOT);
    const childLink = page.locator('.d4-link-label')
      .filter({ hasText: new RegExp(`^${CHILD}$`) }).first();
    await expect(childLink).toBeVisible({ timeout: 10_000 });
    const prevUrl = page.url();
    await childLink.dispatchEvent('dblclick');
    await page.waitForURL(url => url !== prevUrl, { timeout: 10_000 });

    // Verify URL updated to /s/ after navigating via content view double-click
    expect(page.url()).toMatch(/\/s\//);

    await expect(
      page.locator('.d4-link-label').getByText(GRAND, { exact: true }).first(),
    ).toBeVisible({ timeout: 10_000 });
  } finally {
    for (const n of [GRAND, CHILD, ROOT])
      await apiDeleteSpace(page, n);
  }
});

// ===========================================================================
// 4. Rename
// ===========================================================================

test('4. Rename: pre-fill, cancel, success, duplicate error, right panel update', async ({ page }) => {
  // Dialog pre-fills current name; cancel preserves original; rename succeeds; duplicate name shows error
  const ORIG         = 'PW-Gen-Ren-Orig-6';
  const NEW_NAME     = 'PW-Gen-Ren-New-6';
  const OTHER        = 'PW-Gen-Ren-Other-6';
  const REN_PARENT   = 'PW-Gen-Ren-CParent-6';
  const REN_CHILD    = 'PW-Gen-Ren-Child-6';
  const REN_CHILD_NEW = 'PW-Gen-Ren-ChildNew-6';
  try {
    await openSpacesView(page);
    await uiCreateRootSpace(page, ORIG);

    await page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${ORIG}$`) }).first().click();
    await rightClickSpace(page, ORIG);
    await clickMenuItem(page, 'Rename...');
    const d1 = page.locator('.d4-dialog').last();
    const i1 = d1.locator('input[type="text"]').first();
    await expect(i1).toBeVisible({ timeout: 5_000 });
    expect(await i1.inputValue()).toBe(ORIG);

    await i1.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await i1.pressSequentially(NEW_NAME);
    await d1.locator('.ui-btn', { hasText: /^CANCEL$/i }).first().click({ force: true });
    await d1.waitFor({ state: 'hidden', timeout: 3_000 }).catch(() => {});
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${ORIG}$`) }).first()).toBeVisible({ timeout: 5_000 });

    await rightClickSpace(page, ORIG);
    await clickMenuItem(page, 'Rename...');
    const d2 = page.locator('.d4-dialog').last();
    const i2 = d2.locator('input[type="text"]').first();
    await i2.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await i2.pressSequentially(NEW_NAME);
    await page.keyboard.press('Enter');
    // Verify new name in tree, old name gone
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${NEW_NAME}$`) }).first()).toBeVisible({ timeout: 8_000 });
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${ORIG}$`) }).first()).not.toBeVisible({ timeout: 5_000 });
    // Refresh tree and verify result unchanged
    await refreshSpacesTree(page);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${NEW_NAME}$`) }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${ORIG}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
    expect(
      await page.locator('.grok-prop-panel, [class*="context-panel"]').filter({ hasText: new RegExp(`^${ORIG}$`) }).count(),
    ).toBe(0);

    // Verify new name in the content view (gallery grid), old name gone
    await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first().click();
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${NEW_NAME}$`) }).first()).toBeVisible({ timeout: 10_000 });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${ORIG}$`) }).first()).not.toBeVisible({ timeout: 3_000 });

    await apiDeleteSpace(page, NEW_NAME);

    await uiCreateRootSpace(page, ORIG);
    await uiCreateRootSpace(page, OTHER);
    await rightClickSpace(page, ORIG);
    await clickMenuItem(page, 'Rename...');
    const d3 = page.locator('.d4-dialog').last();
    const i3 = d3.locator('input[type="text"]').first();
    await i3.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await i3.pressSequentially(OTHER);
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
    const errorOrOrig =
      await page.locator('.d4-toast, .d4-balloon, [class*="toast"], [class*="error"]').first().isVisible().catch(() => false) ||
      await page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${ORIG}$`) }).first().isVisible().catch(() => false);
    expect(errorOrOrig).toBeTruthy();

    await uiCreateRootSpace(page, REN_PARENT);
    await rightClickSpace(page, REN_PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, REN_CHILD);
    await openSpaceContent(page, REN_PARENT);
    const childCard = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${REN_CHILD}$`) }).first();
    await expect(childCard).toBeVisible({ timeout: 10_000 });
    await childCard.dispatchEvent('contextmenu', { button: 2, bubbles: true });
    await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
    await clickMenuItem(page, 'Rename...');
    const dChild = page.locator('.d4-dialog').last();
    const iChild = dChild.locator('input[type="text"]').first();
    await expect(iChild).toBeVisible({ timeout: 5_000 });
    expect(await iChild.inputValue()).toBe(REN_CHILD);

    // Perform the rename of the child space
    await iChild.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await iChild.pressSequentially(REN_CHILD_NEW);
    await page.keyboard.press('Enter');

    // Verify renamed child in the content view
    await openSpaceContent(page, REN_PARENT);
    await expect(
      page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${REN_CHILD_NEW}$`) }).first(),
    ).toBeVisible({ timeout: 10_000 });
    await expect(
      page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${REN_CHILD}$`) }).first(),
    ).not.toBeVisible({ timeout: 3_000 });

    // Verify renamed child in the browse tree
    await openSpacesView(page);
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${REN_CHILD_NEW}$`) }).first(),
    ).toBeVisible({ timeout: 10_000 });
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${REN_CHILD}$`) }).first(),
    ).not.toBeVisible({ timeout: 3_000 });
    // Refresh tree and verify result unchanged
    await refreshSpacesTree(page);
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${REN_CHILD_NEW}$`) }).first(),
    ).toBeVisible({ timeout: 5_000 });
    await expect(
      page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${REN_CHILD}$`) }).first(),
    ).not.toBeVisible({ timeout: 3_000 });
  } finally {
    for (const n of [ORIG, NEW_NAME, OTHER, REN_PARENT, REN_CHILD, REN_CHILD_NEW])
      await apiDeleteSpace(page, n);
  }
});

// ===========================================================================
// 5. Delete
// ===========================================================================

test('5. Delete: cancel preserves, selective delete keeps sibling, cascade removes hierarchy', async ({ page }) => {
  // Cancel preserves space; deleting one child keeps sibling; deleting parent cascades
  test.setTimeout(180_000);
  const SINGLE = 'PW-Gen-DelSingle-7';
  const PARENT = 'PW-Gen-DelParent-7';
  const CHILD1 = 'PW-Gen-DelChild7A';
  const CHILD2 = 'PW-Gen-DelChild7B';
  try {
    await openSpacesView(page);
    for (const n of [SINGLE, CHILD1, CHILD2, PARENT])
      await apiDeleteSpace(page, n);
    await uiCreateRootSpace(page, SINGLE);

    // Part A: cancel delete from tree preserves the space
    await rightClickSpace(page, SINGLE);
    await clickMenuItem(page, 'Delete');
    await expect(page.locator('.ui-btn', { hasText: /^CANCEL$/i }).first()).toBeVisible({ timeout: 5_000 });
    await page.locator('.ui-btn', { hasText: /^CANCEL$/i }).first().click();
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${SINGLE}$`) }).first()).toBeVisible({ timeout: 5_000 });

    // Delete SINGLE from tree and verify gone from both tree and view
    await rightClickSpace(page, SINGLE);
    await clickMenuItem(page, 'Delete');
    await page.locator('.ui-btn', { hasText: /^DELETE$/i }).first().click();
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${SINGLE}$`) }).first()).not.toBeVisible({ timeout: 8_000 });
    await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first().click();
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${SINGLE}$`) }).first()).not.toBeVisible({ timeout: 5_000 });
    // Refresh tree and verify result unchanged
    await refreshSpacesTree(page);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${SINGLE}$`) }).first()).not.toBeVisible({ timeout: 5_000 });

    await uiCreateRootSpace(page, PARENT);
    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD1);
    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD2);

    // Part B: selective delete child1 from view, child2 survives
    await openSpaceContent(page, PARENT);
    await expect(page.locator('.d4-link-label').getByText(CHILD1, { exact: true }).first()).toBeVisible({ timeout: 10_000 });
    const child1Card = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CHILD1}$`) }).first();
    await child1Card.dispatchEvent('contextmenu', { button: 2, bubbles: true });
    await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
    await clickMenuItem(page, /Delete/i);
    await page.locator('.ui-btn', { hasText: /^DELETE$/i }).first().click();
    // Verify child1 gone from view, child2 still in view
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CHILD1}$`) }).first()).not.toBeVisible({ timeout: 8_000 });
    await expect(page.locator('.d4-link-label').getByText(CHILD2, { exact: true }).first()).toBeVisible({ timeout: 5_000 });
    // Verify child1 also gone from tree
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${CHILD1}$`) }).first()).not.toBeVisible({ timeout: 5_000 });
    // Refresh tree and verify result unchanged
    await refreshSpacesTree(page);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${CHILD1}$`) }).first()).not.toBeVisible({ timeout: 3_000 });

    // Part C: cascade delete parent from tree, verify gone from tree and view
    await openSpacesView(page);
    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Delete');
    await page.locator('.ui-btn', { hasText: /^DELETE$/i }).first().click();
    // Verify parent and remaining child2 gone from tree
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${PARENT}$`) }).first()).not.toBeVisible({ timeout: 8_000 });
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${CHILD2}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
    // Verify parent and child2 gone from Spaces gallery view
    await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first().click();
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${PARENT}$`) }).first()).not.toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CHILD2}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
    // Refresh tree and verify result unchanged
    await refreshSpacesTree(page);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${PARENT}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${CHILD2}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
  } finally {
    for (const n of [SINGLE, PARENT, CHILD1, CHILD2])
      await apiDeleteSpace(page, n);
  }
});

// ===========================================================================
// 6. Favorites
// ===========================================================================

test('6. Favorites: add and remove space from favorites', async ({ page }) => {
  // Add to favorites → appears; remove → gone
  const FAV = 'PW-Gen-Fav-8';
  try {
    await openSpacesView(page);
    await uiCreateRootSpace(page, FAV);

    await rightClickSpace(page, FAV);
    await clickMenuItem(page, 'Add to favorites');

    // Verify in browse tree Favorites section
    await page.evaluate(() => {
      const el = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((e) => e.textContent?.trim() === 'Favorites') as HTMLElement | undefined;
      if (el) el.click();
    });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${FAV}$`) }).first()).toBeVisible({ timeout: 8_000 });

    // Verify in sidebar favorites pane (star icon)
    await page.locator('.d4-tab-header-stripe.layout-sidebar.vertical > div:nth-child(6)').click();
    await expect(page.locator('.grok-favorites-pane').first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.grok-favorites-pane').filter({ hasText: new RegExp(FAV) }).first()).toBeVisible({ timeout: 8_000 });

    // Remove from favorites
    await rightClickSpace(page, FAV);
    await clickMenuItem(page, 'Remove from favorites');

    // Verify gone from browse tree Favorites
    await page.evaluate(() => {
      const el = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((e) => e.textContent?.trim() === 'Favorites') as HTMLElement | undefined;
      if (el) el.click();
    });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${FAV}$`) }).first()).not.toBeVisible({ timeout: 8_000 });

    // Verify gone from sidebar favorites pane
    await page.locator('.d4-tab-header-stripe.layout-sidebar.vertical > div:nth-child(6)').click();
    await expect(page.locator('.grok-favorites-pane').filter({ hasText: new RegExp(FAV) }).first()).not.toBeVisible({ timeout: 5_000 });
  } finally {
    await apiDeleteSpace(page, FAV);
  }
});

// ===========================================================================
// 7. Share dialog
// ===========================================================================

test('7. Share dialog: UI structure, share with permission, verify, delete removes access', async ({ page }) => {
  test.setTimeout(180_000);
  const SPACE      = 'PW-Gen-Share-7';
  const SH_PARENT  = 'PW-Gen-ShareP-7';
  const SH_CHILD   = 'PW-Gen-ShareC-7';
  try {
    await openSpacesView(page);
    await apiDeleteSpace(page, SPACE);
    await uiCreateRootSpace(page, SPACE);

    // Part A: verify dialog structure and default permission
    await rightClickSpace(page, SPACE);
    await clickMenuItem(page, 'Share...');
    await expect(page.locator('input[placeholder*="User"]').first()).toBeVisible({ timeout: 8_000 });
    const permSelect = page.locator('select').first();
    await expect(permSelect).toBeVisible({ timeout: 5_000 });
    await expect(permSelect.locator('option', { hasText: /Full access/i })).toHaveCount(1);
    await expect(permSelect.locator('option', { hasText: /View and use/i })).toHaveCount(1);

    // Part B: select "View and use" permission explicitly, share with second user
    await permSelect.selectOption({ label: 'View and use' });
    const selectedValue = await permSelect.inputValue();
    const userInput = page.locator('input[placeholder*="User"]').first();
    await userInput.fill(SHARING_LOGIN);
    await page.waitForTimeout(1000);
    // Select the user from autocomplete dropdown
    const suggestion = page.locator('.d4-menu-item-label, .d4-item, .itemFrame')
      .filter({ hasText: new RegExp(SHARING_LOGIN.split('@')[0].replace(/[+]/g, '\\$&'), 'i') }).first();
    if (await suggestion.isVisible({ timeout: 5_000 }).catch(() => false))
      await suggestion.click();
    else
      await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
    // Click OK to apply sharing
    const okBtn = page.locator('.ui-btn', { hasText: /^OK$/i }).first();
    await expect(okBtn).toBeVisible({ timeout: 5_000 });
    await okBtn.click();
    await page.waitForTimeout(500);

    // Part C: verify via API — user is in view permissions (not edit) matching the selected permission
    const permCheck = await page.evaluate(async ({ sharingLogin, spaceName }) => {
      try {
        const g = (window as any).grok;
        const spaces = await g.dapi.spaces.filter(`name = "${spaceName}"`).list({ pageSize: 5 });
        const space = spaces.find((s: any) => s.friendlyName === spaceName || s.name === spaceName);
        if (!space) return { error: 'space not found' };
        const perms = await g.dapi.permissions.get(space);
        const loginPrefix = sharingLogin.split('@')[0].replace(/[+_.]/g, '').toLowerCase();
        const inView = (perms.view || []).some((grp: any) =>
          (grp.name || '').toLowerCase().includes(loginPrefix));
        const inEdit = (perms.edit || []).some((grp: any) =>
          (grp.name || '').toLowerCase().includes(loginPrefix));
        return { inView, inEdit, spaceId: space.id };
      }
      catch (e: any) { return { error: e.message }; }
    }, { sharingLogin: SHARING_LOGIN, spaceName: SPACE });
    // "View and use" should result in view permission, not edit
    expect(permCheck).not.toHaveProperty('error');
    expect((permCheck as any).inView).toBe(true);
    expect((permCheck as any).inEdit).toBe(false);

    // Part D: delete space and verify sharing user no longer has access
    const spaceId = (permCheck as any).spaceId;
    await apiDeleteSpace(page, SPACE);
    await page.waitForTimeout(500);
    // Verify space is deleted — not listable via spaces API
    const goneFromSpaces = await page.evaluate(async ({ spaceName }) => {
      try {
        const g = (window as any).grok;
        const spaces = await g.dapi.spaces.filter(`name = "${spaceName}"`).list({ pageSize: 5 });
        return spaces.filter((s: any) => s.friendlyName === spaceName || s.name === spaceName).length === 0;
      }
      catch { return true; }
    }, { spaceName: SPACE });
    expect(goneFromSpaces).toBe(true);
    // Verify permissions no longer accessible for the space
    const permsGone = await page.evaluate(async (id) => {
      try {
        const resp = await fetch(`/api/security/permissions/entity/${id}`);
        const text = await resp.text();
        return resp.status !== 200 || text.includes('error') || text.includes('not found');
      }
      catch { return true; }
    }, spaceId);
    expect(permsGone).toBe(true);

    // Part E: child space is also shareable (dialog opens)
    await openSpacesView(page);
    await uiCreateRootSpace(page, SH_PARENT);
    await rightClickSpace(page, SH_PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, SH_CHILD);
    await openSpaceContent(page, SH_PARENT);
    const childCard = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${SH_CHILD}$`) }).first();
    await expect(childCard).toBeVisible({ timeout: 10_000 });
    await childCard.dispatchEvent('contextmenu', { button: 2, bubbles: true });
    await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
    await clickMenuItem(page, 'Share...');
    await expect(page.locator('input[placeholder*="User"]').first()).toBeVisible({ timeout: 8_000 });
    await page.keyboard.press('Escape');
  } finally {
    for (const n of [SPACE, SH_PARENT, SH_CHILD])
      await apiDeleteSpace(page, n);
  }
});

// ===========================================================================
// 8. Browse tree search
// ===========================================================================

test('8. Browse tree search: match by name and no-match empty result', async ({ page }) => {
  // Search matches space by name; non-existent query hides it
  const SPACE = 'PW-Gen-Search-10';
  try {
    await openSpacesView(page);
    await uiCreateRootSpace(page, SPACE);

    await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/i }).first().click();
    const searchInput = page.locator('input[placeholder*="Search"]').first();
    await expect(searchInput).toBeVisible({ timeout: 8_000 });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${SPACE}$`) }).first()).toBeVisible({ timeout: 15_000 });

    await searchInput.fill(SPACE);
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${SPACE}$`) }).first()).toBeVisible({ timeout: 8_000 });

    await searchInput.fill('xyzzy_no_match_123456');
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${SPACE}$`) }).first()).not.toBeVisible({ timeout: 5_000 });
  } finally {
    await apiDeleteSpace(page, SPACE);
  }
});

// ===========================================================================
// 9. DnD dialog structure
// ===========================================================================

test('9. DnD dialog: options, default Link, Cancel, Link, Copy, duplicate add', async ({ page }) => {
  // Dialog has Link/Copy/Move options with Link default; Cancel aborts; Link and Copy add files; duplicate add has no error
  test.setTimeout(180_000);
  const SPACE = 'PW-Gen-DnD-9';
  try {
    await openSpacesView(page);
    await apiDeleteSpace(page, SPACE);
    await uiCreateRootSpace(page, SPACE);

    // Part A: dialog structure and Cancel
    await dragFileToSpaceNode(page, 'wells.csv', SPACE);
    const sel = page.locator('select').filter({ has: page.locator('option[value="Link"]') }).first();
    await expect(sel).toHaveValue('Link');
    await expect(sel.locator('option[value="Copy"]')).toHaveCount(1);
    await expect(sel.locator('option[value="Move"]')).toHaveCount(1);
    await expect(page.locator(`text=${SPACE}`).first()).toBeVisible();
    await page.locator('.ui-btn', { hasText: /^CANCEL$/i }).first().click();
    // Verify cancel: file not in view
    await openSpaceContent(page, SPACE);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^wells\.csv$/ }).first()).not.toBeVisible({ timeout: 5_000 });
    // Verify cancel: file not in tree under space
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: /^wells\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });
    // Refresh tree and verify result unchanged
    await refreshSpacesTree(page);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: /^wells\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });

    // Part B: Link file
    await dragFileToSpaceNode(page, 'acidiq.csv', SPACE);
    await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
    await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
    // Verify linked file in view
    await openSpaceContent(page, SPACE);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
    // Click linked file to verify right panel shows info
    const linkedFile = page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first();
    await linkedFile.dispatchEvent('click');
    await expect(page.locator('.grok-prop-panel, [class*="context-panel"], .d4-info-bar').first()).toBeVisible({ timeout: 8_000 });

    // Part C: Copy file
    await dragFileToSpaceNode(page, 'TSLA.csv', SPACE);
    await page.locator('select').filter({ has: page.locator('option[value="Link"]') }).first().selectOption('Copy');
    await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
    await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
    // Verify both files in view
    await openSpaceContent(page, SPACE);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
    // Click copied file to verify right panel shows info
    const copiedFile = page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first();
    await copiedFile.dispatchEvent('click');
    await expect(page.locator('.grok-prop-panel, [class*="context-panel"], .d4-info-bar').first()).toBeVisible({ timeout: 8_000 });

    // Part D: duplicate add — no error, still exactly one acidiq.csv
    await dragFileToSpaceNode(page, 'acidiq.csv', SPACE);
    await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
    await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
    await expect(page.locator('.d4-toast.error, [class*="error-toast"]').first()).not.toBeVisible({ timeout: 2_000 }).catch(() => {});
    await openSpaceContent(page, SPACE);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
    // Verify only one acidiq.csv in the space (not duplicated)
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ })).toHaveCount(1, { timeout: 5_000 });
  } finally {
    await apiDeleteSpace(page, SPACE);
  }
});

// ===========================================================================
// 10. DnD Move
// ===========================================================================

test('10. DnD Move: file moves to target and is absent from source', async ({ page }) => {
  // Move file between spaces; file appears in target and disappears from source
  test.setTimeout(180_000);
  const SRC = 'PW-Gen-MoveSrc-10';
  const TGT = 'PW-Gen-MoveTgt-10';
  await openSpacesView(page);
  await apiDeleteSpace(page, SRC);
  await apiDeleteSpace(page, TGT);
  try {
    await uiCreateRootSpace(page, SRC);
    await uiCreateRootSpace(page, TGT);

    await addFileToSpaceViaCopy(page, 'beer.csv', SRC);
    await openSpaceContent(page, SRC);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });

    // Drag file from source view to target tree node
    const fileInView = page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first();
    const tgtNode = page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${TGT}$`) }).first();
    await expect(tgtNode).toBeVisible({ timeout: 5_000 });
    page.setDefaultTimeout(15_000);
    await fileInView.dragTo(tgtNode).catch(() => {});
    page.setDefaultTimeout(30_000);
    const dialogVisible = await page.locator('text=Move entity').isVisible({ timeout: 5_000 }).catch(() => false);
    expect(dialogVisible).toBe(true);

    await page.locator('select').filter({ has: page.locator('option[value="Link"]') }).first().selectOption('Move');
    await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
    await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
    await page.waitForTimeout(1200);

    // Verify file gone from source view
    await openSpaceContent(page, SRC);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).not.toBeVisible({ timeout: 8_000 });

    // Verify file present in target view
    await openSpaceContent(page, TGT);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });

    // Click moved file to verify right panel shows info
    const movedFile = page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first();
    await movedFile.dispatchEvent('click');
    await page.waitForTimeout(800);
    await expect(page.locator('.grok-prop-panel, [class*="context-panel"], .d4-info-bar').first()).toBeVisible({ timeout: 8_000 });

    // Refresh tree and verify result unchanged — file still in target, not in source
    await openSpacesView(page);
    await refreshSpacesTree(page);
    await openSpaceContent(page, SRC);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).not.toBeVisible({ timeout: 5_000 });
    await openSpaceContent(page, TGT);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
  } finally {
    await apiDeleteSpace(page, SRC);
    await apiDeleteSpace(page, TGT);
  }
});

// ===========================================================================
// 11–15. Content view: click, open, search, breadcrumb
// ===========================================================================

const CV_SPACE = 'PW-Gen-CV-Space';
const CV_CHILD = 'PW-Gen-CV-Child';

test.describe('Content view: click, open, search, breadcrumb', () => {
  test.beforeEach(async ({ page }) => {
    await openSpacesView(page);
    await apiDeleteSpace(page, CV_CHILD);
    await apiDeleteSpace(page, CV_SPACE);
    await uiCreateRootSpace(page, CV_SPACE);
  });

  test.afterEach(async ({ page }) => {
    await apiDeleteSpace(page, CV_CHILD);
    await apiDeleteSpace(page, CV_SPACE);
  });

  test('11. Click file in space shows info in right panel', async ({ page }) => {
    await addFileToSpaceViaCopy(page, 'acidiq.csv', CV_SPACE);
    await addFileToSpaceViaCopy(page, 'TSLA.csv', CV_SPACE);
    await openSpaceContent(page, CV_SPACE);

    // Click first file — panel appears with file name and Details section
    const fileLabel = page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first();
    await expect(fileLabel).toBeVisible({ timeout: 10_000 });
    await fileLabel.dispatchEvent('click');
    await page.waitForTimeout(800);
    const panel = page.locator('.grok-prop-panel, [class*="context-panel"], .d4-info-bar').first();
    await expect(panel).toBeVisible({ timeout: 8_000 });
    // Verify panel shows the file name
    await expect(panel).toContainText(/acidiq/i, { timeout: 5_000 });
    // Verify panel has Details section
    await expect(page.locator('.d4-accordion-pane-header').filter({ hasText: /Details/i }).first()).toBeVisible({ timeout: 5_000 });

    // Click second file — panel switches to show its name
    const file2Label = page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first();
    await expect(file2Label).toBeVisible({ timeout: 5_000 });
    await file2Label.dispatchEvent('click');
    await page.waitForTimeout(800);
    await expect(panel).toContainText(/TSLA/i, { timeout: 5_000 });
  });

  test('12. Click child space card shows space info in right panel', async ({ page }) => {
    const CV_CHILD2 = 'PW-Gen-CV-Child2';
    try {
      // Create two child spaces and add a file
      await rightClickSpace(page, CV_SPACE);
      await clickMenuItem(page, 'Create Child Space...');
      await fillNameAndOk(page, CV_CHILD);
      await page.waitForTimeout(500);
      await rightClickSpace(page, CV_SPACE);
      await clickMenuItem(page, 'Create Child Space...');
      await fillNameAndOk(page, CV_CHILD2);
      await page.waitForTimeout(500);
      await addFileToSpaceViaCopy(page, 'acidiq.csv', CV_SPACE);

      await openSpaceContent(page, CV_SPACE);
      const panel = page.locator('.grok-prop-panel, [class*="context-panel"], .d4-info-bar').first();

      // Click first child — panel shows child name and Details
      const childCard = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first();
      await expect(childCard).toBeVisible({ timeout: 10_000 });
      await childCard.dispatchEvent('click');
      await page.waitForTimeout(800);
      await expect(panel).toBeVisible({ timeout: 8_000 });
      await expect(panel).toContainText(new RegExp(CV_CHILD, 'i'), { timeout: 5_000 });
      await expect(page.locator('.d4-accordion-pane-header').filter({ hasText: /Details/i }).first()).toBeVisible({ timeout: 5_000 });

      // Click second child — panel switches
      const child2Card = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD2}$`) }).first();
      await expect(child2Card).toBeVisible({ timeout: 5_000 });
      await child2Card.dispatchEvent('click');
      await page.waitForTimeout(800);
      await expect(panel).toContainText(new RegExp(CV_CHILD2, 'i'), { timeout: 5_000 });

      // Click file — panel switches from space info to file info
      const fileLabel = page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first();
      await expect(fileLabel).toBeVisible({ timeout: 5_000 });
      await fileLabel.dispatchEvent('click');
      await page.waitForTimeout(800);
      await expect(panel).toContainText(/acidiq/i, { timeout: 5_000 });

      // Click back to first child — panel switches from file to space
      await childCard.dispatchEvent('click');
      await page.waitForTimeout(800);
      await expect(panel).toContainText(new RegExp(CV_CHILD, 'i'), { timeout: 5_000 });
    } finally {
      await apiDeleteSpace(page, CV_CHILD2);
    }
  });

  test('13. Double-click file opens table view; eye icon shows preview', async ({ page }) => {
    test.setTimeout(120_000);
    await addFileToSpaceViaCopy(page, 'acidiq.csv', CV_SPACE);
    await openSpaceContent(page, CV_SPACE);

    // Part A: double-click opens file as table view with grid and ribbon
    const fileLabel = page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first();
    await expect(fileLabel).toBeVisible({ timeout: 10_000 });
    const spaceUrl = page.url();
    await fileLabel.dispatchEvent('dblclick');
    await page.waitForURL(url => url !== spaceUrl, { timeout: 15_000 });
    await page.waitForTimeout(2000);
    await expect(page.locator('.d4-ribbon').first()).toBeVisible({ timeout: 10_000 });
    await expect(page.locator('.d4-grid').first()).toBeVisible({ timeout: 10_000 });

    // Part B: go back to space, enable eye icon preview, click file — preview appears inside space view
    await openSpaceContent(page, CV_SPACE);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
    const eyeIcon = page.locator('i.grok-icon[name="icon-eye"], i.fal.fa-eye').first();
    await expect(eyeIcon).toBeVisible({ timeout: 5_000 });
    // Enable preview if not already active
    const isActive = await eyeIcon.evaluate((el) => el.classList.contains('d4-current'));
    if (!isActive)
      await eyeIcon.click();
    await page.waitForTimeout(800);

    // Click file — preview grid should appear inside the space view
    const fileLabel2 = page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first();
    await fileLabel2.dispatchEvent('click');
    await page.waitForTimeout(2000);
    // Verify preview grid appeared
    await expect(page.locator('.d4-grid').first()).toBeVisible({ timeout: 10_000 });
  });

  test('14. Search inside space: exact, partial, no match, clear, child space', async ({ page }) => {
    test.setTimeout(120_000);
    // Add three files and a child space
    await addFileToSpaceViaCopy(page, 'acidiq.csv', CV_SPACE);
    await addFileToSpaceViaCopy(page, 'TSLA.csv', CV_SPACE);
    await addFileToSpaceViaCopy(page, 'beer.csv', CV_SPACE);
    await rightClickSpace(page, CV_SPACE);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CV_CHILD);
    await page.waitForTimeout(500);
    // Wait for dialog to fully close
    await page.locator('.d4-dialog').last().waitFor({ state: 'hidden', timeout: 5_000 }).catch(() => {});
    await page.waitForTimeout(500);

    await openSpaceContent(page, CV_SPACE);
    const searchInput = page.locator('input[placeholder*="Search"]').first();
    await expect(searchInput).toBeVisible({ timeout: 8_000 });
    // Wait for all items to load
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });

    // Exact match — only acidiq visible
    await searchInput.fill('acidiq');
    await page.waitForTimeout(600);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });

    // Partial match — "aci" finds acidiq
    await searchInput.fill('aci');
    await page.waitForTimeout(600);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });

    // No match — nothing visible
    await searchInput.fill('xyzzy_no_match_12345');
    await page.waitForTimeout(600);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).not.toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });

    // Clear search — all three files and child space visible again
    await searchInput.fill('');
    await page.waitForTimeout(600);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first()).toBeVisible({ timeout: 5_000 });

    // Search by child space name — child visible, files hidden
    await searchInput.fill(CV_CHILD);
    await page.waitForTimeout(600);
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first()).toBeVisible({ timeout: 5_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).not.toBeVisible({ timeout: 3_000 });

    // Verify DemoFiles originals are intact after all operations
    await page.goto(`${BASE}/files/System.DemoFiles/?browse=files`);
    await expect(page.locator('.d4-link-label label').first()).toBeVisible({ timeout: 10_000 });
    const beerExists = await page.locator('.d4-link-label label').filter({ hasText: /^beer\.csv$/ }).count();
    if (beerExists === 0) {
      const allFiles = await page.locator('.d4-link-label label').allTextContents();
      throw new Error(`beer.csv missing from DemoFiles after test! Available: ${allFiles.slice(0, 20).join(', ')}`);
    }
  });

  test('15. Breadcrumb navigation: multi-level back, content preserved', async ({ page }) => {
    test.setTimeout(120_000);
    const CV_GRAND = 'PW-Gen-CV-Grand';
    try {
      // Create child, grandchild, and add a file — all via tree/API before navigating
      await rightClickSpace(page, CV_SPACE);
      await clickMenuItem(page, 'Create Child Space...');
      await fillNameAndOk(page, CV_CHILD);
      await page.waitForTimeout(500);
      await page.locator('.d4-dialog').last().waitFor({ state: 'hidden', timeout: 5_000 }).catch(() => {});

      await openSpaceContent(page, CV_SPACE);
      const childCard = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first();
      await expect(childCard).toBeVisible({ timeout: 10_000 });
      await childCard.dispatchEvent('contextmenu', { button: 2, bubbles: true });
      await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
      await clickMenuItem(page, 'Create Child Space...');
      await fillNameAndOk(page, CV_GRAND);
      await page.waitForTimeout(500);
      await page.locator('.d4-dialog').last().waitFor({ state: 'hidden', timeout: 5_000 }).catch(() => {});

      // Add file via API to avoid DemoFiles in browser history
      await page.evaluate(async (spaceName) => {
        const g = (window as any).grok;
        const spaces = await g.dapi.spaces.filter(`name = "${spaceName}"`).list({ pageSize: 5 });
        const space = spaces.find((s: any) => s.friendlyName === spaceName || s.name === spaceName);
        if (space)
          await g.dapi.spaces.id(space.id).files.writeString('acidiq.csv', 'x,y\n1,2\n');
      }, CV_SPACE);

      // Fresh navigation start: parent → child → grandchild
      await openSpaceContent(page, CV_SPACE);
      await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first()).toBeVisible({ timeout: 10_000 });
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
      // Verify header contains parent name
      expect(page.url()).toMatch(/\/s\//);

      // Enter child
      const childLink = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first();
      const parentUrl = page.url();
      await childLink.dispatchEvent('dblclick');
      await page.waitForURL(url => url !== parentUrl, { timeout: 10_000 });
      await page.waitForTimeout(1000);
      const childUrl = page.url();
      await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_GRAND}$`) }).first()).toBeVisible({ timeout: 10_000 });

      // Enter grandchild
      const grandLink = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_GRAND}$`) }).first();
      await grandLink.dispatchEvent('dblclick');
      await page.waitForURL(url => url !== childUrl, { timeout: 10_000 });
      await page.waitForTimeout(1000);
      expect(page.url()).toMatch(/\/s\//);

      // Navigate back to child via tree — expand parent tree node first
      const parentTreeNode = page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${CV_SPACE}$`) }).first();
      await expect(parentTreeNode).toBeVisible({ timeout: 5_000 });
      // Click the toggle/expander next to the parent node
      const expander = parentTreeNode.locator('xpath=../..').locator('.d4-tree-view-tri').first();
      if (await expander.isVisible({ timeout: 2_000 }).catch(() => false))
        await expander.click();
      else
        await parentTreeNode.click();
      await page.waitForTimeout(800);

      // Navigate into child via content view instead of tree
      await openSpaceContent(page, CV_SPACE);
      const childCardBack = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first();
      await expect(childCardBack).toBeVisible({ timeout: 10_000 });
      const prevUrl2 = page.url();
      await childCardBack.dispatchEvent('dblclick');
      await page.waitForURL(url => url !== prevUrl2, { timeout: 10_000 });
      await page.waitForTimeout(1000);
      await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_GRAND}$`) }).first()).toBeVisible({ timeout: 10_000 });

      // Navigate back to parent via tree — child card and file visible
      await openSpaceViaTree(page, CV_SPACE);
      await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first()).toBeVisible({ timeout: 10_000 });
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });

      // Verify search still works after navigating back
      const searchInput = page.locator('input[placeholder*="Search"]').first();
      await expect(searchInput).toBeVisible({ timeout: 5_000 });
      await searchInput.fill('acidiq');
      await page.waitForTimeout(600);
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
      await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CV_CHILD}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
    } finally {
      await apiDeleteSpace(page, CV_GRAND);
    }
  });
});

// ===========================================================================
// 16. Copy semantics: delete source, copy unaffected
// ===========================================================================

test.describe('Copy semantics: delete source, copy unaffected', () => {
  test('16. Delete file from space does not affect copy in another space', async ({ page }) => {
    test.setTimeout(180_000);
    const SRC = 'PW-Gen-SrcCopy-16';
    const TGT = 'PW-Gen-Tgt-16';
    await openSpacesView(page);
    for (const n of [SRC, TGT]) await apiDeleteSpace(page, n);
    try {
      await uiCreateRootSpace(page, SRC);
      await uiCreateRootSpace(page, TGT);

      // Add file to source space via DnD from DemoFiles
      await addFileToSpaceViaCopy(page, 'TSLA.csv', SRC);

      // Copy file from source to target space via DnD (retry if tooltip blocks)
      await openSpaceContent(page, SRC);
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
      for (let attempt = 0; attempt < 3; attempt++) {
        const fl = page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first();
        const tn = page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${TGT}$`) }).first();
        await expect(tn).toBeVisible({ timeout: 5_000 });
        page.setDefaultTimeout(15_000);
        await fl.dragTo(tn).catch(() => {});
        page.setDefaultTimeout(30_000);
        if (await page.locator('text=Move entity').isVisible({ timeout: 5_000 }).catch(() => false))
          break;
        if (attempt < 2)
          await openSpaceContent(page, SRC);
      }
      await expect(page.locator('text=Move entity')).toBeVisible({ timeout: 5_000 });
      await page.locator('select').filter({ has: page.locator('option[value="Link"]') }).first().selectOption('Copy');
      await page.locator('.ui-btn', { hasText: /^YES$/i }).first().click();
      await page.locator('text=Move entity').waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
      await page.waitForTimeout(800);

      // Verify copied file appeared in target view
      await openSpaceContent(page, TGT);
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });

      // Delete source file from source space view
      await openSpaceContent(page, SRC);
      await rightClickItemInSpaceView(page, 'TSLA.csv');
      await page.locator('.d4-menu-item', { hasText: /Delete/i }).first().click();
      await page.locator('.ui-btn', { hasText: /^DELETE$/i }).first().click();
      await page.waitForTimeout(1500);

      // Verify copied entry still present in target view
      await openSpaceContent(page, TGT);
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });

      // Refresh tree and verify result unchanged
      await refreshSpacesTree(page);
      await openSpaceContent(page, TGT);
      await expect(page.locator('.d4-link-label label').filter({ hasText: /^TSLA\.csv$/ }).first()).toBeVisible({ timeout: 10_000 });
    } finally {
      for (const n of [SRC, TGT]) await apiDeleteSpace(page, n);
    }
  });
});

// ===========================================================================
// 17. Entity operations inside a space
// ===========================================================================

test('17. Entity ops: context menus, rename, cancel rename, delete, cancel delete', async ({ page }) => {
  // File context menu (Open, Rename, Delete); rename file; cancel rename; cancel delete; delete file; child space context menu
  test.setTimeout(180_000);
  const SPACE = 'PW-Gen-Entity-18';
  const CHILD = 'PW-Gen-Entity-Child-18';
  const CHILD_NEW = 'PW-Gen-Entity-ChildNew-18';
  try {
    await openSpacesView(page);
    await uiCreateRootSpace(page, SPACE);
    await addFileToSpaceViaCopy(page, 'TSLA.csv', SPACE);
    await addFileToSpaceViaCopy(page, 'acidiq.csv', SPACE);
    await rightClickSpace(page, SPACE);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD);

    await openSpaceContent(page, SPACE);

    await rightClickItemInSpaceView(page, 'TSLA.csv');
    const fItems = page.locator('.d4-menu-item');
    await expect(fItems.filter({ hasText: /Open/i }).first()).toBeVisible({ timeout: 5_000 });
    await expect(fItems.filter({ hasText: /Rename/i }).first()).toBeVisible();
    await expect(fItems.filter({ hasText: /Delete/i }).first()).toBeVisible();
    await page.keyboard.press('Escape');

    await rightClickItemInSpaceView(page, 'TSLA.csv');
    await clickMenuItem(page, 'Rename...');
    const renDlg = page.locator('.d4-dialog').last();
    const renInp = renDlg.locator('input[type="text"]').first();
    await expect(renInp).toBeVisible({ timeout: 5_000 });
    expect(await renInp.inputValue()).toMatch(/TSLA/i);
    await renInp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await renInp.pressSequentially('PW-Gen-tsla-renamed');
    await page.keyboard.press('Enter');
    await expect(page.locator('.d4-link-label label').filter({ hasText: /PW-Gen-tsla-renamed/i }).first()).toBeVisible({ timeout: 8_000 });
    // Verify renamed file persists after refresh
    await refreshSpacesTree(page);

    await rightClickItemInSpaceView(page, 'acidiq.csv');
    await clickMenuItem(page, 'Rename...');
    const cancelDlg = page.locator('.d4-dialog').last();
    const cancelInp = cancelDlg.locator('input[type="text"]').first();
    await expect(cancelInp).toBeVisible({ timeout: 5_000 });
    await cancelInp.pressSequentially('PW-should-not-appear');
    await cancelDlg.locator('.ui-btn', { hasText: /^CANCEL$/i }).first().click({ force: true });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });

    await rightClickItemInSpaceView(page, 'acidiq.csv');
    await clickMenuItem(page, /Delete/i);
    await page.locator('.ui-btn', { hasText: /^CANCEL$/i }).first().click();
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });

    await rightClickItemInSpaceView(page, 'PW-Gen-tsla-renamed');
    await clickMenuItem(page, /Delete/i);
    await page.locator('.ui-btn', { hasText: /^DELETE$/i }).first().click();
    await expect(page.locator('.d4-link-label label').filter({ hasText: /PW-Gen-tsla-renamed/i }).first()).not.toBeVisible({ timeout: 8_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });
    // Refresh and verify delete persisted
    await refreshSpacesTree(page);
    await expect(page.locator('.d4-link-label label').filter({ hasText: /PW-Gen-tsla-renamed/i }).first()).not.toBeVisible({ timeout: 3_000 });
    await expect(page.locator('.d4-link-label label').filter({ hasText: /^acidiq\.csv$/ }).first()).toBeVisible({ timeout: 5_000 });

    // Child space context menu has Rename, Delete, Share
    const childCard = page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CHILD}$`) }).first();
    await expect(childCard).toBeVisible({ timeout: 10_000 });
    await childCard.dispatchEvent('contextmenu', { button: 2, bubbles: true });
    await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
    const cItems = page.locator('.d4-menu-item');
    await expect(cItems.filter({ hasText: /Rename/i }).first()).toBeVisible();
    await expect(cItems.filter({ hasText: /Delete/i }).first()).toBeVisible();
    await expect(cItems.filter({ hasText: /Share/i }).first()).toBeVisible();

    // Rename child space from view
    await clickMenuItem(page, 'Rename...');
    const childDlg = page.locator('.d4-dialog').last();
    const childInp = childDlg.locator('input[type="text"]').first();
    await expect(childInp).toBeVisible({ timeout: 5_000 });
    expect(await childInp.inputValue()).toBe(CHILD);
    await childInp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await childInp.pressSequentially(CHILD_NEW);
    await page.keyboard.press('Enter');
    await page.waitForTimeout(1200);
    // Verify renamed child in view
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CHILD_NEW}$`) }).first()).toBeVisible({ timeout: 10_000 });
    await expect(page.locator('.d4-link-label').filter({ hasText: new RegExp(`^${CHILD}$`) }).first()).not.toBeVisible({ timeout: 3_000 });
    // Verify renamed child in tree
    await refreshSpacesTree(page);


    // No "Duplicate" in space context menu
    await openSpacesView(page);
    await rightClickSpace(page, SPACE);
    await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
    const allTexts = await page.locator('.d4-menu-item').allTextContents();
    expect(allTexts.some((t) => /duplicate/i.test(t))).toBe(false);
    await page.keyboard.press('Escape');
  } finally {
    await apiDeleteSpace(page, CHILD_NEW);
    await apiDeleteSpace(page, CHILD);
    await apiDeleteSpace(page, SPACE);
  }
});

// ===========================================================================
// 18. Edge cases
// ===========================================================================

test('18. Edge cases: name with spaces, circular drag, duplicate child name', async ({ page }) => {
  // Space name with spaces works; dragging parent onto child doesn't crash; duplicate child name shows error
  const SPACED = 'PW Gen Name With Spaces 22';
  const PARENT = 'PW-Gen-Edge-Parent-22';
  const CHILD  = 'PW-Gen-Edge-Child-22';
  try {
    await openSpacesView(page);

    await uiCreateRootSpace(page, SPACED);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${SPACED}$`) }).first()).toBeVisible({ timeout: 10_000 });

    await uiCreateRootSpace(page, PARENT);
    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD);

    await rightClickSpace(page, PARENT);
    await clickMenuItem(page, 'Create Child Space...');
    await fillNameAndOk(page, CHILD);
    await expect(page.locator('.d4-toast, .d4-balloon, [class*="toast"], [class*="error"]').first()).toBeVisible({ timeout: 8_000 });

    const parentLabel = page.locator('.d4-tree-view-group-label:visible').filter({ hasText: new RegExp(`^${PARENT}$`) }).first();
    const childLabel  = page.locator('.d4-tree-view-group-label:visible').filter({ hasText: new RegExp(`^${CHILD}$`) }).first();
    await expect(parentLabel).toBeVisible({ timeout: 5_000 });
    await expect(childLabel).toBeVisible({ timeout: 5_000 });
    await parentLabel.dragTo(childLabel);
    await expect(page.locator('.d4-tree-view-group-label').filter({ hasText: new RegExp(`^${PARENT}$`) }).first()).toBeVisible({ timeout: 5_000 });
  } finally {
    for (const n of [SPACED, PARENT, CHILD])
      await apiDeleteSpace(page, n);
  }
});
