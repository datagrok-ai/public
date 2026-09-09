import { test, expect } from '@playwright/test';
import {
  CONTEXT_MENU,
  TREE_EXPAND_ARROW,
  TREE_EXPAND_ARROW_EXPANDED,
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

const MH_PATH = ['Apps', 'Compute', 'Model-Hub'];
const UNCAT_PATH = [...MH_PATH, 'Uncategorized'];

test.describe('Browse Model Hub (Browse-ModelHub-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-ModelHub-01 — Apps > Compute > Model Hub group is reachable without errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, 'Compute');
    const mh = treeNodeByPath(page, MH_PATH);
    await expect(mh, 'Model Hub node must be present').toBeVisible({ timeout: 10_000 });

    await mh.click();
    await page.waitForTimeout(1500);

    await expectNoErrors(page, sink);
  });

  test('Browse-ModelHub-02 — single click on a model in the tree updates Context Panel without errors', async ({ page }) => {

    test.fail(true, 'Platform regression GROK-19740 — model click throws inside Dart code.');

    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, 'Compute');
    await expandTreeGroup(page, 'Model Hub');
    await expandTreeGroup(page, 'Uncategorized');

    const model = page.locator(`[name^="tree-${UNCAT_PATH.join('---')}---"]`).first();
    await expect(model, 'At least one model must exist under Uncategorized').toBeVisible({ timeout: 10_000 });
    await model.click();
    await page.waitForTimeout(1500);
    await model.hover();
    await page.waitForTimeout(800);

    await expectNoErrors(page, sink);
  });

  test('Browse-ModelHub-03 — double click on a model does not crash, right-click Run is reachable', async ({ page }) => {

    test.fail(true, 'Platform regression GROK-19965 — double-click on a model throws in Dart code.');

    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, 'Compute');
    await expandTreeGroup(page, 'Model Hub');
    await expandTreeGroup(page, 'Uncategorized');

    const model = page.locator(`[name^="tree-${UNCAT_PATH.join('---')}---"]`).first();
    await expect(model).toBeVisible({ timeout: 10_000 });
    await model.dblclick();
    await page.waitForTimeout(2000);
    await expectNoErrors(page, sink);

    const label = model.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();
    await label.click({ button: 'right' });
    await expect(page.locator(CONTEXT_MENU)).toBeVisible({ timeout: 5_000 });
    const runItem = contextMenuItem(page, 'Run');
    expect(await runItem.count(), 'A "Run" item should be present in a model context menu')
      .toBeGreaterThanOrEqual(1);

    await page.keyboard.press('Escape');
    await expectNoErrors(page, sink);
  });

  test('Browse-ModelHub-04 — Uncategorized models group expands without errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, 'Compute');
    await expandTreeGroup(page, 'Model Hub');

    const uncat = treeNodeByPath(page, UNCAT_PATH);
    await uncat.waitFor({ state: 'visible', timeout: 10_000 });
    const tri = uncat.locator(TREE_EXPAND_ARROW).first();
    if (!(await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false))) {
      await tri.click();
      await page.waitForTimeout(1500);
    }
    await expect(tri, 'Uncategorized arrow should be expanded (ref: GROK-19628)')
      .toHaveClass(/d4-tree-view-tri-expanded/);

    await expectNoErrors(page, sink);
  });

});
