import { test, expect } from '@playwright/test';
import { treeGroupByName, treeNodeByPath, viewTabHandle } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

test.describe('Browse Apps (Browse-Apps-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Apps-01 — Apps list loads in reasonable time without errors', async ({ page }) => {
    const sink = watchErrors(page);

    const t0 = Date.now();
    await expandTreeGroup(page, 'Apps');

    await expect(treeNodeByPath(page, ['Apps', 'Tutorials']),
      'A known Apps child should appear after expand').toBeVisible({ timeout: 30_000 });
    const elapsedMs = Date.now() - t0;

    const THRESHOLD_MS = 5_000;
    expect(elapsedMs, `Apps list should load within ${THRESHOLD_MS} ms (took ${elapsedMs} ms)`)
      .toBeLessThan(THRESHOLD_MS);

    await expectNoErrors(page, sink);
  });

  test('Browse-Apps-02 — opening a reference app from the tree works without errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    const tutorials = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await tutorials.waitFor({ state: 'visible', timeout: 10_000 });
    await tutorials.click();

    await expect(viewTabHandle(page, 'Tutorials'), 'Tutorials view should be selected')
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Apps-03 — hover tooltip / details do not throw errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    const tutorials = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await tutorials.waitFor({ state: 'visible', timeout: 10_000 });

    await tutorials.hover();
    await page.waitForTimeout(1500);

    await expectNoErrors(page, sink);
  });
});
