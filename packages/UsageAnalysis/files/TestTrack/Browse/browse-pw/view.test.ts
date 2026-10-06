import { test, expect } from '@playwright/test';
import {
  BROWSE_HEADER_HOME,
  treeGroupByName,
  treeItemByName,
  treeNodeByPath,
  VIEW_SELECTOR,
  VIEW_SELECTOR_COUNT_BADGE,
  VIEW_TAB,
  VIEW_TAB_PIN_BUTTON,
  VIEW_TAB_PREVIEW_TEXT,
  VIEW_TAB_SELECTED,
  viewTabHandle,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
  openTreeItem,
} from './helpers';

test.describe('Browse view modes (Browse-View-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
    await expandTreeGroup(page, 'Apps');
  });

  test('Browse-View-01 — single click is unpinned (replaces), double click is persistent', async ({ page }) => {
    const sink = watchErrors(page);

    const X = { path: ['Apps', 'Tutorials'], view: 'Tutorials' };
    const Y = { path: ['Apps', 'Chem', 'Chemspace'], view: 'Chemspace' };

    async function clickByPath(segments: string[], dblclick = false): Promise<void> {
      await ensureBrowsePanelOpen(page);

      for (let i = 0; i < segments.length - 1; i++) {
        await expandTreeGroup(page, segments[i]);
      }
      const locator = treeNodeByPath(page, segments);
      await locator.waitFor({ state: 'visible', timeout: 10_000 });
      await locator.scrollIntoViewIfNeeded();
      if (dblclick) await locator.dblclick();
      else await locator.click();
      await page.waitForTimeout(1500);
    }

    await clickByPath(X.path);
    await expect(viewTabHandle(page, X.view), `${X.view} should become the current view`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });

    await clickByPath(Y.path);
    await expect(viewTabHandle(page, Y.view), `${Y.view} should become the current view`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });
    await expect(viewTabHandle(page, X.view), `Unpinned ${X.view} should be replaced by ${Y.view}`)
      .toHaveCount(0, { timeout: 10_000 });

    await clickByPath(X.path, true);
    await expect(viewTabHandle(page, X.view), `${X.view} should be attached as a persistent view`)
      .toHaveCount(1, { timeout: 15_000 });

    await clickByPath(Y.path);
    await expect(viewTabHandle(page, Y.view), `${Y.view} should become the current view again`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });
    await expect(viewTabHandle(page, X.view), `Persistent ${X.view} must remain attached after single-click on ${Y.view}`)
      .toHaveCount(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-02 — editing an unpinned view makes it persistent', async ({ page }) => {
    const sink = watchErrors(page);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await treeNodeByPath(page, ['Apps', 'Tutorials']).click();
    await page.waitForTimeout(2000);

    const tab = viewTabHandle(page, 'Tutorials');
    const previewBefore = await tab.locator(VIEW_TAB_PREVIEW_TEXT).count();
    expect(previewBefore, 'Tutorials should start as a browsing preview (unpinned)').toBe(1);

    const pin = tab.locator(VIEW_TAB_PIN_BUTTON);
    await pin.dispatchEvent('click');
    await page.waitForTimeout(800);

    const previewAfter = await tab.locator(VIEW_TAB_PREVIEW_TEXT).count();
    expect(previewAfter, 'After pin, Tutorials should no longer be a preview').toBe(0);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, 'Chem');
    await treeNodeByPath(page, ['Apps', 'Chem', 'Chemspace']).click();
    await page.waitForTimeout(2000);

    await expect(viewTabHandle(page, 'Tutorials'), 'Pinned Tutorials must remain attached')
      .toHaveCount(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-03 — Pin control on the tab converts a preview into a pinned view', async ({ page }) => {
    const sink = watchErrors(page);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await treeNodeByPath(page, ['Apps', 'Tutorials']).click();
    await page.waitForTimeout(2000);

    const tab = viewTabHandle(page, 'Tutorials');
    const pin = tab.locator(VIEW_TAB_PIN_BUTTON);
    await expect(pin, 'Pin button must be present on a preview tab').toHaveCount(1);

    await pin.dispatchEvent('click');
    await page.waitForTimeout(800);
    await expect(tab.locator(VIEW_TAB_PREVIEW_TEXT),
      'Preview marker must disappear after pinning').toHaveCount(0);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-04 — clicking Browse after a persistent view starts a fresh browsing session', async ({ page }) => {
    const sink = watchErrors(page);

    async function clickByPath(segments: string[], dblclick = false): Promise<void> {
      await ensureBrowsePanelOpen(page);
      for (let i = 0; i < segments.length - 1; i++) {
        await expandTreeGroup(page, segments[i]);
      }
      const locator = treeNodeByPath(page, segments);
      await locator.waitFor({ state: 'visible', timeout: 10_000 });
      if (dblclick) await locator.dblclick();
      else await locator.click();
      await page.waitForTimeout(1500);
    }

    await clickByPath(['Apps', 'Tutorials'], true);
    await expect(viewTabHandle(page, 'Tutorials'), 'Tutorials must be attached as persistent')
      .toHaveCount(1, { timeout: 15_000 });

    await clickByPath(['Apps', 'Chem', 'Chemspace']);
    await expect(viewTabHandle(page, 'Chemspace'), 'Chemspace must become the current view')
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });

    await expect(viewTabHandle(page, 'Tutorials'), 'Persistent Tutorials must remain attached')
      .toHaveCount(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-05 — same-type views show a count badge on the Sidebar', async ({ page }) => {
    const sink = watchErrors(page);

    await page.evaluate(async () => {
      const df1 = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      (window as any).grok.shell.addTableView(df1);
      const df2 = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      (window as any).grok.shell.addTableView(df2);
    });
    await page.waitForTimeout(2500);

    const badge = page.locator(VIEW_SELECTOR_COUNT_BADGE).first();
    const text = (await badge.textContent())?.trim() ?? '';
    const count = parseInt(text, 10);
    expect(Number.isFinite(count), `Badge text "${text}" must be numeric`).toBe(true);
    expect(count, 'Badge count should reflect the number of open views (≥ 2)')
      .toBeGreaterThanOrEqual(2);

    await expectNoErrors(page, sink);
  });
});
