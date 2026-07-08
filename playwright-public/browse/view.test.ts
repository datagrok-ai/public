import { test, expect, Page } from '@playwright/test';
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
  resolveNestedApp,
  selectedViewName,
  clickTreePath,
} from './helpers';

/**
 * Nested app under `Apps > Chem`, in preference order. `Chemspace` ships in its own
 * (optional) package; `MPO profiles` ships with Chem itself, so it is the reliable
 * fallback when the Chemspace package isn't deployed. Resolved per-test against the
 * live tree so the suite adapts to whatever is installed instead of hard-failing.
 */
const NESTED_CHEM_APP_CANDIDATES = ['Chemspace', 'MPO profiles'];

/**
 * Open a tree node and return the resulting selected view's tab-handle name, waiting
 * deterministically for the selection to change from whatever was selected before.
 */
async function openAndCapture(page: Page, segments: string[], dblclick = false): Promise<string> {
  const before = await selectedViewName(page);
  await clickTreePath(page, segments, { dblclick });
  await expect
    .poll(async () => selectedViewName(page), { timeout: 15_000, intervals: [200, 400, 800] })
    .not.toBe(before);
  const name = await selectedViewName(page);
  expect(name, `a view must be selected after opening ${segments.join(' / ')}`).toBeTruthy();
  return name!;
}

/**
 * Browsing mode (unpinned) vs persistent view (pinned by double-click) — section 3.
 * Apps tree items are item-leaves: single click opens / replaces, double click pins.
 */
test.describe('Browse view modes (Browse-View-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
    await expandTreeGroup(page, 'Apps');
  });

  test('Browse-View-01 — single click is unpinned (replaces), double click is persistent', async ({ page }) => {
    const sink = watchErrors(page);

    // X = top-level item; Y = nested item under Chem. Y is resolved against the live tree
    // (Chemspace's package is optional) so the test adapts to whatever nested app is deployed.
    const xPath = ['Apps', 'Tutorials'];
    const yPath = await resolveNestedApp(page, 'Chem', NESTED_CHEM_APP_CANDIDATES);
    test.skip(yPath === null, `no nested Chem app installed (${NESTED_CHEM_APP_CANDIDATES.join(', ')})`);

    // 1. Single click X -> unpinned view becomes current. Capture the real view names.
    const xView = await openAndCapture(page, xPath);
    await expect(viewTabHandle(page, xView), `${xView} should become the current view`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });

    // 2. Single click Y -> Y current; unpinned X must be removed.
    const yView = await openAndCapture(page, yPath!);
    await expect(viewTabHandle(page, yView), `${yView} should become the current view`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });
    await expect(viewTabHandle(page, xView), `Unpinned ${xView} should be replaced by ${yView}`)
      .toHaveCount(0, { timeout: 10_000 });

    // 3. Double click X -> X opens as a persistent view.
    await openAndCapture(page, xPath, true);
    await expect(viewTabHandle(page, xView), `${xView} should be attached as a persistent view`)
      .toHaveCount(1, { timeout: 15_000 });

    // 4. Single click Y again -> Y current; persistent X must remain attached.
    await openAndCapture(page, yPath!);
    await expect(viewTabHandle(page, yView), `${yView} should become the current view again`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });
    await expect(viewTabHandle(page, xView), `Persistent ${xView} must remain attached after single-click on ${yView}`)
      .toHaveCount(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-02 — editing an unpinned view makes it persistent', async ({ page }) => {
    const sink = watchErrors(page);

    // Nested app to single-click at the end — resolved against the live tree.
    const yPath = await resolveNestedApp(page, 'Chem', NESTED_CHEM_APP_CANDIDATES);
    test.skip(yPath === null, `no nested Chem app installed (${NESTED_CHEM_APP_CANDIDATES.join(', ')})`);

    // 1) Open Tutorials as unpinned (single click).
    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await treeNodeByPath(page, ['Apps', 'Tutorials']).click();
    await page.waitForTimeout(2000);

    // Verify it's currently a browsing preview (tab text has the .grok-browse-preview class).
    const tab = viewTabHandle(page, 'Tutorials');
    const previewBefore = await tab.locator(VIEW_TAB_PREVIEW_TEXT).count();
    expect(previewBefore, 'Tutorials should start as a browsing preview (unpinned)').toBe(1);

    // 2) "Edit" the view: click the pin button. The button is in the hidden combo-popup
    //    (tabbed mode), so dispatch a click event programmatically instead of UI click.
    const pin = tab.locator(VIEW_TAB_PIN_BUTTON);
    await pin.dispatchEvent('click');
    await page.waitForTimeout(800);

    // After pinning, the .grok-browse-preview class should be gone.
    const previewAfter = await tab.locator(VIEW_TAB_PREVIEW_TEXT).count();
    expect(previewAfter, 'After pin, Tutorials should no longer be a preview').toBe(0);

    // 3) Single-click another item — pinned Tutorials must remain attached.
    await clickTreePath(page, yPath!);
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

    // Click pin via dispatchEvent (hidden in tabbed-mode popup).
    await pin.dispatchEvent('click');
    await page.waitForTimeout(800);
    await expect(tab.locator(VIEW_TAB_PREVIEW_TEXT),
      'Preview marker must disappear after pinning').toHaveCount(0);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-04 — clicking Browse after a persistent view starts a fresh browsing session', async ({ page }) => {
    const sink = watchErrors(page);

    const yPath = await resolveNestedApp(page, 'Chem', NESTED_CHEM_APP_CANDIDATES);
    test.skip(yPath === null, `no nested Chem app installed (${NESTED_CHEM_APP_CANDIDATES.join(', ')})`);

    // Open a persistent view via double click.
    await openAndCapture(page, ['Apps', 'Tutorials'], true);
    await expect(viewTabHandle(page, 'Tutorials'), 'Tutorials must be attached as persistent')
      .toHaveCount(1, { timeout: 15_000 });

    // Return to Browse and single-click another item.
    const yView = await openAndCapture(page, yPath!);
    await expect(viewTabHandle(page, yView), `${yView} must become the current view`)
      .toHaveClass(/tab-handle-selected/, { timeout: 15_000 });

    // The persistent Tutorials must still be there.
    await expect(viewTabHandle(page, 'Tutorials'), 'Persistent Tutorials must remain attached')
      .toHaveCount(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-View-05 — same-type views show a count badge on the Sidebar', async ({ page }) => {
    const sink = watchErrors(page);

    // Open two table views via the JS API (no built-in UI to open multiple identical tables).
    await page.evaluate(async () => {
      const df1 = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      (window as any).grok.shell.addTableView(df1);
      const df2 = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      (window as any).grok.shell.addTableView(df2);
    });
    await page.waitForTimeout(2500);

    // The view-selector count badge should show ≥ 2 (Home + 2 tables = 3 total views).
    const badge = page.locator(VIEW_SELECTOR_COUNT_BADGE).first();
    const text = (await badge.textContent())?.trim() ?? '';
    const count = parseInt(text, 10);
    expect(Number.isFinite(count), `Badge text "${text}" must be numeric`).toBe(true);
    expect(count, 'Badge count should reflect the number of open views (≥ 2)')
      .toBeGreaterThanOrEqual(2);

    await expectNoErrors(page, sink);
  });
});
