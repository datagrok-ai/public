import { test, expect } from '@playwright/test';
import {
  CONTEXT_PANEL,
  CONTEXT_PANEL_COLLAPSE_ALL,
  CONTEXT_PANEL_EXPAND_ALL,
  CONTEXT_PANEL_HEADER,
  CONTEXT_PANEL_INNER,
  treeGroupByName,
  treeItemByName,
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

test.describe('Browse Context Panel (Browse-CtxPanel-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-CtxPanel-01 — F4 toggles the Context Panel', async ({ page }) => {
    const sink = watchErrors(page);

    await expect(page.locator(CONTEXT_PANEL)).toBeVisible();
    await page.keyboard.press('F4');
    await expect(page.locator(CONTEXT_PANEL), 'Context Panel should hide after F4')
      .toBeHidden({ timeout: 5_000 });
    await page.keyboard.press('F4');
    await expect(page.locator(CONTEXT_PANEL), 'Context Panel should reappear after second F4')
      .toBeVisible({ timeout: 5_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-CtxPanel-02 — content updates when switching the focused tree object', async ({ page }) => {
    const sink = watchErrors(page);

    // Click two different tree nodes; assert the Context Panel content changes.
    await expandTreeGroup(page, 'Apps');
    const a = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await a.waitFor({ state: 'visible', timeout: 10_000 });
    await a.click();
    await page.waitForTimeout(1000);
    const contentA = (await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '';

    // Re-open Browse, click a different node.
    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Databases');
    const b = treeNodeByPath(page, ['Databases', 'Postgres']);
    await b.waitFor({ state: 'visible', timeout: 10_000 });
    await b.click();
    await page.waitForTimeout(1500);
    const contentB = (await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '';

    expect(contentB.length, 'Context Panel should have content for the second selection').toBeGreaterThan(0);
    expect(contentB, 'Context Panel content should differ between two distinct objects')
      .not.toEqual(contentA);

    await expectNoErrors(page, sink);
  });

  test('Browse-CtxPanel-03 — Back/Forward navigation through visited objects', async ({ page }) => {
    const sink = watchErrors(page);

    // Visit object A (Tutorials) then object B (Postgres) in Context Panel.
    await expandTreeGroup(page, 'Apps');
    await treeNodeByPath(page, ['Apps', 'Tutorials']).click();
    await page.waitForTimeout(1500);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Databases');
    await treeNodeByPath(page, ['Databases', 'Postgres']).click();
    await page.waitForTimeout(1500);
    const onB = ((await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '');
    expect(onB.toLowerCase(), 'Pre-state: Postgres details should be shown').toContain('postgres');

    // Back → A: Context Panel should show Tutorials details.
    const back = page.locator(`${CONTEXT_PANEL_HEADER} [aria-label="Back"]`);
    await back.click({ force: true });
    await page.waitForTimeout(1000);
    const afterBack = ((await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '');
    expect(afterBack.toLowerCase(), 'After Back, Context Panel must show Tutorials')
      .toContain('tutorials');
    expect(afterBack.toLowerCase(), 'After Back, Postgres details must be gone')
      .not.toContain('postgres');

    // Forward → B: Postgres again.
    const forward = page.locator(`${CONTEXT_PANEL_HEADER} [aria-label="Forward"]`);
    await forward.click({ force: true });
    await page.waitForTimeout(1000);
    const afterForward = ((await page.locator(CONTEXT_PANEL_INNER).innerText().catch(() => '')) ?? '');
    expect(afterForward.toLowerCase(), 'After Forward, Context Panel must show Postgres again')
      .toContain('postgres');

    await expectNoErrors(page, sink);
  });

  test('Browse-CtxPanel-04 — Collapse all / Expand all info panes', async ({ page }) => {
    const sink = watchErrors(page);

    // Focus an entity to get info panes in Context Panel.
    await expandTreeGroup(page, 'Apps');
    await treeNodeByPath(page, ['Apps', 'Tutorials']).click();
    await page.waitForTimeout(2000);

    // Collapse all / Expand all are tiny icons in the Context Panel titlebar that may be
    // CSS-hidden (visibility:hidden) until hover; click via force to bypass visibility checks.
    const collapse = page.locator(CONTEXT_PANEL_COLLAPSE_ALL);
    const expand = page.locator(CONTEXT_PANEL_EXPAND_ALL);
    await expect(collapse, 'Collapse all icon should exist in Context Panel header').toHaveCount(1);
    await expect(expand, 'Expand all icon should exist in Context Panel header').toHaveCount(1);
    await collapse.click({ force: true });
    await page.waitForTimeout(500);
    await expand.click({ force: true });
    await page.waitForTimeout(500);

    // The header buttons remain in the DOM after operating; no errors should appear.
    await expectNoErrors(page, sink);
  });
});
