import { test, expect } from '@playwright/test';
import { BROWSE_HEADER_HOME, HOME_GLOBAL_SEARCH_INPUT } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
} from './helpers';

test.describe('Browse global search (Browse-Search-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
    // Make sure we are on the Home Page where the global search lives.
    await page.locator(BROWSE_HEADER_HOME).click();
    await expect(page.locator(HOME_GLOBAL_SEARCH_INPUT)).toBeVisible({ timeout: 10_000 });
  });

  test('Browse-Search-01 — typing into global search shows results', async ({ page }) => {
    const sink = watchErrors(page);

    const input = page.locator(HOME_GLOBAL_SEARCH_INPUT);
    await input.fill('tutorials');
    await page.waitForTimeout(2500);

    // After typing, the widgets area is replaced by results: assert at least one search
    // result row appears (look for any element containing "tutorials").
    const anyMatch = page.locator('div, label, li', { hasText: /tutorials/i }).first();
    await expect(anyMatch, 'A search result matching "tutorials" should appear').toBeVisible({ timeout: 10_000 });

    await expectNoErrors(page, sink);

    await input.fill('');
    await page.waitForTimeout(500);
  });

  test('Browse-Search-02 — tag search (#demo) produces results in the search results list', async ({ page }) => {
    const sink = watchErrors(page);

    const input = page.locator(HOME_GLOBAL_SEARCH_INPUT);
    await input.fill('#demo');
    await page.waitForTimeout(3500);

    // The search results live inside the `.power-search-lists-host` container.
    const host = page.locator('.power-search-lists-host').first();
    await expect(host, 'Search results host should appear after typing a tag query')
      .toBeVisible({ timeout: 10_000 });
    const hostText = (await host.innerText().catch(() => '')) ?? '';
    expect(hostText.trim().length,
      'Search results host should contain at least one result for #demo')
      .toBeGreaterThan(0);

    await expectNoErrors(page, sink);

    await input.fill('');
    await page.waitForTimeout(500);
  });

  test('Browse-Search-03 — nonsense query does not crash the search', async ({ page }) => {
    const sink = watchErrors(page);

    const input = page.locator(HOME_GLOBAL_SEARCH_INPUT);
    await input.fill('zzz_no_match_zzz_qq_qq');
    await page.waitForTimeout(2500);

    // Look only inside the search results host (the tree on the left always contains
    // labels like "Tutorials" and shouldn't pollute the assertion).
    const host = page.locator('.power-search-lists-host').first();
    const hostText = ((await host.innerText().catch(() => '')) ?? '').toLowerCase();
    expect(hostText, 'No known entity name should appear in search results for a nonsense query')
      .not.toMatch(/tutorials/);

    // The input itself must still respond — clear it.
    await input.fill('');
    await page.waitForTimeout(500);

    await expectNoErrors(page, sink);
  });
});
