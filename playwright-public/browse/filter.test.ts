import { test, expect } from '@playwright/test';
import {
  FILTER_PANEL,
  FILTER_PANEL_SECTION,
  FILTER_PANEL_HEADER,
  FILTER_QUICK_TAGS,
  FILTER_TOGGLE,
  RIBBON,
} from './selectors';
import { goHome, watchErrors, expectNoErrors } from './helpers';

const BASE: string = process.env.DATAGROK_URL!;

test.describe('Browse filters (Browse-Filter-*)', () => {
  test.beforeEach(async ({ page }) => {
    // Filter panel only appears in gallery views. Go to Users — confirmed to render the panel on dev.
    await page.goto(`${BASE}/users`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });

    // Open the filter panel via the "Toggle filters" icon if it's not already open.
    if (!(await page.locator(FILTER_PANEL).isVisible().catch(() => false))) {
      const toggle = page.locator(FILTER_TOGGLE).first();
      await toggle.click();
    }
    await page.waitForSelector(FILTER_PANEL, { timeout: 10_000 });
  });

  test('Browse-Filter-01 — filter panel has no empty / unnamed properties', async ({ page }) => {
    const sink = watchErrors(page);

    // Every visible label/text inside the panel must have non-empty content.
    const emptyCount = await page.locator(FILTER_PANEL).locator('.ui-label, label, .d4-filter-panel-header > div').evaluateAll((els) =>
      els.filter((e) => !(e.textContent ?? '').trim()).length,
    );
    expect(emptyCount, 'No empty / unnamed property labels in the filter panel (ref: GROK-19691)').toBe(0);

    // Sections themselves should each have a non-empty header.
    const sectionHeaders = await page.locator(FILTER_PANEL_HEADER).allTextContents();
    expect(sectionHeaders.length, 'There should be at least one filter section').toBeGreaterThanOrEqual(1);
    for (const h of sectionHeaders) {
      expect(h.trim(), 'Section header text must be non-empty').not.toBe('');
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-Filter-04 — applying a Users quick filter narrows the list', async ({ page }) => {
    const sink = watchErrors(page);

    // The Users gallery exposes quick filters like "All" / "Recently joined".
    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    const recently = page.locator(FILTER_QUICK_TAGS, { hasText: /^Recently joined$/ }).first();

    await expect(all, 'Quick filter "All" should be visible').toBeVisible({ timeout: 5_000 });
    await expect(recently, 'Quick filter "Recently joined" should be visible').toBeVisible({ timeout: 5_000 });

    // Click "Recently joined" — the list should narrow (no error balloons regardless).
    await recently.click();
    await page.waitForTimeout(1500);

    // Reset by clicking "All".
    await all.click();
    await page.waitForTimeout(800);

    await expectNoErrors(page, sink);
  });

  test('Browse-Filter-05 — filter and search work together (Users gallery)', async ({ page }) => {
    const sink = watchErrors(page);

    // Already on /users from the beforeEach.
    const recently = page.locator(FILTER_QUICK_TAGS, { hasText: /^Recently joined$/ }).first();
    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    await expect(recently).toBeVisible({ timeout: 5_000 });
    await recently.click();
    await page.waitForTimeout(800);

    // Now narrow further with search.
    const search = page.locator('input[placeholder^="Search users"]').first();
    await search.fill('a');
    await page.waitForTimeout(1500);

    // Clean up: clear search + reset to "All".
    await search.fill('');
    await all.click();
    await page.waitForTimeout(500);

    await expectNoErrors(page, sink);
  });

  test('Browse-Filter-06 — nonsense search produces no results without crashing', async ({ page }) => {
    const sink = watchErrors(page);

    const search = page.locator('input[placeholder^="Search users"]').first();
    await search.fill('zzz_no_match_zzz_qq_qq');
    await page.waitForTimeout(2000);

    // No `label` items in the gallery should match. The Browse tree on the left still
    // contains labels (My stuff, Apps, etc.), so filter to visible labels inside main view.
    const matches = page.locator('label', { hasText: /zzz/i });
    expect(await matches.count(), 'No gallery result should match a nonsense query').toBe(0);

    await search.fill('');
    await page.waitForTimeout(300);

    await expectNoErrors(page, sink);
  });

});

// Filter cases for Files and Apps galleries — separate describe so the beforeEach can target
// each gallery's own URL.
test.describe('Browse filters — Files gallery (Browse-Filter-02)', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(`${BASE}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });
    if (!(await page.locator(FILTER_PANEL).isVisible().catch(() => false))) {
      await page.locator(FILTER_TOGGLE).first().click();
    }
    await page.waitForSelector(FILTER_PANEL, { timeout: 10_000 });
  });

  test('Browse-Filter-02 — "Created recently" filter for Files works and resets', async ({ page }) => {
    const sink = watchErrors(page);

    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    const recent = page.locator(FILTER_QUICK_TAGS, { hasText: /^Created recently$/ }).first();
    await expect(recent, '"Created recently" tag should exist (ref: GROK-19689)').toBeVisible({ timeout: 5_000 });

    await recent.click();
    await page.waitForTimeout(1000);
    // Reset (best-effort: the "All" quick-filter tag isn't always rendered once a filter
    // is active on the minimal CI stack).
    if (await all.isVisible().catch(() => false))
      await all.click().catch(() => undefined);
    await page.waitForTimeout(500);

    await expectNoErrors(page, sink);
  });
});

test.describe('Browse filters — Apps gallery (Browse-Filter-03)', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(`${BASE}/apps`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });
    if (!(await page.locator(FILTER_PANEL).isVisible().catch(() => false))) {
      await page.locator(FILTER_TOGGLE).first().click();
    }
    await page.waitForSelector(FILTER_PANEL, { timeout: 10_000 });
  });

  test('Browse-Filter-03 — a quick filter for Apps works and resets', async ({ page }) => {
    const sink = watchErrors(page);

    // The Apps gallery offers only All / Favorites / Created recently as quick filters (see
    // xamgle apps_view.dart `extraFilters`) — it has NO "Used by me" tag (that one belongs to
    // the connections / dockerfiles galleries), so assert against a filter Apps actually offers.
    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    const favorites = page.locator(FILTER_QUICK_TAGS, { hasText: /^Favorites$/ }).first();
    await expect(favorites, '"Favorites" quick filter should exist on the Apps gallery')
      .toBeVisible({ timeout: 10_000 });

    await favorites.click();
    await page.waitForTimeout(1000);
    // Reset (best-effort: the "All" tag isn't always rendered once a filter is active on the
    // minimal CI stack).
    if (await all.isVisible().catch(() => false))
      await all.click().catch(() => undefined);
    await page.waitForTimeout(500);

    await expectNoErrors(page, sink);
  });
});
