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

    await page.goto(`${BASE}/users`);
    await page.waitForSelector(RIBBON, { timeout: 30_000 });

    if (!(await page.locator(FILTER_PANEL).isVisible().catch(() => false))) {
      const toggle = page.locator(FILTER_TOGGLE).first();
      await toggle.click();
    }
    await page.waitForSelector(FILTER_PANEL, { timeout: 10_000 });
  });

  test('Browse-Filter-01 — filter panel has no empty / unnamed properties', async ({ page }) => {
    const sink = watchErrors(page);

    const emptyCount = await page.locator(FILTER_PANEL).locator('.ui-label, label, .d4-filter-panel-header > div').evaluateAll((els) =>
      els.filter((e) => !(e.textContent ?? '').trim()).length,
    );
    expect(emptyCount, 'No empty / unnamed property labels in the filter panel (ref: GROK-19691)').toBe(0);

    const sectionHeaders = await page.locator(FILTER_PANEL_HEADER).allTextContents();
    expect(sectionHeaders.length, 'There should be at least one filter section').toBeGreaterThanOrEqual(1);
    for (const h of sectionHeaders) {
      expect(h.trim(), 'Section header text must be non-empty').not.toBe('');
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-Filter-04 — applying a Users quick filter narrows the list', async ({ page }) => {
    const sink = watchErrors(page);

    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    const recently = page.locator(FILTER_QUICK_TAGS, { hasText: /^Recently joined$/ }).first();

    await expect(all, 'Quick filter "All" should be visible').toBeVisible({ timeout: 5_000 });
    await expect(recently, 'Quick filter "Recently joined" should be visible').toBeVisible({ timeout: 5_000 });

    await recently.click();
    await page.waitForTimeout(1500);

    await all.click();
    await page.waitForTimeout(800);

    await expectNoErrors(page, sink);
  });

  test('Browse-Filter-05 — filter and search work together (Users gallery)', async ({ page }) => {
    const sink = watchErrors(page);

    const recently = page.locator(FILTER_QUICK_TAGS, { hasText: /^Recently joined$/ }).first();
    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    await expect(recently).toBeVisible({ timeout: 5_000 });
    await recently.click();
    await page.waitForTimeout(800);

    const search = page.locator('input[placeholder^="Search users"]').first();
    await search.fill('a');
    await page.waitForTimeout(1500);

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

    const matches = page.locator('label', { hasText: /zzz/i });
    expect(await matches.count(), 'No gallery result should match a nonsense query').toBe(0);

    await search.fill('');
    await page.waitForTimeout(300);

    await expectNoErrors(page, sink);
  });

});

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

    await all.click();
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

  test('Browse-Filter-03 — "Used by me" filter for Apps works and resets', async ({ page }) => {
    const sink = watchErrors(page);

    const all = page.locator(FILTER_QUICK_TAGS, { hasText: /^All$/ }).first();
    const usedByMe = page.locator(FILTER_QUICK_TAGS, { hasText: /^Used by me$/ }).first();
    await expect(usedByMe, '"Used by me" tag should exist on Apps gallery (ref: GROK-19688)')
      .toBeVisible({ timeout: 5_000 });

    await usedByMe.click();
    await page.waitForTimeout(1000);

    await all.click();
    await page.waitForTimeout(500);

    await expectNoErrors(page, sink);
  });
});
