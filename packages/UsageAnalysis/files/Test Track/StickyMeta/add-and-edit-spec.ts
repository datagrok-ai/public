import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

// Prerequisite: Create schema and type-spec.ts must have run (TestSchema1 must exist
// and be associated with a Molecule entity type)

test.describe('Sticky Meta / Add and edit', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/files/System.DemoFiles?f=SPGI.csv&browse=files`, {
      waitUntil: 'networkidle', timeout: 30000,
    });
    await page.waitForTimeout(3000);
    // Open SPGI.csv
    await page.evaluate(async () => {
      const spgi = Array.from(document.querySelectorAll('a, .d4-link, label, span'))
        .find(el => el.textContent?.trim() === 'SPGI.csv');
      spgi?.dispatchEvent(new MouseEvent('dblclick', { bubbles: true, cancelable: true }));
      await new Promise(r => setTimeout(r, 5000));
    });
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 2: Select Structure cell in row 1', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv?.dataFrame;
      if (!df) return;
      df.currentRowIdx = 0;
      df.currentCol = df.col('Structure');
    });
    await page.waitForTimeout(500);
    await expect(page.locator('text=Sticky meta')).toBeVisible({ timeout: 5000 });
  });

  test('Step 3: Open Sticky Meta section in Context Panel', async () => {
    await page.evaluate(async () => {
      const header = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Sticky meta');
      header?.click();
      await new Promise(r => setTimeout(r, 1000));
    });
    // Sticky meta pane should expand and show schemas
    const paneContent = await page.evaluate(() => {
      const header = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Sticky meta');
      return header?.nextElementSibling?.textContent?.trim() || '';
    });
    // Some schemas should appear (platform-level schemas like apisamples-int-meta-ntg6l)
    expect(paneContent.length).toBeGreaterThan(0);
  });

  test('Step 1 (sticky column): Add sticky column via Context Panel', async () => {
    // Click + button in sticky meta pane to add a sticky column
    await page.evaluate(async () => {
      const plusBtns = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok .fal.fa-plus'))
        .map(el => el.closest('.ui-btn'));
      const addBtn = plusBtns[plusBtns.length - 1];
      addBtn?.click();
      await new Promise(r => setTimeout(r, 500));
      // Close any combo popup
      document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
      await new Promise(r => setTimeout(r, 300));
    });
    await page.waitForTimeout(1000);
    // Verify a sticky column appeared in the grid
    const hasStickyCol = await page.evaluate(() => {
      const grid = grok.shell.tv?.grid;
      return grid ? Array.from(grok.shell.tv.dataFrame.columns).some(
        (c: any) => c.name === 'int-meta'
      ) : false;
    });
    expect(hasStickyCol).toBe(true);
  });
});
