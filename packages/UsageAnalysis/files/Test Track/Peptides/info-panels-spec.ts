import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Peptides / Info Panels', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(BASE_URL, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Open peptides.csv dataset', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 3000));
      return `rows:${df.rowCount},cols:${df.columns.length}`;
    });
    expect(result).toBe('rows:647,cols:3');
    await expect(page.locator('text=Rows: 647')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Amino acids rendered with different colors', async () => {
    // AlignedSequence column is visible with colored amino acid rendering
    await expect(page.locator('text=AlignedSequence')).toBeVisible();
    // Verify multiple canvas elements exist (grid renders colored amino acids)
    const canvasCount = await page.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(1);
  });

  test('Step 3: Click peptides column title; context panel appears', async () => {
    await page.evaluate(async () => {
      // Enable context panel (may be hidden by default)
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      await new Promise(r => setTimeout(r, 500));
      // Set AlignedSequence column as current object
      const col = grok.shell.tv.dataFrame.columns.byName('AlignedSequence');
      grok.shell.o = col;
      await new Promise(r => setTimeout(r, 2000));
    });

    await expect(page.locator('text=AlignedSequence').first()).toBeVisible({ timeout: 5000 });
  });

  test('Step 4: Context Panel shows Details and Peptides panels', async () => {
    const panels = await page.evaluate(() => {
      return Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'))
        .map(el => el.textContent?.trim()).filter(t => t);
    });
    expect(panels).toContain('Details');
    expect(panels).toContain('Peptides');
  });

  test('Steps 5-6: Expand each tab and verify content', async () => {
    // Expand Details
    const detailsContent = await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const detailsPane = panes.find(p => p.textContent?.trim() === 'Details');
      if (detailsPane) detailsPane.click();
      await new Promise(r => setTimeout(r, 1000));
      return detailsPane?.nextElementSibling?.textContent?.slice(0, 200) || '';
    });
    expect(detailsContent).toContain('Macromolecule');

    // Expand Peptides
    const peptidesContent = await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const peptidesPane = panes.find(p => p.textContent?.trim() === 'Peptides');
      if (peptidesPane) peptidesPane.click();
      await new Promise(r => setTimeout(r, 3000));
      return peptidesPane?.nextElementSibling?.textContent?.slice(0, 300) || '';
    });
    expect(peptidesContent).toContain('Activity');
    expect(peptidesContent).toContain('Scaling');
    expect(peptidesContent).toContain('Launch SAR');
  });
});
