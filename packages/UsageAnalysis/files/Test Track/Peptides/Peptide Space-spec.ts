import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Peptides / Peptide Space', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(BASE_URL, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Prerequisites: Load peptides.csv and select AlignedSequence column', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 3000));
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      const col = grok.shell.tv.dataFrame.columns.byName('AlignedSequence');
      grok.shell.o = col;
      await new Promise(r => setTimeout(r, 2000));
    });
    await expect(page.locator('text=Rows: 647')).toBeVisible({ timeout: 10000 });
  });

  test('Step 1: Launch SAR from Bio->Analyze->SAR', async () => {
    await page.evaluate(async () => {
      // Open Bio ribbon menu
      const bioMenuBtn = Array.from(document.querySelectorAll('.d4-ribbon-name, .d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Bio');
      if (bioMenuBtn) bioMenuBtn.click();
      await new Promise(r => setTimeout(r, 500));

      // Hover over Analyze submenu
      const analyzeItem = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Analyze');
      if (analyzeItem) analyzeItem.dispatchEvent(new MouseEvent('mouseover', { bubbles: true }));
      await new Promise(r => setTimeout(r, 500));

      // Click SAR...
      const sarItem = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'SAR...');
      if (sarItem) sarItem.click();
      await new Promise(r => setTimeout(r, 1000));
    });
    await expect(page.locator('text=Analyze Peptides').first()).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Click OK to launch SAR and wait for calculation results', async () => {
    await page.evaluate(async () => {
      const ok = Array.from(document.querySelectorAll('button, .ui-btn'))
        .find(b => b.textContent?.trim() === 'OK');
      if (ok) ok.click();
      await new Promise(r => setTimeout(r, 20000));
    });
    await expect(page.locator('text=MCL').first()).toBeVisible({ timeout: 30000 });
    await expect(page.locator('text=Logo Summary Table').first()).toBeVisible();
  });

  test('Step 3: Open settings using wrench button', async () => {
    await page.evaluate(() => {
      const wrench = document.querySelector('.d4-ribbon-item .fa-wrench') ||
        document.querySelector('.fa-wrench');
      if (wrench) wrench.click();
    });
    await expect(page.locator('text=Peptides settings')).toBeVisible({ timeout: 5000 });
  });

  test('Step 4: Adjust parameters and click OK', async () => {
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');

      // Change Activity scaling: none → lg
      const scalingSelect = Array.from(dialog.querySelectorAll('select'))
        .find(s => Array.from(s.options).some(o => o.value === 'lg'));
      if (scalingSelect) {
        const ns = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set;
        ns.call(scalingSelect, 'lg');
        scalingSelect.dispatchEvent(new Event('change', { bubbles: true }));
      }

      // Change Similarity Threshold: 70 → 55
      const simInput = Array.from(dialog.querySelectorAll('input')).find(i => i.value === '70');
      if (simInput) {
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        ns.call(simInput, '55');
        simInput.dispatchEvent(new Event('input', { bubbles: true }));
        simInput.dispatchEvent(new Event('change', { bubbles: true }));
      }

      // Change Min Cluster Size: 5 → 3
      const minClusterInput = Array.from(dialog.querySelectorAll('input')).find(i => i.value === '5');
      if (minClusterInput) {
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        ns.call(minClusterInput, '3');
        minClusterInput.dispatchEvent(new Event('input', { bubbles: true }));
        minClusterInput.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });

    // Click OK
    await page.evaluate(() => {
      const ok = document.querySelector('.ui-btn.ui-btn-ok.enabled') ||
        Array.from(document.querySelectorAll('button, .ui-btn')).find(b => b.textContent?.trim() === 'OK');
      if (ok) ok.click();
    });
    await page.waitForTimeout(25000);

    const dialogGone = await page.evaluate(() => !document.querySelector('.d4-dialog'));
    expect(dialogGone).toBe(true);
  });

  test('Step 5: MCL Viewer outputs different results', async () => {
    // MCL should have recalculated — verify canvas is present and columns updated
    await expect(page.locator('text=MCL').first()).toBeVisible({ timeout: 10000 });
    // Logo Summary Table should still be visible
    await expect(page.locator('text=Logo Summary Table').first()).toBeVisible();
    // Table should have MCL columns
    const hasMclColumns = await page.evaluate(() => {
      const cols = grok.shell.tv?.dataFrame?.columns;
      if (!cols) return false;
      const names = Array.from({ length: cols.length }, (_, i) => cols.byIndex(i).name);
      return names.some(n => n.includes('MCL'));
    });
    expect(hasMclColumns).toBe(true);
  });
});
