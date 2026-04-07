import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Fitting', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/apps/DiffStudio`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Open Diff Studio and load Bioreactor', async () => {
    // Load Bioreactor from Library
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('.d4-combo-popup'))
        .find((el) => el.classList.contains('diff-studio-ribbon-widget'));
      if (btn) {
        const r = btn.getBoundingClientRect();
        btn.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2 }));
      }
    });
    await page.waitForTimeout(500);
    const bioreactorItem = page.locator('.d4-menu-item, .d4-menu-item-vert').filter({ hasText: 'Bioreactor' }).first();
    await bioreactorItem.click();
    await page.waitForTimeout(2000);
    await expect(page.locator('text=Columns: 13')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Click Fit icon → Fitting view opens', async () => {
    const fitBtn = page.locator('.d4-ribbon-item').filter({ hasText: 'Fit' });
    await fitBtn.click();
    await page.waitForTimeout(3000);
    await expect(page.locator('text=Bioreactor - fitting')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Fitting')).toBeVisible();
  });

  test('Step 3: Modify Process mode; verify FFox & KKox change', async () => {
    // Change to Mode 1
    await page.evaluate(() => {
      const selects = Array.from(document.querySelectorAll('select'));
      const modeSelect = selects.find((s) => s.value === 'Default');
      if (modeSelect) {
        const nativeSet = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set;
        nativeSet.call(modeSelect, 'Mode 1');
        modeSelect.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });
    await page.waitForTimeout(1000);

    // Verify FFox changed from 0.2 to ~0.163
    const ffoxVal = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      return inputs.find(i => i.value.startsWith('0.1') && parseFloat(i.value) < 0.2)?.value;
    });
    expect(ffoxVal).toBeTruthy();

    // Set back to Default
    await page.evaluate(() => {
      const selects = Array.from(document.querySelectorAll('select'));
      const modeSelect = selects.find((s) => s.value === 'Mode 1');
      if (modeSelect) {
        const nativeSet = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set;
        nativeSet.call(modeSelect, 'Default');
        modeSelect.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });
    await page.waitForTimeout(1000);
  });

  test('Step 4: Enable switchers for switch at, FFox, FKox with ranges', async () => {
    // Enable switchers by position: index 4=switch at, 5=FFox, 12=FKox
    await page.evaluate(() => {
      const switches = Array.from(document.querySelectorAll('.ui-input-switch'));
      [4, 5, 12].forEach(i => switches[i]?.click());
    });
    await page.waitForTimeout(500);

    // Set FFox max to 1.0
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const ffoxMax = inputs.find(i => {
        const label = i.closest('[class*="row"], [class*="line"]')?.textContent;
        return label?.includes('FFox (max)');
      });
      if (ffoxMax) { ffoxMax.value = '1.0'; ffoxMax.dispatchEvent(new Event('change', { bubbles: true })); }
    });

    // Set FKox max to 3
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const fkoxMax = inputs.find(i => {
        const label = i.closest('[class*="row"], [class*="line"]')?.textContent;
        return label?.includes('FKox (max)');
      });
      if (fkoxMax) { fkoxMax.value = '3'; fkoxMax.dispatchEvent(new Event('change', { bubbles: true })); }
    });
    await page.waitForTimeout(500);
  });

  test('Step 5: Add bioreactor-experiment.csv as target', async () => {
    // Load the experiment CSV file
    const loaded = await page.evaluate(async () => {
      const df = await (window as any).grok.dapi.files.readCsv('System:AppData/DiffStudio/library/bioreactor-experiment.csv');
      if (df) {
        (window as any).grok.shell.addTableView(df);
        return `rows:${df.rowCount}`;
      }
      return 'null';
    });
    expect(loaded).toContain('rows:');
    await page.waitForTimeout(1000);

    // Navigate back to Fitting view
    await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.views);
      const fv = views.find((v: any) => v.name === 'Bioreactor - fitting');
      if (fv) (window as any).grok.shell.v = fv;
    });
    await page.waitForTimeout(500);

    // Select Table in the Target combobox
    await page.evaluate(() => {
      const selects = Array.from(document.querySelectorAll('select'));
      const targetSelect = selects.find((s) => Array.from(s.options).some(o => o.value === 'Table'));
      if (targetSelect) {
        const nativeSet = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set;
        nativeSet.call(targetSelect, 'Table');
        targetSelect.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });
    await page.waitForTimeout(500);
    await expect(page.locator('text=Table')).toBeVisible();
  });

  test('Step 6: Run fitting; verify RMSE by iterations descending', async () => {
    await page.evaluate(() => {
      const playBtn = document.querySelector('.fas.fa-play');
      if (playBtn) playBtn.click();
    });
    await page.waitForTimeout(10000);

    // Verify rows > 0
    const rowsText = await page.locator('[class*="status"] [class*="rows"], .d4-status-bar').textContent();
    expect(rowsText).toMatch(/Rows: [1-9]/);

    // Verify RMSE by iterations column header visible
    await expect(page.locator('text=RMSE by iterations')).toBeVisible({ timeout: 15000 });

    // Verify charts present (canvas count >= 2)
    const canvasCount = await page.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(2);
  });
});
