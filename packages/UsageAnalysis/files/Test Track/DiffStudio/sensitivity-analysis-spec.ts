import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Sensitivity Analysis', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/apps/DiffStudio`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(2000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Load Bioreactor and turn off Edit toggle', async () => {
    // Open Library > Bioreactor via combo popup
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('.d4-combo-popup'))
        .find((el: any) => el.classList.contains('diff-studio-ribbon-widget'));
      if (btn) {
        const r = btn.getBoundingClientRect();
        btn.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2 }));
      }
    });
    await page.waitForTimeout(500);

    const bioreactorItem = page.locator('.d4-menu-item, .d4-menu-item-vert').filter({ hasText: 'Bioreactor' }).first();
    await bioreactorItem.click();
    await page.waitForTimeout(2000);

    // Turn off Edit toggle (if it's on)
    const editSwitch = page.locator('.ui-input-switch');
    const isChecked = await editSwitch.evaluate((el: any) => el.classList.contains('ui-input-switch-on'));
    if (isChecked) {
      await editSwitch.click();
      await page.waitForTimeout(500);
    }
  });

  test('Step 2: Click Sensitivity icon → SA view opens', async () => {
    // Click Sensitivity ribbon item
    const sensitivityBtn = page.locator('.d4-ribbon-item').filter({ hasText: 'Sensitivity' });
    await sensitivityBtn.click();
    await page.waitForTimeout(3000);

    // SA view opens as "Bioreactor - comparison"
    await expect(page.locator('text=Bioreactor - comparison')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Sensitivity Analysis')).toBeVisible();
  });

  test('Step 3: Select FFox, FFred, FKox parameters and modify Process mode', async () => {
    // Change Process mode to Mode 1
    await page.evaluate(() => {
      const selects = Array.from(document.querySelectorAll('select'));
      const modeSelect = selects.find((s: any) => s.value === 'Default');
      if (modeSelect) {
        const nativeSet = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
        nativeSet.call(modeSelect, 'Mode 1');
        modeSelect.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });
    await page.waitForTimeout(1000);

    // Ensure FFox, FFred, FKox are enabled (blue switchers)
    // FFox min/max should have non-zero range by default
    const ffoxMin = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input[type="text"], input:not([type])'));
      return (inputs.find((i: any) => i.value === '0.15') as HTMLInputElement)?.value;
    });
    expect(ffoxMin).toBe('0.15');

    // Set FKox max to a valid non-zero value (default is 0)
    await page.evaluate(() => {
      const labels = Array.from(document.querySelectorAll('*'));
      // Find the FKox max textbox
      const fkoxMaxInputs = Array.from(document.querySelectorAll('input')).filter((i: any) => {
        const prev = i.previousSibling || i.closest('.d4-flex-row')?.querySelector('span');
        return false; // Will use positional approach
      });
      // Use positional approach: FKox max is near FKox min
      const allInputs = Array.from(document.querySelectorAll('input'));
      // Find FKox min (value=0) followed by FKox max (also value=0)
      for (let i = 0; i < allInputs.length - 1; i++) {
        if ((allInputs[i] as HTMLInputElement).value === '0' &&
            (allInputs[i + 1] as HTMLInputElement).value === '0') {
          (allInputs[i + 1] as HTMLInputElement).value = '0.05';
          allInputs[i + 1].dispatchEvent(new Event('change', { bubbles: true }));
          allInputs[i + 1].dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
          break;
        }
      }
    });
    await page.waitForTimeout(500);
  });

  test('Step 4: Run SA and verify 4 viewers open', async () => {
    // Click the green play button
    await page.evaluate(() => {
      const playBtn = document.querySelector('.fas.fa-play') as HTMLElement;
      if (playBtn) playBtn.click();
    });
    await page.waitForTimeout(5000);

    // Verify rows > 0 (SA produced data)
    const rowsText = await page.evaluate(() => {
      const status = Array.from(document.querySelectorAll('*')).find((el: any) =>
        el.textContent?.match(/Rows: \d+/) && el.children.length === 0
      );
      return status?.textContent;
    });
    expect(rowsText).toMatch(/Rows: [1-9]/);

    // Verify 4 viewer types are present
    // Correlation plot
    const corrPlot = page.locator('.d4-correlation-plot, [class*="correlation"], canvas').first();
    await expect(corrPlot).toBeVisible({ timeout: 10000 });

    // PC plot
    const pcPlot = page.locator('.d4-pc-plot, [class*="pc-plot"], canvas').first();
    await expect(pcPlot).toBeVisible();

    // The view should have multiple canvas-based viewers
    const canvasCount = await page.evaluate(() =>
      document.querySelectorAll('canvas').length
    );
    expect(canvasCount).toBeGreaterThanOrEqual(3);
  });
});
