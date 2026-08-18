import { test, expect } from './helpers';
import {
  clickDialogPrimary, clickTopMenuLeaf, currentViewerTypes, openDemoCsv,
  resetShell, waitForDialog,
} from './helpers';

test.describe.serial('EDA / Multivariate Analysis', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('MVA on cars.csv produces Grid + 3 Scatter plots + 2 Bar charts', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'cars.csv');

    await clickTopMenuLeaf(page, 'div-ML---Analyze---Multivariate-Analysis...');
    await waitForDialog(page, 'Multivariate Analysis');

    await clickDialogPrimary(page, ['Run', 'RUN', 'OK']);

    await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /Multivariate Analysis/i }))
      .toHaveCount(0, { timeout: 15_000 });

    await page.waitForFunction(() => {
      const tv = (window as any).grok?.shell?.tv;
      if (!tv) return false;
      const lc = Array.from(tv.viewers).map((v: any) => String(v?.type ?? '').toLowerCase());
      return lc.filter((t) => t.includes('scatter')).length >= 3
        && lc.filter((t) => t.includes('bar')).length >= 2
        && lc.some((t) => t === 'grid');
    }, undefined, { timeout: 120_000 });

    const lc = (await currentViewerTypes(page)).map((t) => t.toLowerCase());
    expect(lc.filter((t) => t.includes('scatter')).length, 'scatter viewers').toBeGreaterThanOrEqual(3);
    expect(lc.filter((t) => t.includes('bar')).length, 'bar viewers').toBeGreaterThanOrEqual(2);
    expect(lc.some((t) => t === 'grid'), 'grid viewer').toBe(true);

    const rowCount = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.rowCount);
    expect(rowCount).toBe(30);
  });
});
