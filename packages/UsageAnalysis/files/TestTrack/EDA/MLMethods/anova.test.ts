import { test, expect } from './helpers';
import {
  clickDialogPrimary, clickTopMenuLeaf, currentViewerTypes, openDemoCsv,
  resetShell, visibleTabLabels, waitForDialog,
} from './helpers';

test.describe.serial('EDA / ANOVA', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('ANOVA on demog.csv produces Box plot viewer with Analysis and F-test tabs', async ({ page }) => {
    test.setTimeout(120_000);

    await openDemoCsv(page, 'demog.csv');

    await clickTopMenuLeaf(page, 'div-ML---Analyze---ANOVA...');
    await waitForDialog(page, 'ANOVA');

    await clickDialogPrimary(page, ['Run', 'RUN', 'OK']);

    await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /^ANOVA$/i })).toHaveCount(0, { timeout: 10_000 });

    await page.waitForFunction(() => {
      const tv = (window as any).grok?.shell?.tv;
      return !!tv && Array.from(tv.viewers).some((v: any) => /box ?plot/i.test(String(v?.type ?? '')));
    }, undefined, { timeout: 30_000 });
    expect(await currentViewerTypes(page))
      .toEqual(expect.arrayContaining([expect.stringMatching(/box ?plot/i)]));

    await page.waitForFunction(() => {
      const labels = Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
        .map((el) => (el.textContent ?? '').trim());
      return labels.some((l) => /^Analysis$/i.test(l)) && labels.some((l) => /^F-test$/i.test(l));
    }, undefined, { timeout: 15_000 });
    expect(await visibleTabLabels(page))
      .toEqual(expect.arrayContaining([
        expect.stringMatching(/^Analysis$/i),
        expect.stringMatching(/^F-test$/i),
      ]));
  });
});
