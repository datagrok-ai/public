import { test, expect } from './helpers';
import {
  clickTopMenuLeaf, inputHost, openDemoCsv, resetShell, selectAllColumnsInPicker,
  setInputValue, topDialog, waitForDialog,
} from './helpers';

test.describe.serial('EDA / PLS', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('PLS dialog opens on cars.csv with all expected inputs and accepts Components=3', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'cars.csv');

    await clickTopMenuLeaf(page, 'div-ML---Analyze---PLS...');
    await waitForDialog(page, 'PLS');

    const dlg = topDialog(page);
    for (const host of ['Predict', 'Using', 'Components', 'Quadratic']) {
      await expect(dlg.locator(inputHost(host)).first()).toBeVisible({ timeout: 5_000 });
    }

    await selectAllColumnsInPicker(page, 'Using');
    await expect(dlg.locator(`${inputHost('Using')} .ui-input-editor`)).toContainText(/all|\(\d+\)/i);

    await setInputValue(page, 'Components', '3');

    await page.waitForFunction(() => {
      const d = Array.from(document.querySelectorAll('.d4-dialog')).slice(-1)[0];
      const el = d?.querySelector('[name="input-host-Components"] .ui-input-editor') as HTMLInputElement | null;
      return ((el?.value ?? el?.textContent?.trim() ?? '')) === '3';
    }, undefined, { timeout: 5_000 });
    const componentsValue = await dlg.locator(`${inputHost('Components')} .ui-input-editor`)
      .first().evaluate((el) => (el as HTMLInputElement).value ?? el.textContent?.trim() ?? '');
    expect(componentsValue).toBe('3');

    await dlg.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /^PLS$/i }))
      .toHaveCount(0, { timeout: 10_000 });
  });
});
