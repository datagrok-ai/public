import { test, expect } from './helpers';
import {
  clickDialogPrimary, clickTopMenuLeaf, currentColumnNames, openDemoCsv, resetShell,
  selectAllColumnsInPicker, setBoolInputOn, setInputValue, waitForColumns, waitForDialog,
} from './helpers';

// Test Track scenario: EDA/pca.md
// 1. Open cars.csv from Demo files.
// 2. Top Menu > ML > Analyze > PCA...
// 3. Select all Features, set Components = 3.
// 4. Click OK; expect PC1, PC2, PC3 columns added.
// 5. Repeat with Center and Scale; expect PC1 (2), PC2 (2), PC3 (2) additional columns.

const PC_BASE = ['PC1', 'PC2', 'PC3'];

async function runPcaDialog(
  page: import('@playwright/test').Page,
  opts: { components: string; center?: boolean; scale?: boolean },
): Promise<void> {
  await clickTopMenuLeaf(page, 'div-ML---Analyze---PCA...');
  await waitForDialog(page, 'PCA');

  await selectAllColumnsInPicker(page, 'Features');
  await setInputValue(page, 'Components', opts.components);
  if (opts.center) await setBoolInputOn(page, 'Center');
  if (opts.scale) await setBoolInputOn(page, 'Scale');

  await clickDialogPrimary(page, ['OK', 'Run', 'RUN']);
  await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /^PCA$/i }))
    .toHaveCount(0, { timeout: 15_000 });
}

test.describe.serial('EDA / PCA', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('PCA on cars.csv adds PC1/PC2/PC3, then PC1 (2)/PC2 (2)/PC3 (2) with Center+Scale', async ({ page }) => {
    test.setTimeout(300_000);

    await openDemoCsv(page, 'cars.csv');
    const before = (await currentColumnNames(page)).length;

    await runPcaDialog(page, { components: '3' });

    // Step 4 expectation: three new columns PC1/PC2/PC3 appear on the current DF.
    // Wait in-browser (no per-attempt report errors), then assert once.
    await waitForColumns(page, PC_BASE);
    expect(await currentColumnNames(page)).toEqual(expect.arrayContaining(PC_BASE));
    const afterFirst = await currentColumnNames(page);
    expect(afterFirst.length).toBeGreaterThanOrEqual(before + 3);

    // Step 5: repeat with Center + Scale checked. The platform names duplicate PC
    // columns with a `" (2)"` suffix on collision.
    await runPcaDialog(page, { components: '3', center: true, scale: true });

    await waitForColumns(page, ['PC1 (2)', 'PC2 (2)', 'PC3 (2)']);
    expect(await currentColumnNames(page)).toEqual(expect.arrayContaining(['PC1 (2)', 'PC2 (2)', 'PC3 (2)']));
  });
});
