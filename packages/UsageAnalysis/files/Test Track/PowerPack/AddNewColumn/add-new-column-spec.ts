import { test, expect, Page } from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const BASE_URL = 'https://release-ec2.datagrok.ai/';

async function login(page: Page) {
  await page.goto(BASE_URL);
  await page.fill('input[placeholder="Login"]', 'claude');
  await page.fill('input[placeholder="Password"]', 'grokclaude');
  await page.click('button:has-text("LOGIN")');
  await page.waitForSelector('.grok-app', { timeout: 30000 });
}

async function openDemog(page: Page) {
  await page.evaluate(async () => {
    const df = await (window as any).grok.data.getDemoTable('demog.csv');
    (window as any).grok.shell.addTableView(df);
  });
  await page.waitForTimeout(2000);
}

async function openAddNewColumnDialog(page: Page) {
  await page.evaluate(() => {
    (window as any).grok.shell.topMenu.find('Edit').find('Add New Column...').click();
  });
  await page.waitForSelector('.d4-dialog', { timeout: 5000 });
}

async function addCalculatedColumn(page: Page, name: string, formula: string) {
  await openAddNewColumnDialog(page);
  const nameInput = page.locator('.d4-dialog input[type="text"]').first();
  await nameInput.clear();
  await nameInput.fill(name);
  const editor = page.locator('.d4-dialog .cm-content');
  await editor.click();
  await page.keyboard.type(formula);
  await page.waitForTimeout(1000);
  await page.click('.d4-dialog button:has-text("OK")');
  await page.waitForTimeout(2000);
}

test.describe('PowerPack: Add New Column (Multi-source)', () => {
  test.beforeEach(async ({ page }) => {
    await login(page);
    await page.evaluate(() => (window as any).grok.shell.closeAll());
  });

  test('Steps 1-4: Open demog, add formula column, add dependent column', async ({ page }) => {
    await openDemog(page);

    // Step 3: Add a formula column
    await addCalculatedColumn(page, 'CalcCol1', '${HEIGHT} + ${WEIGHT}');

    const col1Exists = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.dataFrame.columns.contains('CalcCol1');
    });
    expect(col1Exists).toBe(true);

    // Step 4: Add another column referencing the first
    await addCalculatedColumn(page, 'CalcCol2', '${CalcCol1} * 2');

    const col2Exists = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.dataFrame.columns.contains('CalcCol2');
    });
    expect(col2Exists).toBe(true);
  });

  test('Step 5: Change column values — dependent columns update', async ({ page }) => {
    await openDemog(page);
    await addCalculatedColumn(page, 'CalcCol1', '${HEIGHT} + ${WEIGHT}');
    await addCalculatedColumn(page, 'CalcCol2', '${CalcCol1} * 2');

    // Get initial values
    const initial = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      return {
        calc1: df.col('CalcCol1').get(0),
        calc2: df.col('CalcCol2').get(0),
      };
    });

    // Modify HEIGHT value
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.col('HEIGHT').set(0, 200.0);
      df.fireValuesChanged();
    });
    await page.waitForTimeout(2000);

    // Verify propagation
    const updated = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      return {
        calc1: df.col('CalcCol1').get(0),
        calc2: df.col('CalcCol2').get(0),
      };
    });

    expect(updated.calc1).not.toEqual(initial.calc1);
    expect(updated.calc2).toEqual(updated.calc1 * 2);
  });
});
