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
  return await page.evaluate(async () => {
    const df = await (window as any).grok.data.getDemoTable('demog.csv');
    (window as any).grok.shell.addTableView(df);
    return df.rowCount;
  });
}

async function addCalculatedColumn(page: Page, name: string, formula: string) {
  await page.evaluate(() => {
    (window as any).grok.shell.topMenu.find('Edit').find('Add New Column...').click();
  });
  await page.waitForSelector('.d4-dialog', { timeout: 5000 });

  // Set column name
  const nameInput = page.locator('.d4-dialog input[type="text"]').first();
  await nameInput.clear();
  await nameInput.fill(name);

  // Set formula
  const editor = page.locator('.d4-dialog .cm-content');
  await editor.click();
  await page.keyboard.type(formula);
  await page.waitForTimeout(1000);

  // Click OK
  await page.click('.d4-dialog button:has-text("OK")');
  await page.waitForTimeout(2000);
}

test.describe('PowerPack: Formula Refreshing', () => {
  test.beforeEach(async ({ page }) => {
    await login(page);
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await openDemog(page);
  });

  test('Steps 2-4: Create chain of calculated columns', async ({ page }) => {
    // Step 2: Weight2 = ${WEIGHT} + 100
    await addCalculatedColumn(page, 'Weight2', '${WEIGHT} + 100');
    let colExists = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.dataFrame.columns.contains('Weight2');
    });
    expect(colExists).toBe(true);

    // Step 3: Weight3 = ${Weight2} + 100
    await addCalculatedColumn(page, 'Weight3', '${Weight2} + 100');
    colExists = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.dataFrame.columns.contains('Weight3');
    });
    expect(colExists).toBe(true);

    // Step 4: Weight4 = Log10(${Weight3}) - 0.2
    await addCalculatedColumn(page, 'Weight4', 'Log10(${Weight3}) - 0.2');
    colExists = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.dataFrame.columns.contains('Weight4');
    });
    expect(colExists).toBe(true);
  });

  test('Step 5: Verify formula dependency propagation', async ({ page }) => {
    // Create the chain first
    await addCalculatedColumn(page, 'Weight2', '${WEIGHT} + 100');
    await addCalculatedColumn(page, 'Weight3', '${Weight2} + 100');
    await addCalculatedColumn(page, 'Weight4', 'Log10(${Weight3}) - 0.2');

    // Get initial values
    const initialValues = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      return {
        weight2: df.col('Weight2').get(0),
        weight3: df.col('Weight3').get(0),
        weight4: df.col('Weight4').get(0),
      };
    });

    // Modify Weight2 value directly
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.col('Weight2').set(0, 200.0);
      df.fireValuesChanged();
    });
    await page.waitForTimeout(2000);

    // Verify propagation
    const updatedValues = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      return {
        weight3: df.col('Weight3').get(0),
        weight4: df.col('Weight4').get(0),
      };
    });

    // Weight3 should be 200 + 100 = 300
    expect(updatedValues.weight3).toBeCloseTo(300, 0);
    // Weight4 should be Log10(300) - 0.2 ≈ 2.277
    expect(updatedValues.weight4).toBeCloseTo(2.277, 2);

    // Cleanup
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.columns.remove('Weight4');
      df.columns.remove('Weight3');
      df.columns.remove('Weight2');
    });
  });
});
