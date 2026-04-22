import { test, expect, Page } from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const BASE_URL = 'https://public.datagrok.ai/';

async function waitForGrok(page: Page) {
  await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
}

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

async function trainModel(page: Page, opts: {
  file: string;
  predict: string;
  features: string[];
  engine: string;
  components?: number;
  modelName: string;
}) {
  await page.evaluate((f) => (window as any).grok.data.getDemoTable(f), opts.file);
  await page.waitForTimeout(2000);

  await page.click('text=ML');
  await page.click('text=Models');
  await page.click('text=Train Model...');
  await page.waitForTimeout(1000);

  // Engine selection
  const engineInput = page.locator('[data-name="engine"] input, [data-name="modelEngine"] input').first();
  if (await engineInput.isVisible()) {
    await engineInput.click();
    await page.click(`text=${opts.engine}`);
  }

  // Components (for PLS)
  if (opts.components != null) {
    const compInput = page.locator('[data-name="components"] input, input[value="1"]').first();
    if (await compInput.isVisible()) await compInput.fill(String(opts.components));
  }

  // Click RUN
  await page.click('button:has-text("RUN"), button:has-text("Train")');
  await page.waitForTimeout(15000);

  // Save
  await page.click('button:has-text("SAVE"), button:has-text("Save")');
  await page.waitForTimeout(500);
  const nameInput = page.locator('input[type="text"]').first();
  await nameInput.fill(opts.modelName);
  await page.click('button:has-text("OK")');
  await page.waitForTimeout(1000);
}

test.describe('Predictive models — full lifecycle with accelerometer data', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await waitForGrok(page);
    await closeAll(page);
  });

  test('Part 1: Train PLS and Linear Regression on accelerometer.csv', async ({ page }) => {
    // Train PLS model
    await trainModel(page, {
      file: 'sensors/accelerometer.csv',
      predict: 'accel_y',
      features: ['accel_z', 'time_offset'],
      engine: 'EDA: PLS',
      components: 3,
      modelName: 'Accelerometer_model_PLS',
    });

    await page.screenshot({ path: 'accel-pls-dashboard.png' });

    // Switch to Linear Regression and retrain
    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Train Model...');
    await page.waitForTimeout(1000);

    const engineInput = page.locator('[data-name="engine"] input, [data-name="modelEngine"] input').first();
    if (await engineInput.isVisible()) {
      await engineInput.click();
      await page.click('text=Linear Regression, text=Eda:Linear');
    }

    await page.click('button:has-text("RUN"), button:has-text("Train")');
    await page.waitForTimeout(15000);

    await page.click('button:has-text("SAVE"), button:has-text("Save")');
    const nameInput = page.locator('input[type="text"]').first();
    await nameInput.fill('Accelerometer_model_LR');
    await page.click('button:has-text("OK")');
    await page.waitForTimeout(1000);

    await page.screenshot({ path: 'accel-lr-dashboard.png' });
  });

  test('Part 2: Apply PLS and LR models to accelerometer.csv', async ({ page }) => {
    await page.evaluate(() => (window as any).grok.data.getDemoTable('sensors/accelerometer.csv'));
    await page.waitForTimeout(2000);

    const initialColCount = await page.evaluate(() =>
      (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0);

    // Apply PLS
    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Apply Model...');
    await page.waitForTimeout(1000);

    const plsModel = page.locator('text=Accelerometer_model_PLS').first();
    if (await plsModel.isVisible()) await plsModel.click();
    await page.click('button:has-text("OK")');
    await page.waitForTimeout(5000);

    const afterPLS = await page.evaluate(() =>
      (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0);
    expect(afterPLS).toBeGreaterThan(initialColCount);

    // Apply LR
    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Apply Model...');
    await page.waitForTimeout(1000);

    const lrModel = page.locator('text=Accelerometer_model_LR').first();
    if (await lrModel.isVisible()) await lrModel.click();
    await page.click('button:has-text("OK")');
    await page.waitForTimeout(5000);

    const afterLR = await page.evaluate(() =>
      (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0);
    expect(afterLR).toBeGreaterThan(afterPLS);

    await page.screenshot({ path: 'accel-apply-result.png' });
  });

  test('Part 3: Apply LR model to random walk test dataset', async ({ page }) => {
    // Open test dataset via Tools > Dev > Open test dataset
    await page.click('text=Tools');
    await page.click('text=Dev');
    await page.click('text=Open test dataset');
    await page.waitForTimeout(500);

    const rowsInput = page.locator('input[placeholder*="rows"], input').nth(0);
    await rowsInput.fill('1000');
    const colsInput = page.locator('input[placeholder*="cols"], input').nth(1);
    await colsInput.fill('10');

    // Select "random walk" option if available
    const typeSelect = page.locator('select, .d4-combo-box').last();
    if (await typeSelect.isVisible()) {
      await typeSelect.click();
      const randomWalk = page.locator('text=random walk').first();
      if (await randomWalk.isVisible()) await randomWalk.click();
    }

    await page.click('button:has-text("OK")');
    await page.waitForTimeout(2000);

    const initialColCount = await page.evaluate(() =>
      (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0);

    // Apply LR model
    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Apply Model...');
    await page.waitForTimeout(1000);

    const lrModel = page.locator('text=Accelerometer_model_LR').first();
    if (await lrModel.isVisible()) await lrModel.click();
    await page.click('button:has-text("OK")');
    await page.waitForTimeout(5000);

    const afterApply = await page.evaluate(() =>
      (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0);
    expect(afterApply).toBeGreaterThan(initialColCount);

    await page.screenshot({ path: 'accel-apply-random-walk.png' });
  });

  test('Part 4: Delete both Accelerometer models', async ({ page }) => {
    await page.goto(`${BASE_URL}models`);
    await page.waitForTimeout(2000);

    for (const modelName of ['Accelerometer_model_LR', 'Accelerometer_model_PLS']) {
      const modelItem = page.locator(`text=${modelName}`).first();
      if (!(await modelItem.isVisible({ timeout: 3000 }))) continue;

      // Check context panel tabs
      await modelItem.click();
      await page.waitForTimeout(500);

      // Delete
      await modelItem.click({ button: 'right' });
      await page.waitForTimeout(500);
      await page.click('[class*="menu"] text=Delete, .d4-menu-item:has-text("Delete")');
      await page.waitForTimeout(500);
      await page.click('button:has-text("Delete"), button:has-text("OK")');
      await page.waitForTimeout(1500);
    }

    // Verify both models are gone
    for (const modelName of ['Accelerometer_model_LR', 'Accelerometer_model_PLS']) {
      await expect(page.locator(`text=${modelName}`)).not.toBeVisible({ timeout: 3000 });
    }

    await page.screenshot({ path: 'accel-models-deleted.png' });
  });
});
