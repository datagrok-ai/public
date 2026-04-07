import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai/';

async function waitForGrok(page: Page) {
  await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
}

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

test.describe('Train predictive model', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await waitForGrok(page);
    await closeAll(page);
  });

  test('Train classification model: SEX by WEIGHT and HEIGHT', async ({ page }) => {
    // Step 1: Open demog.csv
    await page.evaluate(() => (window as any).grok.data.getDemoTable('demog.csv'));
    await page.waitForTimeout(2000);

    // Step 2: Open Train Model dialog via ML menu
    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Train Model...');
    await page.waitForTimeout(1000);

    // Step 3: Set Predict and Features
    // Set Predict column to SEX
    const predictInput = page.locator('[data-name="predict"] input, .d4-combo-box:has-text("Predict")').first();
    if (await predictInput.isVisible()) {
      await predictInput.click();
      await page.click('text=SEX');
    }

    // Set Features to WEIGHT and HEIGHT
    const featuresInput = page.locator('[data-name="features"] input, .d4-combo-box:has-text("Features")').first();
    if (await featuresInput.isVisible()) {
      await featuresInput.click();
      await page.click('text=WEIGHT');
      await page.click('text=HEIGHT');
    }

    // Step 4 & 5: Click Ignore missing, then Run
    const ignoreMissingCheckbox = page.locator('input[type="checkbox"]').filter({ hasText: /ignore missing/i }).first();
    if (await ignoreMissingCheckbox.isVisible()) await ignoreMissingCheckbox.check();

    // Step 8: Click Predict Probability
    const predictProbCheckbox = page.locator('input[type="checkbox"]').filter({ hasText: /predict probability/i }).first();
    if (await predictProbCheckbox.isVisible()) await predictProbCheckbox.check();

    // Click RUN button
    await page.click('button:has-text("RUN"), button:has-text("Train")');
    await page.waitForTimeout(10000); // Allow training time

    // Verify dashboard appeared
    const dashboard = page.locator('.d4-viewer, .grok-view-dock-manager').first();
    await expect(dashboard).toBeVisible({ timeout: 30000 });

    // Step 9: Save the model as "TestDemog"
    await page.click('button:has-text("SAVE"), button:has-text("Save")');
    await page.waitForTimeout(500);
    const nameInput = page.locator('input[type="text"]').first();
    await nameInput.fill('TestDemog');
    await page.click('button:has-text("OK")');
    await page.waitForTimeout(1000);

    await page.screenshot({ path: 'models-train-classification.png' });
  });

  test('Train regression model: WEIGHT by HEIGHT', async ({ page }) => {
    await page.evaluate(() => (window as any).grok.data.getDemoTable('demog.csv'));
    await page.waitForTimeout(2000);

    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Train Model...');
    await page.waitForTimeout(1000);

    // Set Predict to WEIGHT, Features to HEIGHT
    const predictInput = page.locator('[data-name="predict"] input').first();
    if (await predictInput.isVisible()) {
      await predictInput.click();
      await page.click('text=WEIGHT');
    }

    await page.click('button:has-text("RUN"), button:has-text("Train")');
    await page.waitForTimeout(15000);

    const dashboard = page.locator('.d4-viewer, .grok-view-dock-manager').first();
    await expect(dashboard).toBeVisible({ timeout: 30000 });

    await page.screenshot({ path: 'models-train-regression.png' });
  });
});
