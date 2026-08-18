import { test, expect } from './helpers/diff-studio';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import {
  BASE, openDiffStudio, openModelFromLibrary, toggleRibbonSwitch,
  ribbonSwitchOn, setInputValue, inputEditor, inputHost,
} from './helpers/diff-studio';

test('DiffStudio Scripting — Edit toggle, </> JS view, Run, Save with //tags: model, Model Hub', async ({ page }) => {
  test.setTimeout(300_000);
  const { softStep, assertAllPassed } = createSoftStepCollector();
  const monitor = attachErrorMonitor(page);

  await softStep('Setup: Open Diff Studio + Bioreactor', async () => {
    await openDiffStudio(page);
    await openModelFromLibrary(page, 'Bioreactor');
  });

  await softStep('Step 1: Enable Edit toggle, then click </> to open the JS script view', async () => {
    const wasOn = await ribbonSwitchOn(page, 'Edit');
    if (!wasOn) await toggleRibbonSwitch(page, 'Edit');
    const editorVisible = await page.locator('.diff-studio-eqs-editor').count();
    expect(editorVisible).toBeGreaterThan(0);

    await page.locator('.d4-ribbon-name', { hasText: '</>' }).first().click();

    await page.locator('.CodeMirror').first().waitFor({ timeout: 30_000 });
    await page.waitForTimeout(1500);
  });

  await softStep('Step 2: Add "//tags: model" to the script body and Save', async () => {

    const cm = page.locator('.CodeMirror').first();
    await cm.waitFor({ timeout: 10_000 });
    await cm.click();
    await page.keyboard.press('Control+Home'); 
    await page.keyboard.press('ArrowDown');    
    await page.keyboard.press('End');
    await page.keyboard.press('Enter');
    await page.keyboard.type('//tags: model');
    await page.waitForTimeout(500);

    const tagLine = page.locator('.CodeMirror-line', { hasText: '//tags: model' });
    await expect(tagLine.first()).toBeVisible({ timeout: 5_000 });

    await page.locator('[name="button-Save"]').first().click();
    await page.waitForTimeout(3000);
  });

  await softStep('Step 3: Run the script; move Final-at slider — live update', async () => {

    await page.locator('[name="icon-play"]').first().click();

    await page.locator(inputHost('Final')).waitFor({ timeout: 30_000 });
    await page.waitForTimeout(1500);

    const inp = page.locator(inputEditor('Final'));
    const before = await inp.inputValue();
    await setInputValue(page, 'Final', '500');
    await page.waitForTimeout(2500);
    const after = await inp.inputValue();
    expect(after).toBe('500');
    expect(after).not.toBe(before);

    await expect(page.locator(inputHost('Process-mode'))).toHaveCount(0);
  });

  await softStep('Step 4: Access the saved model in Model Hub (Apps > Compute > Model Hub)', async () => {

    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
    await page.waitForTimeout(1500);

    const clickTreeLabel = async (label: string): Promise<boolean> => {
      return await page.evaluate((text) => {
        const candidates = Array.from(document.querySelectorAll(
          '.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item-label')) as HTMLElement[];
        const el = candidates.find(e => e.textContent?.trim() === text);
        if (!el) return false;
        el.scrollIntoView({ behavior: 'instant', block: 'center' });
        el.click();
        return true;
      }, label);
    };

    expect(await clickTreeLabel('Apps')).toBe(true);
    await page.waitForTimeout(800);
    expect(await clickTreeLabel('Compute')).toBe(true);
    await page.waitForTimeout(800);
    expect(await clickTreeLabel('Model Hub')).toBe(true);
    await page.waitForTimeout(3000);

    await expect(page.getByText('Bioreactor', { exact: true }).first())
      .toBeVisible({ timeout: 20_000 });
  });

  await softStep('Step 5: Interact with the saved model; adjust Final', async () => {

    const items = page.getByText('Bioreactor', { exact: true });
    const count = await items.count();
    expect(count).toBeGreaterThan(0);
    await items.nth(count - 1).dblclick();

    await page.locator(inputHost('Final')).waitFor({ timeout: 30_000 });
    await setInputValue(page, 'Final', '800');
    await page.waitForTimeout(2500);
    const after = await page.locator(inputEditor('Final')).inputValue();
    expect(after).toBe('800');

  });

  assertAllPassed();
  monitor.assertNone();
});
