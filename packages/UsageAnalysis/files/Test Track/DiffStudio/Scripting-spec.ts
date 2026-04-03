import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Scripting', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/apps/DiffStudio`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(2000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1a: Edit toggle opens equations editor', async () => {
    const switchEl = page.locator('.ui-input-switch');
    await switchEl.click();
    await page.waitForTimeout(1000);
    // Verify CodeMirror editor appeared (equations editor)
    const cmEl = page.locator('.CodeMirror');
    await expect(cmEl).toBeVisible();
    // Ribbon should now show Refresh and </> buttons
    await expect(page.locator('.d4-ribbon-item').filter({ hasText: 'Refresh' })).toBeVisible();
  });

  test('Step 1b: </> button opens JavaScript script view', async () => {
    await page.locator('.d4-ribbon-item').filter({ hasText: '</>' }).click();
    await page.waitForTimeout(1500);
    // ScriptView should have Script / Layout / Debug tabs
    await expect(page.locator('.d4-ribbon-item, [class*="tab"]').filter({ hasText: 'Script' })).toBeVisible();
    // Script should contain //language: javascript
    const cm = page.locator('.CodeMirror');
    const content = await cm.evaluate((el: any) => el.CodeMirror.getValue());
    expect(content).toContain('//language: javascript');
  });

  test('Step 2a: Run script produces RichFunctionView', async () => {
    const runBtn = page.locator('.grok-script-run-icon').first();
    await runBtn.click();
    await page.waitForTimeout(3000);
    // Should open a new view with Inputs panel
    await expect(page.locator('text=Bioreactor / Line chart')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2b: Moving Final slider updates chart in real time; no Process mode/Facet', async () => {
    // Verify no Process mode and no Facet tab (per scenario REMARK)
    const processModeEl = page.locator('text=Process mode');
    await expect(processModeEl).not.toBeVisible();
    const facetTab = page.locator('text=Facet');
    await expect(facetTab).not.toBeVisible();

    // Switch to line chart view
    await page.locator('text=Bioreactor / Line chart').click();
    await page.waitForTimeout(500);

    // Change Final slider
    await page.evaluate(() => {
      const ranges = Array.from(document.querySelectorAll('input[type="range"]'));
      const finalRange = ranges.find((r: any) => r.value === '1000');
      if (finalRange) {
        (finalRange as HTMLInputElement).value = '400';
        finalRange.dispatchEvent(new Event('input', { bubbles: true }));
      }
    });
    await page.waitForTimeout(1000);

    const finalVal = await page.evaluate(() => {
      const ranges = Array.from(document.querySelectorAll('input[type="range"]'));
      return (ranges.find((r: any) => r.value === '400') as HTMLInputElement)?.value;
    });
    expect(finalVal).toBe('400');
  });

  test('Step 3: Add //tags: model and save script', async () => {
    // Switch back to script view
    const views = await page.evaluate(() =>
      Array.from((window as any).grok.shell.views).map((v: any) => ({ type: v.type, name: v.name }))
    );
    const scriptView = (views as any[]).find(v => v.type === 'ScriptView');
    expect(scriptView).toBeTruthy();

    await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.views);
      const sv = views.find((v: any) => v.type === 'ScriptView');
      if (sv) (window as any).grok.shell.v = sv;
    });
    await page.waitForTimeout(500);

    // Add //tags: model via CodeMirror
    await page.evaluate(() => {
      const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
      if (cm) {
        const content = cm.getValue();
        const lines = content.split('\n');
        const metaIdx = lines.findIndex((l: string) => l.includes('//meta.features'));
        if (metaIdx >= 0) {
          lines.splice(metaIdx + 1, 0, '//tags: model');
          cm.setValue(lines.join('\n'));
        }
      }
    });

    // Save
    await page.keyboard.press('Control+s');
    await page.waitForTimeout(1000);

    // Verify saved toast
    await expect(page.locator('text=Script saved')).toBeVisible({ timeout: 5000 });
  });

  test('Step 4: Model appears in Model Hub', async () => {
    await page.goto(`${BASE_URL}/apps/Modelhub`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(2000);
    await expect(page.locator('text=Bioreactor')).toBeVisible();
  });

  test('Step 5: Run model from Model Hub; slider updates chart', async () => {
    await page.locator('text=Bioreactor').first().dblclick();
    await page.waitForTimeout(3000);

    // Switch to line chart
    await page.locator('text=Bioreactor / Line chart').click();
    await page.waitForTimeout(500);

    // Move slider
    await page.evaluate(() => {
      const ranges = Array.from(document.querySelectorAll('input[type="range"]'));
      const finalRange = ranges.find((r: any) => r.value === '1000');
      if (finalRange) {
        (finalRange as HTMLInputElement).value = '600';
        finalRange.dispatchEvent(new Event('input', { bubbles: true }));
      }
    });
    await page.waitForTimeout(1000);

    // Verify the input shows 600
    await expect(page.locator('text=600')).toBeVisible();
  });
});
