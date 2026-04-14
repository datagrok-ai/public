import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('DiffStudio Files & Sharing: Open pk.ivp, modify inputs, share URL', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Setup
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;
  });

  // Step 1: Navigate to library and open pk.ivp
  await softStep('Step 1: Navigate to library and open pk.ivp', async () => {
    await page!.goto(`${baseUrl}/files/system.appdata/DiffStudio/library`, {
      waitUntil: 'networkidle', timeout: 30000,
    });
    await page!.waitForTimeout(3000);

    // Click pk.ivp in the file list
    const pkLink = page!.locator('text=pk.ivp').first();
    await pkLink.click();
    await page!.waitForTimeout(5000);

    // Verify DiffStudio opened with PK model
    await expect(page!.locator('text=Dosing')).toBeVisible({timeout: 10000});

    const info = await page!.evaluate(() => {
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      return {rows: df?.rowCount, cols: df?.columns?.length, name: tv?.name};
    });
    expect(info.rows).toBe(1201);
    expect(info.cols).toBe(3);
  });

  // Step 2: Modify Step to 0.1 and Count to 4
  await softStep('Step 2: Set Step to 0.1 and Count to 4', async () => {
    // Change step: find input by current value, use native setter + change event
    await page!.evaluate(async () => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const stepInput = inputs.find(i => i.value === '0.01');
      if (stepInput) {
        const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        nativeSetter.call(stepInput, '0.1');
        stepInput.dispatchEvent(new Event('input', {bubbles: true}));
        stepInput.dispatchEvent(new Event('change', {bubbles: true}));
        stepInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 2000));

      // Change count
      const inputs2 = Array.from(document.querySelectorAll('input'));
      const countInput = inputs2.find(i => i.value === '1');
      if (countInput) {
        const nativeSetter2 = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        nativeSetter2.call(countInput, '4');
        countInput.dispatchEvent(new Event('input', {bubbles: true}));
        countInput.dispatchEvent(new Event('change', {bubbles: true}));
        countInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 3000));
    });

    // Verify URL updated
    const url = page!.url();
    expect(url).toContain('step=0.1');
    expect(url).toContain('count=4');

    // Verify row count changed
    const rows = await page!.evaluate(() => grok.shell.tv?.dataFrame?.rowCount);
    expect(rows).toBe(484);
  });

  // Step 3: Sharing — open URL in new tab, verify same model loads
  await softStep('Step 3: Open shared URL in new tab', async () => {
    const sharedUrl = page!.url();
    expect(sharedUrl).toContain('params:');
    expect(sharedUrl).toContain('step=0.1');
    expect(sharedUrl).toContain('count=4');

    // Open in new tab within the same context
    const newPage = await context.newPage();
    await newPage.goto(sharedUrl, {waitUntil: 'networkidle', timeout: 30000});
    await newPage.waitForTimeout(5000);

    // Verify same model loaded with same inputs
    await expect(newPage.locator('text=Dosing')).toBeVisible({timeout: 10000});

    const info = await newPage.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const stepVal = inputs.find(i => i.value === '0.1' || i.value === '0.10')?.value;
      const countVal = inputs.find(i => i.value === '4')?.value;
      const rows = grok.shell.tv?.dataFrame?.rowCount;
      return {stepVal, countVal, rows};
    });
    expect(info.stepVal).toBeTruthy();
    expect(info.countVal).toBe('4');
    expect(info.rows).toBe(484);

    await newPage.close();
  });

  // Step 4: REMARK — no Multiaxis/Facet for 2-curve model
  await softStep('Step 4 (REMARK): Only linechart for 2-curve model', async () => {
    const hasMultiaxis = await page!.evaluate(() => {
      return Array.from(document.querySelectorAll('*')).some(el =>
        el.textContent?.trim() === 'Multiaxis' && el.children.length === 0
      );
    });
    expect(hasMultiaxis).toBe(false);

    const canvasCount = await page!.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(1);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
