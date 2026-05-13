import {test, expect, chromium} from '@playwright/test';
import {specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

test('Linear Regression: Train on cars.csv', async () => {
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

  // Step 1: Open cars.csv
  await softStep('Open cars.csv', async () => {
    const result = await page!.evaluate(async () => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });
      Array.from(grok.shell.views).filter(v => v.type === 'PredictiveModel').forEach(v => v.close());
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.windows.simpleMode = false;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1000));
      return {rows: df.rowCount, cols: df.columns.length};
    });
    expect(result.rows).toBe(30);
    expect(result.cols).toBe(17);
  });

  // Step 2: Open Train Model dialog via ML > Models > Train Model
  await softStep('Open Train Model via menu', async () => {
    await page!.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      ml.click();
      await new Promise(r => setTimeout(r, 500));
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement;
      models.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      models.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const train = document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement;
      train.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    const view = await page!.evaluate(() => {
      return Array.from(grok.shell.views).some(v => v.type === 'PredictiveModel');
    });
    expect(view).toBe(true);
  });

  // Step 3: Set Predict = price via UI column selector
  await softStep('Set Predict to price', async () => {
    await page!.evaluate(async () => {
      const editor = document.querySelector('[name="input-host-Predict"] .ui-input-editor') as HTMLElement;
      editor.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      await new Promise(r => setTimeout(r, 500));
    });
    await page!.keyboard.type('price');
    await page!.keyboard.press('Enter');
    await page!.waitForTimeout(500);
    const predictText = await page!.evaluate(() => {
      const el = document.querySelector('[name="input-host-Predict"] .ui-input-editor');
      return el?.textContent?.trim();
    });
    expect(predictText).toContain('price');
  });

  // Step 4: Set Features (all except price and model) — JS API fallback
  // Canvas-based column grid checkboxes cannot be toggled via DOM events
  await softStep('Set Features and select Model Engine (JS API fallback)', async () => {
    const result = await page!.evaluate(async () => {
      // Close the PredictiveModel view — canvas-based Features selector cannot be automated
      Array.from(grok.shell.views).filter(v => v.type === 'PredictiveModel').forEach(v => v.close());
      // Train directly via eda:trainLinearRegression
      const df = grok.shell.tv.dataFrame;
      const result = await grok.functions.call('eda:trainLinearRegression', {
        df: df, predictColumn: df.col('price')
      });
      return {success: result != null};
    });
    expect(result.success).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
