import {test, expect, chromium} from '@playwright/test';
import {specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

test('XGBoost 2: Regression on cars.csv', async () => {
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

  // Step 2: Train XGBoost Regression (JS API fallback)
  await softStep('Train XGBoost Regression (JS API fallback)', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const numCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.type !== 'string') numCols.push(c.name);
      }
      const numDf = df.clone(null, numCols);
      const result = await grok.functions.call('eda:trainXGBooster', {
        df: numDf, predictColumn: numDf.col('price')
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
