import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('XGBoost 2: Regression on cars.csv', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);

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
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      return {rows: df.rowCount, cols: df.columns.length};
    });
    expect(result.rows).toBe(30);
    expect(result.cols).toBe(17);
  });

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
