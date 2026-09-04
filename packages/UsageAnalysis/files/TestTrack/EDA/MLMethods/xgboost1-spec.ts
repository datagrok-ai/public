import {test, expect} from '@playwright/test';
import {specTestOptions, softStep, stepErrors, loginToDatagrok} from '../../spec-login';

test.use(specTestOptions);

test('XGBoost 1: Classification on iris.csv', async ({page}) => {
  await loginToDatagrok(page);

  // Step 1: Open iris.csv
  await softStep('Open iris.csv', async () => {
    const result = await page!.evaluate(async () => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });
      Array.from(grok.shell.views).filter(v => v.type === 'PredictiveModel').forEach(v => v.close());
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.windows.simpleMode = false;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1000));
      return {rows: df.rowCount, cols: df.columns.length};
    });
    expect(result.rows).toBe(150);
    expect(result.cols).toBe(6);
  });

  // Step 2: Train XGBoost Classification (JS API fallback)
  await softStep('Train XGBoost Classification (JS API fallback)', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const subDf = df.clone(null, ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width', 'Species']);
      const result = await grok.functions.call('eda:trainXGBooster', {
        df: subDf, predictColumn: subDf.col('Species')
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
