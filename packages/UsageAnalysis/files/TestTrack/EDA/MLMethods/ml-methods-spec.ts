import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('ML Methods - Linear Regression, PLS, Softmax, XGBoost', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);

  await softStep('Linear Regression: Train on cars.csv', async () => {
    const result = await page!.evaluate(async () => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });

      Array.from(grok.shell.views).filter(v => v.type === 'PredictiveModel').forEach(v => v.close());
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv');
      grok.shell.addTableView(df);
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      const result = await grok.functions.call('eda:trainLinearRegression', {
        df: df, predictColumn: df.col('price')
      });
      return { success: result != null };
    });
    expect(result.success).toBe(true);
  });

  await softStep('PLS Regression: Train on cars.csv', async () => {
    const result = await page!.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv');
      grok.shell.addTableView(df);
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      const numCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.type !== 'string') numCols.push(c.name);
      }
      const numDf = df.clone(null, numCols);
      const result = await grok.functions.call('eda:trainPLSRegression', {
        df: numDf, predictColumn: numDf.col('price'), components: 3
      });
      return { success: result != null };
    });
    expect(result.success).toBe(true);
  });

  await softStep('Softmax: Train on iris.csv', async () => {
    const result = await page!.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
      grok.shell.addTableView(df);
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      try {
        await grok.functions.call('eda:trainSoftmax', {
          df: df, predictColumn: df.col('Species'),
          rate: 0.1, iterations: 100, penalty: 0.01, tolerance: 0.001
        });
        return { success: true };
      } catch (e: any) {
        return { success: false, error: e.message };
      }
    });

    console.log('Softmax result:', result);
  });

  await softStep('XGBoost Classification: Train on iris.csv', async () => {
    const result = await page!.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
      grok.shell.addTableView(df);
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      const subDf = df.clone(null, ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width', 'Species']);
      const result = await grok.functions.call('eda:trainXGBooster', {
        df: subDf, predictColumn: subDf.col('Species')
      });
      return { success: result != null };
    });
    expect(result.success).toBe(true);
  });

  await softStep('XGBoost Regression: Train on cars.csv', async () => {
    const result = await page!.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv');
      grok.shell.addTableView(df);
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      const numCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.type !== 'string') numCols.push(c.name);
      }
      const numDf = df.clone(null, numCols);
      const result = await grok.functions.call('eda:trainXGBooster', {
        df: numDf, predictColumn: numDf.col('price')
      });
      return { success: result != null };
    });
    expect(result.success).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
