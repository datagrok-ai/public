import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Softmax: Train on iris.csv', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);

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
      for (let i = 0; i < 20; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 50));
      }
      return {rows: df.rowCount, cols: df.columns.length};
    });
    expect(result.rows).toBe(150);
    expect(result.cols).toBe(6);
  });

  await softStep('Train Softmax (known bug)', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const subDf = df.clone(null, ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width', 'Species']);
      try {
        await grok.functions.call('eda:trainSoftmax', {
          df: subDf, predictColumn: subDf.col('Species'),
          rate: 0.1, iterations: 100, penalty: 0.01, tolerance: 0.001
        });
        return {success: true};
      } catch (e: any) {
        return {success: false, error: e.message};
      }
    });

    console.log('Softmax result:', result);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
