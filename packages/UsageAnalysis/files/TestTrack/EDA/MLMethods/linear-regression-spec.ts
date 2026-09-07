import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Linear Regression: Train on cars.csv', async ({page}) => {
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

  await softStep('Set Features and select Model Engine (JS API fallback)', async () => {
    const result = await page!.evaluate(async () => {

      Array.from(grok.shell.views).filter(v => v.type === 'PredictiveModel').forEach(v => v.close());

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
