import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const MODEL_NAME_CLS = 'TestDemog';
const MODEL_NAME_REG = 'TestDemog2';

test('Train predictive model: classification + regression', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // Pre-cleanup: delete leftover models from previous runs so save dialog doesn't conflict.
  await page.evaluate(async ({names}) => {
    const g: any = (window as any).grok;
    for (const n of names) {
      const list = await g.dapi.models.filter(`friendlyName = "${n}"`).list();
      for (const m of list)
        await g.dapi.models.delete(m);
    }
  }, {names: [MODEL_NAME_CLS, MODEL_NAME_REG]});

  // Setup: open demog.csv, selenium class / Tabs mode, wait for grid
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('Open demog.csv', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBeGreaterThan(0);
  });

  await softStep('ML > Models > Train Model... opens Predictive model view', async () => {
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('Models menu item not found');
      const r = models.getBoundingClientRect();
      const ev = (type: string) => new MouseEvent(type, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Train-Model..."]').click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel', null, {timeout: 15_000});
  });

  await softStep('Set Predict to SEX', async () => {
    await page.locator('[name="input-host-Predict"] [name="div-Predict"]').click();
    await page.keyboard.type('SEX');
    await page.keyboard.press('Enter');
    await expect(page.locator('[name="input-host-Predict"] .d4-column-selector-column')).toHaveText('SEX', {timeout: 10_000});
  });

  await softStep('Set Features to HEIGHT and WEIGHT', async () => {
    await page.locator('[name="input-host-Features"] [name="div-Features"]').click();
    await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
    // Canvas-based column grid — click checkbox column at the right pixel offsets.
    // Order of rows matches DataFrame: USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT(5), WEIGHT(6).
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
      if (!overlay) throw new Error('column-grid overlay not found');
      const click = (x: number, y: number) => {
        const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: x, clientY: y};
        overlay.dispatchEvent(new MouseEvent('mousedown', opts));
        overlay.dispatchEvent(new MouseEvent('mouseup', opts));
        overlay.dispatchEvent(new MouseEvent('click', opts));
      };
      click(826, 397); // HEIGHT
      click(826, 425); // WEIGHT
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('2');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names')).toContainText('HEIGHT', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names')).toContainText('WEIGHT');
  });

  await softStep('Select Impute missing — k-NN Imputation dialog opens', async () => {
    await page.locator('[name="input-host-Impute-missing"] input[type="checkbox"]').click();
    await page.locator('[name="dialog-k-NN-Imputation"]').waitFor({timeout: 10_000});
    await expect(page.locator('[name="dialog-k-NN-Imputation"] .d4-dialog-title')).toHaveText('k-NN Imputation');
  });

  await softStep('Click RUN — model trains and dashboard renders', async () => {
    await page.locator('[name="dialog-k-NN-Imputation"] [name="button-RUN"]').click();
    await page.locator('[name="dialog-k-NN-Imputation"]').waitFor({state: 'detached', timeout: 60_000});
    // Dashboard: confusion matrix + performance metrics
    await page.waitForFunction(() => /Confusions|Performance|Sensitivity/.test(document.body.innerText), null, {timeout: 60_000});
  });

  await softStep('Unselect Impute missing', async () => {
    await page.locator('[name="input-host-Impute-missing"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-Impute-missing"] input[type="checkbox"]')).not.toBeChecked();
  });

  await softStep('Select Ignore missing', async () => {
    await page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]')).toBeChecked();
    await page.waitForFunction(() => /Confusions|Performance/.test(document.body.innerText), null, {timeout: 60_000});
  });

  await softStep('Select Predict Probability — engine becomes regression', async () => {
    await page.locator('[name="input-host-Predict-probability"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-Predict-probability"] input[type="checkbox"]')).toBeChecked();
    // Regression-style ROC curve and Positive class cutoff appear
    await page.waitForFunction(() => /ROC Curve|Positive class cutoff/.test(document.body.innerText), null, {timeout: 60_000});
  });

  await softStep(`Save model as ${MODEL_NAME_CLS}`, async () => {
    await page.locator('button[name="button-Save"], button:has-text("SAVE")').first().click();
    const dlg = page.locator('.d4-dialog [name="input-host-Name"] input');
    await dlg.waitFor({timeout: 10_000});
    await dlg.focus();
    await page.keyboard.type(MODEL_NAME_CLS);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(async (name: string) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length > 0;
    }, MODEL_NAME_CLS, {timeout: 30_000});
  });

  await softStep('Repeat: train regression model — Predict WEIGHT by HEIGHT', async () => {
    // Change Predict to WEIGHT
    await page.locator('[name="input-host-Predict"] [name="div-Predict"]').click();
    await page.keyboard.type('WEIGHT');
    await page.keyboard.press('Enter');
    await expect(page.locator('[name="input-host-Predict"] .d4-column-selector-column')).toHaveText('WEIGHT', {timeout: 10_000});

    // Update Features: leave only HEIGHT (uncheck WEIGHT after dialog reorders)
    await page.locator('[name="input-host-Features"] [name="div-Features"]').click();
    await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
      const click = (x: number, y: number) => {
        const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: x, clientY: y};
        overlay.dispatchEvent(new MouseEvent('mousedown', opts));
        overlay.dispatchEvent(new MouseEvent('mouseup', opts));
        overlay.dispatchEvent(new MouseEvent('click', opts));
      };
      // After reopen, checked rows (HEIGHT, WEIGHT) move to the top: y≈257 (HEIGHT), y≈285 (WEIGHT)
      click(826, 285); // uncheck WEIGHT
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('1');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names')).toHaveText(/^\(1\)\s*HEIGHT/, {timeout: 10_000});

    // Re-enable Ignore missing (cleared on input change) and wait for regression dashboard
    await page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]').click();
    await page.waitForFunction(() => /MSE|RMSE|Predicted WEIGHT/.test(document.body.innerText), null, {timeout: 60_000});

    // Save as TestDemog2
    await page.locator('button[name="button-Save"], button:has-text("SAVE")').first().click();
    const dlg = page.locator('.d4-dialog [name="input-host-Name"] input');
    await dlg.waitFor({timeout: 10_000});
    await dlg.focus();
    await page.keyboard.type(MODEL_NAME_REG);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(async (name: string) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length > 0;
    }, MODEL_NAME_REG, {timeout: 30_000});
  });

  if (stepErrors.length > 0)
    throw new Error(`Soft step failures:\n${stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n')}`);
});
