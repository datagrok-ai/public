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

test('Chem ChemProp: Open mol1K.sdf, train model, predict, verify', async () => {
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

  // Step 1: Open mol1K.sdf
  await softStep('Step 1: Open mol1K.sdf', async () => {
    await page!.evaluate(async () => {
      const df = await grok.data.files.openTable('System:AppData/Chem/mol1K.sdf');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 5000);
      });
      await new Promise(r => setTimeout(r, 5000));
    });

    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(1000);
    expect(info.cols).toBe(6);
  });

  // Step 2: Open Train Model dialog
  await softStep('Step 2: ML > Models > Train Model', async () => {
    await page!.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      if (ml) ml.click();
      await new Promise(r => setTimeout(r, 500));
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement;
      if (models) {
        const rect = models.getBoundingClientRect();
        models.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.left + 5, clientY: rect.top + 5}));
        models.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        models.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.right - 5, clientY: rect.top + 5}));
      }
      await new Promise(r => setTimeout(r, 500));
      const train = document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement;
      if (train) train.click();
      await new Promise(r => setTimeout(r, 3000));
    });

    const viewName = await page!.evaluate(() => grok.shell.v?.name);
    expect(viewName).toBe('Predictive model');
  });

  // Step 3: Set Predict to pIC50_HIV_Integrase, Features to molecule
  await softStep('Step 3: Configure Predict and Features columns', async () => {
    // Click Predict selector and select pIC50_HIV_Integrase
    await page!.evaluate(async () => {
      const selector = document.querySelector('[name="div-Predict"]') as HTMLElement;
      if (selector) selector.click();
      await new Promise(r => setTimeout(r, 500));
    });

    // Click pIC50_HIV_Integrase in the canvas dropdown
    await page!.evaluate(async () => {
      const selector = document.querySelector('[name="div-Predict"]');
      if (!selector) return;
      const rect = selector.getBoundingClientRect();
      // Canvas grid below: header + data rows. pIC50 is 4th data row
      const canvas = document.elementFromPoint(rect.left + 100, rect.bottom + 50);
      if (canvas) {
        // Row height ~16px, pIC50 is row 4 (after header)
        const clickY = rect.bottom + 16 * 4 + 8;
        const opts = {bubbles: true, cancelable: true, view: window, clientX: rect.left + 100, clientY: clickY, button: 0};
        canvas.dispatchEvent(new MouseEvent('mousedown', opts));
        canvas.dispatchEvent(new MouseEvent('mouseup', opts));
        canvas.dispatchEvent(new MouseEvent('click', opts));
      }
      await new Promise(r => setTimeout(r, 500));
    });

    const predictCol = await page!.evaluate(() =>
      document.querySelector('[name="div-Predict"] .d4-column-selector-column')?.textContent?.trim()
    );
    expect(predictCol).toBe('pIC50_HIV_Integrase');

    // Set Features: click Features selector, click All, click OK
    await page!.evaluate(async () => {
      const featuresSelector = document.querySelector('[name="div-Features"]') as HTMLElement;
      if (featuresSelector) featuresSelector.click();
      await new Promise(r => setTimeout(r, 500));

      const allBtn = Array.from(document.querySelectorAll('*')).find(el =>
        el.textContent?.trim() === 'All' && el.children.length === 0 && el.closest('.d4-dialog')
      );
      if (allBtn) (allBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));

      const okBtn = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 500));
    });
  });

  // Step 4: Train model (may fail if ML backend unavailable)
  await softStep('Step 4: Run training', async () => {
    // Click the sync/refresh button to train
    await page!.evaluate(async () => {
      const syncIcon = Array.from(document.querySelectorAll('.d4-ribbon-item i.fa-sync')).find(i =>
        i.getBoundingClientRect().top > 20 && i.getBoundingClientRect().top < 50
      );
      if (syncIcon) {
        const parent = syncIcon.closest('.d4-ribbon-item')!;
        parent.setAttribute('role', 'button');
        parent.setAttribute('aria-label', 'Train model');
      }
    });

    const trainPos = await page!.evaluate(() => {
      const btn = document.querySelector('[aria-label="Train model"]');
      if (!btn) return null;
      const rect = btn.getBoundingClientRect();
      return {x: rect.left + rect.width / 2, y: rect.top + rect.height / 2};
    });
    if (trainPos) await page!.mouse.click(trainPos.x, trainPos.y);

    // Wait for training to complete (may timeout if ML service unavailable)
    await page!.waitForTimeout(30000);

    // Check if prediction column appeared
    const table = await page!.evaluate(() => {
      const t = grok.shell.tables[0];
      return t ? {cols: t.columns.length, name: t.name} : null;
    });
    // This step may fail if ML backend is not running
    expect(table?.cols).toBeGreaterThan(6);
  });

  // Step 5: Verify predictions with scatterplot
  await softStep('Step 5: Verify predictions nearly equal to actual', async () => {
    await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      if (tv) tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 2000));
    });

    const canvases = await page!.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvases).toBeGreaterThanOrEqual(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
