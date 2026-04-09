import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/cars.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Multivariate Analysis', async () => {
  // Connect to existing Chrome via CDP and reuse the Datagrok page
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

  // Close any open dialogs, then open dataset
  await page.evaluate(async (path) => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Test steps

  await softStep('Step 1: Open dataset verified', async () => {
    const rowCount = await page!.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(rowCount).toBe(30);
  });

  await softStep('Step 2: Open Multivariate Analysis and run', async () => {
    // Click ML menu
    await page!.locator('[name="div-ML"]').click();
    await page!.waitForTimeout(300);

    // Hover over Analyze submenu to expand it
    const analyzeItem = page!.locator('[name="div-ML---Analyze"]');
    await analyzeItem.dispatchEvent('mouseenter');
    await analyzeItem.hover();
    await page!.waitForTimeout(500);

    // Click Multivariate Analysis
    await page!.locator('[name="div-ML---Analyze---Multivariate-Analysis..."]').click();

    // Wait for dialog
    await page!.locator('.d4-dialog').waitFor({timeout: 5000});
    await expect(page!.locator('.d4-dialog')).toContainText('Multivariate Analysis');

    // Click RUN
    const runBtn = page!.locator('[name="button-RUN"], [name="button-Run"]');
    await runBtn.click();

    // Wait for viewers to appear
    await page!.locator('[name="viewer-Scatter-plot"]').first().waitFor({timeout: 15000});
    await page!.waitForTimeout(1000);

    // Verify all expected viewers exist
    const viewerInfo = await page!.evaluate(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      return {
        count: viewers.length,
        types: viewers.map(v => v.type)
      };
    });
    expect(viewerInfo.count).toBeGreaterThanOrEqual(5);
  });

  await softStep('Step 3: Check interactivity', async () => {
    // Select a row and verify it propagates
    const result = await page!.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.set(0, true);
      return {
        selCount: df.selection.trueCount,
        colCount: df.columns.length  // should be 24 after MVA (17 original + 7 PLS cols)
      };
    });
    expect(result.selCount).toBe(1);
    expect(result.colCount).toBeGreaterThan(17);

    // Clear selection
    await page!.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
