import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('ANOVA', async () => {
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

  // Close any open dialogs first, then open dataset
  await page.evaluate(async (path) => {
    // Close any existing dialogs
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
    expect(rowCount).toBe(5850);
  });

  await softStep('Step 2: Open ANOVA dialog via ML > Analyze > ANOVA', async () => {
    // Click ML menu
    await page!.locator('[name="div-ML"]').click();
    await page!.waitForTimeout(300);

    // Hover over Analyze submenu to expand it
    const analyzeItem = page!.locator('[name="div-ML---Analyze"]');
    await analyzeItem.dispatchEvent('mouseenter');
    await analyzeItem.hover();
    await page!.waitForTimeout(500);

    // Click ANOVA
    await page!.locator('[name="div-ML---Analyze---ANOVA..."]').click();

    // Wait for ANOVA dialog
    await page!.locator('.d4-dialog').waitFor({timeout: 5000});
    await expect(page!.locator('.d4-dialog')).toContainText('ANOVA');
  });

  await softStep('Step 3: Click RUN and verify results', async () => {
    // Click RUN button (name may be "button-RUN" or "button-Run")
    const runBtn = page!.locator('[name="button-RUN"], [name="button-Run"]');
    await runBtn.click();

    // Wait for results: box plot and Analysis/F-test tabs
    await page!.locator('[name="viewer-Box-plot"]').waitFor({timeout: 15000});

    // Verify Analysis tab exists
    const analysisTab = page!.locator('.d4-tab-header[name="Analysis"]');
    await expect(analysisTab).toBeVisible({timeout: 5000});

    // Verify F-test tab exists
    const ftestTab = page!.locator('.d4-tab-header[name="F-test"]');
    await expect(ftestTab).toBeVisible({timeout: 5000});

    // Verify box plot is visible
    const boxPlot = page!.locator('[name="viewer-Box-plot"]');
    await expect(boxPlot).toBeVisible();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
