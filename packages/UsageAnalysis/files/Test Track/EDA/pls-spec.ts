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

test('PLS', async () => {
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

  // Open dataset
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
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1: Open dataset verified', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.tv.dataFrame.rowCount,
      cols: grok.shell.tv.dataFrame.columns.length
    }));
    expect(info.rows).toBe(30);
    expect(info.cols).toBe(17);
  });

  await softStep('Step 2-3: Open PLS dialog via menu and click RUN', async () => {
    // Click ML menu
    await page!.locator('[name="div-ML"]').click();
    await page!.waitForTimeout(300);

    // Hover Analyze submenu
    const analyzeItem = page!.locator('[name="div-ML---Analyze"]');
    await analyzeItem.dispatchEvent('mouseenter');
    await analyzeItem.hover();
    await page!.waitForTimeout(500);

    // Click PLS
    await page!.locator('[name="div-ML---Analyze---PLS..."]').click();

    // Wait for dialog
    await page!.locator('.d4-dialog').waitFor({timeout: 5000});
    await expect(page!.locator('.d4-dialog')).toContainText('PLS');

    // Click RUN
    const runBtn = page!.locator('[name="button-RUN"], [name="button-Run"]');
    await runBtn.click();
    await page!.waitForTimeout(2000);
  });

  await softStep('Step 4: Verify PLS1, PLS2, PLS3 columns added', async () => {
    const result = await page!.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const colNames: string[] = [];
      for (let i = 0; i < df.columns.length; i++)
        colNames.push(df.columns.byIndex(i).name);
      return {
        totalCols: df.columns.length,
        plsCols: colNames.filter(n => n.startsWith('PLS'))
      };
    });
    expect(result.totalCols).toBe(20);
    expect(result.plsCols).toEqual(['PLS1', 'PLS2', 'PLS3']);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
