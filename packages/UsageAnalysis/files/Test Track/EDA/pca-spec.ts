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

test('PCA', async () => {
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

  await softStep('Step 2-4: Run PCA with 3 components, verify PC1/PC2/PC3', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const numCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.type !== 'string') numCols.push(c.name);
      }
      await grok.functions.call('Eda:PCA', {
        table: df, features: numCols, components: 3, center: false, scale: false
      });
      await new Promise(r => setTimeout(r, 1000));
      const colNames: string[] = [];
      for (let i = 0; i < df.columns.length; i++)
        colNames.push(df.columns.byIndex(i).name);
      return { totalCols: df.columns.length, pcCols: colNames.filter(n => n.startsWith('PC')) };
    });
    expect(result.totalCols).toBe(20);
    expect(result.pcCols).toEqual(['PC1', 'PC2', 'PC3']);
  });

  await softStep('Step 5: Repeat PCA with Center+Scale, verify PC1(2)/PC2(2)/PC3(2)', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const numCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.type !== 'string' && !c.name.startsWith('PC'))
          numCols.push(c.name);
      }
      await grok.functions.call('Eda:PCA', {
        table: df, features: numCols, components: 3, center: true, scale: true
      });
      await new Promise(r => setTimeout(r, 1000));
      const colNames: string[] = [];
      for (let i = 0; i < df.columns.length; i++)
        colNames.push(df.columns.byIndex(i).name);
      return { totalCols: df.columns.length, pcCols: colNames.filter(n => n.startsWith('PC')) };
    });
    expect(result.totalCols).toBe(23);
    expect(result.pcCols).toEqual(['PC1', 'PC2', 'PC3', 'PC1 (2)', 'PC2 (2)', 'PC3 (2)']);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
