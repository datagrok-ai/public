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

test('Tile Viewer: add viewer, calculated column, table switch', async () => {
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

  // Setup: open SPGI.csv with chem wait
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 5000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
  });

  // Step 1: Add Tile viewer
  await softStep('Step 1: Add Tile viewer', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const tileViewer = tv.addViewer('Tile viewer');
      await new Promise(r => setTimeout(r, 3000));
      return {type: tileViewer?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('Tile Viewer');
    expect(result.viewerCount).toBeGreaterThanOrEqual(2);
  });

  // Step 3: Open demog, switch back to SPGI, verify tile viewer persists
  await softStep('Step 3: Table switch — tile viewer persists', async () => {
    const result = await page!.evaluate(async () => {
      const demog = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(demog);
      await new Promise(r => setTimeout(r, 2000));

      // Switch back to SPGI
      const views = Array.from(grok.shell.views);
      const spgiView = views.find(v => v.name === 'Table');
      if (spgiView) grok.shell.v = spgiView;
      await new Promise(r => setTimeout(r, 1000));

      const viewers = Array.from(grok.shell.tv.viewers);
      return {tilePresent: viewers.some(v => v.type === 'Tile Viewer')};
    });
    expect(result.tilePresent).toBe(true);
  });

  // Step 5: Create calculated column
  await softStep('Step 5: Create calculated column cc', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      df.columns.addNewCalculated('cc', '${Average Mass} + 5');
      await new Promise(r => setTimeout(r, 2000));
      return {hasCC: df.col('cc') !== null, value: df.col('cc')?.get(0)};
    });
    expect(result.hasCC).toBe(true);
    expect(result.value).toBeGreaterThan(0);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
