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

test('Grid Viewer: add second grid, modify properties, switch table', async () => {
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

    // Also open SPGI-linked1
    const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked1.csv');
    grok.shell.addTableView(df2);
    await new Promise(r => setTimeout(r, 2000));

    // Switch back to SPGI
    const views = Array.from(grok.shell.views);
    const spgiView = views.find(v => v.name === 'Table');
    if (spgiView) grok.shell.v = spgiView;
  });

  // Steps 1-4: Add second Grid viewer
  await softStep('Steps 1-4: Open datasets, add Grid viewer', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const gv = tv.addViewer('Grid');
      await new Promise(r => setTimeout(r, 2000));
      return {
        type: gv?.type,
        viewerCount: Array.from(tv.viewers).length,
        tableCount: Array.from(grok.shell.tables).length,
      };
    });
    expect(result.type).toBe('Grid');
    expect(result.viewerCount).toBeGreaterThanOrEqual(2);
    expect(result.tableCount).toBe(2);
  });

  // Steps 9-10: Modify properties (row height)
  await softStep('Steps 9-10: Modify row height', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const grids = viewers.filter((v: any) => v.type === 'Grid');
      const secondGrid = grids[grids.length - 1];
      secondGrid.setOptions({rowHeight: 40});
      await new Promise(r => setTimeout(r, 1000));
      const opts = secondGrid.getOptions();
      return {rowHeight: opts.look?.rowHeight};
    });
    expect(result.rowHeight).toBe(40);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
