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

test('Form Viewer: add viewer, multiple tables, verify persistence', async () => {
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

  // Setup: open SPGI with chem wait
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

  // Step 1: Add Form viewer
  await softStep('Step 1: Add Form viewer to SPGI', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const fv = tv.addViewer('Form');
      await new Promise(r => setTimeout(r, 3000));
      return {type: fv?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('Form');
    expect(result.viewerCount).toBeGreaterThanOrEqual(2);
  });

  // Step 1.1: Open additional tables, verify Form viewer persists
  await softStep('Step 1.1: Open SPGI-linked tables, verify persistence', async () => {
    const result = await page!.evaluate(async () => {
      const df1 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked1.csv');
      grok.shell.addTableView(df1);
      await new Promise(r => setTimeout(r, 1000));
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked2.csv');
      grok.shell.addTableView(df2);
      await new Promise(r => setTimeout(r, 1000));

      // Switch back to SPGI
      const views = Array.from(grok.shell.views);
      const spgiView = views.find(v => v.name === 'Table');
      if (spgiView) grok.shell.v = spgiView;
      await new Promise(r => setTimeout(r, 1000));

      const viewers = Array.from(grok.shell.tv.viewers);
      return {
        tableCount: Array.from(grok.shell.tables).length,
        formPresent: viewers.some(v => v.type === 'Form'),
      };
    });
    expect(result.tableCount).toBe(3);
    expect(result.formPresent).toBe(true);
  });

  // Step 2: Verify Form viewer has options
  await softStep('Step 2: Verify Form viewer options', async () => {
    const result = await page!.evaluate(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const fv = viewers.find((v: any) => v.type === 'Form');
      if (!fv) return {error: 'not found'};
      const opts = fv.getOptions();
      return {lookKeys: Object.keys(opts.look || {})};
    });
    expect(result.lookKeys).toContain('sketchState');
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
