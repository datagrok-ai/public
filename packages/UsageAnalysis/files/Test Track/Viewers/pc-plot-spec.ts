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

test('PC Plot: add viewer, set columns, verify properties', async () => {
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
  });

  // Section 1: Add PC Plot
  await softStep('Section 1: Add PC Plot', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const pc = tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 3000));
      return {type: pc?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('PC Plot');
    expect(result.viewerCount).toBeGreaterThanOrEqual(2);
  });

  // Section 3: Set columns and verify
  await softStep('Section 3: Set columns via properties', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const pc = viewers.find((v: any) => v.type === 'PC Plot');
      if (!pc) return {error: 'not found'};
      pc.setOptions({columnNames: ['Average Mass', 'TPSA', 'CLogP', 'CAST Idea ID']});
      await new Promise(r => setTimeout(r, 2000));
      const opts = pc.getOptions();
      return {columns: opts.look?.columnNames};
    });
    expect(result.columns).toEqual(['Average Mass', 'TPSA', 'CLogP', 'CAST Idea ID']);
  });

  // Verify viewer element
  await softStep('Verify PC Plot viewer element', async () => {
    const found = await page!.evaluate(() =>
      !!document.querySelector('[name="viewer-PC-Plot"]')
    );
    expect(found).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
