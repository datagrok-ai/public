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

test('Row Source: all modes on scatter plot with SPGI', async () => {
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

    tv.addViewer('Scatter plot');
    await new Promise(r => setTimeout(r, 2000));
  });

  // Test all row source modes
  await softStep('Test row source modes', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      const viewers = Array.from(grok.shell.tv.viewers);
      const sp = viewers.find((v: any) => v.type === 'Scatter plot');
      if (!sp) return {error: 'not found'};

      const tested: Record<string, string> = {};

      // Filter
      const col = df.col('CAST Idea ID');
      for (let i = 0; i < df.rowCount; i++)
        df.filter.set(i, col.get(i) < 636500, false);
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 500));

      for (const mode of ['All', 'Selected', 'SelectedOrCurrent', 'FilteredSelected', 'Filtered']) {
        sp.setOptions({rowSource: mode});
        await new Promise(r => setTimeout(r, 300));
        tested[mode] = sp.getOptions().look?.rowSource;
      }

      df.filter.setAll(true);
      return tested;
    });
    expect(result['All']).toBe('All');
    expect(result['Selected']).toBe('Selected');
    expect(result['SelectedOrCurrent']).toBe('SelectedOrCurrent');
    expect(result['FilteredSelected']).toBe('FilteredSelected');
    // 'Filtered' is the default — may not appear explicitly in options
    expect(result['Filtered'] === 'Filtered' || result['Filtered'] === undefined).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
