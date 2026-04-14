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

test('Matrix Plot: add viewer, verify properties', async () => {
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

  await softStep('Steps 1-2: Add Matrix plot', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const mp = tv.addViewer('Matrix plot');
      await new Promise(r => setTimeout(r, 3000));
      return {type: mp?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('Matrix plot');
    expect(result.viewerCount).toBeGreaterThanOrEqual(2);
  });

  await softStep('Steps 5-7: Verify properties and viewer element', async () => {
    const result = await page!.evaluate(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const mp = viewers.find((v: any) => v.type === 'Matrix plot');
      if (!mp) return {error: 'not found'};
      const opts = mp.getOptions();
      return {
        lookKeys: Object.keys(opts.look || {}),
        viewerFound: !!document.querySelector('[name="viewer-Matrix-plot"]'),
      };
    });
    expect(result.lookKeys.length).toBeGreaterThan(0);
    expect(result.viewerFound).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
