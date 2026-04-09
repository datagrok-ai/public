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

test('Network Diagram: add viewer, set nodes, verify properties', async () => {
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

    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });

  await softStep('Steps 1-2: Add Network diagram', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const nd = tv.addViewer('Network diagram');
      await new Promise(r => setTimeout(r, 3000));
      return {type: nd?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('Network diagram');
  });

  await softStep('Steps 4-5: Modify node columns', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const nd = viewers.find((v: any) => v.type === 'Network diagram');
      if (!nd) return {error: 'not found'};
      nd.setOptions({node1ColumnName: 'SEX', node2ColumnName: 'RACE'});
      await new Promise(r => setTimeout(r, 2000));
      const opts = nd.getOptions();
      return {
        node1: opts.look?.node1ColumnName,
        node2: opts.look?.node2ColumnName,
        viewerFound: !!document.querySelector('[name="viewer-Network-diagram"]'),
      };
    });
    expect(result.node1).toBe('SEX');
    expect(result.node2).toBe('RACE');
    expect(result.viewerFound).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
