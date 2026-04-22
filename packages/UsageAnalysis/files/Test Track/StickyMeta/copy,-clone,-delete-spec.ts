import {test, expect, chromium} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

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

test('StickyMeta Copy/Clone/Delete: clone table, verify metadata preserved', async () => {
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

  // Setup: open SPGI.csv
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

  // Step 1: Verify SPGI and Sticky meta schema
  await softStep('Step 1: Open SPGI, verify Sticky meta schema', async () => {
    await page!.evaluate(async () => {
      const col = grok.shell.t.col('Structure');
      grok.shell.o = col;
      grok.shell.windows.showProperties = true;
      await new Promise(r => setTimeout(r, 3000));
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const stickyPane = Array.from(panes).find(p =>
        p.textContent?.trim() === 'Sticky meta' && p.getBoundingClientRect().left > 400
      );
      if (stickyPane) (stickyPane as HTMLElement).click();
      await new Promise(r => setTimeout(r, 2000));
    });

    const hasSchema = await page!.evaluate(() => {
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const stickyPane = Array.from(panes).find(p =>
        p.textContent?.trim() === 'Sticky meta' && p.getBoundingClientRect().left > 400
      );
      const content = stickyPane?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      return content?.textContent?.includes('TestSchema1') ?? false;
    });
    expect(hasSchema).toBe(true);
  });

  // Step 2: Clone table and verify metadata preserved
  await softStep('Step 2: Clone table, verify metadata preserved', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      const cloned = df.clone();
      const tv2 = grok.shell.addTableView(cloned);
      await new Promise(r => setTimeout(r, 2000));

      // Check cloned view has Sticky meta
      const col = cloned.col('Structure');
      grok.shell.o = col;
      grok.shell.windows.showProperties = true;
      await new Promise(r => setTimeout(r, 3000));

      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const stickyPane = Array.from(panes).find(p =>
        p.textContent?.trim() === 'Sticky meta' && p.getBoundingClientRect().left > 400
      );
      if (stickyPane) (stickyPane as HTMLElement).click();
      await new Promise(r => setTimeout(r, 2000));
      const content = stickyPane?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');

      return {
        clonedCols: cloned.columns.length,
        clonedRows: cloned.rowCount,
        hasSchema: content?.textContent?.includes('TestSchema1') ?? false,
      };
    });
    expect(result.clonedRows).toBe(3624);
    expect(result.hasSchema).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
