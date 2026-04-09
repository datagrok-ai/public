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

test('StickyMeta Add and Edit: Open SPGI, check schema, add sticky columns', async () => {
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

  const initialCols = await page.evaluate(() => grok.shell.t?.columns?.length);

  // Step 1: Verify dataset
  await softStep('Step 1: Open SPGI.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(3624);
  });

  // Step 2: Select Structure cell, check Sticky meta panel
  await softStep('Step 2: Select Structure cell, find Sticky meta panel', async () => {
    await page!.evaluate(async () => {
      // Select the Structure column to show column properties
      const col = grok.shell.t.col('Structure');
      grok.shell.o = col;
      grok.shell.windows.showProperties = true;
      await new Promise(r => setTimeout(r, 5000));

      // Expand Sticky meta pane
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const stickyPane = Array.from(panes).find(p =>
        p.textContent?.trim() === 'Sticky meta' && p.getBoundingClientRect().left > 400
      );
      if (stickyPane) {
        stickyPane.scrollIntoView();
        (stickyPane as HTMLElement).click();
      }
      await new Promise(r => setTimeout(r, 2000));
    });

    const contentText = await page!.evaluate(() => {
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const stickyPane = Array.from(panes).find(p =>
        p.textContent?.trim() === 'Sticky meta' && p.getBoundingClientRect().left > 400
      );
      const content = stickyPane?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      return content?.textContent?.trim() ?? '';
    });
    expect(contentText).toContain('TestSchema1');
    expect(contentText).toContain('rating');
  });

  // Step 3: Add sticky columns by clicking + buttons
  await softStep('Step 3: Add sticky columns to grid', async () => {
    await page!.evaluate(async () => {
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const stickyPane = Array.from(panes).find(p =>
        p.textContent?.trim() === 'Sticky meta' && p.getBoundingClientRect().left > 400
      );
      const content = stickyPane?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      const plusIcons = Array.from(content?.querySelectorAll('i.fa-plus') ?? [])
        .filter(i => i.getBoundingClientRect().width > 0);
      for (const icon of plusIcons) {
        (icon as HTMLElement).click();
        await new Promise(r => setTimeout(r, 500));
      }
      await new Promise(r => setTimeout(r, 2000));
    });

    const newCols = await page!.evaluate(() => grok.shell.t?.columns?.length);
    expect(newCols).toBeGreaterThan(initialCols!);
  });

  // Step 4: Verify sticky columns were added (columns increased)
  await softStep('Step 4: Verify sticky columns added', async () => {
    const newCols = await page!.evaluate(() => grok.shell.t?.columns?.length);
    expect(newCols).toBeGreaterThan(initialCols!);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
