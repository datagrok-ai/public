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

test('Chem Info Panels: Open smiles.csv, check column and molecule panels', async () => {
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

  // Setup: open smiles.csv with chem wait
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
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

  // Step 1-2: Verify dataset
  await softStep('Step 1-2: Open smiles.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(1000);
    expect(info.cols).toBe(20);
  });

  // Step 3: Click canonical_smiles column header
  await softStep('Step 3: Select canonical_smiles column', async () => {
    await page!.evaluate(async () => {
      grok.shell.windows.showProperties = true;
      const col = grok.shell.t.col('canonical_smiles');
      grok.shell.o = col;
      await new Promise(r => setTimeout(r, 3000));
    });

    const panels = await page!.evaluate(() => {
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      return Array.from(panes).filter(p => p.getBoundingClientRect().left > 400)
        .map(p => p.textContent?.trim());
    });
    expect(panels.length).toBeGreaterThanOrEqual(5);
    expect(panels).toContain('Chemistry');
  });

  // Step 4: Expand all column info panels — no errors
  await softStep('Step 4: Expand all column info panels', async () => {
    const results = await page!.evaluate(async () => {
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const rightPanes = Array.from(panes).filter(p => p.getBoundingClientRect().left > 400);
      const res: {name: string; hasError: boolean}[] = [];
      for (const pane of rightPanes) {
        const name = pane.textContent?.trim() ?? '';
        (pane as HTMLElement).click();
        await new Promise(r => setTimeout(r, 1000));
        const container = pane.closest('.d4-accordion-pane');
        const content = container?.querySelector('.d4-accordion-pane-content');
        const hasError = !!content?.textContent?.toLowerCase()?.includes('error');
        res.push({name, hasError});
      }
      return res;
    });
    const errors = results.filter(r => r.hasError);
    expect(errors.length).toBe(0);
  });

  // Step 5: Click first molecule cell, expand panels — no errors
  await softStep('Step 5: Click molecule cell, expand all panels', async () => {
    await page!.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      grid.currentCell = grid.cell('canonical_smiles', 0);
      await new Promise(r => setTimeout(r, 5000));
    });

    const results = await page!.evaluate(async () => {
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      const rightPanes = Array.from(panes).filter(p => p.getBoundingClientRect().left > 400);
      const res: {name: string; hasError: boolean}[] = [];
      for (const pane of rightPanes) {
        const name = pane.textContent?.trim() ?? '';
        (pane as HTMLElement).click();
        await new Promise(r => setTimeout(r, 500));
        const container = pane.closest('.d4-accordion-pane');
        const content = container?.querySelector('.d4-accordion-pane-content');
        const hasError = !!content?.textContent?.toLowerCase()?.includes('error');
        res.push({name, hasError});
      }
      return res;
    });
    const errors = results.filter(r => r.hasError);
    expect(errors.length).toBe(0);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
