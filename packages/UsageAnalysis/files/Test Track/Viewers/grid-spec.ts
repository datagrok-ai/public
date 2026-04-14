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

test('Grid Viewer: sort, select, context panel, filter', async () => {
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

  // Section 1: Grid basics
  await softStep('Section 1: Sort and select', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      const grid = grok.shell.tv.grid;

      grid.sort(['Average Mass']);
      await new Promise(r => setTimeout(r, 1000));
      grid.sort(['Average Mass'], [false]);
      await new Promise(r => setTimeout(r, 1000));

      df.selection.set(0, true);
      df.selection.set(1, true);
      df.selection.set(2, true);
      const selectedCount = df.selection.trueCount;

      df.selection.setAll(false);
      grid.sort(['Id']);

      return {selectedCount, rows: df.rowCount};
    });
    expect(result.selectedCount).toBe(3);
    expect(result.rows).toBe(3624);
  });

  // Section 4: Context Panel for column
  await softStep('Section 4: Column Context Panel', async () => {
    const panels = await page!.evaluate(async () => {
      const col = grok.shell.t.col('Average Mass');
      grok.shell.o = col;
      grok.shell.windows.showProperties = true;
      await new Promise(r => setTimeout(r, 3000));

      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      return Array.from(panes).filter(p => p.getBoundingClientRect().left > 400)
        .map(p => p.textContent?.trim());
    });
    expect(panels).toContain('Details');
    expect(panels).toContain('Filter');
    expect(panels).toContain('Colors');
    expect(panels.length).toBeGreaterThanOrEqual(10);
  });

  // Section 8: Filtering
  await softStep('Section 8: Filtering', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 3000));

      const totalBefore = df.filter.trueCount;
      const col = df.col('Average Mass');
      for (let i = 0; i < df.rowCount; i++)
        df.filter.set(i, col.get(i) > 300, false);
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 1000));
      const totalAfter = df.filter.trueCount;

      df.filter.setAll(true);
      return {totalBefore, totalAfter};
    });
    expect(result.totalAfter).toBeLessThan(result.totalBefore);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
