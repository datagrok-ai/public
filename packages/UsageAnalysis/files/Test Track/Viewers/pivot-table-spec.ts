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

test('Pivot Table: add viewer, modify properties, verify no errors', async () => {
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

  // Setup: open demog.csv
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

  // Step 1-2: Add Pivot table viewer
  await softStep('Steps 1-2: Open demog, add Pivot table', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const pv = tv.addViewer('Pivot table');
      await new Promise(r => setTimeout(r, 3000));
      return {type: pv?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('Pivot table');
    expect(result.viewerCount).toBeGreaterThanOrEqual(2);
  });

  // Step 4: Modify pivot properties
  await softStep('Step 4: Modify pivot/groupBy/aggregate properties', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const pv = viewers.find((v: any) => v.type === 'Pivot table');
      if (!pv) return {error: 'not found'};

      pv.setOptions({
        pivotColumnNames: ['RACE'],
        groupByColumnNames: ['SEX'],
        aggregateColumnNames: ['HEIGHT', 'WEIGHT'],
      });
      await new Promise(r => setTimeout(r, 3000));

      const opts = pv.getOptions();
      return {
        pivot: opts.look?.pivotColumnNames,
        groupBy: opts.look?.groupByColumnNames,
        aggregate: opts.look?.aggregateColumnNames,
      };
    });
    expect(result.pivot).toEqual(['RACE']);
    expect(result.groupBy).toEqual(['SEX']);
    expect(result.aggregate).toEqual(['HEIGHT', 'WEIGHT']);
  });

  // Step 5: Verify viewer element exists in DOM
  await softStep('Step 5: Verify Pivot table viewer element', async () => {
    const found = await page!.evaluate(() =>
      !!document.querySelector('[name="viewer-Pivot-table"]')
    );
    expect(found).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
