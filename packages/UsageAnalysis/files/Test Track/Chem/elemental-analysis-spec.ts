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

test('Chem Elemental Analysis: Open smiles.csv, run analysis with all options', async () => {
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

  // Step 1: Verify dataset loaded
  await softStep('Step 1: Open smiles.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(1000);
    expect(info.cols).toBe(20);
  });

  // Step 2: Open Chem > Analyze > Elemental Analysis dialog
  await softStep('Step 2: Open Elemental Analysis dialog', async () => {
    await page!.evaluate(async () => {
      const chem = document.querySelector('[name="div-Chem"]') as HTMLElement;
      if (chem) chem.click();
      await new Promise(r => setTimeout(r, 500));

      const analyze = document.querySelector('[name="div-Chem---Analyze"]') as HTMLElement;
      if (analyze) {
        const rect = analyze.getBoundingClientRect();
        analyze.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.left + 5, clientY: rect.top + 5}));
        analyze.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        analyze.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.right - 5, clientY: rect.top + 5}));
      }
      await new Promise(r => setTimeout(r, 800));

      const ea = document.querySelector('[name="div-Chem---Analyze---Elemental-Analysis..."]') as HTMLElement;
      if (ea) ea.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    const dialogOpen = await page!.evaluate(() => !!document.querySelector('.d4-dialog'));
    expect(dialogOpen).toBe(true);
  });

  // Step 3: Turn all checkboxes on
  await softStep('Step 3: Enable all checkboxes', async () => {
    const result = await page!.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      if (!dialog) return {error: 'no dialog'};
      const checkboxes = dialog.querySelectorAll('input[type="checkbox"]');
      for (const cb of checkboxes) {
        if (!(cb as HTMLInputElement).checked) (cb as HTMLElement).click();
      }
      return {checked: checkboxes.length};
    });
    expect(result.checked).toBeGreaterThanOrEqual(2);
  });

  // Step 4: Click OK and verify results
  await softStep('Step 4: Click OK and verify element columns added', async () => {
    await page!.evaluate(async () => {
      const okBtn = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 10000));
    });

    const info = await page!.evaluate(() => {
      const t = grok.shell.t;
      const colNames = t ? Array.from({length: t.columns.length}, (_, i) => t.columns.byIndex(i).name) : [];
      return {colCount: t?.columns?.length, hasH: colNames.includes('H'), hasC: colNames.includes('C')};
    });
    expect(info.colCount).toBeGreaterThan(20);
    expect(info.hasH).toBe(true);
    expect(info.hasC).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
