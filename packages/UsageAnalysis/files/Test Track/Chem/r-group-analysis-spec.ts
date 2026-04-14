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

test('Chem R-Group Analysis: sar_small, MCS, trellis plot, no-core balloon', async () => {
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

  // Setup: open sar_small.csv
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/sar_small.csv');
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

  // Step 1: Verify dataset
  await softStep('Step 1: Open sar_small.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(200);
  });

  // Steps 2-5: Open R-groups dialog, click MCS, click OK, verify trellis
  await softStep('Steps 2-5: R-groups analysis with MCS', async () => {
    // Open dialog
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
      const rg = document.querySelector('[name="div-Chem---Analyze---R-Groups-Analysis..."]') as HTMLElement;
      if (rg) rg.click();
      await new Promise(r => setTimeout(r, 3000));
    });

    // Click MCS
    await page!.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const mcsBtn = Array.from(dialog!.querySelectorAll('*')).find(el =>
        el.textContent?.trim() === 'MCS' && el.children.length === 0 &&
        (el as HTMLElement).getBoundingClientRect().width > 0
      );
      if (mcsBtn) (mcsBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    // Click OK
    await page!.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const okBtn = dialog?.querySelector('[name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 15000));
    });

    // Verify trellis plot and R-group columns
    const result = await page!.evaluate(() => {
      const trellis = document.querySelectorAll('[name*="Trellis"], [name*="trellis"]').length;
      const colNames = Array.from({length: grok.shell.t.columns.length},
        (_, i) => grok.shell.t.columns.byIndex(i).name);
      const rCols = colNames.filter(n => /^R\d+$/.test(n));
      return {trellisCount: trellis, rGroupCols: rCols, totalCols: colNames.length};
    });
    expect(result.trellisCount).toBeGreaterThanOrEqual(1);
    expect(result.rGroupCols.length).toBeGreaterThanOrEqual(2);
  });

  // Step 11: Run without MCS — expect "No core was provided"
  await softStep('Step 11: Run without MCS — no core balloon', async () => {
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
      const rg = document.querySelector('[name="div-Chem---Analyze---R-Groups-Analysis..."]') as HTMLElement;
      if (rg) rg.click();
      await new Promise(r => setTimeout(r, 2000));

      // Click OK without MCS
      const dialog = document.querySelector('.d4-dialog');
      const okBtn = dialog?.querySelector('[name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 3000));
    });

    const balloonFound = await page!.evaluate(() => {
      const balloons = document.querySelectorAll('.d4-balloon-content, .d4-balloon');
      return Array.from(balloons).some(b => b.textContent?.includes('No core'));
    });
    expect(balloonFound).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
