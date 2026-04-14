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

test('Chem MMP: Open mmp_demo, run analysis, check tabs', async () => {
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

  // Setup: open mmp_demo.csv
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/mmp_demo.csv');
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
  await softStep('Step 1: Open mmp_demo.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(20267);
  });

  // Step 2: Open MMP dialog
  await softStep('Step 2: Chem > Analyze > Matched Molecular Pairs', async () => {
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
      const mmp = document.querySelector('[name="div-Chem---Analyze---Matched-Molecular-Pairs..."]') as HTMLElement;
      if (mmp) mmp.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    const dialogOpen = await page!.evaluate(() => !!document.querySelector('.d4-dialog'));
    expect(dialogOpen).toBe(true);
  });

  // Step 3: Select activities and run
  await softStep('Step 3: Select activities, press OK', async () => {
    // Click Activities selector, select All, OK
    await page!.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog')!;
      const labels = dialog.querySelectorAll('label.ui-input-label span, .ui-label');
      for (const label of labels) {
        if (label.textContent?.trim() === 'Activities') {
          const sel = label.closest('label')?.nextElementSibling;
          if (sel) (sel as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
    });

    await page!.evaluate(async () => {
      const dialogs = document.querySelectorAll('.d4-dialog');
      const inner = dialogs[dialogs.length - 1];
      const allBtn = Array.from(inner.querySelectorAll('*')).find(el =>
        el.textContent?.trim() === 'All' && el.children.length === 0
      );
      if (allBtn) (allBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
      const okBtn = inner.querySelector('[name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 1000));
    });

    // Click OK on main dialog
    await page!.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const okBtn = dialog?.querySelector('[name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 30000));
    });

    // Check for new views (MMP creates tabs)
    const views = await page!.evaluate(() =>
      Array.from(grok.shell.views).map(v => v.name)
    );
    // This may fail if MMP has a bug
    expect(views.length).toBeGreaterThan(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
