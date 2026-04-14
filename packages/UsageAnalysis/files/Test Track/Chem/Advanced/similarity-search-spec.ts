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

test('Chem Similarity Search: Open, run search, modify properties', async () => {
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

  // Step 1: Verify dataset
  await softStep('Step 1: Open smiles.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
    }));
    expect(info.rows).toBe(1000);
  });

  // Step 2: Run Similarity Search
  await softStep('Step 2: Chem > Search > Similarity Search', async () => {
    await page!.evaluate(async () => {
      const chem = document.querySelector('[name="div-Chem"]') as HTMLElement;
      if (chem) chem.click();
      await new Promise(r => setTimeout(r, 500));
      const search = document.querySelector('[name="div-Chem---Search"]') as HTMLElement;
      if (search) {
        const rect = search.getBoundingClientRect();
        search.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.left + 5, clientY: rect.top + 5}));
        search.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        search.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.right - 5, clientY: rect.top + 5}));
      }
      await new Promise(r => setTimeout(r, 800));
      const sim = document.querySelector('[name="div-Chem---Search---Similarity-Search..."]') as HTMLElement;
      if (sim) sim.click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const viewerFound = await page!.evaluate(() =>
      !!document.querySelector('[name="viewer-Chem-Similarity-Search"]')
    );
    expect(viewerFound).toBe(true);
  });

  // Step 3-4: Test property modifications
  await softStep('Steps 3-4: Modify properties without errors', async () => {
    const results = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const simViewer = viewers.find((v: any) => v.type?.includes('Similarity'));
      if (!simViewer) return {error: 'viewer not found'};

      const res: Record<string, boolean> = {};

      // Test fingerprint change
      simViewer.setOptions({fingerprint: 'Pattern'});
      await new Promise(r => setTimeout(r, 2000));
      res.fingerprint = simViewer.getOptions().look?.fingerprint === 'Pattern';

      // Test limit change
      simViewer.setOptions({limit: 5});
      await new Promise(r => setTimeout(r, 2000));
      res.limit = simViewer.getOptions().look?.limit === 5;

      // Test distance metric change
      simViewer.setOptions({distanceMetric: 'Dice'});
      await new Promise(r => setTimeout(r, 2000));
      res.metric = simViewer.getOptions().look?.distanceMetric === 'Dice';

      // Test cutoff = 1
      simViewer.setOptions({cutoff: 1.0});
      await new Promise(r => setTimeout(r, 2000));
      res.cutoff = simViewer.getOptions().look?.cutoff === 1;

      // Reset
      simViewer.setOptions({fingerprint: 'Morgan', limit: 12, distanceMetric: 'Tanimoto', cutoff: 0.01});

      return res;
    });
    expect(results.fingerprint).toBe(true);
    expect(results.limit).toBe(true);
    expect(results.metric).toBe(true);
    expect(results.cutoff).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
